"""A small HTTP server for remote RefgetStore tests.

Serves a folder at http://127.0.0.1:<port>/<prefix>/ with byte ranges,
404 (or 403 with --listing-denied) for missing files, fault injection and a
JSON-lines request log. Run it in a subprocess; it prints "port N" first.

    python http_fixture.py --directory D [--prefix s] [--log LOG]
        [--listing-denied] [--fail-first N --fail-status S --match SUBSTR]
"""

import argparse
import json
import os
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer


def make_handler(args):
    lock = threading.Lock()
    state = {"failures_left": args.fail_first}
    root = os.path.abspath(args.directory)
    prefix = "/" + args.prefix.strip("/") + "/"

    class Handler(SimpleHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, *a):
            pass

        def _log(self, status):
            if args.log:
                with lock, open(args.log, "a") as f:
                    f.write(json.dumps({
                        "path": self.path,
                        "range": self.headers.get("Range"),
                        "status": status,
                    }) + "\n")

        def _send_empty(self, status, extra=None):
            self.send_response(status)
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.send_header("Content-Length", "0")
            self.end_headers()
            self._log(status)

        def do_GET(self):
            path = self.path.split("?", 1)[0]
            if args.match and args.match in path:
                with lock:
                    fail = state["failures_left"] > 0
                    if fail:
                        state["failures_left"] -= 1
                if fail:
                    return self._send_empty(args.fail_status)
            missing = 403 if args.listing_denied else 404
            if path.rstrip("/") + "/" == prefix:
                if args.listing_denied:
                    return self._send_empty(403)
                body = b"<html>listing</html>"
                self.send_response(200)
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return self._log(200)
            if not path.startswith(prefix) or ".." in path:
                return self._send_empty(missing)
            full = os.path.join(root, path[len(prefix):])
            if not os.path.isfile(full):
                return self._send_empty(missing)
            with open(full, "rb") as f:
                data = f.read()
            status = 200
            headers = {"Accept-Ranges": "bytes", "ETag": '"%d-%d"' % (len(data), hash(data) & 0xFFFFFFFF)}
            rng = self.headers.get("Range")
            if rng and rng.startswith("bytes="):
                start_s, _, end_s = rng[len("bytes="):].partition("-")
                start = int(start_s or 0)
                end = int(end_s) if end_s else len(data) - 1
                if start >= len(data):
                    return self._send_empty(416, {"Content-Range": "bytes */%d" % len(data)})
                end = min(end, len(data) - 1)
                headers["Content-Range"] = "bytes %d-%d/%d" % (start, end, len(data))
                data = data[start:end + 1]
                status = 206
            self.send_response(status)
            for k, v in headers.items():
                self.send_header(k, v)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            self._log(status)

    return Handler


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--directory", required=True)
    p.add_argument("--prefix", default="s")
    p.add_argument("--log")
    p.add_argument("--listing-denied", action="store_true")
    p.add_argument("--fail-first", type=int, default=0)
    p.add_argument("--fail-status", type=int, default=503)
    p.add_argument("--match")
    args = p.parse_args()
    server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(args))
    print("port %d" % server.server_address[1], flush=True)
    server.serve_forever()


if __name__ == "__main__":
    main()
