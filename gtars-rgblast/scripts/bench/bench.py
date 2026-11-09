import sys, subprocess, json, time, os, shlex
# usage: bench.py <results.jsonl> <reps> <label>::<cleanup cmd>::<cmd> [...]
# Runs each command <reps> times, interleaved (rep 1 of all, then rep 2 ...).
# Each command is wrapped in GNU time; load average recorded before each run.
out, reps = sys.argv[1], int(sys.argv[2])
jobs = [a.split("::") for a in sys.argv[3:]]
TIME = ["/usr/bin/time", "-f", "TIME %e %M"]
for r in range(reps):
    for label, cleanup, cmd in jobs:
        if cleanup: subprocess.run(cleanup, shell=True)
        subprocess.run("sync", shell=True)
        load = open("/proc/loadavg").read().split()[:3]
        t0 = time.time()
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable="/bin/bash")
        wall = time.time() - t0
        tl = [l for l in p.stderr.splitlines() if l.startswith("TIME ")]
        el, rss = (tl[-1].split()[1:3] if tl else (None, None))
        rec = dict(label=label, rep=r, wall=round(wall, 2), time_elapsed=el, maxrss_kb=rss, load=load, rc=p.returncode,
                   tail=(p.stderr.strip().splitlines()[-3:] if p.returncode or "native" in label else []))
        print(json.dumps(rec), flush=True)
        with open(out, "a") as f: f.write(json.dumps(rec) + "\n")
