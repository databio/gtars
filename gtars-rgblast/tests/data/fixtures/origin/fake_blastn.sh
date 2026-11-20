#!/bin/sh
# A stand-in for blastn in the heartbeat, cancel and quiet tests: answers -version, writes one
# stderr line, sleeps $FAKE_BLASTN_SLEEP seconds (default 5) and exits $FAKE_BLASTN_EXIT
# (default 3).
case "$1" in
  -version) echo "blastn: 2.17.0+"; echo " Package: blast 2.17.0, build fake"; exit 0 ;;
esac
echo "fake blastn: sleeping" >&2
sleep "${FAKE_BLASTN_SLEEP:-5}"
exit "${FAKE_BLASTN_EXIT:-3}"
