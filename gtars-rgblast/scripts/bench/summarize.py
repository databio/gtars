import json, sys, statistics, collections
rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
by = collections.defaultdict(list)
for r in rows: by[r['label']].append(r)
for k, v in by.items():
    t = [float(r['time_elapsed']) for r in v if r['time_elapsed']]
    m = [int(r['maxrss_kb']) for r in v if r['maxrss_kb']]
    loads = [r['load'][0] for r in v]
    print(f"{k:22s} median {statistics.median(t):8.2f}s  runs {[round(x,2) for x in t]}  maxRSS median {statistics.median(m)/1024:7.0f} MB  load1 {loads}")
