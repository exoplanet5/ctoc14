"""s15b: priority scheduler for tools/s15_twintail.py jobs under the machine-wide CPU cap.

SCHED.json = {"queues": ["results/s15b/queue_x.json", ...], "cap": 8, "stop": false, "hold": ["TAG", ...]}
Re-read every cycle, so queues / order / holds can be edited while it runs.  A job is PENDING when its OUT/meta.json does
not exist, it is not held, and no process runs (or sits paused on) its job file.  Every cycle the scheduler starts the
first pending jobs (queue order, then job order) while the system-wide number of CPU-bound Python processes
(ps pcpu > 50, 'Python.app') plus the jobs it started in the last 90 s stays below the cap.  "stop": true ends it
(running jobs continue; they resume from their checkpoints if restarted).

usage: nohup nice -n 5 python tools/s15b_sched.py results/s15b/sched.json > results/s15b/sched.out 2>&1 &"""
import os, sys, json, time, pathlib, subprocess
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14')
JDIR = ROOT / 'results/s15/twin_jobs'


def ps_lines():
    try:
        return subprocess.run(['ps', '-Ao', 'pid,stat,pcpu,command'], capture_output=True, text=True).stdout.splitlines()[1:]
    except Exception:
        return []


def busy(lines):
    """CPU-bound Python processes: not stopped (state T) and either pcpu > 50 or one of our compute tools (a process
    resumed from SIGSTOP shows a low decayed pcpu for minutes, so the command decides, not only pcpu).  Batch runners
    and this scheduler sleep and are not counted."""
    n = 0
    for l in lines:
        p = l.split(None, 3)
        if len(p) != 4 or 'Python.app' not in p[3] or p[1].startswith('T'):
            continue
        cmd = p[3]
        if 's15b_sched' in cmd or ' batch ' in cmd or 'ipykernel' in cmd or 'jupyter' in cmd:
            continue
        try:
            hot = float(p[2]) > 50
        except ValueError:
            hot = False
        if hot or '/tools/' in cmd or ' tools/' in cmd:
            n += 1
    return n


def main(sf):
    sf = ROOT / sf; recent = {}; children = []
    log = open(ROOT / 'results/s15b/sched.log', 'a')
    def say(s):
        line = f'[{time.strftime("%H:%M:%S")}] {s}'; print(line, flush=True); log.write(line + '\n'); log.flush()
    say(f'scheduler start ({sf})')
    while True:
        try:
            S = json.load(open(sf))
        except Exception as e:
            say(f'bad sched file: {e}'); time.sleep(30); continue
        if S.get('stop'):
            say('stop flag: exit'); break
        cap = min(int(S.get('cap', 8)), 8); hold = set(S.get('hold', []))
        children = [p for p in children if p.poll() is None]
        lines = ps_lines()
        cmdl = [l.split(None, 3)[3] for l in lines if len(l.split(None, 3)) == 4]
        now = time.time(); recent = {k: t for k, t in recent.items() if now - t < 90}
        b = busy(lines)
        # jobs started < 90 s ago may not show pcpu > 50 yet: count them on top (conservative double count)
        load = b + len(recent)
        pend = []
        for q in S.get('queues', []):
            try:
                jobs = json.load(open(ROOT / q))
            except Exception:
                continue
            for J in jobs:
                tag = J['tag']
                if tag in hold or (ROOT / J['out'] / 'meta.json').exists():
                    continue
                if any(f'twin_jobs/{tag}.json' in c for c in cmdl):
                    continue
                pend.append(J)
        while pend and load < cap:
            J = pend.pop(0); jf = JDIR / f'{J["tag"]}.json'; JDIR.mkdir(parents=True, exist_ok=True)
            json.dump(J, open(jf, 'w')); (ROOT / J['out']).mkdir(parents=True, exist_ok=True)
            p = subprocess.Popen(['nice', '-n', '5', sys.executable, str(ROOT / J.get('runner', 'tools/s15_twintail.py')), 'run', str(jf)],
                                 cwd=ROOT, stdout=open(ROOT / J['out'] / 'stdout.txt', 'a'), stderr=subprocess.STDOUT)
            children.append(p); recent[J['tag']] = time.time(); load += 1
            say(f'started {J["tag"]} (busy {b}, recent {len(recent)}, cap {cap}); {len(pend)} pending')
        time.sleep(15)


if __name__ == '__main__':
    main(sys.argv[1])
