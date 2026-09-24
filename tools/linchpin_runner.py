"""Run the shell commands of a jobs file (one per line) with at most -P at a time; log start/end/rc per job.
Usage: linchpin_runner.py jobs.txt -P 8"""
import sys, time, subprocess, argparse

ap = argparse.ArgumentParser(); ap.add_argument('jobs'); ap.add_argument('-P', type=int, default=8)
a = ap.parse_args()
cmds = [l.strip() for l in open(a.jobs) if l.strip() and not l.startswith('#')]
run = {}; k = 0
log = lambda s: print(f'[{time.strftime("%H:%M:%S")}] {s}', flush=True)
while k < len(cmds) or run:
    while k < len(cmds) and len(run) < a.P:
        run[subprocess.Popen(cmds[k], shell=True)] = (k, time.time()); log(f'start {k}: {cmds[k][-120:]}'); k += 1
    time.sleep(2)
    for p in list(run):
        if p.poll() is not None:
            j, t0 = run.pop(p); log(f'done {j} rc {p.returncode} in {time.time() - t0:.0f} s')
log('all done')
