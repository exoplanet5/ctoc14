"""Stage 15: re-queue the single globally-capped twin runner: kill it (running jobs keep going), put NEW jobs (json
files, in order) in front of the old queue's not-yet-started jobs, restart it on results/s15/twin/queue_<ts>.json.
usage: s15_requeue.py NEW.json [NEW2.json ...] [--front]"""
import os, sys, json, time, glob, pathlib, subprocess
ROOT = pathlib.Path('/Users/mickey/solarsystem/ctoc14'); os.chdir(ROOT)
S = json.load(open('results/s15/stage2_state.json'))
cur = S.get('queue_file', 'results/s15/twin/round4_q2.json')
subprocess.run(['pkill', '-f', f's15_twintail.py batch {cur}'])
time.sleep(1)
old = json.load(open(cur)) if pathlib.Path(cur).exists() else []
new = []
for f in sys.argv[1:]:
    if f.startswith('--'):
        continue
    new += json.load(open(f))
started = lambda j: (ROOT / j['out'] / 'log.txt').exists()
tags = set()
q = []
for j in new + old:
    if started(j) or j['tag'] in tags:
        continue
    tags.add(j['tag']); q.append(j)
qf = f'results/s15/twin/queue_{time.strftime("%H%M%S")}.json'
json.dump(q, open(qf, 'w'), indent=1)
env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1')
subprocess.Popen([sys.executable, 'tools/s15_twintail.py', 'batch', qf, '--procs', '8', '--global-cap', '8'], cwd=ROOT,
                 stdout=open(qf.replace('.json', '.out'), 'w'), stderr=subprocess.STDOUT, env=env, start_new_session=True)
S['queue_file'] = qf; S['updated'] = time.strftime('%m-%d %H:%M')
S.setdefault('queue_history', []).append(qf)
json.dump(S, open('results/s15/stage2_state.json', 'w'), indent=1)
print(f'queue {qf}: {[j["tag"] for j in q]}')
