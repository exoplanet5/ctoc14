"""Does the driver's pool pattern (fork Pool, maxtasksperchild=1, imap_unordered) survive a worker that dies (SIGKILL /
segfault / OOM)?  Same construction as s18_rhfa.run/propose.  A 45 s alarm marks the hang."""
import os, sys, time, signal, multiprocessing as mp
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
import s18_rhfa as R

def job(k):
    if k == 1:
        os.kill(os.getpid(), signal.SIGKILL)      # a crashed worker (segfault / OOM killer look the same to the pool)
    time.sleep(0.5); return k

def alarm(sig, frm):
    print('HANG: imap_unordered did not return within 45 s after one worker died -> the slice would wait forever', flush=True)
    os._exit(3)

if __name__ == '__main__':
    signal.signal(signal.SIGALRM, alarm); signal.alarm(45)
    pool = mp.get_context('fork').Pool(2, initializer=R._worker_init, maxtasksperchild=1)
    tic = time.time(); got = []
    for r in pool.imap_unordered(job, [0, 1, 2, 3]):
        got.append(r); print('got', r, f'{time.time()-tic:.1f} s', flush=True)
    print('all results returned', got)
