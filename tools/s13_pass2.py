"""Stage 13 [P], revision 2: tools/s13_pass.py (imported, unchanged) with the head-route walk fixed.
Found when generation 1 started (2026-09-22 11:16-11:29, results/s13/gen/g01_jitall/):

 * DEEPER CANDIDATE FRONT.  run_member kept only the 8 deepest min-fuel levels.  On route 1 all 8 settled above the
   1050 kg head cap (fronts 48-49 flybys at planner 1026-1035 kg -> twin 1107-1178 kg), the walk found nothing, and
   the pass ended with 0 routes: all six full passes of the first attempt.  The front now keeps 24 levels, so the
   walk (deepest first, batches of 6, at most 18 settles) steps down to the deepest level that settles under the
   cap.  Tail routes are unaffected in practice: their walk takes the deepest candidates across the 5 members first,
   and those are the same.
 * The beam cache key carries the front size (no 8-level front is reused).
 * S13_JITTER_MODE=prem (environment) restricts --jitter to the premium set (targets with a non-zero base lambda);
   the default 'all' is s13_pass's behaviour (U(-j, j) on all 300 lambdas).  Measured on route 1 (scratchpad jit/):
   no jitter, premium-set jitter and all-target jitter give the same front (48 @ 1028-1031 kg planner), so the
   jitter was NOT the cause of the heavy heads; the premium SET was (gc16's 67-target hard set with dv 1.8 legs:
   deepest settled under 1050 kg = 38 flybys; dv 1.2: 43; FREE's N6a used the 45-target U1 set: 48 @ 1048).
Everything else (members, walk, repair, pilot, library, T7, outputs, restart) is s13_pass.

usage: as tools/s13_pass.py."""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, json, hashlib, pathlib
sys.path.insert(0, '/Users/mickey/solarsystem/ctoc14/tools')
import numpy as np
import s13_pass as SP                       # chdirs to ROOT
SB = SP.SB
FRONT2 = 24
JITTER_MODE = os.environ.get('S13_JITTER_MODE', 'all')


def beam_key(job):
    avail = job['avail']
    return SB.Cache.key(dict(v=SB.CODE_TAG + '/p2', front=SP.FRONT, jit=JITTER_MODE, kind=job['kind'], m=job['m'],
                             beam=job['beam'], avail=avail,
                             prize=[round(float(job['prize'][t - 1]), 12) for t in avail],
                             prem=sorted(set(job['prem']) & set(avail)), dvt=job['dvt'],
                             **({'covered': job['covered'], 'k0': job['k0'], 'N': job['N'], 'caps': job['caps']}
                                if job['kind'] == 'rollout' else {})))


class Pass2(SP.Pass):
    def premiums(self):
        a = self.a
        base = np.zeros(300)
        if a.prem:
            d = json.load(open(a.prem if os.path.isabs(a.prem) else pathlib.Path(SP._CWD0) / a.prem))
            d = d.get('lambda', d) if isinstance(d, dict) else d
            if isinstance(d, dict):
                for t, v in d.items():
                    base[int(t) - 1] = float(v)
            else:
                base[:] = np.asarray(d, float)[:300]
        elif a.delta:
            for t in SP.hard_set():
                base[t - 1] = a.delta
        base = np.clip(base, -0.01, 0.02)
        mask = (np.abs(base) > 1e-12) if JITTER_MODE == 'prem' else np.ones(300, bool)
        lam = base.copy()
        if a.jitter > 0:
            rng = np.random.default_rng(a.seed)
            lam = np.clip(lam + rng.uniform(-a.jitter, a.jitter, 300) * mask, -0.01, 0.02)
        self.lam = lam
        self.prize = np.clip(1.0 + lam, 1.00, 1.02)
        self.prem = sorted(t for t in SP.ALL if base[t - 1] > 1e-9)
        ALLi = np.array(SP.ALL) - 1
        self.lam_digest = dict(n_prem=len(self.prem), n_pos=int((lam > 1e-9).sum()), n_neg=int((lam < -1e-9).sum()),
                               n_jittered=int(mask.sum()) if a.jitter > 0 else 0,
                               mean=round(float(lam.mean()), 5), max=round(float(lam.max()), 4),
                               prize_min=round(float(self.prize[ALLi].min()), 4),
                               prize_max=round(float(self.prize[ALLi].max()), 4),
                               sha1=hashlib.sha1(np.round(lam, 9).tobytes()).hexdigest()[:12],
                               src=a.prem or f'delta {a.delta} on gc16 hard set', jitter=a.jitter, seed=a.seed,
                               jitter_mode=JITTER_MODE, front=SP.FRONT, tool='s13_pass2')


SP.FRONT = FRONT2
SP.beam_key = beam_key
SP.Pass = Pass2

if __name__ == '__main__':
    SP.main()
