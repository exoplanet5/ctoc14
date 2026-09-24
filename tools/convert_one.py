"""Convert one tour JSON -> submission fragment (single spacecraft) + info JSON. Usage: convert_one.py tour.json out.txt [m0] [--h days]"""
import sys, json, time, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate
h_days = 1.0
if '--h' in sys.argv:
    k = sys.argv.index('--h'); h_days = float(sys.argv[k + 1]); del sys.argv[k:k + 2]
tour = json.load(open(sys.argv[1])); out = sys.argv[2]
if len(sys.argv) > 3: tour['m0'] = float(sys.argv[3])
eph = Ephemeris(); tic = time.time()
from ctoc14.constants import DAY
traj, info = convert_tour(eph, tour, h=h_days * DAY, verbose=True)
info['runtime_s'] = time.time() - tic; info['m0'] = tour['m0']
write_submission(out, [traj], header=f'tour {sys.argv[1]} m0={tour["m0"]}')
json.dump(info, open(out.replace('.txt', '_info.json'), 'w'), indent=1, default=float)
rep = validate(out)
print(f'DONE {sys.argv[1]}: flybys {info["n_flybys"]} fuel {info["fuel"]:.1f} kg dropped {info["dropped"]} valid={rep.ok} runtime {info["runtime_s"]:.0f}s')
