"""Convert a tour JSON to a submission and validate. Usage: run_convert.py tour.json out.txt [h_days] [tcap]"""
import sys, json, time, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.kepler import Ephemeris
from ctoc14.lowthrust import convert_tour
from ctoc14.submission import write_submission
from ctoc14.validator import validate
from ctoc14.constants import DAY
tour = json.load(open(sys.argv[1])); out = sys.argv[2]
h = float(sys.argv[3]) * DAY if len(sys.argv) > 3 else 1.0 * DAY
tcap = float(sys.argv[4]) if len(sys.argv) > 4 else 0.45
eph = Ephemeris(); tic = time.time()
traj, info = convert_tour(eph, tour, h=h, tcap=tcap, verbose=True)
print(f'conversion done in {time.time()-tic:.0f} s: {info["n_flybys"]} flybys, fuel {info["fuel"]:.1f} kg, dropped {info["dropped"]}')
write_submission(out, [traj], header=f'CTOC14 single spacecraft; tour {sys.argv[1]}')
json.dump(info, open(out.replace('.txt', '_info.json'), 'w'), indent=1)
rep = validate(out)
