"""Build the PDF example (launch -> ballistic flyby of 174 -> 41-sample thrust arc -> coast -> end), write it, validate it."""
import sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ctoc14.submission import SCTrajectory, write_submission
from ctoc14.validator import validate
from ctoc14.constants import DAY
r1 = np.array([-1.9500328647e+07, 1.4581848347e+08, -7.6372048832e+03]); v1 = np.array([-2.9389477847e+01, -4.3166912970e+00, -3.9417153392e+00])
sc = SCTrajectory(0.0, r1, v1, 2000.0)
sc.coast_to(1.0008479152e+07); sc.flyby(174)
t3 = 1.0008489152e+07
ts = t3 + np.arange(41) * DAY
mag = 0.5 * np.sin(np.pi * np.arange(41) / 40.0)
az = np.radians(40.107 + 9.0 * np.arange(41)); el = np.radians(17.189)
Ts = mag[:, None] * np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el) * np.ones(41)], axis=1)
sc.thrust_arc(ts, Ts)
sc.coast(30 * DAY, node=True); sc.coast(100 * DAY); sc.end()
out = pathlib.Path(__file__).resolve().parents[1] / 'results' / 'test_example_submission.txt'
write_submission(out, [sc], header='CTOC14 test: PDF example reconstruction')
print('rows written:', len(sc.rows), '->', out)
# compare with PDF lines 4 and 5
for k, (m_pdf, r_pdf) in enumerate([(1.9999567651e+03, [-1.1238880278e+08, -8.9304194192e+07, -1.6350356961e+07]), (1.9998273580e+03, [-1.1075466879e+08, -9.1364740090e+07, -1.6171502174e+07])]):
    row = sc.rows[3 + k]
    print(f'line {4+k}: mass diff {row[4]-m_pdf:+.6f} kg, pos diff {np.linalg.norm(row[2]-np.array(r_pdf)):.4f} km')
rep = validate(out)
