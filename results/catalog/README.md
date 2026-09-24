# results/catalog -- sorted inventory of every CTOC14 result (written 2026-09-24 by the CATALOG agent)

Nothing here is a contest record; the records are `results/CTOC14_Result_*.txt` and `results/validator2_*.log`, which this
directory only reads. Layout:

- `inventory_raw.json` -- raw measurements: every `results/CTOC14_Result_*.txt` (craft, flybys, covered, missed, per-craft m0 /
  J_i / launch and end day / fuel, md5, size, date, validator verdict) and every `fleet.json` / `route_*.npz` directory under
  `results/` (n, covered, flybys, tanks, sum J_i, J). Produced by `inventory.py` (here); re-run it to refresh
  (`nice -n 5 ~/.venvs/astro313/bin/python results/catalog/inventory.py`, ~1 min).
- `fleets.json` -- the catalog: the same records with stage, provenance, keep/derived/superseded status, the decoded board
  history (`board_history`), duplicate detection (`same_md5_as`), stale-log detection (`stale_validator_log`) and the viewer
  selection summary. Produced by `build_catalog.py` (here) from `inventory_raw.json`; it also writes `viz_selection.json`,
  `README.md` and `docs/results_catalog.md`.
- `viz_selection.json` -- the 9 fleets chosen for the 3D viewer (keys: b300, v4, mix1, t10d, s16E, s16a, p13, p19, g2), each with
  source path, kind (submission | twin), craft count, covered, missed, raw J, date, m0 per craft, title and note. The FLEETS
  agent converts these into `results/viz/data/fleets/<key>.json` (contract: docs/viz_spec.md section 2).
- `validator2_<tag>.log` -- `tools/validator2.py --quiet` output for the 10 submission files that had no contest-time log
  (b300, E1, E3_fleet, milp2, mix1, mix2, swarm_prelim, swarm_v1, swarm_v2, v3), run 2026-09-24.
- Human-readable table: `docs/results_catalog.md`.

Conventions: dates are file mtimes (local, CST); "submitted" times come from the leaderboard snapshots (`leaderboard/<UTC stamp>/`);
raw J = sum J_i + (300 - covered); twin J = impulsive-model estimate, exact J = submission file value.
