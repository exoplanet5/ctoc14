# Identification of the 300 CTOC14 targets (MEA.txt)

**Verdict: all 300 are real, catalogued near-Earth asteroids.** MEA.txt is a verbatim extract of the Minor Planet Center
orbit catalogue (MPCORB `NEA.txt`, standard epoch K2669 = 2026-06-09.0 = MJD 61200.0, the contest epoch). The selection rule
is *potentially hazardous asteroids (PHA: Earth MOID ≤ 0.05 AU, H ≤ 22) with absolute magnitude H ≤ 18.49*, i.e. the ~300 largest
PHAs (diameters ~0.5–5 km). The MPC list today has 299 PHAs with H ≤ 18.49, all in MEA.txt; the 300th member (#280, (288807) 2004 RW164,
H = 18.41) has since lost its PHA flag but was evidently flagged in the snapshot the organisers used.

## Method
1. `tools/identify_nea.py`: downloaded the full JPL SBDB NEO catalogue (`ssd-api.jpl.nasa.gov/sbdb_query.api`, 42k asteroids + comets,
   full precision) and matched each MEA row by nearest neighbour in (a, e, i, Ω, ω); then queried the **JPL Horizons API**
   (`ssd.jpl.nasa.gov/api/horizons.api`, EPHEM_TYPE=ELEMENTS, heliocentric ecliptic J2000, TLIST = JD 2461200.5 TDB) for each
   matched body and compared all six elements with MEA.txt.
2. `tools/identify_nea_mpc.py`: matched all 300 rows against MPCORB `NEA.txt` (`data/MPC_NEA.txt`, 42 253 NEAs) on all six elements.

## Results
- MPC: 265/300 rows agree with the MPC catalogue in all six elements to < 1e-6 (worst |diff| 3.9e-04); every object is identified.
- JPL SBDB: 297/300 have a JPL orbit with essentially identical elements (JPL's standard epoch is also JD 2461200.5). The 3 exceptions are
  poorly-observed objects where the JPL and MPC orbit solutions differ: #144 = 1999 XS35 (JPL a = 17.8174 vs MPC 17.8136 AU),
  #206 = 2015 BO519 (JPL lists a different solution), #278 = 2001 VB (JPL a = 2.392 vs MPC 2.331 AU).
- JPL Horizons osculating elements at the contest epoch: 286/300 objects agree with MEA.txt to better than 5.4e-04
  (relative in a, e; degrees in angles); median |ΔΩ| 8.4e-06°, |ΔM| 4.1e-06°, |Δa/a| 1.5e-08.
  The other 14 objects (IDs [64, 139, 144, 147, 166, 192, 206, 214, 235, 251, 259, 264, 278, 287]) are short-arc/poorly determined orbits (e.g. #144 arc 88 d) where JPL and MPC
  solutions differ mainly in mean anomaly (tens of degrees) and slightly in a — same real objects, different orbit fits. **The contest
  dynamics use the MEA.txt (MPC) elements as truth, so these differences are irrelevant for the competition.**
- H from 14.07 (#1 (3122) Florence) to 18.49; median 17.72; the list is sorted by H (brightest first). 274 numbered, 26 unnumbered, 30 named.
- Well-known members: (3122) Florence, (4183) Cuno, (3200) Phaethon, (1620) Geographos, (1981) Midas, (4179) Toutatis, (2201) Oljato,
  (4486) Mithra, (12923) Zephyr, (52768) 1998 OR2, (161989) Cacus; #131 = 2025 VP (retrograde, i = 134°, a = 8.6 AU); #144 = 1999 XS35 (a = 17.8 AU, P = 75 yr).
- SBDB physical data: diameters known for 114 objects (median 1.15 km, max 7.00 km = Florence); SBDB PHA flag Y for 298
  (the 2 others are the JPL-solution outliers or borderline MOID); SBDB Earth MOID ≤ 0.05 AU for 299.

Full table: `data/nea_catalog_ids.csv` (MEA elements, MPC designation/name/H/arc/flags, SBDB name/spkid/diameter/MOID/PHA/arc, Horizons deviations).
Raw downloads: `data/sbdb_neo_asteroids.json`, `data/sbdb_neo_comets.json`, `data/MPC_NEA.txt`; Horizons log: `results/identify_nea.log`.
