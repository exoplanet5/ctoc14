import json, sys, collections
for tag in ("s16E", "s15a"):
    R = json.load(open(f"results/s16/verify_E/v2_{tag}.json"))
    fb = R["flybys"]
    print(f"==== {tag}: verdict {R['verdict']} J {R['J']:.6f} sumJi {R['sum_Ji']:.6f} Ncov {R['N_covered']} errors {len(R['errors'])} warnings {len(R['warnings'])}")
    for w in R["warnings"]: print("  WARN", w)
    st = collections.Counter(f["status"] for f in fb)
    print("  flyby rows", len(fb), dict(st))
    top = sorted(fb, key=lambda f: -f["d_km"])[:5]
    print("  5 largest flyby distances:")
    for f in top: print("    sc %d line %d ast %d t %.3f d %.6f km vrel %.4f km/s %s" % (f["sc"], f["line"], f["asteroid"], f["t"], f["d_km"], f["vrel_kms"], f["status"]))
    lo = sorted(fb, key=lambda f: f["vrel_kms"])[:5]
    print("  5 smallest vrel:")
    for f in lo: print("    sc %d line %d ast %d t %.3f d %.6f km vrel %.6f km/s %s" % (f["sc"], f["line"], f["asteroid"], f["t"], f["d_km"], f["vrel_kms"], f["status"]))
    cov = sorted({f["asteroid"] for f in fb if f["status"] == "COUNTED"})
    print("  missing:", sorted(set(range(1, 301)) - set(cov)))
    print("  max flyby t %.3f  (T_MISSION 473364000)" % max(f["t"] for f in fb))
    print("  max pos err %.6e km at %s; vel %.6e km/s at %s; mass %.6e kg at %s" % (R["max_pos_err_km"], R["max_pos_err_at"], R["max_vel_err_kms"], R["max_vel_err_at"], R["max_mass_err_kg"], R["max_mass_err_at"]))
    print("  max T sample %.9f, interp %.9f at %s" % (R["max_T_sample"], R["max_T_interp"], R["max_T_interp_at"]))
    sc = R["spacecraft"]
    print("  max vinf %.6f km/s; max launch pos err %.3e km; min row mass %.6f; min int mass %.6f; max t_end %.1f; min t_launch %.1f" % (
        max(s["v_inf_kms"] for s in sc), max(s["launch_pos_err_km"] for s in sc), min(s["min_row_mass"] for s in sc),
        min(s["min_integrated_mass"] for s in sc), max(s["t_end"] for s in sc), min(s["t_launch"] for s in sc)))
    print("  m0:", [round(s["m0"], 6) for s in sc], "max m0", max(s["m0"] for s in sc))
    print("  rk4 xcheck %s km; coast xcheck %.3e km" % (R["crosscheck_rk4_max_km"], R["coast_kepler_xcheck_km"]))
    print("  per-sc covered:", collections.Counter(f["sc"] for f in fb if f["status"] == "COUNTED"))
