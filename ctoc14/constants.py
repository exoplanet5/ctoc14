"""Official constants of CTOC14 Problem A (Table 2/3 of the problem statement). Units: km, s, kg, N, deg."""
import numpy as np

MU = 1.32712440018e11        # km^3/s^2, solar gravitational parameter
AU = 149597870.7             # km
G0 = 9.80665                 # m/s^2
ISP = 4000.0                 # s
TMAX = 0.5                   # N
VE = ISP * G0 / 1000.0       # exhaust velocity, km/s  (39.2266 km/s)
MDOT_MAX = TMAX / (ISP * G0) # kg/s at full thrust (1.2746e-5)
M_DRY = 600.0                # kg
M0_MAX = 2000.0              # kg
FUEL_MAX = M0_MAX - M_DRY    # 1400 kg
VINF_MAX = 4.0               # km/s
D_FLYBY = 1000.0             # km
DT_MIN_THRUST = 8640.0       # s, min spacing of Event=1 samples
DAY = 86400.0
T0_MJD = 62502.0             # 2030-01-01 00:00 UTC, mission t = 0
T_MISSION = 4.73364e8        # s = 15 yr = 5478.75 d
T_EPH_AST_MJD = 61200.0      # 2026-06-09 asteroid element epoch
T_EPH_EARTH_MJD = 60676.0    # 2025-01-01 Earth element epoch
# seconds from each element epoch to t0
T_OFF_AST = (T0_MJD - T_EPH_AST_MJD) * DAY      # 1.124928e8 s
T_OFF_EARTH = (T0_MJD - T_EPH_EARTH_MJD) * DAY  # 1.577664e8 s

# Earth heliocentric ecliptic elements (a AU, e, i deg, Omega deg, omega deg, M0 deg) at MJD 60676.0
EARTH_ELEMENTS = np.array([1.0009175020, 0.017566762041, 0.002976847126, 189.953211282428, 273.196254000254, 357.4135031077])

# validator tolerances
TOL_POS = 1.0      # km
TOL_VEL = 1e-3     # km/s
TOL_MASS = 0.01    # kg
TOL_VINF = 0.01    # km/s

def cost_sc(m0):
    """Per-spacecraft cost J_i = 1 + x + x^2, x = (m0-600)/1400."""
    x = (np.asarray(m0, dtype=float) - M_DRY) / FUEL_MAX
    return 1.0 + x + x * x
