// Kepler propagation for the CTOC14 viewer. Mirrors ctoc14/kepler.py (Body.state) so that the page reproduces the
// contest ephemeris: docs/viz_kepler.md did not exist when this was written, so this is the standard Newton solution
// with the R3(-Om) R1(-i) R3(-om) rotation, the same initial guess, step clip and tolerance as the Python code.
//   M(t) = M0 + n (t - t_epoch),  n = sqrt(mu / a^3)  [a in km, n in rad/s],  E - e sin E = M,
//   x_p = a (cos E - e),  y_p = a sqrt(1 - e^2) sin E,  r = R [x_p, y_p, 0].
// Times are days since t0 (MJD 62502.0); positions are returned in AU (computed in km, divided by AU).

export const DAY = 86400.0;
const DEG = Math.PI / 180.0;
const TWO_PI = 2.0 * Math.PI;

export function solveKepler(M, e) {
  M = M % TWO_PI;
  if (M < 0) M += TWO_PI;
  let E = e < 0.8 ? M + e * Math.sin(M) : Math.PI;
  for (let k = 0; k < 60; k++) {
    const f = E - e * Math.sin(E) - M;
    const fp = 1.0 - e * Math.cos(E);
    let dE = -f / fp;
    if (dE > 1.0) dE = 1.0; else if (dE < -1.0) dE = -1.0;
    E += dE;
    if (Math.abs(dE) < 1e-13) break;
  }
  return E;
}

// Perifocal -> inertial rotation R3(-Om) R1(-i) R3(-w), row-major 3x3.
export function rotationMatrix(OmRad, iRad, wRad) {
  const cO = Math.cos(OmRad), sO = Math.sin(OmRad), ci = Math.cos(iRad), si = Math.sin(iRad), cw = Math.cos(wRad), sw = Math.sin(wRad);
  return [
    cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si,
    sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si,
    sw * si, cw * si, ci,
  ];
}

export class KeplerBody {
  /** el: {a [AU], e, i, Om, om, M0 [deg], epoch_mjd}; t0Mjd: mission t0; muKm: km^3/s^2; auKm: km. */
  constructor(el, t0Mjd, muKm, auKm) {
    this.auKm = auKm;
    this.aKm = el.a * auKm;
    this.e = el.e;
    this.M0 = el.M0 * DEG;
    this.n = Math.sqrt(muKm / (this.aKm * this.aKm * this.aKm));  // rad/s
    this.tOff = (t0Mjd - el.epoch_mjd) * DAY;                       // s from element epoch to t0
    this.R = rotationMatrix(el.Om * DEG, el.i * DEG, el.om * DEG);
    this.sq = Math.sqrt(1.0 - el.e * el.e);
    this.periodDays = TWO_PI / this.n / DAY;
    this.muKm = muKm;
  }
  /** Position in AU at mission time tDays into out[3] (or a new array). */
  positionAU(tDays, out) {
    const M = this.M0 + this.n * (tDays * DAY + this.tOff);
    const E = solveKepler(M, this.e);
    const cE = Math.cos(E), sE = Math.sin(E);
    const x = this.aKm * (cE - this.e), y = this.aKm * this.sq * sE;
    const R = this.R, k = 1.0 / this.auKm;
    out = out || [0, 0, 0];
    out[0] = (R[0] * x + R[1] * y) * k;
    out[1] = (R[3] * x + R[4] * y) * k;
    out[2] = (R[6] * x + R[7] * y) * k;
    return out;
  }
  /** Position [km] and velocity [km/s] at tDays (same arithmetic as ctoc14.kepler.Body.state). */
  stateKm(tDays) {
    const M = this.M0 + this.n * (tDays * DAY + this.tOff);
    const E = solveKepler(M, this.e);
    const cE = Math.cos(E), sE = Math.sin(E);
    const r = this.aKm * (1.0 - this.e * cE);
    const x = this.aKm * (cE - this.e), y = this.aKm * this.sq * sE;
    const fac = Math.sqrt(this.muKm * this.aKm) / r;
    const vx = -fac * sE, vy = fac * this.sq * cE;
    const R = this.R;
    return {
      r: [R[0] * x + R[1] * y, R[3] * x + R[4] * y, R[6] * x + R[7] * y],
      v: [R[0] * vx + R[1] * vy, R[3] * vx + R[4] * vy, R[6] * vx + R[7] * vy],
    };
  }
  /** One full orbit sampled uniformly in eccentric anomaly (closed loop, n points) as a Float32Array (n*3), AU. */
  orbitPoints(n) {
    const arr = new Float32Array(n * 3);
    const R = this.R, k = 1.0 / this.auKm;
    for (let j = 0; j < n; j++) {
      const E = TWO_PI * j / n;
      const x = this.aKm * (Math.cos(E) - this.e), y = this.aKm * this.sq * Math.sin(E);
      arr[3 * j] = (R[0] * x + R[1] * y) * k;
      arr[3 * j + 1] = (R[3] * x + R[4] * y) * k;
      arr[3 * j + 2] = (R[6] * x + R[7] * y) * k;
    }
    return arr;
  }
}

/** Illustrative geocentric mean lunar orbit (bodies.json "moon" block) with linearly precessing Om and om. */
export class MoonOrbit {
  constructor(moon, t0Mjd, auKm) {
    this.m = moon;
    this.t0Mjd = t0Mjd;
    this.auKm = auKm;
    this.n = TWO_PI / moon.period_d;   // rad/day
    this.sq = Math.sqrt(1 - moon.e * moon.e);
  }
  /** Geocentric ecliptic position in AU at mission time tDays. */
  positionGeoAU(tDays, out) {
    const m = this.m, d = (this.t0Mjd + tDays) - m.epoch_mjd;
    const M = m.M0 * DEG + this.n * d;
    const Om = (m.Om + m.Om_rate_deg_per_day * d) * DEG, om = (m.om + m.om_rate_deg_per_day * d) * DEG;
    const E = solveKepler(M, m.e);
    const x = m.a_km * (Math.cos(E) - m.e), y = m.a_km * this.sq * Math.sin(E);
    const R = rotationMatrix(Om, m.i * DEG, om), k = 1.0 / this.auKm;
    out = out || [0, 0, 0];
    out[0] = (R[0] * x + R[1] * y) * k;
    out[1] = (R[3] * x + R[4] * y) * k;
    out[2] = (R[6] * x + R[7] * y) * k;
    return out;
  }
  /** Geocentric orbit loop at epoch tDays (n points, AU) for the ring in Earth-Moon zoom mode. */
  orbitPoints(tDays, n, arr) {
    const m = this.m, d = (this.t0Mjd + tDays) - m.epoch_mjd;
    const Om = (m.Om + m.Om_rate_deg_per_day * d) * DEG, om = (m.om + m.om_rate_deg_per_day * d) * DEG;
    const R = rotationMatrix(Om, m.i * DEG, om), k = 1.0 / this.auKm;
    arr = arr || new Float32Array(n * 3);
    for (let j = 0; j < n; j++) {
      const E = TWO_PI * j / n;
      const x = m.a_km * (Math.cos(E) - m.e), y = m.a_km * this.sq * Math.sin(E);
      arr[3 * j] = (R[0] * x + R[1] * y) * k;
      arr[3 * j + 1] = (R[3] * x + R[4] * y) * k;
      arr[3 * j + 2] = (R[6] * x + R[7] * y) * k;
    }
    return arr;
  }
}

/** MJD -> "YYYY-MM-DD HH:MM UTC". */
export function mjdToUTC(mjd) {
  const ms = (mjd - 40587.0) * 86400000.0;
  const dt = new Date(ms);
  const p = (v) => String(v).padStart(2, '0');
  return `${dt.getUTCFullYear()}-${p(dt.getUTCMonth() + 1)}-${p(dt.getUTCDate())} ${p(dt.getUTCHours())}:${p(dt.getUTCMinutes())} UTC`;
}
