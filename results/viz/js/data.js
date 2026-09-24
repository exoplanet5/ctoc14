// Data loading: real files (data/index.json + data/bodies.json) with fallback to the fixture (data/fixture/...).
// Fleet files are resolved relative to the directory of the index that was loaded.

async function fetchJSON(url) {
  const r = await fetch(url, { cache: 'no-cache' });
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return r.json();
}

export async function loadCatalog() {
  const candidates = [
    { dir: 'data/', label: 'REAL data (results/viz/data/)' },
    { dir: 'data/fixture/', label: 'FIXTURE data (results/viz/data/fixture/) - real index.json/bodies.json missing' },
  ];
  // ?data=<dir> forces a data directory (e.g. ?data=data/fixture/ or a test set); it is tried first
  const forced = new URLSearchParams(location.search).get('data');
  if (forced) candidates.unshift({ dir: forced.endsWith('/') ? forced : forced + '/', label: `FORCED data (?data=${forced})` });
  const errors = [];
  for (const c of candidates) {
    try {
      const [index, bodies] = await Promise.all([fetchJSON(c.dir + 'index.json'), fetchJSON(c.dir + 'bodies.json')]);
      if (!Array.isArray(index.fleets) || !index.fleets.length) throw new Error(c.dir + 'index.json has no fleets');
      if (!Array.isArray(bodies.asteroids) || !bodies.asteroids.length) throw new Error(c.dir + 'bodies.json has no asteroids');
      console.log(`[viewer] loaded ${c.label}: ${index.fleets.length} fleet(s), ${bodies.asteroids.length} asteroid(s)`);
      return { dir: c.dir, label: c.label, index, bodies };
    } catch (e) {
      errors.push(String(e.message || e));
      console.warn(`[viewer] ${c.dir} not usable: ${e.message || e}`);
    }
  }
  throw new Error('no data found: ' + errors.join(' | '));
}

/** Load one fleet file and derive per-craft typed arrays and event tables. */
export async function loadFleet(dir, entry) {
  const raw = await fetchJSON(dir + entry.file);
  const fleet = prepareFleet(raw);
  fleet.entry = entry;
  console.log(`[viewer] fleet ${fleet.key}: ${fleet.craft.length} craft, ${fleet.allFlybys.length} flyby events, ${fleet.nSamples} samples, tmax ${fleet.tmax.toFixed(2)} d`);
  return fleet;
}

export function prepareFleet(raw) {
  let tmax = 0, nSamples = 0;
  const allFlybys = [];
  const craft = raw.craft.map((c, ci) => {
    const n = c.t.length;
    const t = Float64Array.from(c.t);
    const xyz = new Float64Array(n * 3);
    for (let k = 0; k < n; k++) { xyz[3 * k] = c.xyz[k][0]; xyz[3 * k + 1] = c.xyz[k][1]; xyz[3 * k + 2] = c.xyz[k][2]; }
    for (let k = 1; k < n; k++) if (!(t[k] > t[k - 1])) console.error(`[viewer] craft ${c.id}: t not strictly increasing at sample ${k}`);
    const m = c.m && c.m.length === n ? Float64Array.from(c.m) : null;
    const flybys = (c.flybys || []).map(f => {
      const idx = indexOfTime(t, f.t);
      if (idx < 0) console.warn(`[viewer] craft ${c.id}: flyby epoch ${f.t} not an exact sample`);
      return { t: f.t, ast: f.ast, idx, craft: ci };
    }).sort((a, b) => a.t - b.t);
    flybys.forEach(f => allFlybys.push(f));
    const launch = c.launch_d != null ? c.launch_d : t[0];
    const end = c.end_d != null ? c.end_d : t[n - 1];
    tmax = Math.max(tmax, end, t[n - 1]);
    nSamples += n;
    return { id: c.id, index: ci, m0: c.m0, fuel_kg: c.fuel_kg != null ? c.fuel_kg : (c.m0 - 600), launch, end, n, t, xyz, m,
             thrust: (c.thrust || []).map(a => [a[0], a[1]]), flybys, n_flybys: c.n_flybys != null ? c.n_flybys : flybys.length,
             visible: true };
  });
  allFlybys.sort((a, b) => a.t - b.t);
  // first detection per asteroid (only the first flyby of an asteroid counts)
  const firstSeen = new Map();
  for (const f of allFlybys) if (!firstSeen.has(f.ast)) firstSeen.set(f.ast, f.t);
  return { key: raw.key, title: raw.title, kind: raw.kind, source: raw.source, n_craft: raw.n_craft || craft.length,
           covered: raw.covered, missed: raw.missed || [], J_raw: raw.J_raw, sumJi: raw.sumJi, date: raw.date,
           craft, allFlybys, firstSeen, tmax: Math.max(tmax, 1), nSamples };
}

function indexOfTime(t, value) {
  let lo = 0, hi = t.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (t[mid] < value - 1e-9) lo = mid + 1; else if (t[mid] > value + 1e-9) hi = mid - 1; else return mid;
  }
  return -1;
}

/** Largest index k with t[k] <= value (binary search); -1 if value < t[0]. */
export function upperIndex(t, value) {
  let lo = 0, hi = t.length - 1;
  if (value < t[0]) return -1;
  if (value >= t[hi]) return hi;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (t[mid] <= value) lo = mid; else hi = mid; }
  return lo;
}

/** Linear interpolation of the craft position at time v into out[3]; returns the sample index k (t[k] <= v) or -1 / n-1 when outside. */
export function interpPosition(c, v, out) {
  const k = upperIndex(c.t, v);
  if (k < 0) { out[0] = c.xyz[0]; out[1] = c.xyz[1]; out[2] = c.xyz[2]; return -1; }
  if (k >= c.n - 1) { const j = 3 * (c.n - 1); out[0] = c.xyz[j]; out[1] = c.xyz[j + 1]; out[2] = c.xyz[j + 2]; return c.n - 1; }
  const f = (v - c.t[k]) / (c.t[k + 1] - c.t[k]), j = 3 * k;
  out[0] = c.xyz[j] + f * (c.xyz[j + 3] - c.xyz[j]);
  out[1] = c.xyz[j + 1] + f * (c.xyz[j + 4] - c.xyz[j + 1]);
  out[2] = c.xyz[j + 2] + f * (c.xyz[j + 5] - c.xyz[j + 2]);
  return k;
}

export function interpMass(c, v) {
  if (!c.m) return null;
  const k = upperIndex(c.t, v);
  if (k < 0) return c.m[0];
  if (k >= c.n - 1) return c.m[c.n - 1];
  const f = (v - c.t[k]) / (c.t[k + 1] - c.t[k]);
  return c.m[k] + f * (c.m[k + 1] - c.m[k]);
}
