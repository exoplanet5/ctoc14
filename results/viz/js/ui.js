// HUD, control panel, timeline canvas and HTML label overlay.
import * as THREE from 'three';
import { craftColor } from './fleet.js';
import { interpMass } from './data.js';
import { mjdToUTC } from './kepler.js';

const $ = (id) => document.getElementById(id);
const SPEED_MIN = 1, SPEED_MAX = 200;

export class UI {
  constructor(app) {
    this.app = app;
    this.el = {
      fleetSelect: $('fleet-select'), fleetNote: $('fleet-note'), play: $('btn-play'), speed: $('speed'), speedVal: $('speed-val'),
      time: $('time'), dateVal: $('date-val'), tdayVal: $('tday-val'), covered: $('hud-covered'), craft: $('hud-craft'),
      launched: $('hud-launched'), fuel: $('hud-fuel'), fuelLabel: $('hud-fuel-label'), j: $('hud-j'), sumji: $('hud-sumji'),
      table: $('craft-table').querySelector('tbody'), follow: $('follow-select'), optOrbits: $('opt-orbits'), optLabels: $('opt-labels'),
      optAstLabels: $('opt-astlabels'), optMoon: $('opt-moon'), legendMoon: $('legend-moon'), tl: $('tl-canvas'), labels: $('labels'),
      status: $('status'), fps: $('fps'), panel: $('panel'), togglePanel: $('toggle-panel'),
    };
    this.labelNodes = new Map();
    this.rows = [];
    this.tlStatic = document.createElement('canvas');
    this.tlDirty = true;
    this.lastTl = -1;
    this.bind();
  }

  bind() {
    const a = this.app, e = this.el;
    e.play.addEventListener('click', () => a.togglePlay());
    $('btn-step-back').addEventListener('click', () => a.setTime(a.t - 10));
    $('btn-step-fwd').addEventListener('click', () => a.setTime(a.t + 10));
    $('btn-reset').addEventListener('click', () => { a.pause(); a.setTime(0); });
    $('btn-home').addEventListener('click', () => { e.follow.value = ''; a.followIndex = -1; e.optMoon.checked = false; a.setZoom(false); a.scene.homeView(); });
    $('btn-top').addEventListener('click', () => { e.follow.value = ''; a.followIndex = -1; e.optMoon.checked = false; a.setZoom(false); a.scene.topView(); });
    e.speed.addEventListener('input', () => { a.speed = this.speedFromSlider(+e.speed.value); e.speedVal.textContent = `${a.speed.toFixed(a.speed < 10 ? 1 : 0)} d/s`; });
    e.time.addEventListener('input', () => a.setTime(+e.time.value, true));
    e.fleetSelect.addEventListener('change', () => a.selectFleet(e.fleetSelect.value));
    e.follow.addEventListener('change', () => { a.followIndex = e.follow.value === '' ? -1 : +e.follow.value; if (a.followIndex >= 0) { e.optMoon.checked = false; a.setZoom(false); } });
    e.optOrbits.addEventListener('change', () => { a.scene.orbitLines.visible = e.optOrbits.checked; });
    e.optLabels.addEventListener('change', () => this.updateLabelVisibility());
    e.optAstLabels.addEventListener('change', () => this.updateLabelVisibility());
    e.optMoon.addEventListener('change', () => { if (e.optMoon.checked) { e.follow.value = ''; a.followIndex = -1; } a.setZoom(e.optMoon.checked); e.legendMoon.style.display = e.optMoon.checked ? 'block' : 'none'; });
    e.togglePanel.addEventListener('click', () => e.panel.classList.toggle('collapsed'));
    window.addEventListener('keydown', (ev) => {
      if (ev.target && (ev.target.tagName === 'INPUT' || ev.target.tagName === 'SELECT')) return;
      if (ev.code === 'Space') { ev.preventDefault(); a.togglePlay(); }
      else if (ev.code === 'ArrowLeft') a.setTime(a.t - 10);
      else if (ev.code === 'ArrowRight') a.setTime(a.t + 10);
    });
    const scrub = (ev) => {
      const r = e.tl.getBoundingClientRect(); const x = (ev.clientX - r.left) / r.width;
      a.setTime(Math.max(0, Math.min(1, x)) * a.tmax, true);
    };
    let dragging = false;
    e.tl.addEventListener('pointerdown', (ev) => { dragging = true; e.tl.setPointerCapture(ev.pointerId); scrub(ev); });
    e.tl.addEventListener('pointermove', (ev) => { if (dragging) scrub(ev); });
    e.tl.addEventListener('pointerup', () => { dragging = false; });
    e.speed.value = String(this.sliderFromSpeed(30));
    e.speedVal.textContent = '30 d/s';
    if (window.innerWidth <= 760) e.panel.classList.add('collapsed');
  }

  speedFromSlider(v) { return SPEED_MIN * Math.pow(SPEED_MAX / SPEED_MIN, v / 1000); }
  sliderFromSpeed(s) { return Math.round(1000 * Math.log(s / SPEED_MIN) / Math.log(SPEED_MAX / SPEED_MIN)); }

  setStatus(msg, sticky) {
    this.el.status.textContent = msg; this.el.status.classList.toggle('show', !!msg);
    if (msg && !sticky) setTimeout(() => { if (this.el.status.textContent === msg) this.el.status.classList.remove('show'); }, 6000);
  }

  setPlaying(on) { this.el.play.textContent = on ? 'pause' : 'play'; }

  populateFleets(entries, current) {
    this.el.fleetSelect.innerHTML = '';
    for (const f of entries) {
      const o = document.createElement('option'); o.value = f.key;
      o.textContent = f.title || `${f.key}: ${f.n_craft} craft, ${f.covered} covered`;
      this.el.fleetSelect.appendChild(o);
    }
    this.el.fleetSelect.value = current;
  }

  /** Called after a fleet is loaded: craft table, follow menu, timeline statics, labels. */
  setFleet(fleet, scene) {
    const e = this.el, entry = fleet.entry || {};
    e.fleetSelect.value = fleet.key;
    const parts = [];
    if (entry.kind) parts.push(`kind: ${entry.kind}`);
    if (entry.submitted) parts.push(`submitted ${entry.submitted}`);
    if (entry.shown != null) parts.push(`shown ${entry.shown}`);
    if (fleet.source) parts.push(`source: ${fleet.source}`);
    if (entry.note) parts.push(entry.note);
    e.fleetNote.textContent = parts.join(' | ');
    e.craft.textContent = String(fleet.craft.length);
    e.j.textContent = fleet.J_raw != null ? Number(fleet.J_raw).toFixed(4) : '-';
    e.sumji.textContent = fleet.sumJi != null ? Number(fleet.sumJi).toFixed(4) : '-';
    e.time.max = String(fleet.tmax);
    // timeline strip height grows with the craft count (one row per craft), capped for phones
    const phone = window.innerWidth <= 760, n = fleet.craft.length;
    const tlh = Math.max(phone ? 96 : 118, Math.min(phone ? 130 : 230, 60 + 9 * n));
    document.documentElement.style.setProperty('--timeline-h', `${tlh}px`);
    this.tlDirty = true;
    e.follow.innerHTML = '<option value="">none (free camera)</option>';
    e.table.innerHTML = ''; this.rows = [];
    fleet.craft.forEach((c, i) => {
      const o = document.createElement('option'); o.value = String(i); o.textContent = `craft ${c.id}`; e.follow.appendChild(o);
      const tr = document.createElement('tr');
      tr.innerHTML = `<td><input type="checkbox" checked></td><td><span class="sw" style="background:${craftColor(i)}"></span>SC ${c.id}</td>` +
                     `<td class="num fb">0 / ${c.n_flybys}</td><td class="num fuel">0 / ${c.fuel_kg.toFixed(0)}</td><td><button class="small">f</button></td>`;
      const cb = tr.querySelector('input');
      cb.addEventListener('change', () => { c.visible = cb.checked; tr.classList.toggle('off', !cb.checked); this.tlDirty = true; });
      tr.querySelector('button').addEventListener('click', () => { e.follow.value = String(i); e.follow.dispatchEvent(new Event('change')); });
      e.table.appendChild(tr);
      this.rows.push({ tr, fb: tr.querySelector('.fb'), fuel: tr.querySelector('.fuel') });
    });
    // labels
    this.el.labels.innerHTML = ''; this.labelNodes.clear();
    this.addLabel('sun', 'Sun', '');
    this.addLabel('earth', 'Earth', '');
    this.addLabel('moon', 'Moon', '');
    fleet.craft.forEach((c, i) => { const n = this.addLabel(`craft${i}`, `SC ${c.id}`, 'craft'); n.style.color = craftColor(i); });
    scene.astIds.forEach((id) => this.addLabel(`ast${id}`, String(id), 'ast'));
    this.updateLabelVisibility();
    this.tlDirty = true;
  }

  addLabel(key, text, cls) {
    const d = document.createElement('div'); d.className = 'lbl ' + cls; d.textContent = text; this.el.labels.appendChild(d);
    this.labelNodes.set(key, d); return d;
  }
  updateLabelVisibility() {
    const on = this.el.optLabels.checked, ast = on && this.el.optAstLabels.checked;
    for (const [k, n] of this.labelNodes) n.classList.toggle('hidden', !(k.startsWith('ast') ? ast : on));
    this.el.labels.style.display = on ? 'block' : 'none';
  }

  /** Per-frame HUD numbers (cheap) and labels. */
  update(app) {
    const e = this.el, t = app.t, fleet = app.fleet;
    e.time.value = String(t);
    e.dateVal.textContent = mjdToUTC(app.t0Mjd + t);
    e.tdayVal.textContent = `t = ${t.toFixed(1)} d`;
    if (fleet) {
      let covered = 0; for (const [, tf] of fleet.firstSeen) if (tf <= t) covered++;
      e.covered.textContent = `${covered} / ${app.nReachable}`;
      let launched = 0, fuel = 0, interp = false;
      fleet.craft.forEach((c, i) => {
        const on = t >= c.launch; if (on) launched++;
        let fb = 0; for (const f of c.flybys) if (f.t <= t) fb++;
        const row = this.rows[i];
        if (row) { row.fb.textContent = `${fb} / ${c.n_flybys}`; row.tr.classList.toggle('done', fb === c.n_flybys && fb > 0); }
        // per-craft column: "used so far / total"; used = m0 - m(t) when a mass array exists (clamped at 0: the file's
        // m is rounded to 0.01 kg), else the loaded fuel once the craft is launched; 0 before launch
        // a finished craft reports its accounted fuel_kg (twin files: m ends at the 601.5 kg margin, fuel_kg = tank - 600)
        let used = 0;
        if (c.m) interp = true;
        if (on) {
          if (t >= c.end) used = c.fuel_kg;
          else if (c.m) used = Math.max(0, c.m0 - interpMass(c, t));
          else used = c.fuel_kg;
        }
        fuel += used;
        if (row) row.fuel.textContent = `${used.toFixed(0)} / ${c.fuel_kg.toFixed(0)}`;
      });
      e.launched.textContent = `${launched} / ${fleet.craft.length}`;
      e.fuel.textContent = `${Math.max(0, fuel).toFixed(1)} kg`;
      e.fuelLabel.textContent = interp ? 'fuel used (m(t))' : 'fuel loaded (launched)';
    }
    this.updateLabels(app);
    this.drawTimeline(app);
  }

  updateLabels(app) {
    if (!this.el.optLabels.checked) return;
    const cam = app.scene.camera, w = app.scene.container.clientWidth, h = app.scene.container.clientHeight;
    const v = this._v || (this._v = new THREE.Vector3());
    const place = (node, x, y, z) => {
      v.set(x, y, z).project(cam);
      if (v.z > 1 || v.x < -1.1 || v.x > 1.1 || v.y < -1.1 || v.y > 1.1) { node.style.display = 'none'; return; }
      node.style.display = '';
      node.style.transform = `translate(${((v.x + 1) * 0.5 * w).toFixed(0)}px, ${((1 - v.y) * 0.5 * h).toFixed(0)}px) translate(-50%, -140%)`;
    };
    const sc = app.scene;
    place(this.labelNodes.get('sun'), 0, 0, 0);
    place(this.labelNodes.get('earth'), sc.earthPos.x, sc.earthPos.y, sc.earthPos.z);
    const moon = this.labelNodes.get('moon');
    if (sc.zoomMode && sc.moon) place(moon, sc.earthPos.x + sc.moonPosGeo[0], sc.earthPos.y + sc.moonPosGeo[1], sc.earthPos.z + sc.moonPosGeo[2]);
    else moon.style.display = 'none';
    app.fleetLayer.items.forEach((it, i) => {
      const n = this.labelNodes.get(`craft${i}`);
      if (it.marker.visible && it.craft.visible) place(n, it.pos.x, it.pos.y, it.pos.z); else n.style.display = 'none';
    });
    if (this.el.optAstLabels.checked) {
      const a = sc.astPos;
      sc.astIds.forEach((id, k) => place(this.labelNodes.get(`ast${id}`), a[3 * k], a[3 * k + 1], a[3 * k + 2]));
    }
  }

  drawTimeline(app) {
    const c = this.el.tl, dpr = Math.min(window.devicePixelRatio || 1, 2);
    // size from the wrapper, capped: the canvas' CSS box must never depend on its own intrinsic size (feedback loop)
    const wrap = c.parentElement, w = Math.min(wrap.clientWidth, 4096), h = Math.min(wrap.clientHeight, 512);
    if (!w || !h) return;
    if (c.width !== Math.round(w * dpr) || c.height !== Math.round(h * dpr)) { c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); this.tlDirty = true; }
    if (this.tlDirty) this.drawTimelineStatic(app, w, h, dpr);
    if (!this.tlDirty && this.lastTl === app.t) return;
    this.lastTl = app.t; this.tlDirty = false;
    const g = c.getContext('2d');
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.drawImage(this.tlStatic, 0, 0);
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const x = 40 + (w - 48) * (app.t / app.tmax);
    g.strokeStyle = '#4cc9f0'; g.lineWidth = 1; g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke();
  }

  drawTimelineStatic(app, w, h, dpr) {
    const s = this.tlStatic; s.width = Math.round(w * dpr); s.height = Math.round(h * dpr);
    const g = s.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.fillStyle = '#0d111b'; g.fillRect(0, 0, w, h);
    const fleet = app.fleet; if (!fleet) return;
    const x0 = 40, x1 = w - 8, tx = (t) => x0 + (x1 - x0) * (t / app.tmax);
    // year ticks
    g.strokeStyle = '#1c2433'; g.fillStyle = '#8a94a6'; g.font = '9px Menlo, monospace'; g.textAlign = 'center';
    const yearStep = app.tmax > 2000 ? 365.25 : app.tmax > 300 ? 100 : 30;
    for (let t = 0; t <= app.tmax + 1e-6; t += yearStep) {
      const x = tx(t); g.beginPath(); g.moveTo(x, 0); g.lineTo(x, h); g.stroke();
      g.fillText(yearStep === 365.25 ? String(2030 + Math.round(t / 365.25)) : `${Math.round(t)} d`, x, h - 2);
    }
    const n = fleet.craft.length, rowH = Math.max(3, Math.min(12, (h - 14) / n));
    g.textAlign = 'right';
    fleet.craft.forEach((c, i) => {
      const y = 2 + i * rowH + rowH / 2, col = craftColor(i), alpha = c.visible ? 1 : 0.25;
      g.globalAlpha = alpha;
      g.fillStyle = '#8a94a6'; if (rowH >= 7) g.fillText(`SC${c.id}`, x0 - 4, y + 3);
      g.strokeStyle = col; g.lineWidth = Math.max(1, rowH * 0.18); g.beginPath(); g.moveTo(tx(c.launch), y); g.lineTo(tx(c.end), y); g.stroke();
      g.lineWidth = Math.max(2, rowH * 0.5);
      for (const [a, b] of c.thrust) { g.beginPath(); g.moveTo(tx(a), y); g.lineTo(tx(Math.max(b, a + app.tmax / (x1 - x0))), y); g.stroke(); }
      g.fillStyle = '#ffffff'; const th = Math.max(3, rowH * 0.8);
      for (const f of c.flybys) g.fillRect(tx(f.t) - 0.5, y - th / 2, 1, th);
      g.globalAlpha = 1;
    });
  }
}
