// Entry point: load data, build the scene and fleet layer, drive the clock and the render loop.
import * as THREE from 'three';
import { loadCatalog, loadFleet } from './data.js';
import { SolarScene } from './scene.js';
import { FleetLayer } from './fleet.js';
import { UI } from './ui.js';

const MISSION_END_D = 5478.75;

class App {
  constructor() {
    this.t = 0; this.speed = 30; this.playing = false; this.tmax = MISSION_END_D;
    this.fleet = null; this.followIndex = -1; this.zoom = false;
    this.lastFrame = performance.now(); this.frameTimes = []; this.fpsAcc = 0; this.fpsN = 0; this.fpsLast = performance.now();
    this._v = new THREE.Vector3();
  }

  async init() {
    this.ui = new UI(this);
    const cat = await loadCatalog();
    this.catalog = cat;
    this.t0Mjd = cat.bodies.epoch_mjd_t0;
    this.nReachable = cat.bodies.asteroids.length - (cat.bodies.unreachable || []).length;
    this.scene = new SolarScene(document.getElementById('view'), cat.bodies);
    this.fleetLayer = new FleetLayer(this.scene.scene, this.scene.camera);
    this.ui.setStatus(cat.label.startsWith('FIXTURE') ? 'Loaded FIXTURE data (real results/viz/data/index.json not found)' : '', false);
    window.addEventListener('resize', () => { this.scene.resize(); this.ui.tlDirty = true; });
    // requestAnimationFrame is suspended in hidden tabs: say so instead of showing 0 fps
    document.addEventListener('visibilitychange', () => { if (document.hidden) this.ui.el.fps.textContent = 'paused (tab hidden)'; });
    this.ui.populateFleets(cat.index.fleets, cat.index.fleets[0].key);
    await this.selectFleet(cat.index.fleets[0].key);
    this.scene.resize();
    requestAnimationFrame((now) => this.frame(now));
  }

  async selectFleet(key) {
    const entry = this.catalog.index.fleets.find(f => f.key === key);
    if (!entry) return;
    this.pause();
    try {
      const fleet = await loadFleet(this.catalog.dir, entry);
      this.fleet = fleet;
      this.tmax = Math.min(Math.max(fleet.tmax, 10), MISSION_END_D + 1);
      this.fleetLayer.setFleet(fleet);
      this.followIndex = -1;
      this.ui.setFleet(fleet, this.scene);
      this.setTime(0, true);
      this.ui.setStatus('', false);
    } catch (e) {
      console.error('[viewer] fleet load failed', e);
      this.ui.setStatus(`fleet ${key} failed to load: ${e.message}`, true);
    }
  }

  play() { this.playing = true; this.ui.setPlaying(true); if (this.t >= this.tmax) this.setTime(0, true); }
  pause() { this.playing = false; this.ui.setPlaying(false); }
  togglePlay() { this.playing ? this.pause() : this.play(); }

  /** Set the clock; fires flyby flashes only for events crossed forward by a small step (playback), not for jumps. */
  setTime(t, jump = false) {
    const prev = this.t;
    t = Math.max(0, Math.min(this.tmax, t));
    const small = !jump && t > prev && (t - prev) < 60;
    this.t = t;
    this.scene.update(t);
    this.fleetLayer.update(t);
    if (this.fleet) {
      if (small) {
        for (const f of this.fleet.allFlybys) {
          if (f.t > prev && f.t <= t && this.fleet.craft[f.craft].visible) {
            const p = this.scene.asteroidPosition(f.ast, this._v);
            if (p) this.fleetLayer.flash(p, f.craft);
          }
        }
      }
      this.applyAsteroidStates(t);
    }
    if (t >= this.tmax && this.playing) this.pause();
  }

  applyAsteroidStates(t) {
    const sc = this.scene; let changed = false;
    for (let k = 0; k < sc.astIds.length; k++) {
      const id = sc.astIds[k];
      if (sc.unreachable.has(id)) continue;
      const tf = this.fleet.firstSeen.get(id);
      const st = tf !== undefined && tf <= t ? 1 : 0;
      if (sc.astState[k] !== st) { sc.setAsteroidState(k, st, true); changed = true; }
    }
    if (changed) sc.astPoints.geometry.getAttribute('color').needsUpdate = true;
  }

  setZoom(on) { this.zoom = on; this.scene.setZoomMode(on); }

  frame(now) {
    const dt = Math.min(0.1, (now - this.lastFrame) / 1000); this.lastFrame = now;
    if (this.playing) this.setTime(this.t + this.speed * dt);
    else this.fleetLayer.update(this.t);   // keep flash animation alive when paused
    if (this.followIndex >= 0 && this.fleetLayer.items[this.followIndex]) {
      const it = this.fleetLayer.items[this.followIndex];
      if (it.launched) this.scene.follow(it.pos);
    } else if (this.zoom) this.scene.follow(this.scene.earthPos);
    this.scene.render();
    this.ui.update(this);
    // fps readout (1 s average)
    this.fpsN++;
    if (now - this.fpsLast > 1000) { this.fps = this.fpsN * 1000 / (now - this.fpsLast); this.ui.el.fps.textContent = `${this.fps.toFixed(0)} fps`; this.fpsN = 0; this.fpsLast = now; }
    requestAnimationFrame((n) => this.frame(n));
  }
}

const app = new App();
window.__viz = app;   // for console inspection / automated checks
app.init().catch(e => {
  console.error('[viewer] init failed', e);
  const s = document.getElementById('status'); s.textContent = 'init failed: ' + e.message; s.classList.add('show');
});
