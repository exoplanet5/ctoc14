// Fleet layer: per-craft trajectory (past trail bright, future faint), thrust arcs, craft marker, flyby flash rings.
import * as THREE from 'three';
import { interpPosition } from './data.js';

export const CRAFT_PALETTE = ['#4cc9f0', '#f72585', '#ffd166', '#06d6a0', '#f4a261', '#a78bfa', '#ef476f', '#80ed99',
                              '#ff9f1c', '#48bfe3', '#e0aaff', '#b5e48c'];
export function craftColor(i) { return CRAFT_PALETTE[i % CRAFT_PALETTE.length]; }

const FLASH_SECONDS = 1.4;

function markerTexture(hex) {
  const c = document.createElement('canvas'); c.width = c.height = 32;
  const g = c.getContext('2d');
  g.fillStyle = '#000'; g.beginPath(); g.arc(16, 16, 11, 0, 2 * Math.PI); g.fill();
  g.fillStyle = hex; g.beginPath(); g.arc(16, 16, 9, 0, 2 * Math.PI); g.fill();
  g.fillStyle = '#fff'; g.beginPath(); g.arc(16, 16, 3, 0, 2 * Math.PI); g.fill();
  const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace; return t;
}

export class FleetLayer {
  constructor(scene, camera) {
    this.scene = scene; this.camera = camera;
    this.group = new THREE.Group(); scene.add(this.group);
    this.flashes = [];
    this.flashGeom = new THREE.RingGeometry(0.78, 1.0, 48);
    this.fleet = null;
    this.items = [];
    this.tmp = new THREE.Vector3();
    this._p = [0, 0, 0];
  }

  clear() {
    for (const it of this.items) {
      this.group.remove(it.node);
      it.node.traverse(o => { if (o.geometry) o.geometry.dispose(); if (o.material) { if (o.material.map) o.material.map.dispose(); o.material.dispose(); } });
    }
    for (const f of this.flashes) { this.group.remove(f.mesh); f.mesh.material.dispose(); }
    this.items = []; this.flashes = []; this.fleet = null;
  }

  setFleet(fleet) {
    this.clear();
    this.fleet = fleet;
    fleet.craft.forEach((c, ci) => {
      const color = new THREE.Color(craftColor(ci));
      const node = new THREE.Group();
      const posAttr = new THREE.BufferAttribute(Float32Array.from(c.xyz), 3);
      // full path, faint
      const gFull = new THREE.BufferGeometry(); gFull.setAttribute('position', posAttr);
      const future = new THREE.Line(gFull, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.16 }));
      // flown part, bright (drawRange driven)
      const gPast = new THREE.BufferGeometry(); gPast.setAttribute('position', posAttr); gPast.setDrawRange(0, 0);
      const past = new THREE.Line(gPast, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.9 }));
      // current partial segment (last sample -> interpolated position)
      const gSeg = new THREE.BufferGeometry(); gSeg.setAttribute('position', new THREE.BufferAttribute(new Float32Array(6), 3).setUsage(THREE.DynamicDrawUsage));
      const seg = new THREE.Line(gSeg, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.9 }));
      // thrust arcs: segments (p_k, p_{k+1}) for samples inside a thrust interval, ordered by time; drawRange by time
      const tSeg = [], tSegTime = [];
      for (const [a, b] of c.thrust) {
        for (let k = 0; k < c.n - 1; k++) {
          if (c.t[k] >= a - 1e-9 && c.t[k + 1] <= b + 1e-9) {
            tSeg.push(c.xyz[3 * k], c.xyz[3 * k + 1], c.xyz[3 * k + 2], c.xyz[3 * k + 3], c.xyz[3 * k + 4], c.xyz[3 * k + 5]);
            tSegTime.push(c.t[k + 1]);
          }
        }
      }
      const gThr = new THREE.BufferGeometry(); gThr.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(tSeg), 3)); gThr.setDrawRange(0, 0);
      const thrCol = color.clone().lerp(new THREE.Color(0xffffff), 0.55);
      const thrust = new THREE.LineSegments(gThr, new THREE.LineBasicMaterial({ color: thrCol, transparent: true, opacity: 1.0 }));
      const thrustFuture = new THREE.LineSegments(gThr.clone(), new THREE.LineBasicMaterial({ color: thrCol, transparent: true, opacity: 0.3 }));
      thrustFuture.geometry.setDrawRange(0, Infinity);
      // marker
      const marker = new THREE.Sprite(new THREE.SpriteMaterial({ map: markerTexture(craftColor(ci)), sizeAttenuation: false, depthTest: false, transparent: true }));
      marker.scale.set(0.02, 0.02, 1); marker.renderOrder = 10; marker.visible = false;
      [future, past, seg, thrust, thrustFuture, marker].forEach(o => { o.frustumCulled = false; node.add(o); });
      this.group.add(node);
      this.items.push({ craft: c, node, past, seg, thrust, thrustFuture, marker, tSegTime: Float64Array.from(tSegTime), color, pos: new THREE.Vector3(), launched: false });
    });
  }

  /** Update trails/markers to time t (days). Returns nothing; positions readable in items[i].pos. */
  update(t) {
    if (!this.fleet) return;
    const p = this._p;
    for (const it of this.items) {
      const c = it.craft;
      it.node.visible = c.visible;
      if (!c.visible) continue;
      const k = interpPosition(c, t, p);
      it.pos.set(p[0], p[1], p[2]);
      it.launched = t >= c.launch - 1e-9 && k >= 0;
      it.ended = t > c.end + 1e-9;
      it.marker.visible = it.launched && !it.ended;
      it.marker.position.copy(it.pos);
      const cnt = k < 0 ? 0 : Math.min(k + 1, c.n);
      it.past.geometry.setDrawRange(0, cnt);
      const sa = it.seg.geometry.getAttribute('position');
      if (k >= 0 && k < c.n - 1) {
        sa.array[0] = c.xyz[3 * k]; sa.array[1] = c.xyz[3 * k + 1]; sa.array[2] = c.xyz[3 * k + 2];
        sa.array[3] = p[0]; sa.array[4] = p[1]; sa.array[5] = p[2];
        sa.needsUpdate = true; it.seg.visible = true;
      } else it.seg.visible = false;
      // thrust segments with end time <= t
      const ts = it.tSegTime; let n = 0;
      if (ts.length) { let lo = 0, hi = ts.length; while (lo < hi) { const m = (lo + hi) >> 1; if (ts[m] <= t) lo = m + 1; else hi = m; } n = lo; }
      it.thrust.geometry.setDrawRange(0, n * 2);
    }
    // flashes (wall-clock animation)
    const now = performance.now() / 1000;
    for (let i = this.flashes.length - 1; i >= 0; i--) {
      const f = this.flashes[i], age = now - f.start;
      if (age > FLASH_SECONDS) { this.group.remove(f.mesh); f.mesh.material.dispose(); this.flashes.splice(i, 1); continue; }
      const u = age / FLASH_SECONDS;
      const dist = this.camera.position.distanceTo(f.mesh.position);
      const s = dist * (0.012 + 0.05 * u);
      f.mesh.scale.set(s, s, 1);
      f.mesh.material.opacity = 1 - u;
      f.mesh.quaternion.copy(this.camera.quaternion);
    }
  }

  /** Spawn a flash ring at position (Vector3) with the craft colour. */
  flash(position, craftIndex) {
    const mat = new THREE.MeshBasicMaterial({ color: craftColor(craftIndex), transparent: true, opacity: 1, side: THREE.DoubleSide, depthTest: false });
    const mesh = new THREE.Mesh(this.flashGeom, mat);
    mesh.position.copy(position); mesh.renderOrder = 20; mesh.frustumCulled = false;
    this.group.add(mesh);
    this.flashes.push({ mesh, start: performance.now() / 1000 });
    if (this.flashes.length > 40) { const f = this.flashes.shift(); this.group.remove(f.mesh); f.mesh.material.dispose(); }
  }
}
