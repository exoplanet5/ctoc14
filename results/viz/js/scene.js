// Static solar-system scene: Sun, ecliptic grid, Earth (+ orbit), Moon (Earth-Moon zoom), 300 asteroids (orbits + points).
// 1 scene unit = 1 AU, ecliptic J2000: x toward the vernal equinox, z = ecliptic north (camera.up = +z).
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { KeplerBody, MoonOrbit } from './kepler.js';

export const COLORS = {
  astGrey: new THREE.Color(0x9aa3b2), astGreen: new THREE.Color(0x4ade80), astRed: new THREE.Color(0xf87171),
};
const EARTH_R_NORMAL = 0.012, EARTH_R_ZOOM = 0.00030;   // visual radii in AU (real: 4.26e-5 AU)
const MOON_R_ZOOM = 0.00010;                              // visual (real 1.16e-5 AU); Moon distance is true scale
const ORBIT_PTS = 180;

function glowTexture() {
  const c = document.createElement('canvas'); c.width = c.height = 128;
  const g = c.getContext('2d'), grd = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grd.addColorStop(0, 'rgba(255,240,200,1)'); grd.addColorStop(0.25, 'rgba(255,200,90,0.55)');
  grd.addColorStop(0.6, 'rgba(255,150,40,0.12)'); grd.addColorStop(1, 'rgba(255,120,20,0)');
  g.fillStyle = grd; g.fillRect(0, 0, 128, 128);
  const tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace; return tex;
}

function circleGeometry(radius, n = 256) {
  const a = new Float32Array(n * 3);
  for (let j = 0; j < n; j++) { const th = 2 * Math.PI * j / n; a[3 * j] = radius * Math.cos(th); a[3 * j + 1] = radius * Math.sin(th); }
  const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(a, 3)); return g;
}

export class SolarScene {
  constructor(container, bodies) {
    this.container = container;
    this.bodies = bodies;
    this.t0Mjd = bodies.epoch_mjd_t0;
    this.mu = bodies.mu_sun_km3s2; this.au = bodies.au_km;

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false, powerPreference: 'high-performance' });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    this.renderer.setClearColor(0x07090f, 1);
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(45, 1, 0.002, 400);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(2.4, -3.2, 2.2);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true; this.controls.dampingFactor = 0.08;
    this.controls.minDistance = 0.05; this.controls.maxDistance = 60;
    this.controls.target.set(0, 0, 0);

    // lights
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.35));
    const sunLight = new THREE.PointLight(0xfff4e0, 4.0, 0, 2); this.scene.add(sunLight);

    // Sun
    this.sunMesh = new THREE.Mesh(new THREE.SphereGeometry(0.025, 32, 16), new THREE.MeshBasicMaterial({ color: 0xfff1b0 }));
    this.scene.add(this.sunMesh);
    this.sunGlow = new THREE.Sprite(new THREE.SpriteMaterial({ map: glowTexture(), blending: THREE.AdditiveBlending, depthWrite: false, transparent: true }));
    this.sunGlow.scale.set(0.35, 0.35, 1); this.scene.add(this.sunGlow);

    // ecliptic reference circles + vernal-equinox direction
    this.grid = new THREE.Group();
    for (const r of [0.5, 1.0, 1.5, 2.0]) {
      const mat = new THREE.LineBasicMaterial({ color: 0x2a3446, transparent: true, opacity: r === 1.0 ? 0.55 : 0.4 });
      this.grid.add(new THREE.LineLoop(circleGeometry(r), mat));
    }
    const axG = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, 0, 0), new THREE.Vector3(2.2, 0, 0)]);
    this.grid.add(new THREE.Line(axG, new THREE.LineBasicMaterial({ color: 0x2a3446, transparent: true, opacity: 0.5 })));
    this.scene.add(this.grid);

    // Earth + orbit
    this.earthBody = new KeplerBody(bodies.earth, this.t0Mjd, this.mu, this.au);
    const eo = new THREE.BufferGeometry(); eo.setAttribute('position', new THREE.BufferAttribute(this.earthBody.orbitPoints(720), 3));
    this.earthOrbit = new THREE.LineLoop(eo, new THREE.LineBasicMaterial({ color: 0x3b82f6, transparent: true, opacity: 0.7 }));
    this.scene.add(this.earthOrbit);
    this.earthGroup = new THREE.Group(); this.scene.add(this.earthGroup);
    this.earthMesh = new THREE.Mesh(new THREE.SphereGeometry(1, 32, 16),
      new THREE.MeshPhongMaterial({ color: 0x3b82f6, emissive: 0x0b2a5a, shininess: 30 }));
    this.earthMesh.scale.setScalar(EARTH_R_NORMAL); this.earthGroup.add(this.earthMesh);

    // Moon (shown in Earth-Moon zoom only; in the normal view it sits inside the enlarged Earth marker)
    this.moon = bodies.moon ? new MoonOrbit(bodies.moon, this.t0Mjd, this.au) : null;
    this.moonGroup = new THREE.Group(); this.moonGroup.visible = false; this.earthGroup.add(this.moonGroup);
    if (this.moon) {
      this.moonMesh = new THREE.Mesh(new THREE.SphereGeometry(1, 24, 12), new THREE.MeshPhongMaterial({ color: 0xbdbdbd, emissive: 0x222222 }));
      this.moonMesh.scale.setScalar(MOON_R_ZOOM); this.moonGroup.add(this.moonMesh);
      this.moonRingArr = new Float32Array(128 * 3);
      const mg = new THREE.BufferGeometry(); mg.setAttribute('position', new THREE.BufferAttribute(this.moonRingArr, 3));
      this.moonRing = new THREE.LineLoop(mg, new THREE.LineBasicMaterial({ color: 0x8a94a6, transparent: true, opacity: 0.6 }));
      this.moonGroup.add(this.moonRing);
    }

    // asteroids
    this.unreachable = new Set(bodies.unreachable || []);
    this.astBodies = bodies.asteroids.map(a => new KeplerBody(a, this.t0Mjd, this.mu, this.au));
    this.astIds = bodies.asteroids.map(a => a.id);
    this.astIndexById = new Map(this.astIds.map((id, k) => [id, k]));
    const N = this.astBodies.length;
    this.astPos = new Float32Array(N * 3);
    this.astCol = new Float32Array(N * 3);
    this.astState = new Uint8Array(N);   // 0 grey, 1 green, 2 red
    const pg = new THREE.BufferGeometry();
    pg.setAttribute('position', new THREE.BufferAttribute(this.astPos, 3).setUsage(THREE.DynamicDrawUsage));
    pg.setAttribute('color', new THREE.BufferAttribute(this.astCol, 3).setUsage(THREE.DynamicDrawUsage));
    this.astPoints = new THREE.Points(pg, new THREE.PointsMaterial({ size: 3.5, sizeAttenuation: false, vertexColors: true }));
    this.astPoints.frustumCulled = false;
    this.scene.add(this.astPoints);
    this.resetAsteroidColors();

    // orbit lines: reachable in one LineSegments, unreachable as red dashed lines
    const reach = [], unreach = [];
    this.astBodies.forEach((b, k) => (this.unreachable.has(this.astIds[k]) ? unreach : reach).push(b));
    const segArr = new Float32Array(reach.length * ORBIT_PTS * 2 * 3);
    let o = 0;
    for (const b of reach) {
      const p = b.orbitPoints(ORBIT_PTS);
      for (let j = 0; j < ORBIT_PTS; j++) {
        const j2 = (j + 1) % ORBIT_PTS;
        segArr[o++] = p[3 * j]; segArr[o++] = p[3 * j + 1]; segArr[o++] = p[3 * j + 2];
        segArr[o++] = p[3 * j2]; segArr[o++] = p[3 * j2 + 1]; segArr[o++] = p[3 * j2 + 2];
      }
    }
    const sg = new THREE.BufferGeometry(); sg.setAttribute('position', new THREE.BufferAttribute(segArr, 3));
    this.orbitLines = new THREE.Group();
    this.orbitLines.add(new THREE.LineSegments(sg, new THREE.LineBasicMaterial({ color: 0x6b7280, transparent: true, opacity: 0.22 })));
    for (const b of unreach) {
      const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.BufferAttribute(b.orbitPoints(ORBIT_PTS), 3));
      const l = new THREE.LineLoop(g, new THREE.LineDashedMaterial({ color: 0xf87171, dashSize: 0.03, gapSize: 0.02, transparent: true, opacity: 0.8 }));
      l.computeLineDistances(); this.orbitLines.add(l);
    }
    this.scene.add(this.orbitLines);

    this.earthPos = new THREE.Vector3();
    this.moonPosGeo = [0, 0, 0];
    this.zoomMode = false;
    this._tmp = [0, 0, 0];
    this.resize();
  }

  resize() {
    const w = this.container.clientWidth || 1, h = this.container.clientHeight || 1;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h; this.camera.updateProjectionMatrix();
  }

  resetAsteroidColors() {
    for (let k = 0; k < this.astIds.length; k++) this.setAsteroidState(k, this.unreachable.has(this.astIds[k]) ? 2 : 0, true);
    this.astPoints.geometry.getAttribute('color').needsUpdate = true;
  }
  setAsteroidState(k, state, silent) {
    if (this.astState[k] === state && !silent) return;
    this.astState[k] = state;
    const c = state === 1 ? COLORS.astGreen : state === 2 ? COLORS.astRed : COLORS.astGrey;
    this.astCol[3 * k] = c.r; this.astCol[3 * k + 1] = c.g; this.astCol[3 * k + 2] = c.b;
    if (!silent) this.astPoints.geometry.getAttribute('color').needsUpdate = true;
  }

  /** Propagate Earth, Moon and all asteroids to tDays. */
  update(tDays) {
    const p = this._tmp;
    this.lastT = tDays;
    this.earthBody.positionAU(tDays, p);
    this.earthPos.set(p[0], p[1], p[2]);
    this.earthGroup.position.copy(this.earthPos);
    if (this.moon && this.zoomMode) {
      this.moon.positionGeoAU(tDays, this.moonPosGeo);
      this.moonMesh.position.set(this.moonPosGeo[0], this.moonPosGeo[1], this.moonPosGeo[2]);
      this.moon.orbitPoints(tDays, 128, this.moonRingArr);
      this.moonRing.geometry.getAttribute('position').needsUpdate = true;
    }
    const N = this.astBodies.length, a = this.astPos;
    for (let k = 0; k < N; k++) {
      this.astBodies[k].positionAU(tDays, p);
      a[3 * k] = p[0]; a[3 * k + 1] = p[1]; a[3 * k + 2] = p[2];
    }
    this.astPoints.geometry.getAttribute('position').needsUpdate = true;
  }

  asteroidPosition(id, out) {
    const k = this.astIndexById.get(id);
    if (k === undefined) return null;
    out.set(this.astPos[3 * k], this.astPos[3 * k + 1], this.astPos[3 * k + 2]);
    return out;
  }

  setZoomMode(on) {
    if (on === this.zoomMode) return;
    this.zoomMode = on;
    this.moonGroup.visible = on && !!this.moon;
    this.earthMesh.scale.setScalar(on ? EARTH_R_ZOOM : EARTH_R_NORMAL);
    // the Moon is only propagated in zoom mode: place it now, not at the next clock change (else it sits inside Earth while paused)
    if (on && this.lastT != null) this.update(this.lastT);
    if (on) {
      this.camera.near = 2e-5; this.camera.far = 50; this.controls.minDistance = 4e-4; this.controls.maxDistance = 0.2;
      const dir = new THREE.Vector3().subVectors(this.camera.position, this.controls.target).normalize();
      if (!(dir.lengthSq() > 0)) dir.set(0.4, -0.7, 0.5);
      this.controls.target.copy(this.earthPos);
      this.camera.position.copy(this.earthPos).addScaledVector(dir, 0.011);
    } else {
      this.camera.near = 0.002; this.camera.far = 400; this.controls.minDistance = 0.05; this.controls.maxDistance = 60;
      const dir = new THREE.Vector3().subVectors(this.camera.position, this.controls.target).normalize();
      this.camera.position.copy(this.controls.target).addScaledVector(dir, 4.5);
    }
    this.camera.updateProjectionMatrix();
    this.controls.update();
  }

  /** Move controls target to `target` keeping the camera offset (follow mode). */
  follow(target) {
    const d = new THREE.Vector3().subVectors(target, this.controls.target);
    if (d.lengthSq() === 0) return;
    this.controls.target.copy(target);
    this.camera.position.add(d);
  }

  homeView() {
    this.controls.target.set(0, 0, 0);
    this.camera.position.set(2.4, -3.2, 2.2);
    this.controls.update();
  }
  topView() {
    this.controls.target.set(0, 0, 0);
    this.camera.position.set(0, -0.001, 5.5);
    this.controls.update();
  }

  render() {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }
}
