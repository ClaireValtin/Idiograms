// Microscope image engine shared by the Microscope, Abnormalities and FISH pages.
// A scene is a list of chromosomes (band pieces) placed in CSS pixels. The image is built in cached stages,
// so a slider only recomputes what it changes:
//   field (illumination, stain unevenness, cytoplasm, grain)  <- size, seed
//   geometry (chromatids, banding, fibrous texture, halo, FISH spots)  <- scene, condensation, separation, bending, banding
//   absorbance / fluorescence  <- staining, field
//   focus blur (three box passes, cost independent of the blur radius)  <- focus
//   colour (lookup tables)  <- image type
// While a slider is dragged, pages render at half resolution (q = 0.5) and at full resolution on release.
const SCOPE = (() => {
  function rng(seed) { let a = seed >>> 0; return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
  function noise1(r, n = 64) { const v = new Float32Array(n + 1); for (let i = 0; i <= n; i++) v[i] = r() * 2 - 1; return x => { x = Math.max(0, Math.min(0.9999, x)) * n; const i = Math.floor(x), f = x - i, s = f * f * (3 - 2 * f); return v[i] * (1 - s) + v[i + 1] * s; }; }
  // 2D value noise in [0,1], 256-periodic
  function noise2(r) {
    const N = 256, g = new Float32Array(N * N); for (let i = 0; i < g.length; i++) g[i] = r();
    return (x, y) => {
      let xi = Math.floor(x), yi = Math.floor(y); const xf = x - xi, yf = y - yi; xi &= 255; yi &= 255;
      const x1 = (xi + 1) & 255, y1 = (yi + 1) & 255, a = g[yi * N + xi], b = g[yi * N + x1], c = g[y1 * N + xi], d = g[y1 * N + x1];
      const sx = xf * xf * (3 - 2 * xf), sy = yf * yf * (3 - 2 * yf);
      return a + (b - a) * sx + (c - a) * sy + (a - b - c + d) * sx * sy;
    };
  }
  const smooth = (a, b, x) => { const t = Math.max(0, Math.min(1, (x - a) / (b - a))); return t * t * (3 - 2 * t); };
  // ---------------------------------------------------------------- blur: 3 box passes approximate a Gaussian
  function boxes(sigma) {
    const n = 3, wIdeal = Math.sqrt(12 * sigma * sigma / n + 1); let wl = Math.floor(wIdeal); if (wl % 2 === 0) wl--;
    const m = Math.round((12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4));
    return [0, 1, 2].map(i => ((i < m ? wl : wl + 2) - 1) >> 1);
  }
  function boxLine(src, dst, off, stride, n, r) {
    const inv = 1 / (2 * r + 1), first = src[off], last = src[off + (n - 1) * stride];
    let acc = (r + 1) * first; for (let j = 0; j < r; j++) acc += src[off + Math.min(j, n - 1) * stride];
    for (let i = 0; i < n; i++) { const ia = i + r, ib = i - r - 1; acc += (ia < n ? src[off + ia * stride] : last) - (ib >= 0 ? src[off + ib * stride] : first); dst[off + i * stride] = acc * inv; }
  }
  function blur2(buf, tmp, W, H, sigma) {
    if (sigma < 0.4) return;
    for (const r of boxes(sigma)) { if (r < 1) continue; for (let y = 0; y < H; y++) boxLine(buf, tmp, y * W, 1, W, r); for (let x = 0; x < W; x++) boxLine(tmp, buf, x, W, H, r); }
  }
  // Gaussian blur of a copy; large radii are blurred at reduced resolution and upsampled (same look, a fraction of the cost)
  function blurScaled(src, W, H, sigma) {
    const f = sigma >= 8 ? 4 : sigma >= 3 ? 2 : 1;
    if (f === 1) { const out = Float32Array.from(src); blur2(out, new Float32Array(W * H), W, H, sigma); return out; }
    const w = Math.ceil(W / f), h = Math.ceil(H / f), small = new Float32Array(w * h), inv = 1 / (f * f);
    for (let y = 0; y < H; y++) { const sy = ((y / f) | 0) * w, o = y * W; for (let x = 0; x < W; x++) small[sy + ((x / f) | 0)] += src[o + x]; }
    // edge cells may hold fewer pixels; close enough for a blur
    for (let i = 0; i < small.length; i++) small[i] *= inv;
    blur2(small, new Float32Array(w * h), w, h, sigma / f);
    const out = new Float32Array(W * H);
    for (let y = 0; y < H; y++) { const fy = Math.max(0, (y + 0.5) / f - 0.5), j = Math.min(h - 1, fy | 0), j1 = Math.min(h - 1, j + 1), ty = fy - j, o = y * W;
      for (let x = 0; x < W; x++) { const fx = Math.max(0, (x + 0.5) / f - 0.5), i = Math.min(w - 1, fx | 0), i1 = Math.min(w - 1, i + 1), tx = fx - i;
        out[o + x] = (small[j * w + i] * (1 - tx) + small[j * w + i1] * tx) * (1 - ty) + (small[j1 * w + i] * (1 - tx) + small[j1 * w + i1] * tx) * ty; } }
    return out;
  }
  function blur1(src, sigma) { const out = Float32Array.from(src), tmp = new Float32Array(src.length); if (sigma < 0.4) return out; for (const r of boxes(sigma)) { if (r < 1) continue; boxLine(out, tmp, 0, 1, out.length, r); out.set(tmp); } return out; }
  // low-frequency field: evaluate on a coarse grid, bilinear upsample
  function lowField(W, H, step, fn) {
    const gw = Math.ceil(W / step) + 2, gh = Math.ceil(H / step) + 2, g = new Float32Array(gw * gh), out = new Float32Array(W * H);
    for (let j = 0; j < gh; j++) for (let i = 0; i < gw; i++) g[j * gw + i] = fn(i * step, j * step);
    for (let y = 0; y < H; y++) { const fy = y / step, j = Math.floor(fy), ty = fy - j; for (let x = 0; x < W; x++) { const fx = x / step, i = Math.floor(fx), tx = fx - i, o = j * gw + i;
      out[y * W + x] = (g[o] * (1 - tx) + g[o + 1] * tx) * (1 - ty) + (g[o + gw] * (1 - tx) + g[o + gw + 1] * tx) * ty; } }
    return out;
  }

  // ---------------------------------------------------------------- band pieces
  const Lbp = c => DATA[c]["850"].at(-1)[4];
  const cenBp = c => DATA[c]["550"].find(b => b[0].startsWith("q"))[3];
  // segs: {c,s,e} (s>e = inverted) or [c,s,e,dir]; "?" pieces {c:"?",len} are material of unknown origin (len in bp).
  // Each piece carries d0/d1, its span in derivative base pairs, so FISH signals can be placed on it.
  function pieces(segs, lv = "850") {
    const out = []; let seg = 0, d = 0;
    for (const raw of segs) {
      const sg = Array.isArray(raw) ? { c: raw[0], s: raw[3] > 0 ? raw[1] : raw[2], e: raw[3] > 0 ? raw[2] : raw[1] } : raw;
      if (sg.c === "?") { const sc = DATA["1"][lv].at(-1)[2] / Lbp("1"); out.push({ c: "?", name: "?", st: "unknown", len: sg.len * sc, seg, d0: d, d1: d + sg.len }); d += sg.len; seg++; continue; }
      const lo = Math.min(sg.s, sg.e), hi = Math.max(sg.s, sg.e), part = [];
      for (const b of DATA[sg.c][lv]) {
        if (b[4] <= lo || b[3] >= hi) continue;
        const a = Math.max(lo, b[3]), z = Math.min(hi, b[4]), f = v => b[1] + (v - b[3]) / Math.max(1, b[4] - b[3]) * (b[2] - b[1]);
        const len = f(z) - f(a); if (len > 0.01) part.push({ c: sg.c, name: b[0], st: b[5], len, seg, a, z });
      }
      if (sg.s > sg.e) part.reverse();
      for (const p of part) { const span = p.z - p.a; p.d0 = d; p.d1 = d + span; d += span; delete p.a; delete p.z; out.push(p); }
      seg++;
    }
    return out;
  }
  const whole = c => pieces([{ c, s: 0, e: Lbp(c) }]);
  // derivative base-pair position -> ISCN units along the piece list
  function unitAt(pcs, dbp) {
    let u = 0;
    for (const p of pcs) { if (dbp <= p.d1 || p === pcs.at(-1)) return u + p.len * Math.max(0, Math.min(1, (dbp - p.d0) / Math.max(1, p.d1 - p.d0))); u += p.len; }
    return u;
  }
  const totalOf = pcs => pcs.reduce((a, p) => a + p.len, 0);
  const ACRO = new Set(["13", "14", "15", "21", "22"]);

  // ---------------------------------------------------------------- banding
  // darkness of each NCBI stain class in R- and G-banding (G is used for the DAPI pattern too)
  const DARK = {
    R: { gneg: 1, gpos25: 0.66, gpos50: 0.42, gpos75: 0.24, gpos100: 0.08, acen: 0.5, stalk: 0.22 },
    G: { gneg: 0.07, gpos25: 0.3, gpos50: 0.55, gpos75: 0.78, gpos100: 1, acen: 0.4, stalk: 0.2 },
  };
  function hetTone(mode, het, c) {
    if (het === "pale") return 0.12; if (het === "dark") return 0.9;
    if (mode === "R") return 0.14;
    return c === "1" ? 0.9 : c === "16" ? 0.85 : c === "9" ? 0.2 : c === "Y" ? 0.78 : 0.5;
  }
  // heteromorphism: heterochromatin and stalk length differ between homologues; acrocentric short arms are drawn
  // at about half their idiogram length, closer to how they appear on slides
  function vary(pcs, r) {
    return pcs.map(p => { let f = 1; const x = r();
      if (p.st === "gvar") f = 0.72 + x * 0.62; else if (p.st === "stalk") f = 0.7 + x * 0.7;
      if (ACRO.has(p.c) && p.name[0] === "p" && (p.st === "gvar" || p.st === "stalk")) f *= 0.55;
      return { ...p, len0: p.len, len: p.len * f }; });
  }
  function mapU(pcs, u0) { // original units -> varied units
    let a = 0, b = 0;
    for (const p of pcs) { const l0 = p.len0 ?? p.len; if (u0 <= a + l0) return b + (u0 - a) / l0 * p.len; a += l0; b += p.len; }
    return b;
  }
  function profile(pcs, r, mode, het, sigma) {
    const total = Math.max(2, totalOf(pcs)), n = Math.ceil(total) + 1, prof = new Float32Array(n), stalk = new Float32Array(n), D = DARK[mode];
    let pos = 0, unk = 0;
    pcs.forEach((p, i) => {
      const x = r();
      let v = p.st === "gvar" ? hetTone(mode, het, p.c) : p.st === "unknown" ? (unk++ % 2 ? 0.75 : 0.25) : (D[p.st] ?? 0.5);
      if (mode === "R" && p.st === "gneg") v *= 0.88 + x * 0.24;
      if (mode === "G" && p.st.startsWith("gpos")) v *= 0.9 + x * 0.2;
      if (mode === "R" && (i === 0 || i === pcs.length - 1) && p.st === "gneg") v = Math.min(1.15, v * 1.12);   // R-banded ends are almost always positive
      if (p.st === "unknown") { const step = 70; for (let k = Math.floor(pos); k < Math.min(n, Math.ceil(pos + p.len)); k++) prof[k] = ((k / step) | 0) % 2 ? 0.75 : 0.3; }
      else for (let k = Math.floor(pos); k < Math.min(n, Math.ceil(pos + p.len)); k++) { prof[k] = v; if (p.st === "stalk") stalk[k] = 1; }
      pos += p.len;
    });
    const ac = []; pos = 0; for (const p of pcs) { if (p.st === "acen") ac.push([pos, pos + p.len]); pos += p.len; }
    const tc = ac.length ? (ac[0][0] + ac.at(-1)[1]) / 2 / total : 0.5;
    return { prof: blur1(prof, sigma), stalk: blur1(stalk, Math.max(20, sigma * 0.6)), total, n, tc, hasCen: ac.length > 0 };
  }

  // ---------------------------------------------------------------- geometry
  function axis(r, len, bend) {
    const n1 = noise1(r, 6), n2 = noise1(r, 14), th0 = (r() - 0.5) * 0.6 * bend, steps = Math.max(8, Math.ceil(len)), pts = [];
    let x = 0, y = 0;
    for (let i = 0; i <= steps; i++) { const t = i / steps, th = th0 + bend * (0.9 * n1(t) + 0.25 * n2(t)); pts.push({ x, y, th, t, s: t * len }); x += Math.sin(th) * len / steps; y += Math.cos(th) * len / steps; }
    let mn = Infinity, mx = -Infinity; for (const p of pts) { mn = Math.min(mn, p.x); mx = Math.max(mx, p.x); }
    return { pts, cx: (mn + mx) / 2, h: pts.at(-1).y };
  }
  function ringAxis(len) {
    const R = len / (2 * Math.PI), n = Math.max(12, Math.ceil(len)), pts = [];
    for (let j = 0; j <= n; j++) { const a = j / n * 2 * Math.PI; pts.push({ x: R * Math.sin(a), y: R - R * Math.cos(a), th: Math.PI / 2 - a, t: j / n, s: j / n * len, ring: true }); }
    return { pts, cx: 0, h: 2 * R };
  }
  const SIGC = { O: [1, 0.54, 0.12], G: [0.24, 0.88, 0.35], R: [1, 0.23, 0.23], A: [0.22, 0.84, 1], Y: [1, 0.82, 0.25], B: [0.36, 0.49, 1] };

  function geometry(scene, o, W, H, S) {
    const N = W * H, E = new Float32Array(N), ST = new Float32Array(N), AL = new Float32Array(N), AC = new Float32Array(N);
    const fl = o.mode === "FL", pmode = o.mode === "R" ? "R" : "G";
    const sig = fl ? [new Float32Array(N), new Float32Array(N), new Float32Array(N)] : null;
    const sigma = 16000 / o.res, radRef = (scene.rad || 6) * S;
    const labels = [];
    scene.items.forEach((it, ii) => {
      const r = rng(o.seed * 7919 + ii * 104729 + 17);
      if (it.dm) { dms(it, r, S, sig, E, ST, W, H); return; }
      const pcs = vary(it.pieces, r);
      const P = profile(pcs, r, pmode, o.het, sigma * (1 + (r() - 0.5) * 0.15));
      const k = it.k * S, rad = (it.rad ?? scene.rad) * S, len = P.total * k * (1 + (r() - 0.5) * 0.05);
      const ax = it.ring ? ringAxis(len) : axis(r, len, (it.bend ?? 1) * o.bend);
      let ox, oy;
      if (it.anchor === "cen" && !it.ring) { const pc = ax.pts[Math.round(P.tc * (ax.pts.length - 1))]; ox = it.x * S - pc.x; oy = it.y * S - pc.y; }
      else { ox = it.x * S - ax.cx; oy = it.y * S - ax.h / 2; }
      const slant = (r() - 0.5) * 0.3, slantN = noise1(r, 8), wN = noise1(r, 40), splay = 0.3 + 0.9 * r();
      const tw = !it.ring && r() < 0.3 ? { t0: (() => { let t; do t = 0.08 + r() * 0.84; while (Math.abs(t - P.tc) < 0.08); return t; })(), w: 0.015 + 0.03 * r() } : null;
      const sepBase = 0.38 + 1.0 * o.sep, cenOn = P.hasCen && !it.ring, il = 1 / Math.max(1, len);
      const chromatidPos = [];
      for (let c = 0; c < 2; c++) {
        const side = c ? 1 : -1, jit = noise1(r, 30), sN = noise1(r, 24), cw = 1 + (r() - 0.5) * 0.1, pos = [];
        const npt = ax.pts.length, stride = Math.max(1, Math.floor(rad * 0.22 * (len / Math.max(1, npt - 1))));
        for (let ip = 0; ip < npt; ip++) {
          const p = ax.pts[ip], t = p.t, dc = (t - P.tc) / 0.022, cenF = cenOn ? Math.exp(-dc * dc) : 0;
          const endD = it.ring ? 1 : Math.min(t, 1 - t), endF = Math.min(1, endD / 0.012);
          const iu = Math.min(P.n - 1, Math.floor(t * P.total)), st = P.stalk[iu];
          const sideF = tw ? Math.cos(Math.PI * smooth(tw.t0 - tw.w, tw.t0 + tw.w, t)) : 1;
          const sepD = rad * sepBase * (1 - 0.88 * cenF) * (1 - 0.8 * st) * (0.9 + 0.12 * sN(t)) * (1 + splay * (0.25 + o.sep) * Math.exp(-endD / 0.035)) * side * sideF;
          const rr = Math.max(0.6, rad * cw * (1 - 0.42 * cenF) * (1 - 0.68 * st) * (1 + 0.07 * wN(t) + 0.04 * jit(t)) * (0.8 + 0.2 * endF));
          const nx = Math.cos(p.th), ny = -Math.sin(p.th), tx = -ny, ty = nx;
          const px = ox + p.x + nx * sepD, py = oy + p.y + ny * sepD, sl = (slant + 0.15 * slantN(t)) / Math.max(1, len), vj = 0.95 + 0.1 * jit(t * 1.7);
          pos.push([px, py, rr]);
          if (ip % stride && ip !== npt - 1) continue;   // discs overlap heavily: fill every few axis points, stain is sampled per pixel
          const R = rr * 1.3, x0 = Math.max(0, Math.floor(px - R)), x1 = Math.min(W - 1, Math.ceil(px + R)), y0 = Math.max(0, Math.floor(py - R)), y1 = Math.min(H - 1, Math.ceil(py + R)), ir2 = 1 / (rr * rr);
          for (let y = y0; y <= y1; y++) { const dy = y - py; for (let x = x0; x <= x1; x++) {
            const dx = x - px, q2 = (dx * dx + dy * dy) * ir2; if (q2 >= 1.69) continue;
            const e = 1 - Math.sqrt(q2) / 1.3, kk = y * W + x; if (e <= E[kk]) continue;
            const own = dx * nx + dy * ny, across = own + sepD, along = dx * tx + dy * ty, tt = t + along * il + own * sl;
            E[kk] = e; ST[kk] = P.prof[Math.max(0, Math.min(P.n - 1, Math.floor(tt * P.total)))] * vj; AL[kk] = p.s + along; AC[kk] = across + ii * 9.7;
          } }
        }
        chromatidPos.push(pos);
      }
      if (fl && it.signals) for (const s of it.signals) {
        const u = mapU(pcs, s.u), t = Math.max(0, Math.min(1, u / P.total)), idx = Math.round(t * (ax.pts.length - 1));
        const sr = rad * (0.4 + 0.22 * Math.min(1, (s.kb || 200) / 500)), col = SIGC[s.color] || SIGC.R;
        for (let c = 0; c < 2; c++) { const [px, py] = chromatidPos[c][idx]; spot(sig, W, H, px + (r() - 0.5) * rad * 0.35, py + (r() - 0.5) * rad * 0.35, sr, col, 0.85 + 0.3 * r()); }
      }
      if (it.label) labels.push(it.label);
    });
    // texture: fibres run across the chromatid (noise stretched across the axis), granular chromatin, ragged edges
    const r = rng(o.seed * 31 + 7), fib = noise2(r), gra = noise2(r);
    const fa = Math.max(0.7, radRef * 0.09), fc = Math.max(2, radRef * 0.55), gs = Math.max(1.1, radRef * 0.26), gs2 = gs * 0.45;
    const D = new Float32Array(N);
    for (let kk = 0; kk < N; kk++) {
      const e = E[kk]; if (e <= 0) continue;
      const x = kk % W, y = (kk / W) | 0;
      const fz = fib(AL[kk] / fa, AC[kk] / fc), q = (1 - e) * 1.3;
      const body = smooth(0.16, 0.33, e + 0.16 * (fz - 0.5));
      const cyl = q < 1 ? Math.sqrt(1 - q * q) : 0;
      const g = 0.72 + 0.4 * gra(x / gs, y / gs) + 0.16 * gra(x / gs2 + 91, y / gs2 + 37);
      let d = body * (0.42 + 0.58 * Math.pow(cyl, 0.8));
      if (q > 0.85) d = Math.max(d, Math.max(0, fz - 0.58) * 2.2 * (1.3 - q) / 0.45 * 0.5);
      D[kk] = d * g;
    }
    // phase halo: a faint bright rim just outside the chromatin
    const halo = blurScaled(D, W, H, Math.max(1.5, radRef * 0.9));
    for (let kk = 0; kk < N; kk++) halo[kk] = Math.max(0, halo[kk] - D[kk] * 1.2);
    return { D, ST, halo, sig, labels, radRef };
  }
  function spot(sig, W, H, px, py, sr, col, a) {
    const R = sr * 2.2, s2 = 1 / (2 * (sr * 0.55) ** 2), x0 = Math.max(0, Math.floor(px - R)), x1 = Math.min(W - 1, Math.ceil(px + R)), y0 = Math.max(0, Math.floor(py - R)), y1 = Math.min(H - 1, Math.ceil(py + R));
    for (let y = y0; y <= y1; y++) for (let x = x0; x <= x1; x++) { const w = a * Math.exp(-((x - px) ** 2 + (y - py) ** 2) * s2); if (w < 0.01) continue; const k = y * W + x; sig[0][k] += col[0] * w; sig[1][k] += col[1] * w; sig[2][k] += col[2] * w; }
  }
  // double minutes: tiny paired chromatin dots carrying the amplified signal, scattered around the spread
  function dms(it, r, S, sig, E, ST, W, H) {
    const n = it.dm, cx = it.x * S, cy = it.y * S, R = it.spread * S, sz = Math.max(1.2, (it.rad || 4) * S * 0.45);
    for (let i = 0; i < n; i++) {
      const a = r() * Math.PI * 2, d = R * Math.sqrt(r()), px = cx + Math.cos(a) * d, py = cy + Math.sin(a) * d, b = r() * Math.PI;
      for (const side of [-1, 1]) {
        const qx = px + Math.cos(b) * sz * 0.9 * side, qy = py + Math.sin(b) * sz * 0.9 * side;
        const x0 = Math.max(0, Math.floor(qx - sz * 1.3)), x1 = Math.min(W - 1, Math.ceil(qx + sz * 1.3)), y0 = Math.max(0, Math.floor(qy - sz * 1.3)), y1 = Math.min(H - 1, Math.ceil(qy + sz * 1.3));
        for (let y = y0; y <= y1; y++) for (let x = x0; x <= x1; x++) { const e = 1 - Math.hypot(x - qx, y - qy) / (sz * 1.3); const k = y * W + x; if (e > E[k]) { E[k] = e; ST[k] = 0.5; } }
        if (sig) for (const c of it.colors) spot(sig, W, H, qx, qy, sz * 0.9, SIGC[c] || SIGC.R, 0.8);
      }
    }
  }

  function field(W, H, seed, S) {
    const r = rng(seed * 131 + 3), n = noise2(r), N = W * H;
    const cx = W / 2, cy = H / 2, vmax = Math.hypot(cx, cy), tilt = (r() - 0.5) * 0.08, tilt2 = (r() - 0.5) * 0.08;
    const vig = lowField(W, H, 8, (x, y) => (1 - 0.17 * Math.pow(Math.hypot(x - cx, y - cy) / vmax, 2.2)) * (1 + tilt * (x / W - 0.5) + tilt2 * (y / H - 0.5)));
    const sc = Math.max(W, H);
    const sv = lowField(W, H, 8, (x, y) => 0.82 + 0.36 * n(x / (sc * 0.22) + 11, y / (sc * 0.22) + 5));
    const cyl = lowField(W, H, 4, (x, y) => smooth(0.52, 0.78, n(x / (sc * 0.16) + 47, y / (sc * 0.16) + 83)) * (0.55 + 0.45 * n(x / (14 * S) + 3, y / (14 * S) + 9)));
    return { vig, sv, cyl, gx: (r() * 512) | 0, gy: (r() * 512) | 0 };
  }

  // ---------------------------------------------------------------- absorbance (brightfield) or emission (fluorescence)
  function absorb(G, F, o, W, H, S) {
    const N = W * H, A = new Float32Array(N), fl = o.mode === "FL", floor = o.mode === "R" ? 0.16 : 0.12, st = o.stain, fld = o.clean ? 0 : o.field;
    if (fl) for (let k = 0; k < N; k++) A[k] = G.D[k] * (0.5 + 0.5 * G.ST[k]) * st * (0.85 + 0.3 * (F.sv[k] - 0.82) / 0.36);
    else for (let k = 0; k < N; k++) A[k] = (G.D[k] * (floor + (1 - floor) * G.ST[k]) - 0.12 * G.halo[k]) * 1.25 * st * F.sv[k] + F.cyl[k] * fld * 0.11;
    const r = rng(o.seed * 977 + 1), tex = noise2(r);
    // interphase nuclei at the edge of the field
    if (fld > 0.35) {
      const nN = 1 + (r() < (fld - 0.35) * 1.5 ? 1 : 0);
      for (let i = 0; i < nN; i++) {
        const R = G.radRef * (7 + 4 * r()), edge = Math.floor(r() * 4), along = 0.15 + 0.7 * r();
        const cx = edge === 0 ? -R * 0.25 : edge === 1 ? W + R * 0.25 : along * W, cy = edge === 2 ? -R * 0.25 : edge === 3 ? H + R * 0.25 : along * H, ell = 0.85 + 0.25 * r(), rot = r() * Math.PI;
        const ca = Math.cos(rot), sa = Math.sin(rot), ts = Math.max(1.5, G.radRef * 0.35);
        for (let y = Math.max(0, Math.floor(cy - R * 1.2)); y < Math.min(H, cy + R * 1.2); y++) for (let x = Math.max(0, Math.floor(cx - R * 1.2)); x < Math.min(W, cx + R * 1.2); x++) {
          const dx = x - cx, dy = y - cy, u = (dx * ca + dy * sa) / (R * ell), v = (-dx * sa + dy * ca) / R, d = Math.sqrt(u * u + v * v);
          const e = 1 - smooth(0.92, 1.0, d + 0.05 * (tex(x / (ts * 3), y / (ts * 3)) - 0.5)); if (e <= 0) continue;
          A[y * W + x] += e * (fl ? 0.55 : 0.62) * st * (0.75 + 0.45 * tex(x / ts + 40, y / ts + 40)) * (0.85 + 0.15 * Math.sqrt(Math.max(0, 1 - d * d)));
        }
      }
    }
    // debris: small dots, some refractile (bright core)
    const nDeb = Math.round((4 + 22 * fld) * (W * H) / (900 * 900 * S * S) + 3);
    for (let d = 0; d < nDeb; d++) {
      const px = r() * W, py = r() * H, rr = (0.6 + r() * 2.6) * S, a = (0.25 + r() * 0.8) * (r() < 0.25 ? -0.6 : 1);
      for (let y = Math.max(0, Math.floor(py - rr)); y < Math.min(H, py + rr); y++) for (let x = Math.max(0, Math.floor(px - rr)); x < Math.min(W, px + rr); x++) { const u = ((x - px) ** 2 + (y - py) ** 2) / (rr * rr); if (u < 1) A[y * W + x] += (fl ? Math.abs(a) * 0.25 : a) * (1 - u); }
    }
    return A;
  }

  // ---------------------------------------------------------------- colour
  // camera grain: one shared 512x512 tile, offset per cell
  const GRAIN = (() => { const r = rng(4242), g = new Float32Array(512 * 512); for (let i = 0; i < g.length; i++) g[i] = r() + r() - 1; return g; })();
  const LUTN = 2600, LOFF = 0.3, LSC = 400;
  const luts = {};
  function lut(kind) {
    if (luts[kind]) return luts[kind];
    const bg = kind === "gray" ? [226, 226, 226] : [238, 234, 226], K = kind === "gray" ? [1.55, 1.55, 1.55] : [1.05, 1.85, 0.72];
    const L = [0, 1, 2].map(c => { const t = new Float32Array(LUTN); for (let i = 0; i < LUTN; i++) t[i] = bg[c] * Math.exp(-(i / LSC - LOFF) * K[c]); return t; });
    return (luts[kind] = L);
  }
  function colour(ctx, Ab, sigB, F, o, W, H, S) {
    const img = ctx.createImageData(W, H), D = img.data;
    if (o.mode === "FL") {
      for (let y = 0; y < H; y++) { const go = ((y + F.gy) & 511) << 9; for (let x = 0; x < W; x++) {
        const k = y * W + x, a = Ab[k] > 0 ? Ab[k] : 0, gr = GRAIN[go + ((x + F.gx) & 511)] * 3, j = k << 2;
        let R = 6 + a * 60, G = 7 + a * 92, B = 12 + a * 235;
        if (sigB) { R += sigB[0][k] * 255; G += sigB[1][k] * 255; B += sigB[2][k] * 255; }
        D[j] = R + gr; D[j + 1] = G + gr; D[j + 2] = B + gr; D[j + 3] = 255;
      } }
    } else {
      const gray = o.color === "gray", L = lut(gray ? "gray" : "giemsa"), L0 = L[0], L1 = L[1], L2 = L[2], ca = gray ? 0 : Math.max(1, Math.round(S * 0.6)), M = LUTN - 1;
      const idx = a => { const i = ((a + LOFF) * LSC + 0.5) | 0; return i < 0 ? 0 : i > M ? M : i; };
      for (let y = 0; y < H; y++) { const go = ((y + F.gy) & 511) << 9, row = y * W; for (let x = 0; x < W; x++) {
        const k = row + x, j = k << 2, v = F.vig[k], gr = GRAIN[go + ((x + F.gx) & 511)] * 5;
        const aG = Ab[k], aR = ca ? Ab[row + (x + ca < W ? x + ca : W - 1)] : aG, aB = ca ? Ab[row + (x >= ca ? x - ca : 0)] : aG;
        D[j] = L0[idx(aR)] * v + gr; D[j + 1] = L1[idx(aG)] * v + (gray ? gr : gr * 0.7); D[j + 2] = L2[idx(aB)] * v + gr; D[j + 3] = 255;
      } }
    }
    ctx.putImageData(img, 0, 0);
  }

  // ---------------------------------------------------------------- render
  // scene: {key, cssW, cssH, rad, items:[{pieces, x, y, anchor:"cen"|"mid", k, rad?, bend?, ring?, signals?:[{u,color,kb}], dm?, colors?, spread?}], labels:[{t,x,y,bold}]}
  // o: {mode:"R"|"G"|"FL", res, focus, stain, sep, bend, het, color, field, clean, seed, q}
  function render(cv, scene, o, after) {
    const t0 = performance.now();
    const S = Math.min(1.5, window.devicePixelRatio || 1) * (o.q || 1);   // above 1.5x the blurred image gains nothing visible
    const W = Math.max(8, Math.round(scene.cssW * S)), H = Math.max(8, Math.round(scene.cssH * S));
    const CC = cv.__scope || (cv.__scope = {}), C = CC[o.q || 1] || (CC[o.q || 1] = {}), T = {};
    let tt = performance.now(); const lap = n => { const x = performance.now(); T[n] = Math.round(x - tt); tt = x; };
    if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
    cv.style.width = scene.cssW + "px"; cv.style.height = scene.cssH + "px";
    const fk = [W, H, o.seed].join();
    if (C.fk !== fk) { C.F = field(W, H, o.seed, S); C.fk = fk; C.gk = C.ak = C.bk = null; lap("field"); }
    const gk = [fk, scene.key, o.res, o.sep, o.bend, o.het, o.mode].join("|");
    if (C.gk !== gk) { C.G = geometry(scene, o, W, H, S); C.gk = gk; C.ak = null; lap("geometry"); }
    const ak = [gk, o.stain, o.clean ? 0 : o.field].join("|");
    if (C.ak !== ak) { C.A = absorb(C.G, C.F, o, W, H, S); C.ak = ak; C.bk = null; lap("absorb"); }
    const bk = [ak, o.focus].join("|");
    if (C.bk !== bk) {
      const sigma = Math.max(0.3, C.G.radRef * (0.05 + 0.2 * o.focus));
      C.Ab = blurScaled(C.A, W, H, sigma);
      C.sigB = C.G.sig ? C.G.sig.map(ch => blurScaled(ch, W, H, Math.max(0.6, sigma * 0.75))) : null;
      C.bk = bk; lap("blur");
    }
    const ctx = cv.getContext("2d");
    colour(ctx, C.Ab, C.sigB, C.F, o, W, H, S); lap("colour");
    if (scene.labels && scene.labels.length) {
      ctx.textAlign = "center"; ctx.textBaseline = "alphabetic"; ctx.fillStyle = o.mode === "FL" ? "#d6dae0" : "#24282c";
      for (const l of scene.labels) { ctx.font = `${l.bold ? 600 : 400} ${11 * S}px "Geist Mono", ui-monospace, monospace`; ctx.fillText(l.t, l.x * S, l.y * S); }
    }
    if (after) after(ctx, S);
    return { ms: Math.round(performance.now() - t0), S, T };
  }

  // Slider helper: half resolution while dragging, full resolution on release.
  function bindPreview(input, onInput, draw) {
    let raf = 0, pendingQ = 1;
    const go = q => { pendingQ = q; cancelAnimationFrame(raf); raf = requestAnimationFrame(() => draw(pendingQ)); };
    input.addEventListener("input", () => { onInput(); go(0.5); });
    input.addEventListener("change", () => { onInput(); go(1); });
    return go;
  }

  return { rng, pieces, whole, unitAt, totalOf, render, bindPreview, Lbp, cenBp, DARK };
})();
