// FISH signal engine: rearrangements -> chromosomes -> visible signals -> signal code and ISCN.
// Used by fish.html and by tools/check_probes.mjs (node), so it has no DOM code.
const FISH = (() => {
  const ORDER = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"];
  const CI = c => ORDER.indexOf(c);
  const MIN_KB = 40;          // a probe fragment smaller than this is too dim to score
  const MIN_FRAC = 0.12;      // ... or smaller than this fraction of the whole probe
  const JOIN_SAME = 3.0e6;    // pieces of one probe closer than this read as one signal
  const JOIN_DIFF = 1.8e6;    // different probes closer than this read as one fused signal

  // ---------------------------------------------------------------- cell model
  // A chromosome is a list of segments [chr, start, end, dir] read pter -> qter.
  function normalCell(lengths, cen, sex) {
    const cell = [];
    for (const c of ORDER) {
      const n = c === "X" ? (sex === "XY" ? 1 : 2) : c === "Y" ? (sex === "XY" ? 1 : 0) : 2;
      for (let k = 0; k < n; k++) cell.push({ chr: c, label: c, normal: true, segs: [[c, 1, lengths[c], 1]] });
    }
    return { copies: cell, lengths, cen };
  }
  const segLen = s => s[2] - s[1];
  const copyLen = cp => cp.segs.reduce((a, s) => a + segLen(s), 0);
  // position of genomic point (c,p) on a copy, or -1
  function locate(cp, c, p) {
    let off = 0;
    for (let i = 0; i < cp.segs.length; i++) {
      const s = cp.segs[i];
      if (s[0] === c && p >= s[1] && p <= s[2]) return { i, at: off + (s[3] > 0 ? p - s[1] : s[2] - p) };
      off += segLen(s);
    }
    return null;
  }
  // split a copy at derivative coordinate x into [left, right] segment lists
  function cut(cp, x) {
    const L = [], R = []; let off = 0;
    for (const s of cp.segs) {
      const len = segLen(s);
      if (off + len <= x) L.push(s);
      else if (off >= x) R.push(s);
      else {
        const k = x - off;
        if (s[3] > 0) { L.push([s[0], s[1], s[1] + k, 1]); R.push([s[0], s[1] + k, s[2], 1]); }
        else { L.push([s[0], s[2] - k, s[2], -1]); R.push([s[0], s[1], s[2] - k, -1]); }
      }
      off += len;
    }
    return [L, R];
  }
  const flip = segs => segs.slice().reverse().map(s => [s[0], s[1], s[2], -s[3]]);
  // find a copy carrying (c,p): prefer an already abnormal one, then a normal homologue not in `avoid`
  function pick(cell, c, p, avoid = [], preferAbnormal = false) {
    const cands = cell.copies.filter(cp => !avoid.includes(cp) && locate(cp, c, p));
    if (preferAbnormal) { const ab = cands.find(cp => !cp.normal); if (ab) return ab; }
    return cands.find(cp => cp.normal && cp.chr === c) || cands[0] || null;
  }
  const derLabel = c => `der(${c})`;

  function applyEvent(cell, e) {
    const [k, v] = Object.entries(e)[0];
    const cs = cell.copies;
    if (k === "t") {
      const [[ca, pa], [cb, pb]] = v;
      const A = pick(cell, ca, pa), B = pick(cell, cb, pb, [A]);
      if (!A || !B) return;
      const la = locate(A, ca, pa), lb = locate(B, cb, pb);
      const armA = pa >= cell.cen[ca] ? "q" : "p", armB = pb >= cell.cen[cb] ? "q" : "p";
      const [A1, A2] = cut(A, la.at), [B1, B2] = cut(B, lb.at);
      // the centromere-bearing part of each chromosome keeps it; the distal parts swap
      const keepA = armA === "q" ? A1 : A2, distA = armA === "q" ? A2 : A1;
      const keepB = armB === "q" ? B1 : B2, distB = armB === "q" ? B2 : B1;
      const joinA = armA === "q" ? keepA.concat(armB === "q" ? distB : flip(distB)) : (armB === "p" ? distB : flip(distB)).concat(keepA);
      const joinB = armB === "q" ? keepB.concat(armA === "q" ? distA : flip(distA)) : (armA === "p" ? distA : flip(distA)).concat(keepB);
      A.segs = joinA; A.normal = false; A.label = derLabel(A.chr);
      B.segs = joinB; B.normal = false; B.label = derLabel(B.chr);
      return;
    }
    if (k === "der") {          // unbalanced: keep only the der of chromosome v[2]
      const before = cs.slice();
      applyEvent(cell, { t: [v[0], v[1]] });
      const ders = cs.filter(cp => !cp.normal && !before.some(b => b === cp && b.normal));
      for (const d of cs.filter(cp => !cp.normal)) if (d.chr !== v[2]) { d.segs = [[d.chr, 1, cell.lengths[d.chr], 1]]; d.normal = true; d.label = d.chr; }
      return ders;
    }
    if (k === "del" || k === "del2") {
      const [c, s, e2] = v;
      // endpoints nudged inwards so a breakpoint shared with an earlier rearrangement resolves to the right side
      const s1 = s + 1, e1 = e2 - 1;
      const both = cs.filter(cp => locate(cp, c, s1) && locate(cp, c, e1));
      const targets = k === "del2" ? both.slice(0, 2) : [both.find(cp => !cp.normal) || both[0]];
      for (const cp of targets) {
        if (!cp) continue;
        const a = locate(cp, c, s1), b = locate(cp, c, e1);
        if (!a || !b) continue;
        const x0 = Math.min(a.at, b.at), x1 = Math.max(a.at, b.at);
        const [L] = cut(cp, x0), [, R] = cut(cp, x1);
        cp.segs = L.concat(R); cp.normal = false; cp.label = `del(${cp.chr})`;
      }
      return;
    }
    if (k === "delter") {
      const [c, p] = v; const cp = pick(cell, c, p, [], true); if (!cp) return;
      const l = locate(cp, c, p); const [L, R] = cut(cp, l.at);
      cp.segs = p < cell.cen[c] ? R : L; cp.normal = false; cp.label = `del(${cp.chr})`;
      return;
    }
    if (k === "inv") {
      const [c, s, e2] = v; const cp = pick(cell, c, s); if (!cp) return;
      const a = locate(cp, c, s), b = locate(cp, c, e2);
      const [L, M0] = cut(cp, a.at); const [M, R] = cut({ segs: M0 }, b.at - a.at);
      cp.segs = L.concat(flip(M), R); cp.normal = false; cp.label = `inv(${cp.chr})`;
      return;
    }
    if (k === "gain") { const [c, n] = v; for (let i = 0; i < n; i++) cs.push({ chr: c, label: c, normal: true, gained: true, segs: [[c, 1, cell.lengths[c], 1]] }); return; }
    if (k === "loss") { const i = cs.findIndex(cp => cp.chr === v[0] && cp.normal); if (i >= 0) cs.splice(i, 1); return; }
    if (k === "armloss") {
      const [c, arm] = v; const cp = cs.find(x => x.chr === c && x.normal); if (!cp) return;
      const m = cell.cen[c];
      cp.segs = arm === "p" ? [[c, m, cell.lengths[c], 1]] : [[c, 1, m, 1]]; cp.normal = false; cp.label = `del(${c})`;
      return;
    }
    if (k === "dup_der") { const d = cs.find(cp => !cp.normal && cp.chr === v[0]); if (d) cs.push({ ...d, segs: d.segs.slice(), label: `+der(${d.chr})` }); return; }
    if (k === "gainreg") {
      const [[c, s, e2], n] = v; const pad = 1.5e6;
      for (let i = 0; i < n; i++) cs.push({ chr: c, label: `+${c}`, normal: false, extra: true, segs: [[c, Math.max(1, s - pad), Math.min(cell.lengths[c], e2 + pad), 1]] });
      return;
    }
    if (k === "amp") { const [[c, s, e2], n] = v; cs.push({ chr: c, label: "amp", normal: false, amp: n, segs: [[c, s, e2, 1]] }); return; }
    if (k === "ins") {
      const [tc, tp, sc, ss, se] = v; const cp = pick(cell, tc, tp); if (!cp) return;
      const l = locate(cp, tc, tp); const [L, R] = cut(cp, l.at);
      cp.segs = L.concat([[sc, ss, se, 1]], R); cp.normal = false; cp.label = `der(${cp.chr})`;
      return;
    }
  }

  function buildCell(ref, events, sex = "XX") {
    const cell = normalCell(ref.lengths, ref.cen, sex);
    for (const e of events) applyEvent(cell, e);
    return cell;
  }

  // ---------------------------------------------------------------- signals
  const compKb = comp => comp.segs.filter(s => s[0] !== "Y" || comp.segs.every(t => t[0] === "Y")).reduce((a, s) => a + (s[2] - s[1]), 0) / 1000;

  // pieces of each component on one copy, in derivative coordinates
  function pieces(cp, comp) {
    const out = []; let off = 0;
    for (const s of cp.segs) {
      for (const t of comp.segs) {
        if (t[0] !== s[0]) continue;
        const lo = Math.max(s[1], t[1]), hi = Math.min(s[2], t[2]);
        if (hi <= lo) continue;
        const a = s[3] > 0 ? off + (lo - s[1]) : off + (s[2] - hi), b = a + (hi - lo);
        out.push([a, b]);
      }
      off += segLen(s);
    }
    out.sort((x, y) => x[0] - y[0]);
    const blobs = [];
    for (const p of out) {
      const last = blobs[blobs.length - 1];
      if (last && p[0] - last.b <= JOIN_SAME) { last.b = Math.max(last.b, p[1]); last.bp += p[1] - p[0]; }
      else blobs.push({ a: p[0], b: p[1], bp: p[1] - p[0] });
    }
    return blobs;
  }

  // Visible signals for a probe (optionally one probe set) in a cell.
  function signals(probe, cell, set = null) {
    const comps = probe.comps.map((c, i) => ({ ...c, i })).filter(c => set === null || c.set === set);
    const spots = [];
    cell.copies.forEach((cp, ci) => {
      let blobs = [];
      for (const c of comps) {
        const tot = compKb(c);
        for (const b of pieces(cp, c)) {
          const kb = b.bp / 1000;
          if (kb < Math.min(MIN_KB, tot * 0.5) || kb < tot * MIN_FRAC) continue;
          blobs.push({ ...b, comp: c.i, color: c.color, kb });
        }
      }
      if (cp.amp) { // amplified region: one cluster per component present
        for (const b of blobs) spots.push({ copy: ci, comps: [b.comp], colors: [b.color], kb: b.kb, at: b.a, amp: cp.amp, len: copyLen(cp) });
        return;
      }
      blobs.sort((x, y) => x.a - y.a);
      const groups = [];
      for (const b of blobs) {
        const g = groups.find(g => g.some(x => Math.max(0, Math.max(x.a, b.a) - Math.min(x.b, b.b)) <= JOIN_DIFF));
        if (g) g.push(b); else groups.push([b]);
      }
      // merge groups that became connected
      for (let changed = true; changed;) {
        changed = false;
        outer: for (let i = 0; i < groups.length; i++) for (let j = i + 1; j < groups.length; j++) {
          if (groups[i].some(x => groups[j].some(y => Math.max(0, Math.max(x.a, y.a) - Math.min(x.b, y.b)) <= JOIN_DIFF))) {
            groups[i] = groups[i].concat(groups[j]); groups.splice(j, 1); changed = true; break outer;
          }
        }
      }
      for (const g of groups) {
        const compsIn = [...new Set(g.map(x => x.comp))].sort((a, b) => a - b);
        spots.push({ copy: ci, comps: compsIn, colors: compsIn.map(i => probe.comps[i].color), kb: g.reduce((a, x) => a + x.kb, 0),
          parts: g.map(x => ({ comp: x.comp, color: x.color, kb: x.kb, at: (x.a + x.b) / 2 })), at: g.reduce((a, x) => a + (x.a + x.b) / 2 * x.kb, 0) / g.reduce((a, x) => a + x.kb, 0), len: copyLen(cp) });
      }
    });
    return spots;
  }

  const sets = probe => [...new Set(probe.comps.map(c => c.set))].sort();

  // Signal code such as 1O1G2F, written per probe set
  function code(probe, cell) {
    return sets(probe).map(s => {
      const comps = probe.comps.map((c, i) => ({ ...c, i })).filter(c => c.set === s);
      const colors = [...new Set(comps.map(c => c.color))];
      const isCtrl = c => /^D\d+Z\d|^D[XY]Z\d/.test(c.locus || "");
      const pair = [...new Set(comps.filter(c => !isCtrl(c)).map(c => c.color))];
      const sp = signals(probe, cell, s);
      const counts = new Map();
      for (const x of sp) {
        let lab;
        if (x.amp) lab = "amp " + x.colors.join("");
        else if (pair.length === 2 && x.comps.length >= 2 && new Set(x.colors).size === 2 && x.colors.every(c => pair.includes(c))) lab = "F";
        else lab = [...new Set(comps.filter(c => x.comps.includes(c.i)).map(c => c.color))].join("");
        counts.set(lab, (counts.get(lab) || 0) + 1);
      }
      const keys = [...counts.keys()].sort((a, b) => {
        const ra = a.startsWith("amp") ? 3 : a === "F" ? 2 : a.length > 1 ? 1.5 : 0, rb = b.startsWith("amp") ? 3 : b === "F" ? 2 : b.length > 1 ? 1.5 : 0;
        if (ra !== rb) return ra - rb;
        const oa = colors.indexOf(a[0]), ob = colors.indexOf(b[0]);
        return a.length - b.length || oa - ob;
      });
      return keys.map(k => k.startsWith("amp") ? k : counts.get(k) + k).join("").replace(/([A-Z])amp/, "$1 amp");
    }).join(" | ");
  }

  // parse a code such as "1O1G2F" or "2AO2G | 2O" into a multiset for comparison
  function parseCode(str) {
    return (str || "").split("|").map(part => {
      const m = new Map(); const s = part.trim();
      const amp = s.match(/amp\s*([A-Z]+)/); let rest = s.replace(/amp\s*[A-Z]+/, "");
      if (amp) m.set("amp " + amp[1], 1);
      for (const [, n, lab] of rest.matchAll(/(\d+)([A-Z]+)/g)) m.set([...lab].sort().join(""), (m.get([...lab].sort().join("")) || 0) + +n);
      return m;
    });
  }
  function sameCode(a, b) {
    const A = parseCode(a), B = parseCode(b);
    if (A.length !== B.length) return false;
    return A.every((m, i) => m.size === B[i].size && [...m].every(([k, v]) => B[i].get(k) === v));
  }

  // ---------------------------------------------------------------- ISCN (nuc ish)
  function iscn(probe, cell, normalCell, nCells = 200) {
    const all = signals(probe, cell), norm = signals(probe, normalCell);
    const comps = probe.comps;
    const loc = i => comps[i].locus;
    const pos = i => { const s = comps[i].segs[0]; return [CI(s[0]), (s[1] + s[2]) / 2]; };
    const lociIdx = [...new Map(comps.map((c, i) => [c.locus, i])).values()].sort((a, b) => { const A = pos(a), B = pos(b); return A[0] - B[0] || A[1] - B[1]; });
    const idxOfLocus = new Map(comps.map((c, i) => [c.locus, i]));
    const count = l => all.filter(s => s.comps.some(i => loc(i) === l)).length;
    const ampOf = l => all.some(s => s.amp && s.comps.some(i => loc(i) === l));
    // loci pairs that sit together on a normal chromosome
    const adj = new Set();
    for (const s of norm) for (const a of s.comps) for (const b of s.comps) if (loc(a) !== loc(b)) adj.add(loc(a) + "|" + loc(b));
    const geneOf = l => l.replace(/^[53]'/, "");
    const isBAPair = (a, b) => /^[53]'/.test(a) && /^[53]'/.test(b) && geneOf(a) === geneOf(b);
    // counts list, collapsing 5'/3' halves of one gene when they always travel together
    const loci = lociIdx.map(i => loc(i));
    const parts = []; const used = new Set(); const extra = [];
    for (const l of loci) {
      if (used.has(l)) continue;
      const partner = loci.find(m => m !== l && isBAPair(l, m));
      if (partner) {
        used.add(partner);
        const both = all.filter(s => s.comps.some(i => loc(i) === l) && s.comps.some(i => loc(i) === partner)).length;
        const onlyL = count(l) - both, onlyP = count(partner) - both;
        const five = l.startsWith("5'") ? l : partner, three = l.startsWith("5'") ? partner : l;
        if (onlyL === 0 && onlyP === 0) parts.push([geneOf(l), both]);
        else if (onlyL === onlyP) { parts.push([geneOf(l), both + onlyL]); extra.push(`${five} sep ${three}x${onlyL}`); }
        else { parts.push([five, count(five)]); parts.push([three, count(three)]); if (both) extra.push(`${five} con ${three}x${both}`); }
        continue;
      }
      parts.push([l, ampOf(l) ? "amp" : count(l)]);
    }
    // separations of loci that are adjacent on normal chromosomes (tri-colour break-apart designs)
    for (const key of adj) {
      const [a, b] = key.split("|"); if (a > b || isBAPair(a, b)) continue;
      const onlyA = all.filter(s => s.comps.some(i => loc(i) === a) && !s.comps.some(i => loc(i) === b)).length;
      const onlyB = all.filter(s => s.comps.some(i => loc(i) === b) && !s.comps.some(i => loc(i) === a)).length;
      const k = Math.min(onlyA, onlyB);
      if (k > 0) extra.push(`${a} sep ${b}x${k}`);
    }
    // fusions of loci that are apart on normal chromosomes
    const fus = new Map();
    for (const s of all) {
      if (s.amp) continue;
      const ls = [...new Set(s.comps.map(loc))];
      const inv = new Set();
      for (const a of ls) for (const b of ls) if (a !== b && !adj.has(a + "|" + b) && !isBAPair(a, b)) { inv.add(a); inv.add(b); }
      if (inv.size < 2) continue;
      const key = loci.filter(l => inv.has(l)).join(" con ");
      fus.set(key, (fus.get(key) || 0) + 1);
    }
    for (const [k, n] of fus) extra.push(`${k}x${n}`);
    const nums = parts.map(p => p[1]);
    let head;
    if (nums.every(n => n === nums[0]) && nums[0] !== "amp" && parts.length > 1) head = `(${parts.map(p => p[0]).join(",")})x${nums[0]}`;
    else head = "(" + parts.map(([l, n]) => n === "amp" ? `${l}amp` : `${l}x${n}`).join(",") + ")";
    return `nuc ish${head}${extra.length ? "(" + extra.join(",") + ")" : ""}[${nCells}]`;
  }

  return { ORDER, normalCell, buildCell, applyEvent, signals, code, parseCode, sameCode, iscn, sets, copyLen, compKb };
})();
if (typeof module !== "undefined") module.exports = FISH;
