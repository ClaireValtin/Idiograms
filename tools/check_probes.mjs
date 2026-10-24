// Check every FISH pattern in data/probes.json: the signal code derived from the probe positions
// and the rearrangement must equal the code the manufacturer gives (`expect`).
//   node tools/check_probes.mjs          # summary, non-zero exit on any mismatch
//   node tools/check_probes.mjs -v       # also print every pattern with its ISCN string
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
const require = createRequire(import.meta.url);
const FISH = require("../src/_fish.js");
const root = new URL("..", import.meta.url);
const probes = JSON.parse(readFileSync(new URL("data/probes.json", root))).probes;
const bands = JSON.parse(readFileSync(new URL("data/bands.json", root))).chromosomes;
const lengths = {}, cen = {};
for (const [c, v] of Object.entries(bands)) {
  lengths[c] = v["850"].at(-1)[4];
  const ac = v["850"].filter(b => b[5] === "acen");
  cen[c] = ac[0][4];
}
const ref = { lengths, cen };
const verbose = process.argv.includes("-v");
let bad = 0, checked = 0;
for (const p of probes) {
  for (const pat of p.patterns) {
    const sex = pat.events.some(e => e.gain && e.gain[0] === "Y") ? "XX" : "XX";
    const cell = FISH.buildCell(ref, pat.events, sex);
    const normal = FISH.buildCell(ref, [], sex);
    const got = FISH.code(p, cell);
    const iscn = FISH.iscn(p, cell, normal);
    if (pat.expect !== null && pat.expect !== undefined) {
      checked++;
      const ok = FISH.sameCode(got, pat.expect);
      if (!ok) { bad++; console.log(`MISMATCH ${p.id} / ${pat.name}: expected ${pat.expect}, got ${got}`); }
      else if (verbose) console.log(`ok  ${p.id} / ${pat.name}: ${got}   ${iscn}`);
    } else if (verbose) console.log(`--  ${p.id} / ${pat.name}: ${got}   ${iscn}`);
  }
}
console.log(`${checked} patterns checked, ${bad} mismatches`);
process.exit(bad ? 1 : 0);
