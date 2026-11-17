import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";

import { applyAhoNumbering } from "../aho-scheme.js";
import { parseAlignment } from "../alignment-states.js";
import { applyChothiaNumberingHeavy, applyChothiaNumberingLight } from "../chothia-scheme.js";
import { applyImgtNumbering } from "../imgt-scheme.js";
import { applyKabatNumberingHeavy, applyKabatNumberingLight } from "../kabat-scheme.js";
import { NumberingError, symmetricCdr } from "../numbering-core.js";

// Reference cases: 507 PDB antibody/TCR domains, one per distinct loop-length and insertion
// pattern, with expected numbering from the Python package (cross-checked against ANARCI).
// Set REFERENCE=path/to/other.json to run a larger set made with test/make_reference.py.
const referencePath = process.env.REFERENCE ?? new URL("./fixtures/reference.json", import.meta.url);
const reference = JSON.parse(readFileSync(referencePath, "utf8"));

function run(scheme, rec) {
	const heavy = rec.chain === "H";
	switch (scheme) {
		case "imgt":
			return applyImgtNumbering(rec.annotations, rec.sequence);
		case "kabat":
			return (heavy ? applyKabatNumberingHeavy : applyKabatNumberingLight)(rec.annotations, rec.sequence);
		case "chothia":
			return (heavy ? applyChothiaNumberingHeavy : applyChothiaNumberingLight)(rec.annotations, rec.sequence);
		case "aho":
			return applyAhoNumbering(rec.annotations, rec.sequence, rec.chain);
	}
}

for (const scheme of ["imgt", "kabat", "chothia", "aho"]) {
	test(`${scheme}: matches the reference numbering`, () => {
		let checked = 0;
		for (const rec of reference) {
			const expected = rec.expected[scheme];
			if (!expected) continue;
			if (expected.error) {
				assert.throws(() => run(scheme, rec), NumberingError, `${rec.name} should be refused`);
			} else {
				assert.deepEqual(run(scheme, rec), expected, `${rec.name} (${rec.chain})`);
			}
			checked++;
		}
		assert.ok(checked > 100, `only ${checked} cases checked`);
	});
}

// Trastuzumab VH aligned to the IMGT HMM (gaps at IMGT 10, 31-34, 60-61, 73).
const TRASTUZUMAB_RF = "x".repeat(128);
const TRASTUZUMAB_VH = "EVQLVESGG-GLVQPGGSLRLSCAASGFNI----KDTYIHWVRQAPGKGLEWVARIYPT--NGYTRYADSVK-GRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";

function positionsOf(result) {
	return Object.fromEntries(result.numbering.map((p, i) => [String(p), result.concatSequence[i]]));
}

test("trastuzumab VH: published Kabat numbering", () => {
	const kabat = positionsOf(applyKabatNumberingHeavy(TRASTUZUMAB_RF, TRASTUZUMAB_VH));
	assert.equal(["31", "32", "33", "34", "35"].map((p) => kabat[p]).join(""), "DTYIH");
	assert.equal(kabat["52A"], "P");
	assert.deepEqual(["82A", "82B", "82C"].map((p) => kabat[p]), ["N", "S", "L"]);
	assert.deepEqual(["100", "100A", "100B", "100C", "101", "102"].map((p) => kabat[p]), ["F", "Y", "A", "M", "D", "Y"]);
	assert.equal(kabat["113"], "S");
});

test("trastuzumab VH: IMGT gaps and conserved residues", () => {
	const result = applyImgtNumbering(TRASTUZUMAB_RF, TRASTUZUMAB_VH);
	const imgt = positionsOf(result);
	assert.equal(result.concatSequence.length, result.numbering.length);
	assert.equal(imgt["23"], "C");
	assert.equal(imgt["104"], "C");
	assert.equal(imgt["118"], "W");
	for (const gap of ["10", "31", "32", "33", "34", "60", "61", "73"]) assert.equal(imgt[gap], undefined);
});

test("one number per residue, deletions get no number", () => {
	const r = applyKabatNumberingHeavy(TRASTUZUMAB_RF, TRASTUZUMAB_VH);
	assert.equal(r.concatSequence, TRASTUZUMAB_VH.replaceAll("-", ""));
	assert.equal(r.numbering.length, r.concatSequence.length);
});

test("parseAlignment: insert columns, MSA gaps and flanks", () => {
	// 2 flank residues, states 1-2, an insertion plus an MSA gap column, states 3-4, 2 flank residues.
	const { residues, states } = parseAlignment("..xx..xx..", "mkABc.CDef");
	assert.equal(residues, "MKABCCDEF");
	assert.deepEqual(states, [
		[[1, "m"], 2],
		[[2, "m"], 3],
		[[3, "i"], 4], // insertions carry the following state number
		[[3, "m"], 5],
		[[4, "m"], 6],
	]);
});

test("IMGT CDR3 insertions are symmetric on 111 and 112", () => {
	const labels = (n) => symmetricCdr(n, 13, 105, 118).filter(Boolean).map(([p, ins]) => `${p}${ins.trim()}`);
	assert.deepEqual(labels(14).slice(6, 9), ["111", "112A", "112"]);
	assert.deepEqual(labels(16).slice(6, 11), ["111", "111A", "112B", "112A", "112"]);
	assert.deepEqual(labels(5), ["105", "106", "107", "116", "117"]);
});

test("AHo requires a chain type", () => {
	assert.throws(() => applyAhoNumbering(TRASTUZUMAB_RF, TRASTUZUMAB_VH), NumberingError);
});
