import { parseAlignment } from "./alignment-states.js";
import { LETTERS, numberRegions, range, relabel, toResult } from "./numbering-core.js";
import {
	HEAVY_RELS,
	HEAVY_STATE_STRING,
	LIGHT_REGION_STRING,
	LIGHT_RELS,
	LIGHT_STATE_STRING,
	cdr3,
	fw1Insertions,
	heavyCdr2,
	lightCdr2,
	regionIndex,
} from "./kabat-scheme.js";

const chothiaChar = LETTERS;

// Chothia region definitions (Chothia numbers) and where each region takes insertions.
// Chothia differs from Kabat in CDR1: heavy insertions on 31, light insertions on 30.
const heavyRegions = [
	{ start: 1, end: 25, insertionPos: 6 }, // FR1: insertions at 6
	{ start: 26, end: 32, insertionPos: 31 }, // CDRH1: insertions at 31 (31A, 31B ...)
	{ start: 33, end: 51, insertionPos: null }, // FR2
	{ start: 52, end: 56, insertionPos: 52 }, // CDRH2: insertions at 52
	{ start: 57, end: 94, insertionPos: 82 }, // FR3: 82A-82C, placed by the alignment
	{ start: 95, end: 102, insertionPos: 100 }, // CDRH3: insertions at 100
	{ start: 103, end: 113, insertionPos: null }, // FR4
];
const lightRegions = [
	{ start: 1, end: 23, insertionPos: null }, // FR1
	{ start: 24, end: 34, insertionPos: 30 }, // CDRL1: insertions at 30 (30A ...)
	{ start: 35, end: 49, insertionPos: null }, // FR2
	{ start: 50, end: 56, insertionPos: 52 }, // CDRL2: insertions at 52
	{ start: 57, end: 88, insertionPos: 68 }, // FR3: insertions (and a single deletion) at 68
	{ start: 89, end: 97, insertionPos: 95 }, // CDRL3: insertions at 95
	{ start: 98, end: 107, insertionPos: null }, // FR4
];

// Chothia's heavy CDR1 stretch is shorter than Kabat's (it ends at IMGT state 38, not 40).
const HEAVY_REGION_STRING = "11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";

/** Light FR3: insertions on 68; a single deletion is also put on 68, otherwise follow the alignment. */
function lightFw3(region, insertionPos) {
	const length = region.length;
	const insertions = Math.max(length - 34, 0);
	if (insertions) {
		return relabel(region, [
			...range(55, insertionPos + 1),
			...chothiaChar.slice(0, insertions).map((c) => [insertionPos, c]),
			...range(insertionPos + 1, 89),
		]);
	}
	if (length === 33) return relabel(region, [...range(55, insertionPos), ...range(insertionPos + 1, 89)]);
	return region;
}

/**
 * @param {string} annotations - The `#=GC RF` line: 'x' denotes a countable position (IMGT HMM state), '.' an insertion column.
 * @param {string} sequence - The aligned sequence row: uppercase = match, lowercase = insertion, '-' or '.' = gap.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their
 *   corresponding Chothia numbering, one entry per residue (e.g. 31, "31A", 32).
 */
function applyChothiaNumberingHeavy(annotations, sequence) {
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions } = numberRegions(residues, states, HEAVY_STATE_STRING, HEAVY_REGION_STRING, regionIndex(8), HEAVY_RELS, 8, [0, 2, 4, 6]);

	// CDRH1: insertions on 31.
	const cdr1 = regions[2];
	const length = cdr1.length;
	const cdr1Pos = heavyRegions[1].insertionPos;
	const insertions = Math.max(length - 11, 0);
	const cdr1Annotations = insertions
		? [...range(23, cdr1Pos + 1), ...chothiaChar.slice(0, insertions).map((c) => [cdr1Pos, c]), [32, " "], [33, " "]]
		: // negative slice ends are intentional: they match the reference rules for very short loops
			[...range(23, cdr1Pos + 1).slice(0, length - 2), ...[[32, " "], [33, " "]].slice(0, length)];
	return toResult([
		...fw1Insertions(regions[0], heavyRegions[0].insertionPos),
		...regions[1],
		...relabel(cdr1, cdr1Annotations),
		...regions[3],
		...heavyCdr2(regions[4], heavyRegions[3].insertionPos),
		...regions[5],
		...cdr3(regions[6], "H", 36),
		...regions[7],
	]);
}

/**
 * @param {string} annotations - The `#=GC RF` line: 'x' denotes a countable position (IMGT HMM state), '.' an insertion column.
 * @param {string} sequence - The aligned sequence row: uppercase = match, lowercase = insertion, '-' or '.' = gap.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their
 *   corresponding Chothia numbering (kappa or lambda light chain), one entry per residue.
 */
function applyChothiaNumberingLight(annotations, sequence) {
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions } = numberRegions(residues, states, LIGHT_STATE_STRING, LIGHT_REGION_STRING, regionIndex(7), LIGHT_RELS, 7, [1, 3, 4, 5]);

	// CDRL1: insertions on 30, deleting forwards from 31.
	const cdr1 = regions[1];
	const length = cdr1.length;
	const cdr1Pos = lightRegions[1].insertionPos;
	const cdr1Annotations = [
		...range(24, cdr1Pos + 1).slice(0, Math.max(0, length)),
		...chothiaChar.slice(0, Math.max(length - 11, 0)).map((c) => [cdr1Pos, c]),
		...range(cdr1Pos + 1, 35).slice(Math.abs(Math.min(0, length - 11))),
	];
	return toResult([
		...regions[0],
		...relabel(cdr1, cdr1Annotations),
		...regions[2],
		...lightCdr2(regions[3], lightRegions[3].insertionPos),
		...lightFw3(regions[4], lightRegions[4].insertionPos),
		...cdr3(regions[5], "L", 35),
		...regions[6],
	]);
}

export { applyChothiaNumberingLight };

export { applyChothiaNumberingHeavy };
