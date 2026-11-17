import { parseAlignment } from "./alignment-states.js";
import { LETTERS, NumberingError, kabatChothiaCdr3, numberRegions, range, relabel, toResult } from "./numbering-core.js";

const kabatChar = LETTERS;

// Kabat region definitions (Kabat numbers) and where each region takes insertions.
const heavyRegions = [
	{ start: 1, end: 30, insertionPos: 6 }, // FR1: insertions at 6
	{ start: 31, end: 35, insertionPos: 35 }, // CDRH1: insertions at 35 (35A, 35B ...)
	{ start: 36, end: 49, insertionPos: null }, // FR2
	{ start: 50, end: 65, insertionPos: 52 }, // CDRH2: insertions at 52
	{ start: 66, end: 94, insertionPos: 82 }, // FR3: 82A-82C, placed by the alignment
	{ start: 95, end: 102, insertionPos: 100 }, // CDRH3: insertions at 100
	{ start: 103, end: 113, insertionPos: null }, // FR4
];
const lightRegions = [
	{ start: 1, end: 23, insertionPos: null }, // FR1
	{ start: 24, end: 34, insertionPos: 27 }, // CDRL1: insertions at 27 (27A ...)
	{ start: 35, end: 49, insertionPos: null }, // FR2
	{ start: 50, end: 56, insertionPos: 52 }, // CDRL2: insertions at 52
	{ start: 57, end: 88, insertionPos: null }, // FR3: placed by the alignment
	{ start: 89, end: 97, insertionPos: 95 }, // CDRL3: insertions at 95
	{ start: 98, end: 107, insertionPos: null }, // FR4
];

// Mapping of the 128 IMGT HMM states onto Kabat. 'X' = the state has its own Kabat position,
// 'I' = the state is an insertion in Kabat. The region strings split the states into the
// stretches that are renumbered separately; rels is the Kabat-IMGT offset at each start.
export const HEAVY_STATE_STRING = "XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
const HEAVY_REGION_STRING = "11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";
export const HEAVY_RELS = { 0: 0, 1: -1, 2: -1, 3: -5, 4: -5, 5: -8, 6: -12, 7: -15 };
export const LIGHT_STATE_STRING = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";
export const LIGHT_REGION_STRING = "11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";
export const LIGHT_RELS = { 0: 0, 1: 0, 2: -6, 3: -6, 4: -13, 5: -16, 6: -20 };

export const regionIndex = (n) => Object.fromEntries(Array.from({ length: n }, (_, i) => [String(i + 1), i]));

/** Heavy FR1 (Kabat and Chothia): put any insertions on `insertionPos` (6). */
export function fw1Insertions(region, insertionPos) {
	const insertions = region.filter(([[, ins]]) => ins !== " ").length;
	if (!insertions) return region;
	const start = region[0][0][0];
	return relabel(region, [
		...range(start, insertionPos + 1),
		...kabatChar.slice(0, insertions).map((c) => [insertionPos, c]),
		...range(insertionPos + 1, insertionPos + 4),
	]);
}

/** Heavy CDR2 (Kabat and Chothia): 50-57, insertions on `insertionPos` (52), deletions from 52 backwards. */
export function heavyCdr2(region, insertionPos) {
	const length = region.length;
	const insertions = Math.max(length - 8, 0);
	return relabel(region, [
		...range(50, insertionPos + 1).slice(0, Math.max(0, length - 5)),
		...kabatChar.slice(0, insertions).map((c) => [insertionPos, c]),
		...range(insertionPos + 1, 58).slice(Math.abs(Math.min(0, length - 5))),
	]);
}

/** Light CDR2 (Kabat and Chothia): insertions on `insertionPos` (52); shorter loops follow the alignment. */
export function lightCdr2(region, insertionPos) {
	const insertions = Math.max(region.length - 4, 0);
	if (!insertions) return region;
	return relabel(region, [
		...range(51, insertionPos + 1),
		...kabatChar.slice(0, insertions).map((c) => [insertionPos, c]),
		...range(insertionPos + 1, 55),
	]);
}

/** CDR3 (Kabat and Chothia): insertions on 100 (heavy) or 95 (light), deletions in a fixed order. */
export function cdr3(region, chain, maxLength) {
	if (region.length > maxLength) throw new NumberingError("CDR3 too long to number");
	return relabel(region, kabatChothiaCdr3(region.length, chain));
}

/**
 * @param {string} annotations - The `#=GC RF` line: 'x' denotes a countable position (IMGT HMM state), '.' an insertion column.
 * @param {string} sequence - The aligned sequence row: uppercase = match, lowercase = insertion, '-' or '.' = gap.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their
 *   corresponding Kabat numbering, one entry per residue (e.g. 52, "52A", 53).
 */
function applyKabatNumberingHeavy(annotations, sequence) {
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions } = numberRegions(residues, states, HEAVY_STATE_STRING, HEAVY_REGION_STRING, regionIndex(8), HEAVY_RELS, 8, [2, 4, 6]);

	// CDRH1: positions up to 35, insertions on 35, deleting from 35 backwards.
	const cdr1 = regions[2];
	const cdr1Pos = heavyRegions[1].insertionPos;
	const cdr1Annotations = [
		...range(23, cdr1Pos + 1).slice(0, cdr1.length),
		...kabatChar.slice(0, Math.max(0, cdr1.length - 13)).map((c) => [cdr1Pos, c]),
	];
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
 *   corresponding Kabat numbering (kappa or lambda light chain), one entry per residue.
 */
function applyKabatNumberingLight(annotations, sequence) {
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions } = numberRegions(residues, states, LIGHT_STATE_STRING, LIGHT_REGION_STRING, regionIndex(7), LIGHT_RELS, 7, [1, 3, 5]);

	// CDRL1: insertions on 27, deleting forwards from 28.
	const cdr1 = regions[1];
	const length = cdr1.length;
	const cdr1Pos = lightRegions[1].insertionPos;
	const cdr1Annotations = [
		...range(24, cdr1Pos + 1).slice(0, Math.max(0, length)),
		...kabatChar.slice(0, Math.max(length - 11, 0)).map((c) => [cdr1Pos, c]),
		...range(cdr1Pos + 1, 35).slice(Math.abs(Math.min(0, length - 11))),
	];
	return toResult([
		...regions[0],
		...relabel(cdr1, cdr1Annotations),
		...regions[2],
		...lightCdr2(regions[3], lightRegions[3].insertionPos),
		...regions[4],
		...cdr3(regions[5], "L", 35),
		...regions[6],
	]);
}

export { applyKabatNumberingLight };

export { applyKabatNumberingHeavy };
