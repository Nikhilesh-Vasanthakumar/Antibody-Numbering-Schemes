import { parseAlignment } from "./alignment-states.js";
import { LETTERS, NumberingError, compareAnnotations, numberRegions, range, relabel, toResult } from "./numbering-core.js";

/**
 * @typedef {'L' | 'K' | 'H' | 'A' | 'B' | 'D' | 'G'} ChainType
 */

const ahoChar = LETTERS;

// AHo region definitions (AHo numbers). AHo reserves a position for every structurally
// observed residue, so loops are mostly *gapped*: positions are removed in a fixed order
// around the gap centre. insertionPos is only used for loops longer than AHo allows.
const regions = [
	{ start: 1, end: 10, insertionPos: 8 }, // Region B: indel at 8
	{ start: 11, end: 24, insertionPos: null }, // Region C
	{ start: 25, end: 42, insertionPos: 36 }, // Region D (CDR1): gaps around 28 and 36, see orderedDeletions
	{ start: 43, end: 57, insertionPos: null }, // Region E
	{ start: 58, end: 77, insertionPos: 63 }, // Region F (CDR2): gaps around 63
	{ start: 78, end: 93, insertionPos: 85 }, // Region H: gaps on 86 then 85, insertions on 85
	{ start: 94, end: 106, insertionPos: null }, // Region I
	{ start: 107, end: 138, insertionPos: 123 }, // Region J (CDR3): gaps around 123
	{ start: 139, end: 149, insertionPos: null }, // Region K
];

// CDR1 (AHo 25-42): order in which positions are deleted for short loops, per chain type.
const orderedDeletions = {
	L: [28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42],
	K: [28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	H: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	A: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	B: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	D: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	G: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
};
const cdr2Deletions = [63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77];
const cdr2DeletionsAlpha = [74, 73, 63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 75, 76, 77];
const fw3Deletions = [86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78];
const cdr3Deletions = [
	123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130, 116, 131,
	115, 132, 114, 133, 113, 134, 112, 135, 111, 136, 110, 137, 109, 138, 108, 107,
];

// Mapping of the 128 IMGT HMM states onto AHo regions B-K (A and G are empty) and the
// AHo-IMGT offset at the start of each region.
const STATE_STRING = "X".repeat(128);
const REGION_STRING = "BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK";
const REGION_INDEX = Object.fromEntries("ABCDEFGHIJK".split("").map((c, i) => [c, i]));
const RELS = { 0: 0, 1: 0, 2: 0, 3: 0, 4: 2, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 10: 21 };

/** Keep the positions not yet deleted for this loop length; extra residues become insertions after insertionPos. */
function gapped(region, deletions, insertionPos) {
	const n = deletions.length;
	const length = region.length;
	let annotations = deletions
		.slice(Math.max(n - length, 0))
		.map((p) => [p, " "])
		.sort(compareAnnotations);
	const insertions = Math.max(length - n, 0);
	if (insertions > 26) throw new NumberingError("too many insertions for AHo numbering");
	if (insertions) {
		const at = annotations.findIndex(([p, ins]) => p === insertionPos && ins === " ") + 1;
		annotations = [...annotations.slice(0, at), ...ahoChar.slice(0, insertions).map((c) => [insertionPos, c]), ...annotations.slice(at)];
	}
	return relabel(region, annotations);
}

/** Region B (AHo 1-10): indels on 8, taking N-terminal truncation into account. */
function regionB(region, insertionPos) {
	const length = region.length;
	if (!length) return region;
	const start = region[0][0][0];
	const stretch = 10 - (start - 1);
	let annotations;
	if (length > stretch) {
		annotations = [...range(start, insertionPos + 1), ...ahoChar.slice(0, length - stretch).map((c) => [insertionPos, c]), ...range(insertionPos + 1, 11)];
	} else {
		const ordered = [[insertionPos, " "], ...range(start, 11).filter(([p]) => p !== insertionPos)];
		annotations = ordered.slice(Math.max(stretch - length, 0)).sort(compareAnnotations);
	}
	return relabel(region, annotations);
}

/**
 * @param {string} annotations - The `#=GC RF` line: 'x' denotes a countable position (IMGT HMM state), '.' an insertion column.
 * @param {string} sequence - The aligned sequence row: uppercase = match, lowercase = insertion, '-' or '.' = gap.
 * @param {ChainType} chainType - The type of chain ('L', 'K', 'H', 'A', 'B', 'D', 'G').
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding AHo numbering.
 */
function applyAhoNumbering(annotations, sequence, chainType) {
	if (!(chainType in orderedDeletions)) throw new NumberingError(`AHo numbering needs a chain type (L, K, H, A, B, D, G), got ${chainType}`);
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions: r, end } = numberRegions(residues, states, STATE_STRING, REGION_STRING, REGION_INDEX, RELS, 11, [1, 3, 4, 5, 7, 9]);
	const numbered = [
		...r[0],
		...regionB(r[1], regions[0].insertionPos),
		...r[2],
		...gapped(r[3], orderedDeletions[chainType], regions[2].insertionPos),
		...r[4],
		...gapped(r[5], chainType === "A" ? cdr2DeletionsAlpha : cdr2Deletions, regions[4].insertionPos),
		...r[6],
		...gapped(r[7], fw3Deletions, regions[5].insertionPos),
		...r[8],
		...gapped(r[9], cdr3Deletions, regions[7].insertionPos),
		...r[10],
	];
	// AHo's domain extends one residue beyond IMGT for light chains (position 149).
	const last = numbered.at(-1);
	if (last && last[0][0] === 148 && last[0][1] === " " && end !== null && end + 1 < residues.length) {
		numbered.push([[149, " "], residues[end + 1]]);
	}
	return toResult(numbered);
}

export { applyAhoNumbering };
