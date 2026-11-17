import { parseAlignment } from "./alignment-states.js";
import { NumberingError, numberRegions, symmetricCdr, toResult } from "./numbering-core.js";

// IMGT region definitions. HMM match states are IMGT positions, so these regions
// are used directly. In the CDRs, gaps and insertions are placed symmetrically
// around the centre positions (insertionPos): 111.1, 111.2 ... and ... 112.2, 112.1
// in IMGT notation, returned here as "111A", "111B" ... "112B", "112A".
const regions = [
	{ start: 1, end: 26, insertionPos: null }, // FR1
	{ start: 27, end: 38, insertionPos: [32, 33] }, // CDR1
	{ start: 39, end: 55, insertionPos: null }, // FR2
	{ start: 56, end: 65, insertionPos: [60, 61] }, // CDR2
	{ start: 66, end: 104, insertionPos: null }, // FR3
	{ start: 105, end: 117, insertionPos: [111, 112] }, // CDR3
	{ start: 118, end: 128, insertionPos: null }, // FR4
];

const stateString = "X".repeat(128); // every IMGT state is an IMGT position
const regionString = regions.map((r, i) => String(i + 1).repeat(r.end - r.start + 1)).join("");
const regionIndex = Object.fromEntries(regions.map((_, i) => [String(i + 1), i]));
const rels = Object.fromEntries(regions.map((_, i) => [i, 0]));
const cdrRegions = [1, 3, 5];

/** Renumber a CDR symmetrically between its region start and end. */
function renumberCdr(region, { start, end }) {
	const positions = symmetricCdr(region.length, end - start + 1, start, end + 1).filter((a) => a !== null);
	return region.map(([, aa], i) => [positions[i], aa]);
}

/**
 * Applies IMGT numbering to an amino acid sequence based on Stockholm annotations.
 *
 * @param {string} annotations - The `#=GC RF` line: 'x' denotes a countable position (IMGT HMM state), '.' an insert column.
 * @param {string} sequence - The aligned sequence row: uppercase = match, lowercase = insertion, '-' or '.' = gap.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their
 *   corresponding IMGT numbering, one entry per residue (e.g. 110, 111, "111A", "112A", 112).
 */
function applyImgtNumbering(annotations, sequence) {
	const { residues, states } = parseAlignment(annotations, sequence);
	const { regions: numbered } = numberRegions(residues, states, stateString, regionString, regionIndex, rels, regions.length, cdrRegions);
	// 13 positions + insertions on 111 and 112 (up to "ZZ" each) cover any realistic CDR3.
	if (numbered[5].length > 117) throw new NumberingError("CDR3 too long for IMGT numbering");
	return toResult(
		numbered.flatMap((region, i) => (cdrRegions.includes(i) ? renumberCdr(region, regions[i]) : region)),
	);
}

export { applyImgtNumbering };
