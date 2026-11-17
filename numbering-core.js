/**
 * Shared machinery for converting an IMGT-state alignment into a numbering scheme.
 *
 * Each scheme is described by:
 *  - stateString:  128 chars, one per IMGT state. "X" = has its own scheme position,
 *                  "I" = treated as an insertion on the previous scheme position.
 *  - regionString: 128 chars assigning each IMGT state to a region.
 *  - rels:         offset (scheme number - IMGT state) at the start of each region.
 *  - excludeDeletions: regions (CDRs) the scheme renumbers itself.
 *
 * numberRegions() applies the mapping region by region, then each scheme re-places
 * the insertions/deletions in its CDRs. The rules follow ANARCI (BSD 3-Clause).
 * Internally an insertion code of " " means "no insertion" (it sorts before "A").
 */

const UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
export const ALPHABET = [...UPPER, ...UPPER.map((c) => c + c), " "];
export const LETTERS = ALPHABET.slice(0, -1);

export class NumberingError extends Error {
	constructor(message) {
		super(message);
		this.name = "NumberingError";
	}
}

/** Compare two [number, insertion] annotations like Python tuples. */
export function compareAnnotations(a, b) {
	if (a[0] !== b[0]) return a[0] - b[0];
	return a[1] < b[1] ? -1 : a[1] > b[1] ? 1 : 0;
}

const ENFORCED_PATTERNS = [
	[[25, "m"], [26, "m"], [27, "m"], [28, "i"]],
	[[38, "i"], [38, "m"], [39, "m"], [40, "m"]],
	[[54, "m"], [55, "m"], [56, "m"], [57, "i"]],
	[[65, "i"], [65, "m"], [66, "m"], [67, "m"]],
	[[103, "m"], [104, "m"], [105, "m"], [106, "i"]],
	[[117, "i"], [117, "m"], [118, "m"], [119, "m"]],
];

function bufferRegion(s) {
	if (s < 23) return -1;
	if (s >= 25 && s < 28) return 0;
	if (s > 37 && s <= 40) return 1;
	if (s >= 54 && s < 57) return 2;
	if (s > 64 && s <= 67) return 3;
	if (s >= 103 && s < 106) return 4;
	if (s > 116 && s <= 119) return 5;
	return null;
}

/**
 * Move insertions HMMER may have misplaced: insertions on the edges of the frameworks
 * next to each CDR are pushed into the CDR; FW1 insertions become N-terminal states
 * when at least as many N-terminal states are missing.
 */
export function smoothInsertions(stateVector) {
	let buffer = [];
	const out = [];
	let region = null;
	for (const entry of stateVector) {
		const [[stateId], ] = entry;
		const reg = bufferRegion(stateId);
		if (reg !== null) {
			buffer.push(entry);
			region = reg;
			continue;
		}
		if (buffer.length) {
			const nIns = buffer.filter(([[, k]]) => k === "i").length;
			if (nIns === 0) {
				out.push(...buffer);
			} else if (region === -1) {
				let ntDeletions = buffer[0][0][0] - 1;
				for (const [[, k], idx] of buffer) {
					if (k === "d" || idx === null) ntDeletions++;
					else break;
				}
				if (ntDeletions >= nIns) {
					let newStates = buffer.filter(([[, k]]) => k === "m").map(([s]) => s);
					const first = newStates[0][0];
					const residues = buffer.filter(([[, k]]) => k !== "d");
					const nAdd = residues.length - newStates.length;
					const prefix = [];
					for (let x = first - nAdd; x < first; x++) prefix.push([x, "m"]);
					newStates = [...prefix, ...newStates];
					residues.forEach((r, i) => out.push([newStates[i], r[1]]));
				} else {
					out.push(...buffer);
				}
			} else {
				const residues = buffer.filter(([[, k]]) => k !== "d");
				const pattern = ENFORCED_PATTERNS[region];
				let newStates;
				if (region % 2) {
					newStates = [
						...Array(Math.max(0, residues.length - 3)).fill(pattern[0]),
						...pattern.slice(Math.max(4 - residues.length, 1)),
					];
				} else {
					newStates = [...pattern.slice(0, 3), ...Array(Math.max(0, residues.length - 3)).fill(pattern[2])];
				}
				residues.forEach((r, i) => out.push([newStates[i], r[1]]));
			}
			buffer = [];
		}
		out.push(entry);
	}
	out.push(...buffer); // an alignment ending inside a buffered stretch is kept as aligned
	return out;
}

/**
 * Number the aligned residues region by region using the scheme tables.
 * @returns {{ regions: Array<Array<[[number, string], string]>>, start: number|null, end: number|null }}
 */
export function numberRegions(sequence, stateVector, stateString, regionString, regionIndex, rels, nRegions, excludeDeletions) {
	rels = { ...rels };
	stateVector = smoothInsertions(stateVector);
	const regions = Array.from({ length: nRegions }, () => []);
	let insertion = -1;
	let previousStateId = 1;
	let previousKind = "d";
	let start = null;
	let end = null;
	let region = null;
	const letter = (i) => ALPHABET.at(i);

	for (const [[stateId, kind], si] of stateVector) {
		if (kind !== "i" || region === null) region = regionIndex[regionString[stateId - 1]];

		if (kind === "m") {
			if (stateString[stateId - 1] === "I") {
				if (previousKind !== "d") insertion++;
				rels[region]--;
			} else {
				insertion = -1;
			}
			regions[region].push([[stateId + rels[region], letter(insertion)], sequence[si]]);
			previousStateId = stateId;
			if (start === null) start = si;
			end = si;
			previousKind = kind;
		} else if (kind === "i") {
			insertion++;
			regions[region].push([[previousStateId + rels[region], letter(insertion)], sequence[si]]);
			if (start === null) start = si;
			end = si;
			previousKind = kind;
		} else {
			previousKind = kind;
			if (stateString[stateId - 1] === "I") {
				rels[region]--;
				continue;
			}
			insertion = -1;
			previousStateId = stateId;
		}

		if (insertion >= 25 && excludeDeletions.includes(region)) insertion = 0;
		if (insertion >= 25) throw new NumberingError("too many insertions for the numbering scheme to handle");
	}
	return { regions, start, end };
}

/** Give the residues of `region` the positions in `annotations` (in order). */
export function relabel(region, annotations) {
	return region.map(([, aa], i) => [annotations[i], aa]);
}

/** Positions [start, end) as annotations without insertion codes. */
export function range(start, end) {
	const out = [];
	for (let p = start; p < end; p++) out.push([p, " "]);
	return out;
}

/**
 * IMGT-style symmetric numbering of a CDR between `start` and `end` (exclusive).
 * Positions fill from both ends towards the centre; residues beyond `maxLength` become
 * insertions on the two central positions. null marks an unoccupied position.
 */
export function symmetricCdr(length, maxLength, start, end) {
	const n = Math.max(length, maxLength);
	const ann = new Array(n).fill(null);
	if (length === 0) return ann;
	if (length === 1) {
		ann[0] = [start, " "];
		return ann;
	}
	let front = 0;
	let back = -1;
	for (let i = 0; i < Math.min(length, maxLength); i++) {
		if (i % 2) {
			ann[n + back] = [end + back, " "];
			back--;
		} else {
			ann[front] = [start + front, " "];
			front++;
		}
	}
	const centre = [];
	ann.forEach((v, i) => v === null && centre.push(i));
	if (!centre.length) return ann;
	const centreLeft = ann[Math.min(...centre) - 1][0];
	const centreRight = ann[Math.max(...centre) + 1][0];
	const half = Math.floor(maxLength / 2);
	const frontFactor = maxLength % 2 === 0 ? half : half + 1;
	const backFactor = half;
	const reverse = [...LETTERS].reverse();
	for (let i = 0; i < Math.max(0, length - maxLength); i++) {
		if (i % 2 === 0) {
			ann[n + back] = [centreRight, reverse.at(back + backFactor)];
			back--;
		} else {
			ann[front] = [centreLeft, LETTERS.at(front - frontFactor)];
			front++;
		}
	}
	return ann;
}

/**
 * CDR3 positions for Kabat and Chothia (identical in both schemes).
 * Heavy: 93-102, insertions on 100. Light: 89-97, insertions on 95.
 */
export function kabatChothiaCdr3(length, chain) {
	const heavy = chain === "H";
	const nPositions = heavy ? 10 : 9;
	const insertAt = heavy ? 100 : 95;
	const deletionOrder = heavy ? [100, 99, 98, 97, 96, 95, 101, 102, 94, 93] : [95, 94, 93, 92, 91, 96, 97, 90, 89];
	const insertions = Math.max(length - nPositions, 0);
	if (insertions > 26) throw new NumberingError("too many insertions for the numbering scheme to handle");
	const kept = deletionOrder.slice(Math.max(0, nPositions - length)).map((p) => [p, " "]);
	const ins = LETTERS.slice(0, insertions).map((c) => [insertAt, c]);
	return [...kept, ...ins].sort(compareAnnotations);
}

/**
 * Format the final result in the original API shape: `numbering` holds a plain number
 * for positions without an insertion code and a string such as "52A" otherwise.
 */
export function toResult(numbered) {
	return {
		concatSequence: numbered.map(([, aa]) => aa).join(""),
		numbering: numbered.map(([[n, ins]]) => (ins.trim() ? `${n}${ins.trim()}` : n)),
	};
}
