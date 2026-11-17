/**
 * Convert one row of a Stockholm alignment against the 128-state IMGT HMM into a state vector.
 *
 * @param {string} annotations - The `#=GC RF` line: 'x' = match column (IMGT state), '.' = insert column.
 * @param {string} sequence - The aligned row: uppercase = residue on a match state, lowercase = inserted
 *   residue, '-' = deletion, '.' = gap in an insert column.
 * @returns {{ residues: string, states: Array<[[number, "m"|"i"|"d"], number|null]> }}
 *   `residues` is every residue of the row in order (uppercase). Each state vector entry is
 *   `[[imgtState, kind], residueIndex]` with `residueIndex` pointing into `residues` (null for a
 *   deletion). An insertion ("i") carries the state number of the match state that follows it. Residues in insert columns before the first or after the last match column are
 *   flanking sequence outside the domain (hmmalign shows them that way) and get no state;
 *   deletions at either end of the domain are dropped.
 */
export function parseAlignment(annotations, sequence) {
	if (annotations.length !== sequence.length) {
		throw new Error(`annotation line (${annotations.length}) and sequence row (${sequence.length}) differ in length`);
	}
	const firstMatchColumn = [...annotations].findIndex((a) => a !== ".");
	let lastMatchColumn = annotations.length - 1;
	while (lastMatchColumn >= 0 && annotations[lastMatchColumn] === ".") lastMatchColumn--;

	let residues = "";
	const states = [];
	let state = 0;
	for (let i = 0; i < annotations.length; i++) {
		const isMatch = annotations[i] !== ".";
		const ch = sequence[i];
		const isGap = ch === "-" || ch === ".";
		if (isMatch) {
			state++;
			if (isGap) {
				states.push([[state, "d"], null]);
			} else {
				states.push([[state, "m"], residues.length]);
				residues += ch.toUpperCase();
			}
		} else if (!isGap) {
			// As in HMMER/ANARCI state vectors, an insertion is labelled with the *next* match state.
			if (i > firstMatchColumn && i < lastMatchColumn) states.push([[state + 1, "i"], residues.length]);
			residues += ch.toUpperCase();
		}
	}
	if (state > 128) throw new Error(`alignment has ${state} match columns; expected at most 128 (IMGT HMM)`);

	let first = 0;
	while (first < states.length && states[first][0][1] === "d") first++;
	let last = states.length - 1;
	while (last >= first && states[last][0][1] === "d") last--;
	return { residues, states: states.slice(first, last + 1) };
}
