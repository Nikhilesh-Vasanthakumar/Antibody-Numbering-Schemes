/**
 * @param {string} annotations - A string where 'x' denotes a countable position and '.' denotes an insertion point.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding Chothia numbering.
 */
function applyChothiaNumbering(annotations, sequence) {
	const chothiaLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	let chothiaNumber = 1;
	const chothiaChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	let numbering = [];

	for (let i = 0; i < chothiaLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(chothiaNumber);
					chothiaNumber++;
				} else if (annotations[i] === ".") {
					let insertionCounter = 0;
					while (annotations[i] === "." && sequenceIndex < sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(
							`${chothiaNumber - 1}${chothiaChar[insertionCounter]}`,
						);
						insertionCounter++;
						sequenceIndex++;
						i++;
					}
					i--;
					continue;
				}

				sequenceIndex++;
			} else {
				annotatedSequence += ".";
				chothiaNumber++;
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
			chothiaNumber++;
		}
	}

	let concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyChothiaNumbering };
