/**
 * @param {string} annotations - A string where 'x' denotes a countable position
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding aho numbering.
 */

function applyAhoNumbering(annotations, sequence) {
	const ahoLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	let ahoNumber = 1;
	const ahoChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	let numbering = [];

	for (let i = 0; i < ahoLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(ahoNumber);
					ahoNumber++;
				} else if (annotations[i] === ".") {
					let insertionCounter = 0;
					while (annotations[i] === "." && sequenceIndex < sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(`${ahoNumber - 1}${ahoChar[insertionCounter]}`);
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
				ahoNumber++;
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
			ahoNumber++;
		}
	}

	let concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyAhoNumbering };

// Example usage
const annotations = "xxxx.xxxx....xx";
const sequence = "QVQLVQSGAAAEE";
const result = applyAhoNumbering(annotations, sequence);
console.log(result);
