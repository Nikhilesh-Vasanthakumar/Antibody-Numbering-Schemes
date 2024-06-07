/**
 * @param {string} annotations - A string where 'x' denotes a countable position
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding kabat numbering.
 */
function applyKabatNumberingHeavy(annotations, sequence) {
	const kabatLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	let kabatNumber = 1;
	const kabatChar = [
		"A",
		"B",
		"C",
		"D",
		"E",
		"F",
		"G",
		"H",
		"I",
		"J",
		"K",
        "L",
        "M",
        "N",
        "O",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "U",
        "V",
        "W",
        "X",
        "Y",
        "Z"
	];
	let annotatedSequence = "";
	let numbering = [];

	for (let i = 0; i < kabatLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];
				if (annotations[i + 1] === ".") {
					numbering.push(kabatNumber);
					let insertionCounter = 0;
					sequenceIndex++;
					i++;
					while (annotations[i] === "." && sequenceIndex < sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(String(kabatNumber) + kabatChar[insertionCounter]);
						insertionCounter++;
						sequenceIndex++;
						i++;
					}
					kabatNumber++;
					console.log(JSON.stringify(numbering));
				}
				numbering.push(kabatNumber);
				kabatNumber++;
			} else {
				annotatedSequence += ".";
				kabatNumber++;
			}
			sequenceIndex++;
		} else if (annotations[i] === ".") {
			annotatedSequence += ".";
		} else {
			annotatedSequence += ".";
			if (sequence[sequenceIndex] !== "." && sequence[sequenceIndex] !== "-") {
				sequenceIndex++;
			}
		}
	}

	let concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyKabatNumberingHeavy };
