/**
 * @param {string} annotations - A string where 'x' denotes a countable position
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: number[] }} - Returns both the sequence and their corresponding imgt numbering.
 */
function applyImgtNumbering(annotations, sequence) {
	let imgtLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	let imgtNumber = 1;
	let annotatedSequence = "";
	let numbering = [];

	for (let i = 0; i < imgtLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];
				if (
					annotations[i + 1] === "." &&
					(imgtNumber < 111 || imgtNumber > 117)
				) {
					numbering.push(imgtNumber);
					console.log("general loop", imgtNumber);
					let insertionCounter = 1;
					sequenceIndex++;
					i++;
					while (annotations[i] === "." && sequenceIndex <= sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(imgtNumber + 0.1 * insertionCounter);
						console.log(imgtNumber + 0.1 * insertionCounter);
						insertionCounter++;
						sequenceIndex++;
						i++;
					}
					imgtNumber++;
				} else if (
					imgtNumber >= 111 &&
					imgtNumber <= 117 &&
					annotations[i + 1] === "."
				) {
					if (annotations[i + 1] === ".") {
						numbering.push(imgtNumber);
						console.log("specific loop", imgtNumber);
						let insertions = [];
						while (i < imgtLength && annotations[i + 1] === ".") {
							insertions.push(sequence[i]);
							i++;
						}
						console.log("insertions", insertions);
						if (insertions.length % 2 === 0) {
							let temp = imgtNumber;
							for (let j = 0; j < insertions.length / 2; j++) {
								numbering.push(Number((temp + 0.1).toFixed(1)));
								temp = Number((temp + 0.1).toFixed(1));
								console.log("specific loop even", imgtNumber + 0.1);
							}
							imgtNumber++;
							for (let j = insertions.length / 2; j < insertions.length; j++) {
								let forwardcounter = insertions.length - j;
								numbering.push(imgtNumber + 0.1 * forwardcounter);
							}
						} else if (insertions.length % 2 === 1) {
							console.log("specific loop insert", imgtNumber);
							console.log(Number.isInteger(imgtNumber));
							let temp = imgtNumber;
							for (let j = 0; j < Math.floor(insertions.length / 2) + 1; j++) {
								numbering.push(Number((temp + 0.1).toFixed(1)));
								temp = Number((temp + 0.1).toFixed(1));
								console.log("temp", temp);
							}
							imgtNumber++;
							for (
								let j = Math.floor(insertions.length / 2) + 1;
								j < insertions.length;
								j++
							) {
								let reversecounter = insertions.length - j;
								console.log("reversecounter", reversecounter);
								numbering.push(imgtNumber + 0.1 * reversecounter);
							}
						}
					}
					imgtNumber++;
				}
			}
			numbering.push(imgtNumber);
			imgtNumber++;
		} else {
			annotatedSequence += ".";
			imgtNumber++;
		}
		sequenceIndex++;
	}
	let concatSequence = "";
	for (const seqchar of annotatedSequence) {
		if (seqchar === "." || seqchar === "-") {
			continue;
		} else {
			concatSequence += seqchar;
		}
	}
	return { concatSequence, numbering };
}

export { applyImgtNumbering };

//Maybe we can use different functions for different numbering positions
//CDR1
//CDR1 has a range from 27 (inc.) to 39 (exc.) and has a theoretical maximum length of 12.
//CDR2
//CDR2 has a range from 56 (inc.) to 66 (exc.) and has a theoretical length of 10.
//CDR3
//CDR3 has a range from 105 (inc.) to 118 (exc.). Insertions are placed on 112 and 111 symetrically.
