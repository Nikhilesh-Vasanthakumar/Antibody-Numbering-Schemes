
# VH and VL Sequence Analysis Tool with HMM-Based Numbering System for Antibody Research

This project provides implementations of various antibody numbering schemes including IMGT, Chothia, Aho, Kabat, and alignment using the Stockholm format. These schemes are crucial for the standardized numbering and alignment of antibody sequences, facilitating consistent analysis and comparison.
## Data Used

Sequences used in this project were collected from ANARCI (Antigen Receptor Numbering and Receptor Classification). The collected sequences were used to create Hidden Markov Models (HMMs). These HMMs were then employed to perform the alignment of antibody sequences according to the specified numbering schemes.
Steps Involved:

    - Sequence Collection: Antibody sequences were collected from the ANARCI database.
    - HMM Creation: The collected sequences were used to create HMMs using HMMER.
    - Alignment: The created HMMs were used to align the sequences in accordance with the various numbering schemes implemented in this project.
## Installation

1.To install the necessary packages and set up the project, follow these steps:
```{bash}
git clone https://github.com/Nikhilesh-Vasanthakumar/Antibody-Numbering-Schemes.git
```
2.Navigate to the project directory:
```{bash}
cd Antibody-Numbering-Schemes
```

**Note:This does not contain the full implementation of the Numbering Schemes.It only contains the algorithm behind the Numbering schemes.Please contact BionamicAB for more information. [Email](morten.krogh@bionamic.io)**
## IMGT ( ImMunoGeneTics ) Numbering Algorithm


```{JavaScript}
/**
 * Applies IMGT numbering to an amino acid sequence based on Stockholm annotations.
 *
 * @param {string} annotations - A string where 'x' denotes a countable position.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: number[] }} - Returns both the sequence and their corresponding IMGT numbering.
 */
function applyImgtNumbering(annotations, sequence) {
    // Determine the length to process based on the shorter of annotations or sequence
    const imgtLength = Math.min(annotations.length, sequence.length);
    let sequenceIndex = 0;
    let imgtNumber = 1;
    let annotatedSequence = "";
    const numbering = [];

    for (let i = 0; i < imgtLength; i++) {
        // Handle lowercase sequence characters with '.' annotations
        if (
            sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase() &&
            annotations[i] === "."
        ) {
            annotatedSequence += sequence[sequenceIndex];
            sequenceIndex++;
            imgtNumber++;
            continue;
        }

        // Handle non '-' sequence characters
        if (sequence[sequenceIndex] !== "-") {
            annotatedSequence += sequence[sequenceIndex];

            // Handle insertion points outside of CDH3 region
            if (annotations[i] === "." && (imgtNumber < 111 || imgtNumber > 117)) {
                let insertionCounter = 1;
                imgtNumber--;
                while (annotations[i] === "." && sequenceIndex <= sequence.length) {
                    annotatedSequence += sequence[sequenceIndex];
                    numbering.push(imgtNumber + 0.1 * insertionCounter);
                    insertionCounter++;
                    sequenceIndex++;
                    i++;
                }
                imgtNumber++;
                // Handle insertion points within CDH3 region
            } else if (
                imgtNumber >= 111 &&
                imgtNumber <= 117 &&
                annotations[i] === "."
            ) {
                if (annotations[i + 1] === ".") {
                    numbering.push(imgtNumber);
                    const insertions = [];

                    // Collect insertions
                    while (i < imgtLength && annotations[i + 1] === ".") {
                        insertions.push(sequence[i]);
                        i++;
                    }

                    // Distribute insertions for even count
                    if (insertions.length % 2 === 0) {
                        let temp = imgtNumber;

                        for (let j = 0; j < insertions.length / 2; j++) {
                            numbering.push(Number((temp + 0.1).toFixed(1)));
                            temp = Number((temp + 0.1).toFixed(1));
                        }
                        imgtNumber++;

                        for (let j = insertions.length / 2; j < insertions.length; j++) {
                            const forwardCounter = insertions.length - j;
                            numbering.push(imgtNumber + 0.1 * forwardCounter);
                        }
                    } else if (insertions.length % 2 === 1) {
                        // Distribute insertions for odd count
                        let temp = imgtNumber;

                        for (let j = 0; j < Math.floor(insertions.length / 2) + 1; j++) {
                            numbering.push(Number((temp + 0.1).toFixed(1)));
                            temp = Number((temp + 0.1).toFixed(1));
                        }
                        imgtNumber++;

                        for (
                            let j = Math.floor(insertions.length / 2) + 1;
                            j < insertions.length;
                            j++
                        ) {
                            const reverseCounter = insertions.length - j;
                            numbering.push(imgtNumber + 0.1 * reverseCounter);
                        }
                    }
                }
                imgtNumber++;
            }
        }
        // Add the current IMGT number to the numbering list
        numbering.push(imgtNumber);
        imgtNumber++;
        sequenceIndex++;
    }

    // Remove '.' and '-' characters from the annotated sequence
    let concatSequence = "";
    for (const seqChar of annotatedSequence) {
        if (seqChar !== "." && seqChar !== "-") {
            concatSequence += seqChar;
        }
    }

    return { concatSequence, numbering };
}

export { applyImgtNumbering };
```
## Kabat Numbering Algorithm

```{JavaScript}

/**
 * @param {string} annotations - A string where 'x' denotes a countable position and '.' denotes an insertion point.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding Kabat numbering.
 */
function applyKabatNumberingHeavy(annotations, sequence) {
	const kabatLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	const kabatChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	const numbering = [];

	// Kabat region definitions for heavy chain
	const regions = [
		{ start: 1, end: 29, insertionPos: 6 }, // Region 1: Insertions at position 6
		{ start: 30, end: 35, insertionPos: 35 }, // Region 2 (CDRH1): Insertions at position 35
		{ start: 36, end: 49, insertionPos: null }, // Region 3: No specific insertion position
		{ start: 50, end: 65, insertionPos: 52 }, // Region 4 (CDRH2): Insertions at position 52
		{ start: 66, end: 94, insertionPos: null }, // Region 5: No specific insertion position
		{ start: 95, end: 102, insertionPos: 100 }, // Region 6 (CDRH3): Insertions at position 100
		{ start: 103, end: 113, insertionPos: null }, // Region 7: No specific insertion position
	];

	let currentRegion = 0;
	let kabatNumber = regions[currentRegion].start;
	let insertionCounter = 0;

	for (let i = 0; i < kabatLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase()) {
				annotatedSequence += sequence[sequenceIndex];
				sequenceIndex++;
				continue;
			}
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(kabatNumber);
					kabatNumber++;
				} else if (annotations[i] === ".") {
					while (annotations[i] === "." && sequenceIndex < sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(`${kabatNumber - 1}${kabatChar[insertionCounter]}`);
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
				kabatNumber++;
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
			kabatNumber++;
		}

		// Check if we need to move to the next region
		if (
			kabatNumber > regions[currentRegion].end &&
			currentRegion < regions.length - 1
		) {
			currentRegion++;
			kabatNumber = regions[currentRegion].start;
			insertionCounter = 0; // Reset insertion counter for the new region
		}
	}

	const concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

/**
 * @param {string} annotations - A string where 'x' denotes a countable position and '.' denotes an insertion point.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding Kabat numbering.
 */
function applyKabatNumberingLight(annotations, sequence) {
	const kabatLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	const kabatChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	const numbering = [];

	// Kabat region definitions for light chain
	const regions = [
		{ start: 1, end: 10, insertionPos: 6 }, // Region 1: Insertions at position 6
		{ start: 11, end: 30, insertionPos: null }, // Region 2: No specific insertion position
		{ start: 31, end: 35, insertionPos: 35 }, // Region 3 (CDRH1): Insertions at position 35
		{ start: 36, end: 51, insertionPos: null }, // Region 4: No specific insertion position
		{ start: 52, end: 58, insertionPos: 52 }, // Region 5 (CDRH2): Insertions at position 52
		{ start: 59, end: 93, insertionPos: null }, // Region 6: No specific insertion position
		{ start: 94, end: 103, insertionPos: 100 }, // Region 7 (CDRH3): Insertions at position 100
		{ start: 104, end: 128, insertionPos: null }, // Region 8: No specific insertion position
	];

	let currentRegion = 0;
	let kabatNumber = regions[currentRegion].start;
	let insertionCounter = 0;

	for (let i = 0; i < kabatLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase()) {
				annotatedSequence += sequence[sequenceIndex];
				sequenceIndex++;
				continue;
			}
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(kabatNumber);
					kabatNumber++;
				} else if (annotations[i] === ".") {
					while (annotations[i] === "." && sequenceIndex < sequence.length) {
						annotatedSequence += sequence[sequenceIndex];
						numbering.push(`${kabatNumber - 1}${kabatChar[insertionCounter]}`);
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
				kabatNumber++;
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
			kabatNumber++;
		}

		// Check if we need to move to the next region
		if (
			kabatNumber > regions[currentRegion].end &&
			currentRegion < regions.length - 1
		) {
			currentRegion++;
			kabatNumber = regions[currentRegion].start;
			insertionCounter = 0; // Reset insertion counter for the new region
		}
	}

	const concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyKabatNumberingLight };

export { applyKabatNumberingHeavy };
```
## Chothia Numbering Algorithm

```{JavaScript}
/**
 * @param {string} annotations - A string where 'x' denotes a countable position and '.' denotes an insertion point.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding Chothia numbering.
 */
function applyChothiaNumberingHeavy(annotations, sequence) {
	const chothiaLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	const chothiaChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	const numbering = [];

	// Chothia region definitions for heavy chain
	const regions = [
		{ start: 1, end: 26, insertionPos: 6 },
		{ start: 27, end: 38, insertionPos: 31 },
		{ start: 39, end: 55, insertionPos: null },
		{ start: 56, end: 65, insertionPos: 52 },
		{ start: 66, end: 93, insertionPos: null },
		{ start: 94, end: 103, insertionPos: 100 },
		{ start: 104, end: 118, insertionPos: null },
	];

	let currentRegion = 0;
	let chothiaNumber = regions[currentRegion].start;
	let insertionCounter = 0;

	for (let i = 0; i < chothiaLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase()) {
				annotatedSequence += sequence[sequenceIndex];
				sequenceIndex++;
				continue;
			}
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(chothiaNumber);
					chothiaNumber++;
				} else if (annotations[i] === ".") {
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
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
		}

		// Check if we need to move to the next region
		if (
			chothiaNumber > regions[currentRegion].end &&
			currentRegion < regions.length - 1
		) {
			currentRegion++;
			chothiaNumber = regions[currentRegion].start;
			insertionCounter = 0; // Reset insertion counter for the new region
		}
	}

	const concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

/**
 * @param {string} annotations - A string where 'x' denotes a countable position and '.' denotes an insertion point.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding Chothia numbering.
 */
function applyChothiaNumberingLight(annotations, sequence) {
	const chothiaLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	const chothiaChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	const numbering = [];

	// Chothia region definitions for light chain
	const regions = [
		{ start: 1, end: 23, insertionPos: null },
		{ start: 24, end: 34, insertionPos: 30 },
		{ start: 35, end: 49, insertionPos: null },
		{ start: 50, end: 56, insertionPos: 52 },
		{ start: 57, end: 88, insertionPos: null },
		{ start: 89, end: 98, insertionPos: 95 },
		{ start: 99, end: 107, insertionPos: null },
	];

	let currentRegion = 0;
	let chothiaNumber = regions[currentRegion].start;
	let insertionCounter = 0;

	for (let i = 0; i < chothiaLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase()) {
				annotatedSequence += sequence[sequenceIndex];
				sequenceIndex++;
				continue;
			}
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					numbering.push(chothiaNumber);
					chothiaNumber++;
				} else if (annotations[i] === ".") {
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
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
		}

		// Check if we need to move to the next region
		if (
			chothiaNumber > regions[currentRegion].end &&
			currentRegion < regions.length - 1
		) {
			currentRegion++;
			chothiaNumber = regions[currentRegion].start;
			insertionCounter = 0; // Reset insertion counter for the new region
		}
	}

	const concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyChothiaNumberingLight };

export { applyChothiaNumberingHeavy };
```
## Aho Numbering Algorithm
```{JavaScript}
/**
 * @typedef {'L' | 'K' | 'H' | 'A' | 'B' | 'D' | 'G'} ChainType
 */

/**
 * @param {string} annotations - A string where 'x' denotes a countable position.
 * @param {string} sequence - The amino acid sequence to be numbered.
 * @param {ChainType} chainType - The type of chain ('L', 'K', 'H', 'A', 'B', 'D', 'G').
 * @returns {{ concatSequence: string, numbering: (number|string)[] }} - Returns both the sequence and their corresponding AHo numbering.
 */
function applyAhoNumbering(annotations, sequence, chainType) {
	const ahoLength = Math.min(annotations.length, sequence.length);
	let sequenceIndex = 0;
	const ahoChar = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".split("");
	let annotatedSequence = "";
	const numbering = [];

	// AHo region definitions
	const regions = [
		{ start: 1, end: 10, insertionPos: 8 }, // Region B
		{ start: 11, end: 24, insertionPos: null }, // Region C
		{ start: 25, end: 42, insertionPos: 28 }, // Region D
		{ start: 43, end: 57, insertionPos: null }, // Region E
		{ start: 58, end: 77, insertionPos: 63 }, // Region F
		{ start: 78, end: 93, insertionPos: 86 }, // Region H
		{ start: 94, end: 106, insertionPos: null }, // Region I
		{ start: 107, end: 138, insertionPos: 123 }, // Region J
		{ start: 139, end: 149, insertionPos: null }, // Region K
	];

	const orderedDeletions = {
		L: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		K: [28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		H: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		A: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		B: [28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		D: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
		G: [28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
	};

	let currentRegion = 0;
	let ahoNumber = regions[currentRegion].start;
	let insertionCounter = 0;

	for (let i = 0; i < ahoLength; i++) {
		if (annotations[i] === "x" || annotations[i] === ".") {
			if (sequence[sequenceIndex] === sequence[sequenceIndex].toLowerCase()) {
				annotatedSequence += sequence[sequenceIndex];
				sequenceIndex++;
				continue;
			}
			if (sequence[sequenceIndex] !== "-") {
				annotatedSequence += sequence[sequenceIndex];

				if (annotations[i] === "x") {
					// Check for deletions
					if (orderedDeletions[chainType].includes(ahoNumber)) {
						ahoNumber++;
					}
					numbering.push(ahoNumber);
					ahoNumber++;
				} else if (annotations[i] === ".") {
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
				sequenceIndex++;
			}
		} else {
			annotatedSequence += ".";
		}

		// Check for region change
		if (
			currentRegion < regions.length - 1 &&
			ahoNumber > regions[currentRegion].end
		) {
			currentRegion++;
			ahoNumber = regions[currentRegion].start;
			insertionCounter = 0; // Reset insertion counter for the new region
		}
	}

	const concatSequence = annotatedSequence.replace(/[.-]/g, "");
	return { concatSequence, numbering };
}

export { applyAhoNumbering };
```
## Conclusion

In conclusion, this study explored the role of plasmacytoid dendritic cells (pDCs) in the context of Pseudomonas infection and investigated potential sex differences in pDC function and gene expression. The findings revealed distinct gene expression patterns and pathway activations between males and females at different timepoints during Pseudomonas treatment. Additionally, the study highlighted the integral role of pDCs in antitumor immune responses, autoimmune responses, and their potential as therapeutic targets. Further investigations into sex differences and the involvement of pDCs in infectious diseases and autoimmune disorders, such as lupus, can provide valuable insights for understanding immune responses and developing targeted therapies.