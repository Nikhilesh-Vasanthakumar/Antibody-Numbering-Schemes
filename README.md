
# VH and VL Sequence Analysis Tool with HMM-Based Numbering System for Antibody Research

This project provides implementations of the IMGT, Kabat, Chothia and AHo antibody
numbering schemes in JavaScript. Each function takes a sequence aligned in Stockholm
format to an IMGT-state HMM and returns its numbering. These schemes are the standard
way to number and compare antibody sequences.

## Data Used

The HMMs are built from germline sequences collected from ANARCI (Antigen Receptor
Numbering and Receptor Classification), which come from IMGT.

Steps involved:

    - Sequence collection: IMGT-aligned germline V and J genes from ANARCI.
    - HMM creation: one HMM per species and chain type, built with HMMER (`hmmbuild --hand`),
      so that its 128 match states are exactly IMGT positions 1-128.
    - Alignment: sequences are aligned to the HMMs (e.g. `hmmalign`), and the numbering
      functions below turn the alignment into IMGT, Kabat, Chothia or AHo numbers.

A companion Python package, `antibody-numbering`, runs the whole pipeline from sequence to
numbering. It is the reference these functions are tested against.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Nikhilesh-Vasanthakumar/Antibody-Numbering-Schemes.git
```
2. Navigate to the project directory and run the tests (Node.js 18 or newer, no dependencies):
```bash
cd Antibody-Numbering-Schemes
npm test
```

## Input and output

Every function takes the two lines of a Stockholm alignment against the 128-state IMGT HMM:

* `annotations`: the `#=GC RF` line. `x` = match column (IMGT state), `.` = insert column.
* `sequence`: the aligned row. Uppercase = residue on a match state, lowercase = inserted
  residue, `-` = deletion, `.` = gap in an insert column. Flanking residues before the first or
  after the last match column (as `hmmalign` shows them) are not numbered.

Every function returns `{ concatSequence, numbering }`, with exactly one entry in `numbering`
per residue of `concatSequence`. A position without an insertion code is a number (`52`).
A position with one is a string (`"52A"`). IMGT insertions such as 111.1 / 112.1 are returned
as `"111A"` / `"112A"`. A domain the scheme cannot number (for example an ultralong bovine CDR3
in Kabat) throws a `NumberingError`.

```js
import { applyKabatNumberingHeavy } from "./kabat-scheme.js";

const { concatSequence, numbering } = applyKabatNumberingHeavy(rfLine, alignedRow);
// numbering: [1, 2, ..., 52, "52A", 53, ..., 82, "82A", "82B", "82C", 83, ..., 100, "100A", ...]
```

## IMGT (ImMunoGeneTics) Numbering Algorithm

`applyImgtNumbering(annotations, sequence)` (`imgt-scheme.js`)

The HMM states are IMGT positions, so the frameworks are numbered directly. The CDRs are
renumbered with the IMGT gap rules. CDR1 (27-38), CDR2 (56-65) and CDR3 (105-117) fill from
both ends towards the centre. Short loops therefore lose their central positions, and long
loops get insertions on the central pair: 32/33, 60/61 and 111/112 (111A, 111B ... 112B, 112A).

## Kabat Numbering Algorithm

`applyKabatNumberingHeavy(annotations, sequence)` and `applyKabatNumberingLight(annotations, sequence)` (`kabat-scheme.js`)

The IMGT states are mapped onto Kabat positions region by region. Kabat insertion points:
heavy 6, 35 (CDRH1), 52 (CDRH2), 82 (FR3) and 100 (CDRH3); light 27 (CDRL1), 52 (CDRL2) and
95 (CDRL3). Short CDR3s lose positions in Kabat's fixed deletion order.

## Chothia Numbering Algorithm

`applyChothiaNumberingHeavy(annotations, sequence)` and `applyChothiaNumberingLight(annotations, sequence)` (`chothia-scheme.js`)

Same as Kabat except CDR1, which uses the structural insertion point: 31 for heavy and 30 for
light chains. Light-chain FR3 insertions go on 68.

## Aho Numbering Algorithm

`applyAhoNumbering(annotations, sequence, chainType)` (`aho-scheme.js`), where `chainType` is
`'H'`, `'K'`, `'L'` (or TCR `'A'`, `'B'`, `'G'`, `'D'`).

AHo reserves a position for every structurally observed residue (1-149), so loops are gapped,
not given insertion codes. Positions are removed in a fixed order around 8 (FR1), 28 and 36
(CDR1, with a chain-type specific order), 63 (CDR2), 86/85 (FR3) and 123 (CDR3).

## Files

| File | Contents |
|---|---|
| `imgt-scheme.js`, `kabat-scheme.js`, `chothia-scheme.js`, `aho-scheme.js` | The numbering functions and each scheme's region tables |
| `alignment-states.js` | `parseAlignment`: turns an RF line and an aligned row into IMGT HMM states |
| `numbering-core.js` | Shared region mapping, CDR gap rules and insertion helpers |
| `test/numbering.test.js` | Tests (`npm test`) |
| `test/fixtures/reference.json` | 507 PDB domains, one per distinct loop-length/insertion pattern, with the expected numbering |
| `test/make_reference.py` | Regenerates the fixture from the Python package |

## Testing

`npm test` checks all four schemes against the 507 reference domains, plus published
numbering for trastuzumab. Running the same tests on all 9,930 single-domain antibody and TCR
chains in ANARCI's PDB benchmark set (`REFERENCE=... npm test`) gives identical results. The
Python reference itself agrees with ANARCI on every domain of that set.

## Changes in version 2

Version 1 counted alignment columns, which gave wrong numbering for almost all real
antibodies. In particular:

* deleted (`-`) positions shifted every later number;
* inserted residues were dropped, duplicated or left unnumbered;
* Kabat, Chothia and AHo insertions (52A, 82A-C, 100A, 31A ...) could not be produced, because
  they fall on IMGT match states, not on HMM insert columns;
* the CDR gap and insertion rules of each scheme were not applied.

Version 2 maps the alignment onto each scheme with per-region tables and applies each
scheme's CDR rules. It keeps the same function names, inputs and `{ concatSequence, numbering }`
output.
