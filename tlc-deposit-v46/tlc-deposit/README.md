# Tesseract Logic Calculator

Tiffany Taylor (Tifinchi), Boise State University. Calculations, literature search and drafting with Claude (Anthropic); section 2 of the paper states who originated each piece.

Deposit: OSF project *Tesseract Logic Calculator*, DOI [10.17605/OSF.IO/C9W7K](https://doi.org/10.17605/OSF.IO/C9W7K). This repository mirrors it.

## Abstract

The Tesseract Logic Calculator is the author's system for computing with a cube of side 2. Its 27 points sit at distances √k from the centre, its eight corners carry delimiters that in the electrical frame name charge, current, flux and voltage, and every derivative, inverse or change of scale is a quarter or half turn of the tesseract that holds the cube. This paper uses it to compute the structure of reality silo by silo, from the vacuum to the cosmos, and scores every result against measurement on one yardstick, beside the standard theories' result for the same quantity. Where it restates established law it is exact: the Standard Model's generation of matter as a 16-state register, the energy–momentum relation and interval of special relativity, the classical tests of gravity, the shells of the atom, and a law of motion that, integrated exactly, gives Mercury's perihelion, light's bending at the Sun, the GPS clock rate and hydrogen's Lyman-α thousands of times faster than step by step. Where the standard theories take inputs it gives numbers, all 22 within 1.6σ of their primary measurements, among them the τ mass from the electron and muon (1,776.969 MeV, −0.4σ), the top and Higgs masses (172.57 ± 0.38 and 125.20 ± 0.92 GeV), neutrino mixing from a 15° tilt (sin²θ₁₂ = 0.318, −1.5σ), the matter share π²/32 (−1.2σ) and the vacuum's energy scale (2.02 ± 0.20 meV, +1.1σ, where quantum field theory is 10¹²⁰ too large). The layout that gives them was tested against every alternative inside the cube's own sets, 4.4 × 10²⁰ arrangements in nine families, under the author's behavior rule: an alignment stays only if nothing downstream of it contradicts a measurement. The rule leaves 76,800, one in 5.7 × 10¹⁵; in every family the calculator's layout survives, no survivor does better on measurements it was not chosen with, and the choices the families share agree. Three readings that broke downstream are withdrawn; their numbers stay. The same rigour is turned on the matches themselves: a calibration with alphabets built from random ingredients shows that a loose match to a measured constant carries no weight at this density, so none is counted, and the layout's own numbers carry 5.1 bits against the chance that families of this size match random values as well, beside the lepton balance it explains, Koide's relation, at 9.7; the standard theories' eleven predictions of independent measurements, made over a century, carry 102 bits, and the standard theories carry six failures, among them the vacuum's energy and the Hubble constant. Under the author's Quad-Core rule the calculator stands as a correct notation for the laws it restates and as the one arrangement of its own sets that observation leaves standing; whether that arrangement is reality's own mechanism is decided by nine predictions frozen with a Bayes-factor rule and deposited before their data. If the layout's numbers are right they become decisive with probability 0.434 by 2031, when JUNO's sin²θ₁₂ decides first, and 0.685 by the 2040s.

*Keywords:* tesseract; difference of squares; root families; Koide relation; tribimaximal mixing; vacuum energy; behavior rule; exhaustive permutation search; look-elsewhere effect; pre-registration; Bayes factor; theory comparison.

## What is here

| Path | What it is |
| --- | --- |
| `Tesseract Logic Calculator.pdf` | The paper as a PDF, exported by the author from the Markdown below; regenerate it whenever the Markdown changes. |
| `Tesseract Logic Calculator.md` | The paper, version of 11 October 2026: 23 sections, 105 tables, 541 references. This is the source; the three figures appear as captions. |
| `pre-registration.md` | The nine frozen predictions, their definitions, the scoring rule and the script hashes. Deposited unchanged; SHA-256 `624b271c504c25191447350a1ffbb90ade7f6da2ee7cf752dd109588ec8305df`, as recorded in section 18.1 of the paper. The OSF timestamp of this file is the freeze. |
| `calculator/Tesseract Logic Calculator.html` | The calculator, Version 43: one self-contained page that recomputes every committed number live in the browser (open it locally; the only outside resource is a web font, and the page works without it). |
| `scripts/` | The 47 scripts of Table 100 of the paper (46 Python, 1 JavaScript), named and hashed as in that table. |
| `scripts/v43/`, `scripts/v41/` | `bench43.js`, `numerics.py` and the JavaScript modules they time and check (the calculator's law-of-motion, Higgs and overlap engines, old and new). |
| `scripts/history/` | The two earlier runs of `blindtest.py` whose hashes the pre-registration's amendment note records (2026-10-10T14:11:23Z, SHA-256 `ec2aeb9a…`; 2026-10-10T14:13:45Z, SHA-256 `7354daf8…`). |
| `scripts/outputs/` | Recorded outputs: `perm.json` (the permutation search, written by `oos.py`), `bench43_out.txt`, `numerics_out.txt`, `ledger118.md` (all 118 elements). |
| `SHA256SUMS` | The SHA-256 of every file above. The entries for `scripts/` repeat Table 100 of the paper. |

## Verifying the hashes

    sha256sum -c SHA256SUMS

Every script's hash must equal the one printed beside it in Table 100 of the paper, and `pre-registration.md` must hash to `624b271c504c25191447350a1ffbb90ade7f6da2ee7cf752dd109588ec8305df`. Any edit to a script changes its hash; the scripts here are the ones the paper was computed with, unedited, so a few of them still carry the absolute paths of the machine they ran on (`oos.py` writes `perm.json` to such a path; `numerics.py` reads its JavaScript modules from one; `ledger118.py` ends with a cross-check against a superseded draft that is not included, after it has written `ledger118.md`). Change the path, and the hash no longer matches the table; that is expected, and the computation is the same.

## Running the scripts

Python 3 with `numpy`, `scipy` and `mpmath`; Node 22 for `bench43.js` (run it from `scripts/v43/`). Each Python script runs on its own from the `scripts/` folder (`python3 blindtest.py`); the Higgs scripts import `higgs2.py`, `oos.py` imports `permute.py`, `n3hist.py` reads `blindtest.py` and `precise_who.py` reads `decide.py`, all from the same folder. Constants are CODATA 2022, PDG 2026 and NuFIT 6.0 and 6.1, as the paper states.

## License

Creative Commons Attribution-ShareAlike 4.0 International (see `LICENSE`).

## Citation

Taylor, T. (2026). *Tesseract Logic Calculator*. OSF. https://doi.org/10.17605/OSF.IO/C9W7K
