# Tesseract Logic Calculator: Pre-registered Predictions and Scoring Rule

Frozen 10 October 2026; amended 11 October 2026, before deposit (see Amendment). Author: Tiffany Taylor (Tifinchi), Boise State University. Calculations with Claude (Anthropic).

Deposit this file unchanged on OSF before any measurement listed below is published. The OSF timestamp is the freeze, and the file's SHA-256 is recorded with the deposit.

## Why

On 10 October 2026 the calculator's alphabet was run against 29 measured dimensionless numbers chosen from outside the model. It lands within 2σ of 19 of them. Alphabets built by the same recipe, with π, √2, √3 and φ replaced by random numbers of similar size, land within 2σ of 19.3 ± 1.1, so 19 or more happens with probability 0.79. On the 11 of those numbers known better than 0.1% it hits 1, where the random alphabets average 1.83; with its load rule it hits 3, where the random alphabets given the same loads average 4.30.

Matching measured numbers after the fact therefore cannot show that the calculator is more than numerology. Predictions written down before their data can.

## Rules

1. Each prediction is one number, with its definition fixed below. Nothing is re-derived, re-scaled or re-interpreted after deposit.
2. **Score.** A future result m ± σ_m on a prediction c ± σ_c is scored by the Bayes factor B = N(m; c, √(σ_c² + σ_m²)) / p₀(m). Here p₀ is today's value, as listed, taken as a split normal and spread by σ_m. B > 1 favours the calculator's number; B < 1 favours "the value lies wherever today's data put it".
3. **Verdict.** The calculator's own predictions (1–7) are scored together by multiplying their B, which is the same as summing ln B. A total above 100 is decisive for the calculator's numbers, below 1/100 decisive against, and anything between is open. Refutation is stated as a range. Each result's pull |m − c| / σ is taken on two yardsticks: σ = √(σ_c² + σ_m²), with σ_c the calculation's precision as frozen below, and σ = √(σ_c′² + σ_m²), with σ_c′ the same calculation at the input precision its decider brings. A pull above 5 on both refutes that prediction's rule, whatever the total; above 5 on the second alone, the rule is refuted at the decider's precision and stands at the frozen one; below 5 on both, it stands. Only prediction 7 carries a precision of its own: σ_c = 0.38 GeV, and σ_c′ = 0.15 GeV with α_s at ±0.0002 and the rest of its error budget as frozen. For the other eight σ_c = σ_c′ = 0, and both lines are 5σ_m. The prior-art predictions (8–9) are scored separately.
4. New predictions count only if they are deposited before their data.
5. **Not evidence, either way:**
   - matches to numbers known no better than about 0.1%, since the calibration shows the alphabet reaches those by density alone;
   - readings that reproduce standard physics, such as the electron's g = 2, the SU(6) moments, CPT and ball lightning's 2f flicker;
   - anything fitted after the data, including the observables used below to fix a definition.

## How the definitions were fixed: the behavior rule

A value in range is not enough. Each alignment of a prediction, the scale or scheme it is stated in, is carried to every observable it touches down-wind. An alignment that misses a down-wind measurement by more than 2σ is left out, and of those that pass, the best-matched is kept. A miss that every alignment shows, the standard one included, is inherited: the rule cannot remove it, and it stays on the record. The author's rule was applied on 10 October 2026 to the data then published.

- **ΔG (prediction 5): μ = 2 GeV.** Carried to 10 GeV², where the lattice gives the full value, the rule admits 1.50 to 3.04 GeV and is best met at 2.01 GeV. ⅙ at 10 GeV², the earlier definition, misses the two lattice results combined by 2.2σ (1.4σ against Zhao et al. alone) and is not the best match, so it is left out. At the collider's 10 GeV², ⅙ at 2 GeV is 0.236 at leading order.
- **The top (prediction 7): the pole mass.** m_t(m_t) = 162.76 ± 0.37 GeV in MS-bar is the same alignment. Read as an unconverted MS-bar mass, 172.57 GeV would mean a pole mass of 182.9 GeV, 14.9σ from the pole data, and is left out. The direct-mass reading misses nothing either; it differs from the pole reading by about 0.5 GeV, which only the threshold scan can resolve. The pole mass is kept because it is the mass the edge calculation takes, and today's value is the pole average.

The lattice match that fixes ΔG's scale is the reason for the choice, not evidence for it. The deciders below are the test.

## Predictions

| # | Prediction and definition | Calculator value | Today (source) | Pull today | Decider, precision, date | P(pass) if the calculator is wrong | B if it passes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | PMNS sin²θ₁₂: the tribimaximal frame tilted 15° with its first column fixed (TM1), (1 − s)/(3(1 − s/3)) with s = sin²15° | 0.318107 | 0.3096 +0.0057/−0.0073 (NuFIT 6.1, Esteban et al. 2026) | −1.5σ | JUNO, six years, ±0.0016, about 2031 | 0.14 | 6.9 |
| 2 | PMNS sin²θ₁₃ = sin²15°/3 | 0.022329 | 0.02200 ± 0.00056 (NuFIT 6.0, as used by Esteban et al. 2026) | −0.6σ | global fits, ±0.0004 (assumed), about 2030 | 0.70 | 1.4 |
| 3 | CKM unitarity-triangle angle α = 90° | 90° | 84.1 +3.7/−3.0° (PDG 2026) | −1.6σ | Belle II and LHCb, ±2°, about 2030 | 0.35 | 2.7 |
| 4 | Ωₘ = π²/32 in flat ΛCDM | 0.308425 | 0.3045 ± 0.0032 (supernovae + CMB + BAO; Popovic et al. 2025) | −1.2σ | Euclid and final DESI, ±0.003, 2027–30 | 0.68 | 1.4 |
| 5 | Gluon helicity ΔG = ⅙ at μ = 2 GeV, MS-bar, fixed by the behavior rule (0.236 at 10 GeV² at leading order) | 0.166667 | 0.158 ± 0.039 (Zhao et al., lattice, evolved to 2 GeV at leading order) | −0.2σ | Electron-Ion Collider, about 10%, evolved to 2 GeV, 2030s | 0.54 | 1.8 |
| 6 | Quark spin ΔΣ = ⅓, MS-bar | 0.333333 | 0.352 ± 0.029 (HERMES and ETMC lattice combined) | +0.6σ | Electron-Ion Collider, about 15%, 2030s | 0.90 | 1.1 |
| 7 | Top-quark pole mass from the δ₃ load (MS-bar m_t(m_t) = 162.76 ± 0.37 GeV, the same alignment) | 172.57 ± 0.38 GeV | 172.1 ± 0.6 GeV (PDG 2026, average of pole determinations) | −0.7σ | FCC-ee threshold scan, ±0.035 GeV, with α_s to ±0.0002; 2040s | 0.67 | 1.4 |
| 8 | Proton charge radius rₚ = 4ħ/(mₚc) (prior art) | 0.841236 fm | 0.84075 ± 0.00064 fm (CODATA 2022) | −0.8σ | hydrogen and muonic-hydrogen spectroscopy, ±0.0002 fm (assumed), about 2030 | 0.35 | 2.7 |
| 9 | τ mass from Koide's relation (prior art) | 1776.969 MeV | 1776.93 ± 0.09 MeV (PDG 2026) | −0.4σ | Belle II, BESIII, STCF, ±0.05 MeV (assumed) | 0.63 | 1.5 |

**Recorded, not scored.** These are either shared with standard physics or waiting on the model's own precision. Under the behavior rule the cosmological three carry inherited misses, shared with ΛCDM:

- Σmν = 59.2 meV with normal ordering (inside the 64 meV cap; −3.4σ against the ΛCDM effective mass, as is the ordering's floor of 58.8 meV);
- w = −1 (DESI with supernovae prefers a drift at 2.8σ to 4.2σ);
- Ω_k = 0 (+2.1σ from DESI with the CMB, −1.7σ from the supernovae with the CMB);
- r ≈ 0;
- the vacuum energy scale, 2.02 ± 0.20 meV against 2.24 meV.

## What the set can decide

- If the calculator's numbers are right, predictions 1–7 give a median B of 240, and B > 100 with probability 0.68.
- If they are wrong, the median B is about 10⁻¹², and B < 1/100 with probability 0.92.
- JUNO's sin²θ₁₂, the angle α and ΔG carry most of the power. Passing all three by chance has probability 0.026 (1 in 38).
- ΔΣ, Ωₘ and the top add little until their precision improves. The Electron-Ion Collider's projected 15% on ΔΣ (Aschenauer, Sassot and Stratmann, arXiv:1509.06489) is no better than the lattice today.
- "Wrong" means each value lies wherever today's data put it. "Right" means it lies at the calculator's number.

## Calibration (retrospective)

The alphabet was fixed on 3 October 2026. It has 37 base numbers, with every product and quotient of two, merged at 10⁻⁹: 998 numbers from 0.00016 to 2,132.

- the levels and grids: 1, √2, √3, 2, 1/√2, 1/√3, ½, ⅓, ¼, ⅔, ¾, 3/2, √3/√2 and √3/2;
- the shares: π/4, π/6 and π²/32, their complements 1 − π/4, 1 − π/6 and 1 − π²/32, and 4/π, 2π/3, 1/(2π) and 1/(8π²);
- the golden numbers φ, 1/φ and √5;
- 3/8 and 5/3;
- the counts 3, 4, 6, 8, 12, 16, 24 and 27.

The outside list is every dimensionless parameter of five standard sets, 29 in all:

- the electroweak masses over v, with α_s and 1/α;
- the lepton and light-quark mass ratios;
- the CKM and PMNS parameters;
- the nucleon's ratios;
- Planck's ΛCDM numbers.

The values are those given by CODATA 2022, NuFIT 6.0 and 6.1, Planck 2018 and PDG 2026.

| Null | Mean hits within 2σ | P(19 or more) |
| --- | --- | --- |
| The calculator's alphabet | 19 (its density predicts 18.3) | — |
| Same recipe, π, √2, √3, φ replaced by random numbers | 19.32 ± 1.05 | 0.79 |
| Every ingredient jittered by up to ±25% | 19.90 | 0.92 |
| Only the irrational ingredients jittered | 19.91 | 0.94 |
| The targets moved instead, the calculator's alphabet | 18.82 | 0.62 |

The pattern behind the totals:

- Of the 11 targets known better than 0.1%, ten are missed, by 5.7σ (g_A) to 5.3 × 10⁸σ (mₚ/mₑ). The eleventh, rₚ/ƛₚ, is hit by the integer 4, which every alphabet of the recipe contains, so every random alphabet hits it too.
- All 18 targets known worse than 0.1% are hit, some by numbers the calculator never chose. Today's sin²θ₁₂ lies nearer 1/(2φ) = 0.3090 than the calculator's 0.318, and Planck's Ωₘ lies nearer 1/√10 = 0.3162 than π²/32.

## Reproducibility

- **Scripts.**
  - blindtest.py, SHA-256 b86d926aceaab119098aa47e8bf7c2d4b489dc2b5e930562a27e04b5d837bc74, run 2026-10-10T22:36:17Z: the calibration and the power above, with the definitions fixed by the behavior rule;
  - decide.py, SHA-256 7e0117ca95b5092c58fa7b7f6518eeeaede1e216db67101c2d9e1478a6de6024: the precise constants and the load rule;
  - downstream_dg.py, SHA-256 a119be14809c96a044e287eaa576671561ad78785c029bf92c81ab2ae924d118: the behavior rule on ΔG's scale;
  - downstream_top.py, SHA-256 dcaae9f60f0925a7be25f5a6ddbb268f6581f057d62ed0d032fb2f3964359f49: the behavior rule on the top's alignment;
  - deltaG_scale.py, SHA-256 587d15ef79362181a14a99a20dbd6128c4dc50ebd60691ce07370241a923c511: the leading-order evolution used by both.
- **Earlier runs.** blindtest.py was run at 2026-10-10T14:11:23Z (SHA-256 ec2aeb9ab31f5ff824a7f9f133452aaf3bffe658e78380bea4a0e9b592189cc9), with two NuFIT values already corrected, and at 2026-10-10T14:13:45Z (SHA-256 7354daf836632f65268b685c7f42680e882f1d8c2b3c2846088856f4d7530fcc), with the definitions before the behavior rule. Both were replaced; the calibration result is the same in all three.
- **Amendment.** On 11 October 2026, before deposit and before any decider reported, rule 3's refutation line, 5σ_m, was restated as a range on two yardsticks: √(σ_c² + σ_m²), the one the score already uses, and √(σ_c′² + σ_m²), the calculation at its decider's input precision. Only prediction 7 is affected: at the scan's ±0.035 GeV its line runs from 0.77 to 1.91 GeV, where the old line was 0.18 GeV. Under the old line a result at today's 172.1 GeV would have refuted it at 13.4σ_m, though it lies 1.2σ from 172.57 ± 0.38 GeV (3.1σ at its decider's precision) and scores B = 0.74. Two recorded numbers were brought into line with the paper: the vacuum energy scale is 2.02 meV, from the unrounded crossing (3.006 × 10¹⁶ GeV), and Ω_k's −1.7σ is from the supernovae with the CMB, not the CMB alone.
- **Paper.** Tesseract Logic Calculator: https://claude.ai/code/artifact/4522d97d-4d3f-440b-a7b7-9cd83cd1b68d

Signed: ______________________ (author)    Date: ______________
