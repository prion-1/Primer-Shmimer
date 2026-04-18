# Primer Shmimer

A client-side PCR primer analysis tool that evaluates thermodynamic stability and predicts secondary structures. Runs entirely in the browser with no server, no accounts, no data leaving your machine.

The structural prediction engine (hairpins, dimers) is Primer3's `thal.c` compiled to WebAssembly. Melting temperature is computed in JavaScript using the SantaLucia 1998 nearest-neighbor model with Owczarzy 2008 salt corrections. These are the same algorithms used by Primer-BLAST, IDT OligoAnalyzer, and other established primer design tools.

## How to Use

> [!tip]
You can launch the app directly:
<a href="https://prion-1.github.io/Primer-Shmimer/">
  <img src="https://img.shields.io/badge/>-LAUNCH-2ea44f?style=for-the-badge" alt="Launch">
</a>

Alternatively:

Clone repo and open `index.html` in a browser. No server required — the app works from `file://`.

Enter one or two primer sequences (5'→3', characters A/T/C/G only). Adjust salt concentrations and annealing temperature if your PCR conditions differ from the defaults.

The tool reports, per primer: melting temperature, GC content, a composite quality score with transparent penalty breakdown, self-dimer analysis (global and 3' extensible), hairpin analysis, 3' and 5' end stability, GC clamp, and run/repeat detection.

When two primers are provided, the tool additionally reports cross-primer dimers (global and 3' extensible) and the melting temperature difference between the pair.

## What It Computes

### Melting Temperature (Tm)

Computed in JavaScript using the nearest-neighbor (NN) thermodynamic model.

Parameters: SantaLucia 1998 unified set — 10 unique dinucleotide ΔH and ΔS values derived from 108 oligonucleotide duplexes across seven laboratories. Initiation corrections are applied for terminal A·T vs G·C base pairs. A symmetry correction (ΔS −1.4 cal/(mol·K)) is applied to self-complementary (palindromic) sequences, which also changes the concentration term from C/4 to C in the Tm denominator.

Salt correction: the Owczarzy 2008 algorithm is used by default. This implements a ratio-based decision tree:
- When the ratio R = √[Mg²⁺_free] / [monovalent] is below 0.22, monovalent cations dominate and the correction uses monovalent-only terms.
- When R is between 0.22 and 6.0, a transition formula adjusts the Owczarzy coefficients (a, d, g) based on monovalent concentration.
- When R ≥ 6.0 or monovalent = 0, the full Mg²⁺-dependent correction applies.

Free Mg²⁺ is computed as [Mg²⁺] − [dNTP], since dNTPs chelate magnesium 1:1. If dNTP ≥ Mg²⁺, free Mg²⁺ is treated as effectively zero.

The SantaLucia 1998 Na-equivalent formula (monovalent + 120×√Mg) is available as a fallback via the Advanced Settings dropdown.

The Tm equation for non-self-complementary sequences:

```
Tm = ΔH / (ΔS + R × ln(C/4)) − 273.15
```

For self-complementary sequences, `C/4` is replaced by `C`. When the Owczarzy correction is active, the formula operates on 1/Tm with an additive correction term.

Validated against primer3-py `calc_tm` on a 25-sequence panel at matched parameters: MAE 0.026°C, max deviation 0.64°C.

### Hairpin Prediction

Computed by Primer3's `thal()` function compiled to WebAssembly, called with `thal_hairpin` mode.

Primer3 uses a dynamic programming algorithm that considers Watson-Crick matches, all 8 mismatch types (with context-dependent NN parameters from Allawi, SantaLucia, and Peyret 1997–1999), single-nucleotide bulges, internal loops, dangling ends, terminal mismatches at the loop-stem junction, empirical Turner 2004 loop initiation parameters, special triloop and tetraloop sequences with anomalous stability, and all-cytosine loop penalties. The algorithm finds the most stable (most negative ΔG) hairpin structure the sequence can form.

Output: ΔG (kcal/mol) at the annealing temperature, structure Tm (°C), and an ASCII structure diagram showing the stem, loop, and pairing.

Validated against primer3-py `calc_hairpin`: MAE < 10⁻¹⁵ kcal/mol (numerically identical to floating-point precision).

### Dimer Prediction

Computed by Primer3's `thal()` function compiled to WebAssembly.

Two analyses are performed per primer pair:

**Global (ANY):** Finds the most thermodynamically stable dimer structure at any alignment. Called with `thal_any` mode. This predicts the worst-case thermodynamic sink — how much primer is lost to dimerization.

**3' Extensible:** Determines whether the 3' end of the primer is engaged in a dimer that DNA polymerase could extend, producing primer-dimer artifacts. This is biologically more consequential than the global score because even a weak 3'-anchored dimer can produce visible artifacts through exponential amplification.

The 3' extensible analysis works as follows:
1. Primer3's `thal_end1` result is obtained (forces the 3' terminal base of strand 1 to participate in pairing).
2. The returned structure is parsed to verify that the 3' terminal base is part of a contiguous run of ≥2 paired bases (not an isolated pair bridged by a large internal loop).
3. The opposing strand must extend ≥2 unpaired bases beyond the paired region (providing a template overhang for polymerase extension).
4. If the `thal_end1` structure fails these checks (e.g., it contains a biologically impossible isolated terminal pair), the tool checks whether the global (ANY) structure happens to satisfy the criteria. If so, that structure is reported as 3' extensible.
5. If neither structure qualifies, "No 3' extensible dimer found" is reported.

This post-processing produces more biologically meaningful results than raw Primer3 `thal_end1` output, which can report structures where a single base pair is forced at the 3' terminus across a 20-nucleotide gap.

For self-dimers, both copies of the sequence are evaluated against each other. For cross-dimers, primer 1 is evaluated against primer 2.

Validated against primer3-py `calc_homodimer` / `calc_heterodimer`: MAE < 10⁻¹⁵ kcal/mol.

### Quality Score

A penalty-based composite score from 0 to 100. Starts at 100; each metric that falls outside ideal ranges subtracts points. The score is not from Primer3 or any published formula — it is a heuristic based on community-consensus thresholds from IDT, Primer3 defaults, and published PCR optimization guidelines.

| Metric | No penalty | Moderate | Severe |
|---|---|---|---|
| Tm | 57–63°C | 52–57 or 63–68°C: −5 | <52 or >68°C: −15 |
| GC% | 40–60% | 30–40 or 60–70%: −5 | <30 or >70%: −15 |
| Hairpin Tm | ≤ Ta − 15°C or no structure Tm | Ta − 15°C < Tm ≤ Ta − 5°C: −10 | > Ta − 5°C: −20 |
| Self-dimer global ΔG | > −5 kcal/mol | −5 to −8: −10 | < −8: −20 |
| Self-dimer 3' extensible ΔG | > −5 kcal/mol | −5 to −8: −15 | < −8: −25 |
| 3' end stability (last 5 bp) | > −9 kcal/mol | −9 to −10: −5 | < −10: −10 |
| GC clamp (G/C in last 5) | 1–3 | 0 or 4: −5/−10 | 5: −10 |
| Runs (≥4 identical bases) | None | Detected: −10 | — |
| Length | 18–25 nt | 15–17 or 26–30: −5 | <15 or >30: −10 |

Hairpin penalties use the structure Tm relative to the user-specified annealing temperature (Ta), matching the structure Tm traffic light: near-Ta hairpins are penalized because they can remain folded during annealing. The 3' extensible dimer carries the heaviest penalty (−25) because it directly produces artifacts. The score card is expandable — clicking "Penalties" shows exactly which metrics deducted points and by how much.

Score interpretation: 80–100 (green) — good primer. 60–79 (yellow) — usable but has issues. Below 60 (red) — redesign recommended.

### Structure Visualization

Dimer structures are displayed as alignments (5'→3' top strand, match symbols, 3'→5' bottom strand). Watson-Crick pairs are shown as `|` (teal), wobble pairs as `:` (yellow), mismatches as `.` (red). Internal loops within the dimer are indicated by lowercase loop bases (looping strand) and bridge characters (`─`, contiguous strand).

Hairpin structures are displayed with the stem oriented horizontally and the loop shown at the fold point with `┐`/`┘` turn characters.

### Structure Tm Traffic Light

Each structure's Tm is colored relative to the user-specified annealing temperature (Ta):
- Green: structure Tm < Ta − 15°C (fully melted at annealing temperature, irrelevant)
- Yellow: structure Tm between Ta − 15°C and Ta − 5°C (may partially form during cooling)
- Red: structure Tm > Ta − 5°C (stable at annealing temperature, will compete with template binding)

### Additional Sequence Checks

- **3' End Stability:** ΔG of the terminal 5 bases at the 3' end. Indicates how strongly the primer's 3' end binds. Too stable (< −9 kcal/mol): increased risk of mispriming at partially complementary sites. Too weak: may not extend efficiently.
- **5' End Stability:** ΔG of the first 5 bases. A stable 5' end (GC clamp) helps anchor the primer to the template.
- **GC Clamp:** Count of G and C bases among the last 5 bases at the 3' end. Ideal: 1–3. A GC clamp stabilizes the 3' end for efficient extension without promoting mispriming.
- **Runs:** Stretches of ≥4 identical bases (e.g., AAAA). Increases risk of polymerase slippage and mispriming.
- **Repeats:** Dinucleotide repeats of ≥4 units (e.g., ATATATAT). Increases mispriming risk.

### Target ΔG

The ΔG of the full primer sequence evaluated at the annealing temperature. This represents the free energy of the primer-template duplex under actual reaction conditions. More negative values indicate stronger template binding.

## Default Parameters

| Parameter | Default | Unit | Notes |
|---|---|---|---|
| Na⁺ | 50 | mM | Monovalent cation concentration |
| Mg²⁺ | 1.5 | mM | Divalent cation concentration |
| dNTP | 0.8 | mM | Total dNTP (4 × 0.2 mM each) |
| Primer concentration | 250 | nM | Total primer in reaction |
| Annealing temperature | 60 | °C | Used for structure ΔG evaluation, Tm traffic lights, and hairpin score penalties |
| Evaluation temperature | 37 | °C | Temperature for Tm-independent ΔG calculations |
| Max hairpin loop size | 30 | nt | Maximum loop length considered in hairpin search |
| Salt correction | Owczarzy 2008 | — | Can be switched to SantaLucia 1998 + Na equivalent |

## Architecture

```
index.html          UI layout, input forms, result templates
style.css           Styling, traffic-light colors, monospace alignment
thermodynamics.js   NN parameters for JS Tm calculation
thal_bridge.js      JS ↔ WASM bridge for Primer3 thal engine
algorithms.js       Tm calculation (JS), scoring, dimer/hairpin wrappers
app.js              DOM rendering, structure visualization, event handling
primer3_wasm/
  dist/
    thal.js         Emscripten-generated loader with inlined WASM
```

The Tm calculation runs in pure JavaScript. Hairpin and dimer predictions are delegated to Primer3's `thal.c` compiled to WebAssembly via Emscripten. The WASM binary is inlined as base64 in `thal.js` (using Emscripten's `SINGLE_FILE=1` option), so the app works from `file://` without a server.

No external dependencies. No frameworks. No build step for the frontend. No data transmitted anywhere.

## Limitations

1. **No template interaction.** The tool analyzes primers in isolation. It does not predict mispriming against a target genome, amplicon viability, or primer specificity. For specificity checking, use Primer-BLAST (NCBI) or similar tools.

2. **No multiplex analysis.** Cross-dimer analysis is limited to one primer pair. For multiplex PCR with 3+ primers, all pairwise interactions should be checked — this tool does not automate that.

3. **Tm accuracy.** The Owczarzy salt correction omits the length-dependent `(1/(2*(len-1))) * (e + f·ln[Mg] + g·(ln[Mg])²)` term. Empirical validation showed that including it degraded agreement with Primer3's compiled `oligotm` from MAE 0.026°C to 2.064°C at standard conditions. The term is therefore intentionally excluded. This may affect accuracy for sequences significantly longer than typical primers (~20–30 nt) or at unusual Mg²⁺ concentrations.

4. **Score is heuristic.** The 0–100 quality score is not from Primer3 or any peer-reviewed publication. It is a weighted penalty system based on community-standard thresholds. Different PCR applications (genotyping, cloning, qPCR, sequencing) may warrant different thresholds. The penalty weights reflect general-purpose PCR.

5. **Sequence length.** Primer3's NN two-state model is considered reliable for sequences up to 60 nt. The tool does not enforce a hard length limit, but thermodynamic predictions for very long sequences should be treated with caution.

6. **Input validation.** Only the characters A, T, C, G (auto-uppercased) are accepted. Ambiguous bases (N, R, Y, etc.), and modified bases are not supported.

## Computational Methods — References

- SantaLucia J. "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95:1460–1465 (1998).
- Owczarzy R et al. "Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations." Biochemistry 47:5336–5353 (2008).
- Allawi HT, SantaLucia J. "Thermodynamics and NMR of internal G·T mismatches in DNA." Biochemistry 36:10581–10594 (1997).
- Allawi HT, SantaLucia J. "Nearest neighbor thermodynamic parameters for internal G·A mismatches in DNA." Biochemistry 37:2170–2179 (1998).
- Allawi HT, SantaLucia J. "Nearest-neighbor thermodynamics of internal A·C mismatches in DNA." Biochemistry 37:9435–9444 (1998).
- Allawi HT, SantaLucia J. "Thermodynamics of internal C·T mismatches in DNA." Nucleic Acids Res 26:2694–2701 (1998).
- Peyret N et al. "Nearest-neighbor thermodynamics and NMR of DNA sequences with internal A·A, C·C, G·G, and T·T mismatches." Biochemistry 38:3468–3477 (1999).
- Bommarito S et al. "Thermodynamic parameters for DNA sequences with dangling ends." Nucleic Acids Res 28:1929–1934 (2000).
- SantaLucia J, Hicks D. "The thermodynamics of DNA structural motifs." Annu Rev Biophys Biomol Struct 33:415–440 (2004).
- Untergasser A et al. "Primer3 — new capabilities and interfaces." Nucleic Acids Res 40:e115 (2012).
- Rentzeperis D et al. "Thermodynamics of DNA hairpins: contribution of loop size to hairpin stability and ethidium binding." Nucleic Acids Res 21:2683-2689 (1993).

## License

This application uses Primer3's thermodynamic alignment engine (`thal.c`), compiled to WebAssembly. Primer3 is Copyright (c) 1996–2024 Whitehead Institute for Biomedical Research, Steve Rozen, Andreas Untergasser, and others. Released under the GNU General Public License, version 2 (GPL-2.0). Source: https://github.com/primer3-org/primer3
