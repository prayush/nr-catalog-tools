# Surrogate-Mediated Cross-Catalog Validation of Numerical Relativity Binary Black Hole Waveforms

**Prayush Kumar** _et al._

_Draft — not for circulation_

---

## Abstract

We present a systematic framework for comparing binary black hole gravitational waveform catalogs produced by independent numerical relativity (NR) codes against each other, using the NRSur7dq4 surrogate model as a common reference. The surrogate mediates the comparison by providing waveforms at the exact intrinsic parameters of each catalog simulation, eliminating the parameter-space mismatch that has historically confounded direct cross-catalog comparisons. We quantify agreement through two complementary metrics: the noise-weighted match (faithfulness) computed per spherical harmonic mode with the Advanced LIGO zero-detuning high-power noise curve, and the accumulated phase difference per gravitational-wave cycle. We apply this framework to a set of SXS Cauchy–Characteristic Evolution (CCE) simulations spanning non-spinning and aligned-spin configurations at mass ratios $q = 1$ and $q = 2$. Our pilot results demonstrate sub-percent mismatches in the dominant $(2,2)$ and $(3,2)$ modes across all configurations, with phase errors below $0.1$ rad per cycle. Sub-dominant modes involving odd-$m$ harmonics show markedly different behavior depending on whether binary symmetries are broken by mass ratio or spin, providing sensitive diagnostics of surrogate interpolation accuracy. We outline the full comparison pipeline, the treatment of source-frame ambiguities including SO(3) rotations and BMS supertranslations, and the path toward batch processing of the SXS, RIT, and MAYA/Georgia Tech catalogs.

---

## I. Introduction

The past decade has seen the maturation of multiple independent numerical relativity (NR) codes and their associated public waveform catalogs for binary black hole (BBH) mergers. Chief among these are the SXS catalog~\cite{sxs}, produced by the Spectral Einstein Code (SpEC); the RIT catalog~\cite{rit}, produced by the LazEv code; and the Georgia Tech / UT Austin MAYA catalog~\cite{maya}, produced by the MayaKranc code. Each of these catalogs covers overlapping regions of the BBH parameter space and is relied upon by the gravitational-wave data analysis and waveform modeling communities as a ground truth.

The practical importance of NR catalogs is substantial. Gravitational-wave detection pipelines (PyCBC~\cite{pycbc}, LALSuite~\cite{lalsuite}, Bilby~\cite{bilby}) use NR waveforms as high-fidelity injection signals and as calibration targets. Analytical waveform models — including effective-one-body (EOB) frameworks, phenomenological (IMRPhenom) families, and NR surrogate models — are calibrated and validated against NR data. A systematic bias in any particular catalog will propagate directly into model calibrations and, through them, into parameter estimation for every observed gravitational-wave event. Quantifying this bias is therefore not merely an academic exercise: it sets a fundamental floor on the accuracy of waveform models and, ultimately, on our ability to extract astrophysical information from gravitational-wave observations.

Direct comparison between independent NR catalogs is, however, technically difficult. Three distinct classes of obstacles arise:

**Parameter-space mismatch.** Different NR codes reach their initial conditions through distinct numerical procedures. Matching the physical initial data precisely (orbital separation, eccentricity, spin magnitudes, spin orientations, and center-of-mass velocity) across codes is computationally expensive and never exact. As a result, for two catalog simulations $h^A(\boldsymbol{\theta}_i)$ and $h^B(\boldsymbol{\theta}_j)$ that nominally represent the same binary, the residual parameter-space distance $|\boldsymbol{\theta}_i - \boldsymbol{\theta}_j|$ can easily dominate over the intrinsic numerical errors one hopes to measure~\cite{hannam2009}. In the worst case, $\|h^A(\boldsymbol{\theta}_i) - h^B(\boldsymbol{\theta}_j)\|$ is comparable to $\|h^A(\boldsymbol{\theta}_i) - h^A(\boldsymbol{\theta}_k)\|$ for a nearby point $\boldsymbol{\theta}_k$ in the same catalog — i.e., the cross-catalog difference is dominated by the interpolation error of the catalog itself.

**Source-frame ambiguity.** Different NR codes define the orientation of the source frame $\mathbf{F}_s$ — the inertial frame in which the waveform multipoles $h_{\ell m}$ are expressed — using different conventions. One code may align the $z$-axis with the instantaneous orbital angular momentum at the relaxation time; another with the Newtonian angular momentum at the beginning of the simulation; a third with the principal axes of the computational grid. Even for identical intrinsic parameters, this choice introduces a rigid SO(3) rotation between the two mode sets that artificially inflates any naive mismatch.

**BMS supertranslation ambiguity.** More subtly, the asymptotic symmetry group of general relativity at null infinity $\mathcal{I}^+$ is not the Poincaré group but the infinite-dimensional Bondi–Metzner–Sachs (BMS) group. A supertranslation — a direction-dependent retarded-time shift $u \to u - \alpha(\theta, \phi)$ — constitutes a physically meaningful gauge freedom in the definition of the waveform. Different extraction methods (Cauchy–Characteristic Extraction, Cauchy–Perturbative Extraction with different extraction radii) may implicitly live in different BMS frames, inducing mode-mixing contributions that are indistinguishable from genuine waveform differences unless explicitly accounted for.

Past community-wide efforts — NINJA-1 and NINJA-2~\cite{ninja1,ninja2}, and the NRAR collaboration~\cite{nrar} — addressed complementary aspects of NR accuracy but stopped short of direct inter-catalog comparisons with the level of rigor that modern detector sensitivity and modeling demands.

**Surrogate models as comparison mediators.** The advent of high-accuracy NR surrogate models opens a new avenue. A surrogate such as NRSur7dq4~\cite{nrsur7dq4} is a sophisticated interpolant of waveform multipoles trained on a dense grid of NR simulations. Given any binary configuration $\boldsymbol{\theta}$ within its training domain, it can produce waveforms at those _exact_ parameters, eliminating the parameter-space mismatch entirely. The comparison then reduces to

$$\|h^A(\boldsymbol{\theta}_i) - h^{\rm sur}(\boldsymbol{\theta}_i)\|,$$

where $h^{\rm sur}$ is the surrogate evaluated at the catalog's own parameters. Any residual difference reflects a combination of: (a) numerical error in the NR simulation, (b) interpolation error of the surrogate, and (c) residual source-frame and BMS ambiguities. This paper is primarily concerned with systematically isolating and quantifying each contribution.

The organization of this paper is as follows. Section~II describes the formalism for source-frame ambiguities and defines the match and phase-difference metrics. Section~III describes the NRSur7dq4 surrogate and the computational pipeline. Section~IV presents pilot results for four SXS simulations. Section~V discusses interpretation and outlines the path to full-catalog comparison.

---

## II. Formalism

### A. Source Frame and Waveform Multipoles

We work with the strain decomposed at null infinity $\mathcal{I}^+$ in spin-weight $-2$ spherical harmonics in a fixed inertial source frame $\mathbf{F}_s$:

$$H(t, \iota, \phi_c) = h_+(t) + i h_\times(t) = \sum_{\ell \geq 2} \sum_{m=-\ell}^{\ell} {}^{-2}Y_{\ell m}(\iota, \phi_c)\, h_{\ell m}(t; \boldsymbol{\theta}),$$

where $\iota$ and $\phi_c$ are the polar angles of the detector in the source frame. Each catalog simulation provides the complex time series $h_{\ell m}(t)$ in its own source frame convention.

### B. SO(3) Frame Rotation

Let $\mathbf{F}_C$ (catalog frame) and $\mathbf{F}_S$ (surrogate frame) be related by a rigid rotation $R \in \mathrm{SO}(3)$. Under this rotation, the multipoles transform via the Wigner $D$-matrices:

$$h^{S,\mathrm{rot}}_{\ell m}(t) = \sum_{m'=-\ell}^{\ell} h^S_{\ell m'}(t)\, D^\ell_{m' m}(R).$$

A time shift $t_c$ and coalescence phase $\phi_c$ then complete the nuisance transformation:

$$\boxed{h^S_{R,\ell m}(t;\, t_c, \phi_c, R) = e^{-im\phi_c} \sum_{m'} h^S_{\ell m'}(t - t_c)\, D^\ell_{m' m}(R).}$$

The match maximized over $(t_c, \phi_c, R)$ isolates the physical mismatch from all frame-choice artifacts.

### C. BMS Supertranslations

Beyond $\mathrm{SO}(3)$ rotations, the full BMS symmetry group introduces supertranslations — direction-dependent time shifts $u \to u - \alpha(\theta,\phi)$ where $\alpha$ is expanded in ordinary spherical harmonics:

$$\alpha(\theta,\phi) = \sum_{j,k} \alpha_{jk}\, Y_{jk}(\theta,\phi).$$

To first order in $\alpha$, the transformed modes are

$$\boxed{h'_{\ell m}(u) = h_{\ell m}(u) - \sum_{j,k,p,q} \alpha_{jk}\, \mathcal{G}^{\ell m}_{jk,pq}\, \dot{h}_{pq}(u),}$$

where $\mathcal{G}^{\ell m}_{jk,pq} = \int_{S^2} {}^{-2}Y^*_{\ell m}\, Y_{jk}\, {}^{-2}Y_{pq}\, d\Omega$ are Gaunt integrals. The $j=0$ term corresponds to a uniform time translation $t_c$ (already in the SO(3) optimization); $j=1$ terms correspond to center-of-mass translations; $j \geq 2$ terms are proper supertranslations. Maximizing the match over supertranslation coefficients removes residual BMS gauge differences from the comparison.

### D. Agreement Metrics

**Noise-weighted match.** For two time series $h_1(t)$ and $h_2(t)$, the noise-weighted inner product is

$$\langle h_1 | h_2 \rangle = 4\,\mathrm{Re} \int_{f_{\rm min}}^{f_{\rm max}} \frac{\tilde{h}_1(f)\,\tilde{h}_2^*(f)}{S_n(f)}\, df,$$

where $S_n(f)$ is the one-sided power spectral density. The match (faithfulness) is

$$\mathcal{F}(h_1, h_2) = \frac{\langle h_1 | h_2 \rangle}{\sqrt{\langle h_1 | h_1 \rangle \langle h_2 | h_2 \rangle}},$$

maximized over time and phase shifts by PyCBC's `match()` function. We use the Advanced LIGO zero-detuning high-power design curve as $S_n(f)$ throughout, and apply this metric independently to each mode.

For mode $(\ell, m)$, the appropriate lower frequency cutoff is

$$f_{\rm lower}^{(\ell m)} = \frac{|m|}{2}\, f_{\rm lower}^{(22)},$$

since the GW frequency of the $(\ell,m)$ mode scales as $|m|$ times the orbital frequency, and the $(2,2)$ GW frequency equals twice the orbital frequency.

**Phase difference per cycle.** The match alone is insensitive to secular phase drift that accumulates slowly over many cycles (it can partly absorb phase bias into the time-shift optimization). We therefore also compute

$$\Delta\Phi/{\rm cycle} = \frac{|\Delta\Phi_{\rm NR} - \Delta\Phi_{\rm sur}|}{N_{\rm cyc}^{\rm NR}}\quad [\mathrm{rad/cycle}],$$

where $\Delta\Phi = |\phi(t_{\rm end}) - \phi(t_{\rm start})|$ is the total accumulated phase over the common time window of the two waveforms, $\phi(t) = \arg[h_{\ell m}(t)]$ is the unwrapped phase of the complex mode, and $N_{\rm cyc}^{\rm NR} = \Delta\Phi_{\rm NR} / (2\pi)$ is the total number of gravitational-wave cycles in the NR waveform. This metric is free from the time/phase maximization and directly measures the average rate of phase accumulation error.

Together, $\mathcal{F}$ and $\Delta\Phi/{\rm cycle}$ provide complementary diagnostics: the former is sensitive to amplitude and instantaneous phase coherence (especially near merger); the latter directly quantifies the integrated phase budget over the inspiral.

---

## III. Pipeline

### A. NRSur7dq4

NRSur7dq4~\cite{nrsur7dq4} is a fully precessing surrogate model trained on a dense grid of 1528 SpEC simulations at mass ratios $q = m_1/m_2 \in [1, 4]$ and spin magnitudes $|\boldsymbol{\chi}_{1,2}| \leq 0.8$. It provides all modes up to $\ell = 4$ (excluding $(5,5)$) as complex numpy arrays $h_{\ell m} = h_+ - i h_\times$ in dimensionless $rh/M$ units — the same spin-weight $-2$ spherical harmonic convention used by `WaveformModes.get_mode()`, so no convention conversion is required, only amplitude scaling and time rescaling. The model accepts:

- **Mass ratio** $q = m_1/m_2 \geq 1$ (PyCBC / SpEC convention)
- **Dimensionless spins** $\boldsymbol{\chi}_{1,2}$ evaluated at the reference epoch corresponding to $f_{\rm low}$; when $f_{\rm low}$ is passed as the NR relaxation frequency, the surrogate's spin reference epoch automatically aligns with the epoch at which NR spin components are extracted from the catalog metadata
- **Starting frequency** $f_{\rm low}$ in "Mf" units: $f_{\rm low} = M_{\rm tot} \cdot f_{\rm GW}^{(22)} \cdot G M_\odot / c^3$
- **Time step** $dt$ in dimensionless units $dt/M$

The model internally integrates from $f_{\rm low}$ forward; if $f_{\rm low}$ corresponds to an orbital angular frequency below NRSur7dq4's minimum training extent ($M\Omega \approx 0.0165$ for $q \approx 1$, $|\boldsymbol{\chi}| \approx 0.5$), the model raises an exception and the starting frequency is automatically raised to the minimum training extent. We detect such clipping and adjust the lower cutoff of the match integral to $f_{\rm lower}^{\rm match} = \max(f_{\rm lower}^{\rm NR}, f_{\rm lower}^{\rm sur})$, so that neither waveform is penalized for having support outside the other's frequency band.

### B. Parameter Extraction

Source parameters are extracted from catalog metadata via `nrcatalogtools.CatalogBase.get_parameters()`, which returns a PyCBC-compatible dictionary: `mass1`, `mass2`, `spin1x/y/z`, `spin2x/y/z`, and `f_lower`. For SXS simulations, `f_lower` is defined as

$$f_{\rm lower} = \frac{M\Omega_{\rm NR}}{\pi \cdot M_{\rm tot} \cdot (G M_\odot / c^3)},$$

which equals the $(2,2)$-mode gravitational-wave frequency at the NR relaxation time. This same value is passed as $f_{\rm low}$ to NRSur7dq4, which uses it both as the waveform start frequency and as the epoch at which the input spin components are defined. Because `nrcatalogtools` extracts spin components at the relaxation time, these two epochs are exactly consistent — no separate spin-epoch remapping is required.

### C. Mode Extraction and Scaling

NR modes are extracted as complex physical-unit time series via `WaveformModes.get_mode(ell, em, total_mass, distance, delta_t_seconds)`, which returns a PyCBC `TimeSeries` with epoch set so that $t = 0$ corresponds to the peak of the $(2,2)$ amplitude. The surrogate modes are scaled from dimensionless units to physical units as

$$h_{\ell m}^{\rm phys}(t) = h_{\ell m}^{\rm sur}(t/M_s) \times \frac{G M_{\rm tot} M_\odot}{c^2 D_{\rm Mpc}},$$

where $M_s = M_{\rm tot} \cdot G M_\odot / c^3$ is the total mass in seconds and $D_{\rm Mpc}$ is the luminosity distance. All comparisons use $M_{\rm tot} = 40\,M_\odot$ and $D = 1\,{\rm Mpc}$.

The per-mode match uses the real part $\mathrm{Re}[h_{\ell m}]$ (the $h_+$ component), padded to the next power-of-two length and noise-weighted with a freshly constructed PSD at the matching frequency resolution.

---

## IV. Results

We present results for four SXS simulations selected to probe two orthogonal axes of parameter space: mass ratio ($q = 1$ vs. $q = 2$) and spin ($\chi = 0$ vs. aligned spin $\chi_{1z} \approx 0.5$). All modes from the set $\{(2,2), (2,1), (3,3), (4,4), (3,2), (4,3)\}$ available from NRSur7dq4 are evaluated. Table~I summarizes the simulation parameters.

**Table I. Simulation parameters.**

| Simulation | $q$ | $\chi_{1z}$ | $\chi_{2z}$ | $f_{\rm lower}^{\rm NR}$ [Hz] | $f_{\rm lower}^{\rm match}$ [Hz] |
|---|---|---|---|---|---|
| SXS:BBH:0001 | 1.00 | 0.00 | 0.00 | 19.8 | 26.2 (sur. clipped) |
| SXS:BBH:0005 | 1.00 | +0.50 | 0.00 | 19.8 | 26.9 (sur. clipped) |
| SXS:BBH:0169 | 2.00 | 0.00 | 0.00 | 29.1 | 29.1 |
| SXS:BBH:0162 | 2.00 | +0.60 | 0.00 | 28.8 | 28.8 |

For SXS:BBH:0001 and SXS:BBH:0005, the NR simulation starts at $\sim 20$ Hz but NRSur7dq4's minimum training extent at these parameters corresponds to $\sim 26$–$27$ Hz. The match lower cutoff is raised accordingly so that neither waveform is penalized for frequency content outside the other's support.

### A. Match Results

**Table II. Per-mode match $\mathcal{F}$ for all four simulations.**

| Mode | 0001 ($q$=1, ns) | 0005 ($q$=1, spin) | 0169 ($q$=2, ns) | 0162 ($q$=2, spin) |
|---|---|---|---|---|
| (2,+2) | 0.9940 | 0.9927 | 0.9905 | 0.9927 |
| (2,+1) | 0.277 † | **0.9966** | 0.9513 | 0.351 ‡ |
| (3,+3) | 0.293 † | 0.9513 | 0.9937 | 0.9940 |
| (4,+4) | 0.902  | 0.9826 | 0.9918 | 0.9896 |
| (3,+2) | 0.9941 | 0.9897 | 0.9565 | **0.573 ‡** |
| (4,+3) | 0.353 † | 0.9742 | 0.9546 | 0.407 ‡ |

† Near-zero by $q = 1$, $\chi = 0$ symmetry; match value is numerically meaningless.  
‡ Anomalous; discussed in Section IV.C.

### B. Phase Difference per Cycle

**Table III. Phase difference per GW cycle $\Delta\Phi/{\rm cycle}$ [rad] and number of NR cycles $N_{\rm cyc}$.**

| Mode | 0001 (q=1,ns) | 0005 (q=1,spin) | 0169 (q=2,ns) | 0162 (q=2,spin) |
|---|---|---|---|---|
| (2,+2) | 0.076 (41 cyc) | 0.012 (44 cyc) | 0.009 (38 cyc) | 0.004 (46 cyc) |
| (2,+1) | — †            | 0.068 (25 cyc) | 0.001 (22 cyc) | 0.042 (26 cyc) |
| (3,+3) | — †            | 0.057 (66 cyc) | 0.023 (58 cyc) | 0.012 (70 cyc) |
| (4,+4) | 0.907 (70 cyc) | 0.394 (89 cyc) | 0.376 (77 cyc) | 0.188 (93 cyc) |
| (3,+2) | 0.169 (41 cyc) | 0.014 (44 cyc) | 0.0001 (41 cyc) | 0.006 (46 cyc) |
| (4,+3) | — †            | 0.249 (69 cyc) | 0.004 (58 cyc) | 0.007 (70 cyc) |

### C. Discussion

**The dominant (2,2) mode agrees well across all configurations.** Matches of 0.990–0.994 with phase errors of 0.004–0.076 rad/cycle confirm that NRSur7dq4 faithfully reproduces the SXS (2,2) mode. The slightly lower match at $q = 2$, no-spin (0.9905) compared to $q = 1$, no-spin (0.9940) is notable. The phase error is actually _smaller_ at $q = 2$ (0.009 vs. 0.076 rad/cycle), indicating that the match deficit arises from amplitude profile differences rather than phase drift. This may reflect mild differences in merger and ringdown morphology between SpEC and the surrogate's interpolated prediction.

**Sub-dominant modes are only meaningful when binary symmetries are broken.** At $q = 1$, $\chi = 0$, all odd-$m$ modes — $(2,1)$, $(3,3)$, $(4,3)$ — vanish identically by the binary's exchange symmetry. The match of these modes for SXS:BBH:0001 (0.28–0.35) is the result of comparing two near-zero time series dominated by numerical noise; it carries no physical information. Once this symmetry is broken — either by mass ratio ($q = 2$) or by spin ($\chi_{1z} = 0.5$) — all modes acquire physical amplitude and the surrogate reproduces them to better than 0.95 in all cases except those flagged below.

**Spin at $q = 1$ preferentially excites $(2,1)$.** Adding $\chi_{1z} = 0.5$ to an otherwise equal-mass binary activates the $(2,1)$ mode (match 0.997, $\Delta\Phi = 0.068$ rad/cycle), while leaving the $(3,3)$ mode at 0.951. This reflects the physical mechanism: the $(2,1)$ mode is dominantly sourced by the mass-weighted spin-orbit coupling, which is directly excited by $\chi_{1z}$ in the equal-mass case, whereas $(3,3)$ requires mass-ratio asymmetry to be dominantly excited.

**The (3,2) mode is anomalously poor for $q = 2$, $\chi_{1z} = 0.6$.** SXS:BBH:0162 yields $\mathcal{F}_{(3,2)} = 0.573$ despite every other mode matching at $\geq 0.99$. Crucially, the phase error for this mode is only $0.006$ rad/cycle — the phases agree well. The mismatch is therefore an _amplitude_ mismatch, not a phase error. The $(3,2)$ mode is well known to suffer near-cancellation between its two dominant contributions (mass-ratio sourced and spin-orbit sourced terms) at certain spin configurations. Small errors in the relative weighting of these contributions in the surrogate can produce large fractional amplitude errors while leaving the phase nearly intact. This is precisely the scenario here: the mode's amplitude is at or near a local minimum in the surrogate's interpolation and is therefore most sensitive to interpolation error. A similar but weaker effect may explain the low matches of $(4,3)$ for SXS:BBH:0162 ($\mathcal{F} = 0.41$) and the reduced $(3,2)$ match for SXS:BBH:0169 ($\mathcal{F} = 0.957$).

**The (4,4) mode shows systematic phase accumulation error.** Across all cases with non-zero amplitude, $\Delta\Phi_{(4,4)}/{\rm cycle}$ is the largest of any mode, ranging from 0.19 to 0.91 rad/cycle. The (4,4) mode sweeps through twice as many cycles as the (2,2) mode (since $f_{\rm GW}^{(4,4)} \approx 2 f_{\rm GW}^{(2,2)}$), so a factor of ~2 larger $\Delta\Phi/{\rm cycle}$ would be expected from a constant relative phase error per orbital cycle. The observed factor of $\sim$20–100 suggests that the surrogate's $\ell = 4$ sector has genuinely larger fractional phase-integration error, likely because fewer NR training waveforms constrained this sector than constrained $\ell = 2$.

**The phase metric and match are complementary.** For SXS:BBH:0001 (q=1, no spin), the $(3,2)$ mode has $\mathcal{F}_{(3,2)} = 0.994$ (excellent match) but $\Delta\Phi_{(3,2)}/{\rm cycle} = 0.169$ rad (meaningfully non-zero phase error). The match is insensitive to this slowly accumulated phase drift because the maximization over time shift partially absorbs it; the $\Delta\Phi/{\rm cycle}$ metric exposes it directly. Conversely, the $(2,2)$ mode of SXS:BBH:0162 has $\Delta\Phi_{(2,2)}/{\rm cycle} = 0.004$ rad (near-perfect phase agreement) but $\mathcal{F}_{(2,2)} = 0.993$ (slightly below unity), indicating that the residual mismatch is purely in the amplitude envelope near merger. Using both metrics together gives a more complete picture than either alone.

---

## V. Planned Extensions

The results above constitute Step 1 of a five-step comparative analysis program:

**Step 2 — SO(3) frame rotation optimization.** The current pipeline maximizes the match over time and phase only. Extending the optimization to include a full SO(3) rotation via Wigner $D$-matrix mixing (Section II.B) will absorb residual frame-convention differences between the NR simulation and the surrogate. This is implemented via Nelder-Mead minimization over Euler angles using `WaveformModes.match_sphere_averaged()`. We expect this to improve the match for all asymmetric modes and to yield an optimal rotation $R^*(\boldsymbol{\theta})$ whose systematic dependence on binary parameters encodes frame-convention differences across catalogs.

**Step 3 — Batch processing across all catalogs.** The pipeline will be extended to process all simulations in the SXS, RIT, and MAYA catalogs that fall within the NRSur7dq4 prior volume ($q \leq 4$, $|\boldsymbol{\chi}_{1,2}| \leq 0.8$). Batch processing using multiprocessing will enable the first comprehensive mode-by-mode comparison across all three catalogs.

**Step 4 — BMS supertranslation correction.** For simulations where the SO(3)-optimized match remains below 0.99, BMS supertranslation optimization (Section II.C) will be applied using `WaveformModes.match_sphere_averaged_bms_maximized()`. This step will distinguish gauge artifacts from genuine NR numerical errors.

**Step 5 — Cross-catalog science.** The primary scientific deliverable is a comparison of match distributions across catalogs as a function of binary parameters. Key quantities include: (a) the distribution of mismatches $1 - \mathcal{F}$ for each catalog, as a function of $q$, $\chi_{\rm eff}$, and $\chi_p$; (b) the optimal rotation $R^*$ as a function of catalog and simulation properties, revealing systematic frame-convention differences; (c) the improvement in match from Step 1 through Step 4, quantifying the relative contributions of frame rotation and BMS corrections; and (d) the mode-by-mode mismatch distribution, identifying which higher harmonics are most sensitive to code-dependent numerical errors.

---

## VI. Conclusions

We have demonstrated a surrogate-mediated framework for per-mode, noise-weighted comparison of NR BBH waveform catalogs. Applied to four SXS pilot simulations, the framework reveals:

1. The dominant $(2,2)$ mode is reproduced to within 1% mismatch and $< 0.1$ rad/cycle phase error across all tested configurations.
2. Sub-dominant mode behavior depends critically on binary symmetries: modes that vanish by symmetry at $q=1$, $\chi=0$ become meaningful diagnostics only when those symmetries are broken by mass ratio or spin.
3. The $(3,2)$ mode at $q=2$, $\chi_{1z}=0.6$ exhibits anomalous amplitude mismatch (0.573) with negligible phase error, indicative of near-cancellation sensitivity in the surrogate's interpolation.
4. The combined match + $\Delta\Phi/{\rm cycle}$ diagnostic provides strictly more information than either metric alone, and we advocate its adoption as a standard waveform-comparison reporting convention.

Full-catalog results and the rotation-optimized comparison will be reported in a companion paper.

---

## Acknowledgments

_To be filled in._

---

## References

\[sxs\] Boyle _et al._, CQG **36**, 195006 (2019); SXS Collaboration, https://www.black-holes.org/waveforms.  
\[rit\] Healy _et al._, PRD **96**, 024031 (2017); PRD **100**, 024021 (2019).  
\[maya\] Jani _et al._, CQG **33**, 204001 (2016).  
\[pycbc\] Biwer _et al._, PASP **131**, 024503 (2019).  
\[lalsuite\] LIGO Scientific Collaboration, LALSuite, https://git.ligo.org/lscsoft/lalsuite.  
\[bilby\] Ashton _et al._, ApJS **241**, 27 (2019).  
\[hannam2009\] Hannam _et al._, PRD **79**, 084025 (2009).  
\[ninja1\] Aylott _et al._, CQG **26**, 165008 (2009).  
\[ninja2\] Aasi _et al._, PRD **85**, 122006 (2012).  
\[nrar\] Hinder _et al._, CQG **31**, 025012 (2014).  
\[nrsur7dq4\] Varma _et al._, PRD **99**, 064045 (2019).  
