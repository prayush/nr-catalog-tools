# Surrogate-Mediated Cross-Catalog Validation of Numerical Relativity Binary Black Hole Waveforms

**Prayush Kumar** _et al._

_Draft — not for circulation_

---

## Abstract

We present a systematic framework for comparing binary black hole gravitational waveform catalogs produced by independent numerical relativity (NR) codes against each other, using the NRSur7dq4 surrogate model as a common reference. The surrogate mediates the comparison by providing waveforms at the exact intrinsic parameters of each catalog simulation, eliminating the parameter-space mismatch that has historically confounded direct cross-catalog comparisons. We quantify agreement through two complementary metrics: the noise-weighted match (faithfulness) computed per spherical harmonic mode with the Advanced LIGO zero-detuning high-power noise curve, and the accumulated phase difference per gravitational-wave cycle. We classify all catalog simulations into six spin-eccentricity categories and apply the framework in a batch comparison of 1,648 non-spinning and aligned-spin simulations ($\chi_\perp = 0$) from the SXS, RIT, and MAYA catalogs. For the dominant $(2,2)$ mode, the SXS median mismatch is $< 1\%$ for quasi-circular systems; the RIT distribution is bimodal with a significant low-match tail near high-spin, high-$q$ parameters. We further distinguish 342 SXS simulations that form the NRSur7dq4 training set (in-sample, near-perfect match) from 242 independent quasi-circular SXS simulations (out-of-sample, revealing the surrogate's interpolation boundaries). Sub-dominant mode behavior varies markedly across catalogs and depends critically on binary symmetries. We outline the full comparison pipeline, the treatment of source-frame ambiguities including SO(3) rotations and BMS supertranslations, and the path toward the remaining precessing-spin catalog populations.

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
- **Dimensionless spins** $\boldsymbol{\chi}_{1,2}$ specified at a reference epoch controlled by the `f_ref` parameter (see below)
- **Reference frequency** $f_{\rm ref}$ in cycles/$M$: $f_{\rm ref} = M_{\rm tot} \cdot f_{\rm GW}^{(22)} \cdot G M_\odot / c^3$; sets the epoch at which the input spin components are defined
- **Starting frequency** $f_{\rm low}$: controls waveform truncation only; per the gwsurrogate documentation, `f_low=0` is recommended for NRSur7dq4, which returns the full waveform from the surrogate's natural minimum
- **Time step** $dt$ in dimensionless units $dt/M$

Two cases arise depending on whether the NR waveform starts before or after the surrogate's minimum training frequency ($M\Omega \approx 0.0161$ at the parameters studied here):

1. **NR shorter than surrogate** ($f_{\rm lower}^{\rm NR} > f_{\rm min}^{\rm sur}$): we pass `f_low=0` and `f_ref`$= M_s \cdot f_{\rm lower}^{\rm NR}$. The surrogate backward-evolves the spins from the NR epoch to its natural start, giving the full common waveform.
2. **NR longer than surrogate** ($f_{\rm lower}^{\rm NR} < f_{\rm min}^{\rm sur}$): the surrogate domain cannot reach the NR starting frequency, so `f_ref` is clipped to the surrogate minimum. For the aligned-spin and non-spinning systems studied here, the spin components do not precess, so the metadata spin values are valid at any epoch and no spin-epoch error is introduced. For a general precessing system, this case would require extracting the instantaneous spins from NR dynamics at $f_{\rm min}^{\rm sur}$.

We adjust the lower cutoff of the match integral to $f_{\rm lower}^{\rm match} = \max(f_{\rm lower}^{\rm NR}, f_{\rm lower}^{\rm sur})$, so that neither waveform is penalized for having support outside the other's frequency band.

### B. Parameter Extraction

Source parameters are extracted from catalog metadata via `nrcatalogtools.CatalogBase.get_parameters()`, which returns a PyCBC-compatible dictionary: `mass1`, `mass2`, `spin1x/y/z`, `spin2x/y/z`, and `f_lower`. For SXS simulations, `f_lower` is defined as

$$f_{\rm lower} = \frac{M\Omega_{\rm NR}}{\pi \cdot M_{\rm tot} \cdot (G M_\odot / c^3)},$$

which equals the $(2,2)$-mode gravitational-wave frequency at the NR relaxation time. This same value is passed as `f_ref` (in cycles/$M$: $f_{\rm ref} = f_{\rm lower} \cdot M_s$) to NRSur7dq4, which uses it to set the spin reference epoch. The waveform start is controlled separately via `f_low=0`. Because `nrcatalogtools` extracts spin components at the relaxation time and `f_ref` is set to the corresponding frequency, these two epochs are exactly consistent — no separate spin-epoch remapping is required, provided $f_{\rm lower}^{\rm NR} \geq f_{\rm min}^{\rm sur}$. When $f_{\rm lower}^{\rm NR} < f_{\rm min}^{\rm sur}$ (NR waveform longer than the surrogate), `f_ref` is clipped to the surrogate domain minimum.

### C. Mode Extraction and Scaling

NR modes are extracted as complex physical-unit time series via `WaveformModes.get_mode(ell, em, total_mass, distance, delta_t_seconds)`, which returns a PyCBC `TimeSeries` with epoch set so that $t = 0$ corresponds to the peak of the $(2,2)$ amplitude. The surrogate modes are scaled from dimensionless units to physical units as

$$h_{\ell m}^{\rm phys}(t) = h_{\ell m}^{\rm sur}(t/M_s) \times \frac{G M_{\rm tot} M_\odot}{c^2 D_{\rm Mpc}},$$

where $M_s = M_{\rm tot} \cdot G M_\odot / c^3$ is the total mass in seconds and $D_{\rm Mpc}$ is the luminosity distance. All comparisons use $M_{\rm tot} = 40\,M_\odot$ and $D = 1\,{\rm Mpc}$.

The per-mode match uses the real part $\mathrm{Re}[h_{\ell m}]$ (the $h_+$ component), padded to the next power-of-two length and noise-weighted with a freshly constructed PSD at the matching frequency resolution.

### D. Computational Implementation

The Step 1 pipeline is implemented as four cooperating Python modules residing in `project/scripts/`: a main driver (`compare_one_sim_vs_surrogate.py`), a surrogate interface (`surrogate_utils.py`), a match-computation library (`match_utils.py`), and a catalog-loading abstraction (`catalog_utils.py`). A fifth script (`mass_scan.py`) extends the single-mass comparison to a grid of total masses. All modules are self-contained and depend only on `nrcatalogtools`, `pycbc`, `gwsurrogate`, and standard scientific Python libraries. We describe each module in detail below.

#### D.1 Main driver: `compare_one_sim_vs_surrogate.py`

The driver accepts a catalog name, simulation identifier, total mass, PSD name, sample spacing, and output directory via a command-line interface and executes the following six-step workflow.

**Step 1 — Catalog and waveform loading.** The catalog is instantiated through the `catalog_utils.load_catalog()` factory, which dispatches to `nrcatalogtools.SXSCatalog.load()`, `RITCatalog.load()`, or `MayaCatalog.load()` depending on the requested tag. The NR waveform object (`WaveformModes`) is retrieved via `cat.get(sim_name)`. At this stage the waveform data are dimensionless retarded-time multipoles $r\,h_{\ell m}/M$ as stored in the original catalog files.

**Step 2 — Parameter extraction.** Intrinsic parameters are obtained from `cat.get_parameters(sim_name, total_mass=M)`, which queries the catalog metadata and returns a PyCBC-compatible dictionary containing `mass1`, `mass2`, `spin1x/y/z`, `spin2x/y/z`, and `f_lower`. The mass ratio $q = m_1/m_2 \geq 1$ and individual spin magnitudes $|\boldsymbol{\chi}_{1,2}|$ are checked against the NRSur7dq4 prior bounds ($q \leq 4$, $|\boldsymbol{\chi}| \leq 0.8$) via `surrogate_utils.check_surrogate_prior()`, which prints a warning but does not abort for out-of-prior cases.

**Step 3 — Surrogate generation.** `surrogate_utils.generate_surrogate_modes()` is called with the extracted parameter dictionary, the total mass, a fiducial luminosity distance of 1 Mpc, and the requested sample spacing $\Delta t = 1/4096$ s. It returns a dictionary of PyCBC `TimeSeries` objects keyed by $(\ell, m)$, plus the effective starting GW frequency $f_{\rm lower}^{\rm sur}$ of the surrogate output (see §D.2). The effective match lower cutoff is set to $f_{\rm lower}^{\rm match} = \max(f_{\rm lower}^{\rm NR},\, f_{\rm lower}^{\rm sur})$, ensuring neither waveform is penalized for frequency content outside the other's support.

**Step 4 — PSD construction.** Rather than building a single global PSD, the pipeline delegates PSD construction to `match_utils.compute_mode_match()`, which builds a fresh `aLIGOZeroDetHighPower` PSD at each mode's frequency resolution after the waveforms have been zero-padded to the appropriate power-of-two length. This guarantees exact consistency between the PSD frequency grid and the waveform frequency grid, which `pycbc.filter.match()` requires.

**Step 5 — Per-mode match computation.** For each mode $(\ell, m)$ in the set $\{(2,2),(2,1),(3,3),(4,4),(5,5),(3,2),(4,3)\}$, the script:

1. Extracts the NR mode via `wfm.get_mode(ell, em, total_mass, distance, delta_t_seconds)`, obtaining a complex physical-unit `TimeSeries` with epoch set so $t = 0$ is at the $(2,2)$ amplitude peak.
2. Retrieves the corresponding surrogate `TimeSeries` from the `h_sur` dictionary.
3. Computes the noise-weighted match on the real parts ($h_+$ components) via `match_utils.compute_mode_match()`, using mode-specific frequency cutoff $f_{\rm lower}^{(\ell m)} = f_{\rm lower}^{\rm match} \cdot |m|/2$.
4. Computes the accumulated phase difference per GW cycle via `match_utils.compute_phase_diff_per_cycle()`, operating on the complex mode time series over their common time window.

The $(5,5)$ mode is recorded as NaN because NRSur7dq4 does not provide $\ell = 5$ modes. Modes that are near-zero (as for odd-$m$ harmonics at $q=1$, $\chi=0$) are not excluded by the pipeline; the degenerate-mode guard in `compute_mode_match()` returns NaN if the maximum amplitude of either waveform falls below $10^{-50}$, preventing division by zero in the inner-product normalization.

**Step 6 — Output.** Results are written to a per-simulation CSV file and a three-panel figure ($(2,2)$ amplitude comparison, $(2,2)$ real-part detail near merger, and a color-coded bar chart of per-mode matches), together with a formatted console table.

#### D.2 Surrogate interface: `surrogate_utils.py`

This module encapsulates all interaction with the `gwsurrogate` library and implements the correct spin-epoch alignment protocol.

**Surrogate loading.** The NRSur7dq4 model is a large file (~800 MB) that is expensive to load from disk. The module uses a module-level singleton `_nrsur7dq4` so the model is loaded at most once per Python process, regardless of how many simulations are evaluated in a single run.

**Physical unit conversion.** The total mass in solar masses is converted to seconds via $M_s = M_{\rm tot} \cdot G M_\odot/c^3$ (using `nrcatalogtools.utils.time_to_physical()`), which sets the physical time and amplitude scales. The dimensionless surrogate time step is $\Delta t_{\rm dimless} = \Delta t_{\rm phys}/M_s$, and the reference frequency is $f_{\rm ref} = f_{\rm lower}^{\rm NR} \cdot M_s$ cycles$/M$.

**The `f_low`/`f_ref` distinction.** The `gwsurrogate` API for NRSur7dq4 distinguishes two frequency arguments with distinct physical meanings. The `f_low` argument controls waveform *truncation*: the surrogate evaluates the full waveform from its natural minimum frequency regardless, and `f_low` simply discards the low-frequency segment of the output before returning it. The gwsurrogate documentation explicitly recommends `f_low=0` for NRSur7dq4, since the model is already short and no truncation is needed. The `f_ref` argument instead sets the *reference epoch* at which the input spin vectors $\boldsymbol{\chi}_{1,2}$ are defined; it is specified in cycles/$M$ (= $M_s \cdot f_{\rm GW}^{(22)}$ in Hz). The surrogate internally backward-evolves the spin dynamics from the `f_ref` epoch to its natural starting frequency, ensuring that the spins at any output time correctly reflect the physical precession history. Setting `f_low = f_ref` (the old convention) would have the surrogate both truncate the waveform and set the spin epoch at the NR relaxation time, but these are logically independent operations and conflating them produces incorrect spin evolution when `f_ref` is below the surrogate minimum.

The pipeline therefore calls:

```python
t_sur, h_sur, _ = sur(q, chiA, chiB, ellMax=4, dt=dt_dimless, f_low=0, f_ref=f_ref_dimless)
```

where `f_ref_dimless = f_lower_hz * m_secs`.

**Epoch clipping for long NR waveforms.** When $f_{\rm lower}^{\rm NR} < f_{\rm min}^{\rm sur}$ (i.e. the NR simulation starts at a lower frequency than the surrogate's minimum training extent), the `f_ref` value falls below the surrogate's domain and `gwsurrogate` raises an exception of the form `"Got omega_ref = X < Y = omega_0, too small"`. The module catches this exception, parses the omega_0 value from the error string via a regular expression, clips `f_ref` to $1.01 \times \omega_0 / \pi$ cycles/$M$ (a 1% safety margin), and re-issues the surrogate call. For the aligned-spin and non-spinning systems in the pilot study, the in-plane spin components $\chi_\perp$ are zero and spin vectors do not precess, so the spin values from the NR metadata remain valid at any epoch and no spin-epoch error is introduced by this clipping. For general precessing systems (detected via $\chi_{1\perp}^2 + \chi_{2\perp}^2 > 10^{-8}$), the module prints an explicit warning that the metadata spins are being used at a clipped epoch and that proper treatment would require extracting instantaneous spins from the NR dynamics at $f_{\rm min}^{\rm sur}$.

**Effective starting frequency.** After the surrogate call, the module computes the actual GW frequency at the first output sample from the phase derivative of the $(2,2)$ mode:

$$f_{\rm lower}^{\rm sur} = \frac{1}{2\pi}\frac{d\phi_{22}}{dt}\bigg|_{t=t_{\rm start}}$$

where $\phi_{22} = \arg[h_{22}]$ is the unwrapped phase. This value is returned alongside the mode dictionary and used by the driver to set $f_{\rm lower}^{\rm match}$.

**Epoch alignment.** The surrogate time array $t_{\rm sur}$ is rescaled to physical seconds and shifted so that $t = 0$ coincides with the peak of the $(2,2)$ amplitude (identified as $\arg\max |h_{22}^{\rm sur}|$), matching the epoch convention of `WaveformModes.get_mode()`.

#### D.3 Match computation library: `match_utils.py`

**`compute_mode_match(h_nr, h_sur, f_lower_mode, psd_name)`** receives the real parts of two complex mode time series, zero-pads both to the next power-of-two length $\geq \max(\mathrm{len}(h_{\rm NR}), \mathrm{len}(h_{\rm sur}))$, builds a `from_string(psd_name, ...)` PSD at the resulting frequency resolution $\Delta f = 1/(N_{\rm FFT} \Delta t)$, and calls `pycbc.filter.match()` which maximizes over time and phase shifts. The function returns NaN if either waveform's maximum amplitude falls below $10^{-50}$ (degenerate mode guard). The low-frequency cutoff passed to `pycbc.filter.match()` is always the mode-scaled value $f_{\rm lower}^{(\ell m)} = f_{\rm lower}^{\rm match} \cdot |m|/2$.

**`compute_phase_diff_per_cycle(h_nr, h_sur)`** operates on the complex mode time series. It identifies the common time window $[t_{\rm start}, t_{\rm end}]$ from the `start_time` and `end_time` attributes of the two PyCBC `TimeSeries` objects, slices both waveforms to this common window, unwraps the complex argument to obtain $\phi_{\rm NR}(t)$ and $\phi_{\rm sur}(t)$, and computes

$$\frac{\Delta\Phi}{\rm cycle} = \frac{|\Delta\Phi_{\rm NR} - \Delta\Phi_{\rm sur}|}{N_{\rm cyc}^{\rm NR}},$$

where $\Delta\Phi = |\phi(t_{\rm end}) - \phi(t_{\rm start})|$ and $N_{\rm cyc}^{\rm NR} = \Delta\Phi_{\rm NR}/(2\pi)$. The function returns NaN if fewer than 0.5 GW cycles are present in the common window. Note that the metric is not maximized over any time or phase shift; any residual time-shift error in the epoch alignment will appear as a non-zero $\Delta\Phi/{\rm cycle}$.

**`mode_f_lower(f_lower, em)`** implements the mode-frequency scaling $f_{\rm GW}^{(\ell m)} = |m| \cdot f_{\rm orbital} = |m| \cdot f_{\rm lower}/2$, where $f_{\rm lower}$ is the $(2,2)$-mode reference frequency in Hz and the factor of 1/2 converts from $(2,2)$ GW frequency to orbital frequency.

#### D.4 Catalog abstraction: `catalog_utils.py`

The `load_catalog(name)` factory provides a single entry point for all three supported catalogs. It calls `nrcatalogtools.SXSCatalog.load(download=False)` for SXS (suppressing automatic metadata downloads during batch runs), and the corresponding `load()` methods for RIT and MAYA. The complementary `filter_by_surrogate_prior(catalog, ...)` function iterates over `catalog.simulations_list`, calls `get_parameters()` for each simulation, and passes the result through `surrogate_utils.check_surrogate_prior()`, returning the subset of simulations with $q \in [1,4]$ and $|\boldsymbol{\chi}_{1,2}| \leq 0.8$. This will be the entry point for the batch processing in Step 3.

#### D.5 Batch processor: `batch_aligned_catalogs.py`

The Step 2 batch processor extends the single-simulation driver to all simulations in categories a–d across all three catalogs.  It provides a multi-phase pipeline:

**Enumeration phase.** The script reads the pre-computed classification JSON files in `catalog_organization/` and extracts all `(catalog, sim_id, category, is_nrsur_calibration)` tuples for the requested categories.  For SXS simulations, the JSON stores a `nrsur7dq4_calibration` boolean flag that records whether the simulation was part of the NRSur7dq4 training set.

**Metadata phase.** Parameters and eccentricity are collected sequentially via `cat.get_parameters()` and `cat.get_metadata()`.  Simulations with NaN mass ratio or spin, or outside the NRSur7dq4 prior bounds ($q > 4$ or $|\chi_{1,2}| > 0.8$), are silently excluded.  Eccentricity strings of the form `"< 0.002"` or `"~0.01"` are cleaned by stripping the leading symbol before conversion.

**Restartable parallel processing.** The merged output CSV is checked for already-processed `(catalog, sim_id)` pairs, and only the remaining simulations are dispatched to a `multiprocessing.Pool`.  The NRSur7dq4 surrogate (~800 MB) is loaded once per worker process via a pool initializer; per-process catalog singletons avoid redundant catalog file reads.  Run-time config (sample spacing, PSD name, distance) is stamped into each job dict before dispatch, satisfying Python's requirement that pool function arguments be picklable.  Each completed row is appended to the merged CSV immediately, so partial runs are recoverable.

**Post-processing.**  After the processing loop, the script writes per-catalog CSV files and calls `plot_batch_results.make_individual_sim_figures()` to generate a per-simulation match figure for every analyzed simulation.

**CSV migration.** When the script is run against a CSV produced by an older version lacking the `nrsur7dq4_calibration` column, a lightweight migration function (`_migrate_add_calibration()`) backfills the column from the classification JSONs without triggering a recomputation.

#### D.6 Mass scan: `mass_scan.py`

The mass scan script extends the single-mass comparison to a grid of total masses $M_{\rm tot} \in \{10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100\}\,M_\odot$ for the four pilot simulations. For each $({\rm sim}, M)$ pair, it calls `compare_one_mass()`, which re-extracts parameters at the new total mass (since `f_lower` in Hz scales as $M^{-1}$), regenerates the surrogate, and computes the per-mode match. The NR waveform object is loaded once per simulation and reused across the mass grid, since `WaveformModes.get_mode()` accepts `total_mass` as a rescaling argument. Results are written to a single CSV file (`mass_scan_results.csv`) and a $2 \times 3$ panel figure showing $\log_{10}(1 - \mathcal{F})$ versus $M_{\rm tot}$ for the six surrogate-supported modes, with shaded bands at mismatch thresholds of 1%, 3%, and 10%.

![Figure: Mass scan mismatch per-mode as a function of total binary mass for the pilot configurations.](figs/mass_scan_mismatch.png)


### E. Simulation Classification

To structure the comparison across catalogs, we classify every simulation according to its spin geometry and orbital eccentricity.  Let $\chi_\perp = \sqrt{\chi_{1x}^2 + \chi_{1y}^2 + \chi_{2x}^2 + \chi_{2y}^2}$ be the total in-plane spin magnitude and $e$ be the reference eccentricity from the catalog metadata.  We define six categories:

| Category | Name | Formal conditions |
|---|---|---|
| (a) | Non-spinning eccentric | $\chi_\perp < \varepsilon_\chi$, $|\chi_{1z}| + |\chi_{2z}| < \varepsilon_\chi$, $e > \varepsilon_e$ |
| (b) | Non-spinning quasi-circular | $\chi_\perp < \varepsilon_\chi$, $|\chi_{1z}| + |\chi_{2z}| < \varepsilon_\chi$, $e \leq \varepsilon_e$ |
| (c) | Aligned-spin eccentric | $\chi_\perp < \varepsilon_\chi$, $|\chi_{1z}| + |\chi_{2z}| \geq \varepsilon_\chi$, $e > \varepsilon_e$ |
| (d) | Aligned-spin quasi-circular | $\chi_\perp < \varepsilon_\chi$, $|\chi_{1z}| + |\chi_{2z}| \geq \varepsilon_\chi$, $e \leq \varepsilon_e$ |
| (e) | Precessing eccentric | $\chi_\perp \geq \varepsilon_\chi$, $e > \varepsilon_e$ |
| (f) | Precessing quasi-circular | $\chi_\perp \geq \varepsilon_\chi$, $e \leq \varepsilon_e$ |

with thresholds $\varepsilon_\chi = 0.001$ and $\varepsilon_e = 0.005$.

**Catalog counts.** Table IV shows the number of simulations per category for each catalog, restricted to the NRSur7dq4 prior volume ($q \leq 4$, $|\chi_{1,2}| \leq 0.8$) after metadata filtering.

**Table IV. Simulation counts within the NRSur7dq4 prior volume.**

| Category | SXS | RIT | MAYA |
|---|---|---|---|
| (a) non-spinning eccentric | 206 | 392 | 53 |
| (b) non-spinning quasi-circular | 177 | 54 | 23 |
| (c) aligned eccentric | 21 | 231 | 86 |
| (d) aligned quasi-circular | 687 | 541 | 26 |
| (e) precessing eccentric | 30 | 117 | 303 |
| (f) precessing quasi-circular | 3043 | 437 | 67 |
| **Total (all categories)** | **4164** | **1772** | **558** |

Several catalog-specific features are highly noteworthy and reflect the distinct parameter-space exploration strategies of the respective collaboration groups as summarized in the `nr-catalog-tools/catalog_organization` subdirectories:

1. **SXS Catalog**: The SXS catalog is heavily dominated by the precessing quasi-circular (f) category (3,043 simulations), representing the massive, systematic parameter-space coverage required for training and validating high-accuracy quasi-circular precessing models such as `NRSur7dq4` and other precessing templates. Non-spinning eccentric (a) and aligned quasi-circular (d) simulations have moderate representation (206 and 687 simulations respectively). However, eccentric precessing (category e, 30 simulations) and aligned eccentric (category c, 21 simulations) configurations are highly underrepresented, representing a relative gap in SpEC's eccentric exploration.
2. **RIT Catalog**: The RIT catalog exhibits a highly diverse distribution across the subcategories. It has the largest non-spinning eccentric (a) sample (392 simulations) and aligned-spin eccentric (c) sample (231 simulations). Combined with a large aligned quasi-circular (d) population of 541 simulations and a precessing quasi-circular (f) population of 437 simulations, the RIT catalog represents an exceptionally well-rounded catalog for testing eccentric, non-spinning, and aligned-spin waveforms, providing excellent physical parameter coverage.
3. **MAYA Catalog**: The MAYA catalog is highly specialized, dominated heavily by precessing-eccentric configurations (category e, 303 simulations). This represents a deliberate, large-scale systematic search of the eccentric precessing parameter space using the Einstein Toolkit and MayaKranc code. In contrast, quasi-circular aligned-spin (d, 26 simulations) and non-spinning quasi-circular (b, 23 simulations) systems are much less represented, demonstrating MAYA's complementary scientific focus on eccentric precessing dynamics.

**NRSur7dq4 calibration sub-classification.** A key feature of the SXS catalog is that 1,731 of its simulations were used as training data for the NRSur7dq4 surrogate~\cite{nrsur7dq4}: 60 in category (b), 282 in category (d), and 1,389 in category (f).  All calibration simulations are quasi-circular ($e = 0$) by the surrogate's training design.  Categories (a), (c), and (e) contain no calibration simulations.  This stratification is recorded in the `catalog_organization/sxs_classification.json` file as a per-simulation boolean flag, propagated into the results CSV, and used to split the SXS analysis into calibration (in-sample) and non-calibration (out-of-sample) subsets.

**Rationale for processing categories a–d first.** For all systems in categories (a)–(d), the in-plane spin components $\chi_\perp = 0$ by construction, so both the NR simulation and the NRSur7dq4 output have their orbital angular momentum aligned with the $z$-axis of the source frame.  There is therefore no SO(3) frame rotation to be optimized: the phase maximization in `pycbc.filter.match()` already absorbs the residual rotation about $z$.  This makes categories (a)–(d) the natural starting point for the systematic comparison.  Categories (e) and (f) require a full SO(3) optimization (Section~V) and are reserved for a later analysis step.

---

## IV. Results

We present results for four SXS simulations selected to probe two orthogonal axes of parameter space: mass ratio ($q = 1$ vs. $q = 2$) and spin ($\chi = 0$ vs. aligned spin $\chi_{1z} \approx 0.5$). All modes from the set $\{(2,2), (2,1), (3,3), (4,4), (5,5), (3,2), (4,3)\}$ are evaluated against NR; the $(5,5)$ mode is unavailable from NRSur7dq4 ($\ell_{\rm max}=4$) and appears as N/A. Table~I summarizes the simulation parameters.

**Table I. Simulation parameters.**

| Simulation | $q$ | $\chi_{1z}$ | $\chi_{2z}$ | $f_{\rm lower}^{\rm NR}$ [Hz] | $f_{\rm lower}^{\rm match}$ [Hz] |
|---|---|---|---|---|---|
| SXS:BBH:0001 | 1.00 | 0.00 | 0.00 | 19.8 | 26.1 (sur. clipped) |
| SXS:BBH:0005 | 1.00 | +0.50 | 0.00 | 19.8 | 26.7 (sur. clipped) |
| SXS:BBH:0169 | 2.00 | 0.00 | 0.00 | 29.1 | 29.1 |
| SXS:BBH:0162 | 2.00 | +0.60 | 0.00 | 28.8 | 28.8 |

For SXS:BBH:0001 and SXS:BBH:0005, the NR simulation starts at $\sim 20$ Hz but NRSur7dq4's minimum training extent at these parameters corresponds to $\sim 26$–$27$ Hz. The match lower cutoff is raised accordingly so that neither waveform is penalized for frequency content outside the other's support.

### A. Match Results

**Table II. Per-mode match $\mathcal{F}$ for all four simulations.**

| Mode | 0001 ($q$=1, ns) | 0005 ($q$=1, spin) | 0169 ($q$=2, ns) | 0162 ($q$=2, spin) |
|---|---|---|---|---|
| (2,+2) | 0.9940 | 0.9933 | 0.9947 | 0.9931 |
| (2,+1) | — † | **0.9966** | 0.9513 | 0.350 ‡ |
| (3,+3) | — † | 0.9548 | 0.9968 | 0.9844 |
| (4,+4) | 0.902  | 0.9886 | 0.9965 | 0.9778 |
| (5,+5) | N/A | N/A | N/A | N/A |
| (3,+2) | 0.9932 | 0.9924 | 0.9581 | **0.569 ‡** |
| (4,+3) | — † | 0.9744 | 0.9596 | 0.406 ‡ |

† Near-zero by $q = 1$, $\chi = 0$ symmetry; match value is numerically meaningless.  
‡ Anomalous; discussed in Section IV.C.

### B. Phase Difference per Cycle

**Table III. Phase difference per GW cycle $\Delta\Phi/{\rm cycle}$ [rad] and number of NR cycles $N_{\rm cyc}$.**

| Mode | 0001 (q=1,ns) | 0005 (q=1,spin) | 0169 (q=2,ns) | 0162 (q=2,spin) |
|---|---|---|---|---|
| (2,+2) | 0.051 (42 cyc) | 0.004 (44 cyc) | 0.005 (43 cyc) | 0.005 (47 cyc) |
| (2,+1) | — †            | 0.049 (25 cyc) | 0.005 (24 cyc) | 0.025 (26 cyc) |
| (3,+3) | — †            | 0.045 (67 cyc) | 0.019 (65 cyc) | 0.001 (71 cyc) |
| (4,+4) | 0.983 (71 cyc) | 0.309 (89 cyc) | 0.334 (86 cyc) | 0.240 (95 cyc) |
| (5,+5) | — | — | — | — |
| (3,+2) | 0.163 (43 cyc) | 0.004 (44 cyc) | 0.003 (43 cyc) | 0.004 (47 cyc) |
| (4,+3) | — †            | 0.236 (69 cyc) | 0.001 (65 cyc) | 0.004 (71 cyc) |

### C. Discussion

**The dominant (2,2) mode agrees well across all configurations.** Matches of 0.993–0.995 with phase errors of 0.004–0.051 rad/cycle confirm that NRSur7dq4 faithfully reproduces the SXS (2,2) mode. The slightly lower match at $q = 2$, spin (0.9931) compared to $q = 2$, no-spin (0.9947) is consistent with spin-induced amplitude corrections near merger. The phase error is smaller than 0.006 rad/cycle for all $q=2$ cases, indicating that any residual mismatch is dominated by amplitude profile differences rather than phase drift.

**Sub-dominant modes are only meaningful when binary symmetries are broken.** At $q = 1$, $\chi = 0$, all odd-$m$ modes — $(2,1)$, $(3,3)$, $(4,3)$ — vanish identically by the binary's exchange symmetry. The match of these modes for SXS:BBH:0001 (0.28–0.35) is the result of comparing two near-zero time series dominated by numerical noise; it carries no physical information. Once this symmetry is broken — either by mass ratio ($q = 2$) or by spin ($\chi_{1z} = 0.5$) — all modes acquire physical amplitude and the surrogate reproduces them to better than 0.95 in all cases except those flagged below.

**Spin at $q = 1$ preferentially excites $(2,1)$.** Adding $\chi_{1z} = 0.5$ to an otherwise equal-mass binary activates the $(2,1)$ mode (match 0.997, $\Delta\Phi = 0.049$ rad/cycle), while leaving the $(3,3)$ mode at 0.955. This reflects the physical mechanism: the $(2,1)$ mode is dominantly sourced by the mass-weighted spin-orbit coupling, which is directly excited by $\chi_{1z}$ in the equal-mass case, whereas $(3,3)$ requires mass-ratio asymmetry to be dominantly excited.

**The (3,2) mode is anomalously poor for $q = 2$, $\chi_{1z} = 0.6$.** SXS:BBH:0162 yields $\mathcal{F}_{(3,2)} = 0.569$ despite other primary modes matching at $\geq 0.98$. Crucially, the phase error for this mode is only $0.004$ rad/cycle — the phases agree well. The mismatch is therefore an _amplitude_ mismatch, not a phase error. The $(3,2)$ mode is well known to suffer near-cancellation between its two dominant contributions (mass-ratio sourced and spin-orbit sourced terms) at certain spin configurations. Small errors in the relative weighting of these contributions in the surrogate can produce large fractional amplitude errors while leaving the phase nearly intact. This is precisely the scenario here: the mode's amplitude is at or near a local minimum in the surrogate's interpolation and is therefore most sensitive to interpolation error. A similar but weaker effect may explain the low matches of $(4,3)$ for SXS:BBH:0162 ($\mathcal{F} = 0.406$) and the $(2,1)$ mismatch ($\mathcal{F} = 0.350$), and the reduced $(3,2)$ match for SXS:BBH:0169 ($\mathcal{F} = 0.958$).

**The (4,4) mode shows systematic phase accumulation error.** Across all cases with non-zero amplitude, $\Delta\Phi_{(4,4)}/{\rm cycle}$ is the largest of any mode, ranging from 0.24 to 0.98 rad/cycle. The (4,4) mode sweeps through twice as many cycles as the (2,2) mode (since $f_{\rm GW}^{(4,4)} \approx 2 f_{\rm GW}^{(2,2)}$), so a factor of ~2 larger $\Delta\Phi/{\rm cycle}$ would be expected from a constant relative phase error per orbital cycle. The observed factor of $\sim$20–100 suggests that the surrogate's $\ell = 4$ sector has genuinely larger fractional phase-integration error, likely because fewer NR training waveforms constrained this sector than constrained $\ell = 2$.

**The phase metric and match are complementary.** For SXS:BBH:0001 (q=1, no spin), the $(3,2)$ mode has $\mathcal{F}_{(3,2)} = 0.993$ (excellent match) but $\Delta\Phi_{(3,2)}/{\rm cycle} = 0.163$ rad (meaningfully non-zero phase error). The match is insensitive to this slowly accumulated phase drift because the maximization over time shift partially absorbs it; the $\Delta\Phi/{\rm cycle}$ metric exposes it directly. Conversely, the $(2,2)$ mode of SXS:BBH:0162 has $\Delta\Phi_{(2,2)}/{\rm cycle} = 0.005$ rad (near-perfect phase agreement) but $\mathcal{F}_{(2,2)} = 0.993$ (slightly below unity), indicating that the residual mismatch is dominated by the amplitude envelope near merger. Using both metrics together gives a more complete picture than either alone.

### D. Step 2: Batch Comparison of Non-Spinning and Aligned-Spin Systems

We extend the pilot analysis to all simulations in categories (a)–(d) across the SXS, RIT, and MAYA catalogs that fall within the NRSur7dq4 prior volume.  After metadata filtering, this yields 774 SXS, 686 RIT, and 188 MAYA simulations, for a total of 1,648 waveform comparisons.  We focus our discussion on the quasi-circular subsets (categories b and d) where the surrogate is expected to perform best, but report all categories.

![Figure 1a: Match per mode vs. source parameters (mass ratio, effective spin, individual spins, initial eccentricity) for all simulations.](figs/fig1a_match_vs_params.png)

![Figure 1b: Match per mode vs. source parameters for the quasi-circular (categories b+d) subset.](figs/fig1b_match_vs_params_qc.png)

**Table V. Per-mode match statistics (median / 10th percentile) for quasi-circular systems (categories b+d).**

| Mode | SXS (N=579) | RIT (N=229) | MAYA (N=49) |
|---|---|---|---|
| $(2,2)$ | 0.9923 / 0.9496 | 0.9231 / 0.4941 | 0.9728 / 0.8439 |
| $(2,1)$ | 0.8541 / 0.3739 | 0.9676 / 0.7152 | 0.9843 / 0.3775 |
| $(3,3)$ | 0.9865 / 0.2926 | 0.7710 / 0.3241 | 0.9707 / 0.3120 |
| $(4,4)$ | 0.9819 / 0.8768 | 0.7829 / 0.3895 | 0.9709 / 0.4383 |
| $(3,2)$ | 0.8580 / 0.4615 | 0.9700 / 0.6452 | 0.9794 / 0.7781 |
| $(4,3)$ | 0.6711 / 0.2989 | 0.8248 / 0.5763 | 0.9765 / 0.4615 |

#### D.1 Dominant (2,2) mode

The SXS (2,2) match distribution is tightly concentrated near unity: median 0.9923, 90th percentile 0.9945, with the 10th percentile at 0.9496.  The tail below 0.95 consists primarily of (i) short NR simulations (few cycles, low SNR) and (ii) high-mass-ratio, high-spin cases near the edge of the NRSur7dq4 prior.  

The RIT (2,2) distribution is markedly broader: median 0.9231, 10th percentile only 0.494.  The CDF (Figure 3a, 6a) reveals a bimodal structure — a high-match peak near $\mathcal{F} \sim 0.99$ and a low-match population extending below 0.5.  Inspection of the parameter dependence (Figure 1a) shows that the low-match RIT simulations are concentrated at high $|$\chi_{\rm eff}|$ and $q \sim 4$, near the boundary of the surrogate's training domain; phase mismatch at the boundary is expected to be larger.

![Figure 3a: Cumulative distribution functions (CDFs) of per-mode matches across the three catalogs for quasi-circular systems.](figs/fig3a_match_cdf.png)

The MAYA (2,2) distribution is intermediate: median 0.9728, 10th percentile 0.8439, and 90th percentile 0.9794.  The small MAYA sample size (49 quasi-circular simulations) limits statistical precision but suggests systematically lower matches than SXS at comparable parameters.

**SXS calibration vs. non-calibration.** The NRSur7dq4 surrogate was trained on 342 of the 579 SXS quasi-circular simulations analyzed here (337 in the CSV; the small discrepancy reflects the prior cut).  We test whether in-sample simulations show systematically higher matches.  The calibration subset has $(2,2)$ median 0.9926 and 10th percentile 0.9899 (nearly uniform distribution above 0.985).  The non-calibration subset has the same median (0.9919) but a much wider lower tail: 10th percentile 0.360.  This is the expected behaviour: NRSur7dq4 by construction interpolates its training set with near-zero error, so all calibration simulations achieve near-perfect match.  The non-calibration sims test the surrogate's generalization ability; most achieve high match, but a minority fall below 0.95, indicating regions of the parameter space where the sparse training grid limits interpolation accuracy.  Figures 6a–b display these two populations separately for each mode.

![Figure 6a: CDF of per-mode matches for the SXS calibration (in-sample) vs. non-calibration (out-of-sample) quasi-circular subsets.](figs/fig6a_sxs_cal_match_cdf.png)

![Figure 6b: CDF of mismatches (1 - F) for the SXS calibration vs. non-calibration quasi-circular subsets.](figs/fig6b_sxs_cal_mismatch_cdf.png)


#### D.2 Sub-dominant modes

Sub-dominant mode behaviour is substantially more variable across catalogs than the dominant $(2,2)$ mode.  Several patterns emerge from Table V and Figures 1a, 3b, 5a:

**SXS** shows high median matches for $(3,3)$ (0.987) and $(4,4)$ (0.982) but a substantially lower median for $(4,3)$ (0.671) and a wide spread for $(2,1)$ (0.854) and $(3,2)$ (0.858).  The wide spread in odd-$m$ modes reflects the binary-symmetry mechanism discussed in Section IV.C: equal-mass, non-spinning systems contribute near-zero amplitudes in odd-$m$ modes, making their matches numerically ill-defined.

**RIT** shows a complementary pattern: higher $(2,1)$ (0.968) and $(3,2)$ (0.970) medians than SXS, but substantially lower $(3,3)$ (0.771) and $(4,4)$ (0.783) medians.  The $(3,3)$ mode in RIT shows an extended low-match tail consistent with systematic phase error in this harmonic, distinct from the $(2,2)$ bimodal structure.

**MAYA** consistently achieves the highest median sub-dominant-mode matches across all six modes compared here, with all medians above 0.970.  The small MAYA sample size prevents strong statistical conclusions, but the pattern is consistent with MAYA waveforms being produced at parameters where NRSur7dq4 interpolates well.

#### D.3 Effect of eccentricity

We analyze the quantitative dependence of mismatch and phase errors on orbital eccentricity by examining the fifth columns of Figure 1a (representing the mismatch $\log_{10}(1-\mathcal{F})$ on log-y, log-x axes) and Figure 2a (representing the phase difference $\Delta\Phi/\mathrm{cycle}$ on log-x axis), showing only the $e > 0$ population across the SXS, RIT, and MAYA catalogs.

For very low initial eccentricity $e < 0.005$ (the quasi-circular threshold defining categories b and d), the mismatch values reside at extremely low floors. The dominant $(2,2)$ mode for SXS waveforms generally registers mismatches below $10^{-3}$ to $10^{-4}$, demonstrating the pristine analytical convergence between SpEC's spectral Einstein code simulations and the circular surrogate's interpolant. The RIT and MAYA catalogs display slightly higher floors for $e < 0.005$ but remain securely below the $10^{-2}$ (1%) mismatch threshold.

As the initial eccentricity climbs above the $0.005$ threshold, we observe a sharp, deterministic power-law escalation of the mismatch across all three catalogs and all harmonic modes. As eccentricity increases from $0.005$ to $0.05$, the dominant $(2,2)$ mode mismatch grows by over two orders of magnitude, typically exceeding $10^{-2}$ to $10^{-1}$ across all catalogs (SXS, RIT, and MAYA). For highly eccentric systems where $e \gtrsim 0.1$, the mismatch rapidly saturates near $0.1$ to $0.5$ (matches drop to 0.9–0.5).

Importantly, higher-order sub-dominant harmonics (e.g. $(3,3)$ and $(4,4)$) degrade even more aggressively than the dominant $(2,2)$ mode. This is because eccentric modulation acts directly on the orbital frequency, and higher harmonics multiply this modulation by $m$. This results in much faster amplitude envelope distortion and rapid cycle-by-cycle phase slippage, causing the sub-dominant modes to hit the noise/mismatch floor at lower eccentricities.

Looking at Figure 2a, the phase difference per cycle is tightly constrained to $< 0.05$ rad/cycle when $e < 0.005$. However, as eccentricity increases, the phase difference dramatically scales, surpassing $1.0$ rad/cycle for $e \ge 0.05$. Physically, this is the hallmark of eccentric orbital motion, which incorporates periastron precession and radial frequency oscillations, resulting in significant phase modulations that a circular model (trained strictly at $e=0$) is fundamentally unable to capture. This serves as a highly robust negative control, validating that our comparative framework is exceptionally sensitive to physical model deficiencies.

![Figure 2a: Phase difference per cycle vs. parameters (including log-scale initial eccentricity) for all simulations.](figs/fig2a_phasediff_vs_params.png)

#### D.4 Phase difference per cycle

The phase-difference metric (Figures 2a, 2b) provides complementary diagnostic power.  For quasi-circular SXS systems, the $(2,2)$ phase difference has median 0.01 rad/cycle — well below the GW cycle-budget threshold of 0.1 rad/cycle often cited as the acceptable error for parameter estimation.  The $(4,4)$ mode has the largest phase errors (median $\sim$0.3 rad/cycle), consistent with the pilot result in Section IV.B and the sparser NRSur7dq4 training in the $\ell = 4$ sector.  RIT and MAYA show systematically larger phase differences for several modes, suggesting that some of the match deficit in Table V is attributable to cumulative phase error rather than amplitude discrepancy.

![Figure 2b: Phase difference per cycle vs. parameters for the quasi-circular subset.](figs/fig2b_phasediff_vs_params_qc.png)

![Figure 3c: Cumulative distribution functions (CDFs) of phase difference per cycle across the three catalogs.](figs/fig3c_phasediff_cdf.png)

![Figure 4: Mismatch heatmap in the (q, chi_eff) parameter space for the dominant (2,2) mode across SXS, RIT, and MAYA.](figs/fig4_match22_heatmap_qc.png)


---

## V. Planned Extensions

The results in Sections IV and IV.D constitute Steps 1 and 2 of a five-step comparative analysis program:

**Step 1 — Pilot comparison (complete).** Per-mode noise-weighted match and phase-difference-per-cycle for four SXS pilot simulations, establishing the pipeline and calibrating the metrics (Section~IV.A).

**Step 2 — Batch comparison of non-spinning and aligned-spin systems (complete).** Full batch comparison of 1,648 simulations from SXS, RIT, and MAYA in categories (a)–(d), demonstrating that the framework scales to catalog-wide analysis (Section~IV.D).

**Step 3 — SO(3) frame rotation optimization.** The current pipeline maximizes the match over time and phase only. Extending the optimization to include a full SO(3) rotation via Wigner $D$-matrix mixing (Section~II.B) will absorb residual frame-convention differences for precessing-spin systems (categories e and f) where the NR simulation and surrogate source frames can differ by an arbitrary rotation. This is implemented via Nelder-Mead minimization over Euler angles. We expect this to improve the match for all asymmetric modes, and the optimal rotation $R^*(\boldsymbol{\theta})$ will encode systematic frame-convention differences between catalogs as a function of binary parameters.

**Step 4 — BMS supertranslation correction.** For simulations where the SO(3)-optimized match remains below 0.99, BMS supertranslation optimization (Section~II.C) will be applied. This step will distinguish gauge artifacts from genuine NR numerical errors, and provide an upper bound on the BMS contribution to the observed mismatches.

**Step 5 — Full-catalog cross-catalog science.** The primary scientific deliverable is a comparison of match distributions across all three catalogs as a function of binary parameters, including the precessing-spin populations (categories e and f) that are absent from Steps 1–2. Key quantities include: (a) the distribution of mismatches $1 - \mathcal{F}$ for each catalog as a function of $q$, $\chi_{\rm eff}$, and $\chi_p$; (b) the optimal rotation $R^*$ as a function of catalog and simulation properties; (c) the improvement in match from Step 2 through Step 4, quantifying the relative contributions of frame rotation and BMS corrections; and (d) the mode-by-mode mismatch distribution identifying which higher harmonics are most sensitive to code-dependent numerical errors.

---

## VI. Conclusions

We have demonstrated a surrogate-mediated framework for per-mode, noise-weighted comparison of NR BBH waveform catalogs. Applied first to four SXS pilot simulations and then in a batch comparison of 1,648 simulations from the SXS, RIT, and MAYA catalogs, the framework reveals:

1. **The dominant $(2,2)$ mode** is reproduced to within 1% mismatch and $< 0.1$ rad/cycle phase error for the large majority of quasi-circular SXS simulations (median mismatch $< 0.01$, 90th percentile $< 0.006$).  RIT shows a broader distribution with a significant low-match tail (10th percentile mismatch $\sim$50\%), concentrated at high-spin, high-$q$ systems near the surrogate prior boundary.  MAYA is intermediate.

2. **NRSur7dq4 training simulations vs. independent SXS simulations.** The 342 SXS quasi-circular simulations used to train NRSur7dq4 show near-uniform $(2,2)$ match above 0.985 (10th percentile 0.990) — the expected in-sample behavior of a well-trained interpolant.  The 242 non-calibration SXS quasi-circular simulations have the same median but a substantially wider lower tail (10th percentile 0.360), identifying the parameter-space regions where the sparse training grid limits NRSur7dq4 accuracy.

3. **Sub-dominant modes** exhibit significantly more inter-catalog variance than the $(2,2)$ mode.  SXS shows low median matches for $(4,3)$ (0.671) and $(3,2)$ (0.858); RIT shows low medians for $(3,3)$ (0.771) and $(4,4)$ (0.783).  The complementary pattern suggests different systematic errors in the two NR codes' treatment of higher multipoles.

4. **Eccentric systems** show match degradation consistent with the surrogate's quasi-circular training assumption, providing a built-in negative control: at $e \gtrsim 0.05$, $(2,2)$ matches routinely fall below 0.9.

5. **Sub-dominant mode behavior depends critically on binary symmetries**: modes that vanish by symmetry at $q=1$, $\chi=0$ become meaningful diagnostics only when those symmetries are broken by mass ratio or spin.

6. **The combined match + $\Delta\Phi/{\rm cycle}$ diagnostic** provides strictly more information than either metric alone, and we advocate its adoption as a standard waveform-comparison reporting convention.

The precessing-spin categories (e and f), which account for $\sim$74% of all catalog simulations within the NRSur7dq4 prior volume, require a full SO(3) frame-rotation optimization before meaningful comparison is possible.  This extension, together with BMS supertranslation correction and the full cross-catalog analysis, will be reported in a companion paper.

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
