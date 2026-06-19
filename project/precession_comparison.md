cat > /home/prayush/src/nr-catalog-tools/project/precessing_comparison.md << 'MDEOF'
# Pipeline walkthrough: precessing SXS simulation vs NRSur7dq4

This document traces the full execution path of
`project/scripts/compare_one_sim_vs_surrogate.py` for a precessing SXS binary,
with mathematical detail and critical evaluation of each step.

The CLI script is a thin shim; all logic lives in
`nrcatalogtools/comparisons.py:compare_sim_vs_surrogate`.

---

## Step 1 — Load catalog and waveform

```python
cat = load_catalog("SXS")
wfm = cat.get(sim_name)
```

`wfm` is a `WaveformModes` object wrapping the SXS `rhOverM` strain data.
Internally the SXS python library delivers modes $h_{\ell m}(t) \cdot r / M$ in
**dimensionless NR units** ($t$ in units of $M$, $r$ in $M$).  No physical
scaling happens at this point.

---

## Step 2 — Parameter extraction

Source: `nrcatalogtools/metadata.py:get_source_parameters_from_metadata`

For SXS the code reads the `reference_time` epoch, not `relaxation_time`.
The SXS `reference_time` is chosen to correspond to a specific GW frequency,
making it the more physically meaningful spin epoch.

```python
q     = metadata["reference_mass_ratio"]
spin1 = metadata["reference_dimensionless_spin1"]   # 3-vector, inertial frame
spin2 = metadata["reference_dimensionless_spin2"]
Momega = norm(metadata["reference_orbital_frequency"])   # |M Ω_orb|, dimensionless
f_lower = Momega / π / (M_total × MTSUN_SI)             # [Hz]
```

`reference_orbital_frequency` is the orbital angular velocity 3-vector whose
magnitude is the dimensionless $M\Omega_\text{orb}$ at the reference epoch.
Dividing by $\pi$ converts to the (2,2) GW frequency:

$$f_\text{lower} = \frac{M\Omega_\text{orb}}{\pi \cdot M_\text{seconds}}
                 = \frac{\Omega_\text{orb}}{\pi}\ \text{[Hz]}
                 = f_{22}$$

since $f_{22} = 2 f_\text{orb} = \Omega_\text{orb}/\pi$.
**`f_lower` is therefore the (2,2) GW frequency at the reference epoch,
not the orbital frequency.**

---

## Step 3 — Epoch-aligned spin extraction (precessing systems only)

Source: `nrcatalogtools/surrogate.py:_epoch_align_spins`

For aligned-spin or non-spinning binaries the spin vectors are constant, so
the metadata values are valid at any epoch.  For precessing systems the
instantaneous spin direction rotates with the orbital plane and must be
extracted from the NR horizon data at the epoch that matches the surrogate's
`f_ref`.

### A. Find the reference time

```python
t_start_M = abs(sur._sur_dimless.t_0)        # queried from surrogate at runtime
t_target  = strain.max_norm_time() - t_start_M   # NR coordinate time
t_ref     = clip(t_target, t_h[0], t_h[-1])      # clamped to horizon range
```

`sur._sur_dimless.t_0 ≈ −4300 M` is NRSur7dq4's actual training window start.
For NR simulations shorter than ~4300 M, `t_target` falls before the
simulation start and `t_ref = t_h[0]` (first horizon sample).

### B. Identify heavier / lighter body

The sxs library labels horizons `A` and `B` but does not guarantee any
mass ordering.  gwsurrogate requires
$\hat{n} = \text{lighter} \to \text{heavier}$ (body2→body1).
The code reads the masses at the reference epoch and orders accordingly:

```python
mA = h.A.mass[idx_h];  mB = h.B.mass[idx_h]
if mA >= mB:
    chi_primary, r_primary     = h.A.chi_inertial, h.A.coord_center_inertial
    chi_secondary, r_secondary = h.B.chi_inertial, h.B.coord_center_inertial
else:
    chi_primary, r_primary     = h.B.chi_inertial, h.B.coord_center_inertial
    chi_secondary, r_secondary = h.A.chi_inertial, h.A.coord_center_inertial
```

### C. Build the coprecessing frame

$$\hat{n} = \frac{\vec{r}_\text{primary} - \vec{r}_\text{secondary}}
                 {|\vec{r}_\text{primary} - \vec{r}_\text{secondary}|}$$

$$\dot{\vec{r}}_\text{sep} \approx
  \frac{\vec{r}_\text{sep}(t+\delta t) - \vec{r}_\text{sep}(t-\delta t)}{2\delta t}
  \quad\text{(central difference)}$$

$$\hat{L} = \frac{\vec{r}_\text{sep} \times \dot{\vec{r}}_\text{sep}}
                 {|\vec{r}_\text{sep} \times \dot{\vec{r}}_\text{sep}|}$$

$$\hat{\lambda} = \hat{L} \times \hat{n}$$

The rotation matrix $R$ maps from the inertial frame to the coprecessing frame:

$$R = \begin{pmatrix} \hat{n} \\ \hat{\lambda} \\ \hat{L} \end{pmatrix}
\quad\Rightarrow\quad
\vec{\chi}_\text{coprecessing} = R\,\vec{\chi}_\text{inertial}$$

giving $(\chi_x, \chi_y, \chi_z) = (\chi\cdot\hat{n},\; \chi\cdot(\hat{L}\times\hat{n}),\; \chi\cdot\hat{L})$.

This matches gwsurrogate's documented convention exactly
(`DynamicsSurrogate` docstring: `χ_y = χ · (L̂ × n̂)`).

### D. GW frequency at the epoch

```python
omega22 = d(unwrapped_phase_22)/dt   # central difference in NR time
f_gw    = |omega22| / (2π)           # cycles/M  (dimensionless, = M f_22)
```

This `f_ref_dimless` is passed to the surrogate as the spin-definition epoch.

---

## Step 4 — Surrogate call

```python
t_sur, h_sur, _ = sur(q, chiA, chiB, ellMax=4,
                      dt=dt_dimless, f_low=0, f_ref=f_ref_dimless)
```

- `f_low=0` — return the full waveform from $t \approx -4300\,M$ to ringdown
- `f_ref=f_ref_dimless` — spin-definition epoch as $Mf_{22}$ (dimensionless)
- `dt=dt_dimless = \Delta t_\text{s} / M_\text{s}` — dimensionless time step
- Output `h_sur[(l,m)]` in geometric units ($h \cdot r / M$)

**Amplitude scaling:**

$$h^\text{phys}_{lm}(t) = h^\text{NR}_{lm}(t) \cdot \frac{GM_\text{tot}/c^2}{D_L}$$

**Epoch alignment:** both surrogate and NR waveforms have their `epoch` set so
that `t = 0` corresponds to the peak of $|h_{22}|$.  This makes their absolute
time stamps directly comparable in the matching step.

**`f_lower_effective`:** computed from the actual first-sample frequency of the
surrogate (2,2) mode — the real starting GW frequency of the surrogate output,
which may be higher than the NR `f_lower` for high-mass or high-q systems.

---

## Step 5 — Per-mode match

Source: `nrcatalogtools/waveform/matching.py:compute_mode_match`

For each $(l, m) \in$ `NR_MODES`:

### A. Mode frequency cutoff

```python
f_lower_match  = max(f_lower_NR, f_lower_effective)   # common lower bound
f_lower_mode   = f_lower_match × |m| / 2              # per-mode cutoff
```

Since `f_lower_match` is $f_{22}$ and GW frequency scales as $|m| f_\text{orb}$:

| mode | $f_\text{mode}$ |
|------|----------------|
| $(2,2)$ | $f_{22}$ |
| $(2,1)$ | $f_{22}/2$ |
| $(3,3)$ | $3f_{22}/2$ |
| $(4,4)$ | $2f_{22}$ |

### B. Intersect time windows and taper

```python
t_start = max(h_nr.start_time, h_sur.start_time)
t_end   = min(h_nr.end_time,   h_sur.end_time)
```

Both waveforms are tapered (Tukey window $\alpha=0.2$, i.e. 10% taper on each
end) then zero-padded to the next power of two.  The PSD is built at the
matching frequency resolution $\Delta f = 1/(N_\text{FFT}\,\Delta t)$.

**Only the real part** ($h_+$ component) is used:
```python
h_nr  = h_nr_complex.real()
h_sur_mode = h_sur[(l,m)].real()
```
This is correct: `pycbc.filter.match()` maximises over an overall phase
$e^{i\phi}$, so comparing $\text{Re}(h_{lm})$ is equivalent to the full complex
match — the phase maximisation recovers the $\pi/2$ degree of freedom.

### C. Noise-weighted match

$$\mathcal{F}(h_\text{NR}, h_\text{sur}) =
  \max_{\Delta t,\,\phi}
  \frac{\langle h_\text{NR},\; h_\text{sur}\,e^{i\phi} \rangle(\Delta t)}
       {\sqrt{\langle h_\text{NR}, h_\text{NR}\rangle
              \langle h_\text{sur}, h_\text{sur}\rangle}}$$

where
$\langle a, b\rangle = 4\,\text{Re}\!\int_{f_\text{low}}^{f_\text{Nyq}} \tilde{a}^*(f)\,\tilde{b}(f)/S_n(f)\,df$
using `aLIGOZeroDetHighPower` as $S_n(f)$.

**For precessing systems:** each mode $h_{lm}$ is frame-dependent.  In the
standard (non-`--rotate`) path the comparison is done in a fixed NR frame
without optimising over SO(3) rotations between the two frames.  Any residual
frame mismatch will suppress per-mode matches.  Re-run with `--rotate` to apply
the full SO(3)-optimised match.

---

## Step 6 — Phase drift metric

Source: `nrcatalogtools/waveform/matching.py:compute_phase_diff_per_cycle`

```python
phi_nr  = unwrap(angle(h_nr_complex))    # over common window
phi_sur = unwrap(angle(h_sur[(l,m)]))

ΔΦ_NR  = |phi_nr[-1]  - phi_nr[0]|
ΔΦ_sur = |phi_sur[-1] - phi_sur[0]|
N_cyc  = ΔΦ_NR / (2π)

phase_diff_per_cycle = |ΔΦ_NR − ΔΦ_sur| / N_cyc   [rad/cycle]
```

Taking differences *within* each waveform removes any constant initial-phase
offset, so the result is independent of coalescence-phase convention.  The
metric measures the **cycle-count error** over the common window: how much do
the two waveforms disagree in the total number of GW cycles they accumulate?

This is distinct from the match metric, which uses the match()-optimal time
shift.  The phase drift uses the fixed, peak-aligned time axis (both waveforms
have $t=0$ at the (2,2) peak), which is the physically natural alignment for
comparing inspiral phasing.

---

## Known limitations and open issues

| Item | Status | Details |
|------|--------|---------|
| Frame optimisation for precessing modes | Addressed by `--rotate` | Standard path gives lower bound on match; SO(3)-optimised path is correct |
| `_epoch_align_spins` clamps to horizon start for short NR runs | Documented | For NR waveforms shorter than $\sim4300\,M$ the epoch defaults to the first horizon sample |
| Central-difference accuracy for $\hat{L}$ | Minor | Depends on horizon time-step density; fine for well-sampled horizon data |
| m=0 modes | Not in `NR_MODES` | `mode_f_lower` returns `f_lower` as a conservative bound; physically undefined for memory modes |
| Taper may be short for very high-mass HOM modes | Low risk | Tukey $\alpha=0.2$ (10% each end); revisit if common windows are $< 0.1$ s |
MDEOF