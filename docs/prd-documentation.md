# PRD: Documentation Audit & Improvement

**Status**: Draft  
**Author**: Prayush Kumar  
**Date**: 2026-05-26  
**Scope**: All code in `nrcatalogtools/` and `docs/`

---

## 1. Background and Motivation

`nr-catalog-tools` is used across gravitational-wave research pipelines and is publicly hosted at
`https://github.com/gwnrtools/nr-catalog-tools` with a MkDocs site at
`https://gwnrtools.github.io/nr-catalog-tools/`. The package has grown through several
architectural changes (waveform sub-package split, pyproject.toml migration, catalog plugin
registry, per-catalog YAML schemas) and the documentation has not kept pace uniformly.

### Current State Assessment

| Area | Coverage | Issue |
|------|----------|-------|
| `waveform/modes.py` | Excellent | No issues |
| `waveform/loaders.py` | Excellent | No issues |
| `waveform/matching.py` | Excellent | No issues |
| `waveform/units.py` | Minimal | No module docstring; `ELL_MIN`/`ELL_MAX` undocumented |
| `catalog.py` | Good | Abstract methods lack parameter docs |
| `registry.py` | Excellent | No issues |
| `metadata.py` | Good | Module-level constants undocumented |
| `sxs.py` | Good | `_add_paths_to_metadata()` undocumented |
| `rit.py` | Moderate | `RITCatalogHelper` methods almost all undocumented |
| `maya.py` | Sparse | Most methods beyond `load()` undocumented |
| `lvc.py` | Moderate | Top-level functions documented; ~20 helper functions undocumented |
| `utils.py` | Good | Most public functions documented |
| `docs/index.md` | Good | Module tree stale (references `waveform.py` not `waveform/`) |
| `docs/catalogs.md` | Good | No issues |
| `docs/waveform.md` | Good | No issues |
| `docs/architecture.md` | Good | No issues |
| `docs/package.md` | Good | Mentions `waveform.py`; should reference sub-package |
| `docs/goal.md` | Excellent | No issues |
| `mkdocs.yml` | Good | No auto-generated API reference configured |
| README.md | Good | No issues |

### Problems

1. **Stale module references**: `index.md` and `package.md` reference `waveform.py` (old
   monolithic file) instead of the `waveform/` sub-package introduced in a recent commit.

2. **Undocumented helpers in `lvc.py`**: ~20 internal helpers (e.g., `get_ref_vals()`,
   `check_interp_req()`, many `_get_*` helpers) have no docstrings. They are not private in
   name-mangling terms but are called indirectly from public functions, making them hard to
   understand.

3. **Sparse `RITCatalogHelper` docstrings**: `rit.py` has 998 lines. `RITCatalogHelper` methods
   (`_scrape_metadata()`, `_build_url()`, `_cache_metadata()`, etc.) are undocumented, making
   the scraping and caching logic opaque.

4. **Sparse `MayaCatalog` internals**: `maya.py` lacks docstrings on `reload()`, `clear_cache()`,
   `waveform_filepath_from_simname()`, `download_waveform_data()`, and several private helpers.

5. **No type annotations**: No function signatures carry type hints anywhere in the codebase,
   making editor tooling and static analysis less useful and auto-generated API docs less
   informative.

6. **No auto-generated API reference**: MkDocs is configured but does not use `mkdocstrings`.
   Docstrings are duplicated manually into `.md` files. Drift between code and docs is
   inevitable.

7. **`waveform/units.py` undocumented**: Module docstring absent; two module-level constants
   (`ELL_MIN`, `ELL_MAX`) have no explanation.

8. **No Tutorials section**: The `docs/` directory has a `tutorials/` directory but no content
   and no MkDocs nav entry.

9. **`README.md` links to docs site but no quick local preview instructions**.

10. **GitHub rendering**: Some pages use MathJax which doesn't render on raw GitHub markdown;
    equations should have plain-text fallbacks or be in a format that degrades gracefully.

---

## 2. Goals

1. Every public function, class, and method has a complete docstring (summary line, parameters,
   returns, raises, and an example where useful).
2. Every module has a module-level docstring explaining its purpose and public API.
3. Type annotations added to all public function signatures.
4. MkDocs site uses `mkdocstrings` to auto-generate API reference pages from docstrings;
   manual API docs in `.md` files are removed or reduced to prose-only summaries.
5. All stale references to the old `waveform.py` monolith are replaced with the `waveform/`
   sub-package.
6. A minimal Tutorials section (at least two end-to-end examples) is added to the docs.
7. All documentation pages render correctly both on the GitHub web UI and via MkDocs.

### Non-Goals

- Restructuring the Python package layout or changing the public API.
- Adding tutorials for optional dependencies (scri, gwsurrogate) beyond brief notes.
- Porting docs to Sphinx.

---

## 3. Requirements

### REQ-DOC-1: Module Docstrings

Every Python file in `nrcatalogtools/` must have a module-level docstring describing:
- What the module contains
- Key public classes or functions it exports

**Files needing new or improved module docstrings:**

| File | Action |
|------|--------|
| `nrcatalogtools/__init__.py` | Add module docstring listing all public exports |
| `nrcatalogtools/waveform/__init__.py` | Add module docstring for sub-package |
| `nrcatalogtools/waveform/units.py` | Add module docstring; document `ELL_MIN`, `ELL_MAX` |
| `nrcatalogtools/lvc.py` | Strengthen existing module docstring with section overview |
| `nrcatalogtools/rit.py` | Add module docstring |
| `nrcatalogtools/maya.py` | Add module docstring |

---

### REQ-DOC-2: Function and Method Docstrings

All public methods and functions must have Google-style docstrings (matching the style
already used in `waveform/modes.py` and `registry.py`) with:
- One-line summary
- `Args:` section listing each parameter with type and description
- `Returns:` section
- `Raises:` section for expected exceptions
- `Example:` block for non-trivial functions

**Specific targets (functions with no or stub docstrings):**

#### `lvc.py` — approximately 20 functions

| Function | Priority |
|----------|----------|
| `get_ref_vals()` | High |
| `check_interp_req()` | High |
| `get_nr_to_lal_rotation_angles()` | High |
| All `_get_*` helpers | Medium |
| `rotate_modes()` | High |
| `compute_mismatch_*()` | High |

#### `rit.py` — `RITCatalogHelper` methods

| Method | Priority |
|--------|----------|
| `_scrape_metadata()` | High |
| `_build_url()` / `_build_waveform_url()` | Medium |
| `_cache_metadata()` | Medium |
| `_load_cached_metadata()` | Medium |
| All filename-generation methods | Medium |

#### `maya.py` — `MayaCatalog` methods

| Method | Priority |
|--------|----------|
| `reload()` | High |
| `clear_cache()` | High |
| `waveform_filepath_from_simname()` | High |
| `download_waveform_data()` | High |
| `_build_metadata_dict()` | Medium |

#### `catalog.py` — Abstract base methods

| Method | Priority |
|--------|----------|
| All `@abstractmethod` declarations | High |
| `set_attribute_in_waveform_data_file()` | Medium |

---

### REQ-DOC-3: Type Annotations

Add PEP 484 type annotations to all public function signatures across the package.
Internal/private helpers may use type comments if full annotation is too verbose.

**Priority order:**

1. `catalog.py` — public interface (highest priority; users inherit from this)
2. `waveform/modes.py` — central user-facing class
3. `sxs.py`, `rit.py`, `maya.py` — public `load()`, `get()`, `get_metadata()`, `get_parameters()`
4. `metadata.py` — `get_source_parameters_from_metadata()`
5. `registry.py` — already well-documented; add type hints
6. `utils.py` — public helpers
7. `lvc.py` — public functions first, helpers after
8. `waveform/loaders.py`, `waveform/matching.py`, `waveform/units.py`

**Style**: Use `from __future__ import annotations` for forward references. Do not import
types only needed at runtime. Use `Optional[X]` rather than `X | None` for Python 3.8
compatibility (the package supports Python ≥ 3.8).

---

### REQ-DOC-4: Auto-Generated API Reference with mkdocstrings

Install and configure `mkdocstrings[python]` to generate API reference pages directly from
docstrings.

#### 4a: Install dependency

Add to `pyproject.toml` optional dependencies (docs group):

```toml
[project.optional-dependencies]
docs = [
    "mkdocs-material",
    "mkdocstrings[python]",
    "mkdocs-gen-files",
    "mkdocs-literate-nav",
]
```

#### 4b: Configure `mkdocs.yml`

```yaml
plugins:
  - search
  - mkdocstrings:
      handlers:
        python:
          options:
            docstring_style: google
            show_source: true
            show_root_heading: true
            members_order: source
```

#### 4c: Add API reference pages

Create `docs/api/` directory with one page per module:

```
docs/api/
├── index.md         # Overview linking to each module page
├── catalog.md       # ::: nrcatalogtools.catalog
├── rit.md           # ::: nrcatalogtools.rit
├── sxs.md           # ::: nrcatalogtools.sxs
├── maya.md          # ::: nrcatalogtools.maya
├── waveform.md      # ::: nrcatalogtools.waveform (sub-package)
├── metadata.md      # ::: nrcatalogtools.metadata
├── registry.md      # ::: nrcatalogtools.registry
├── lvc.md           # ::: nrcatalogtools.lvc
└── utils.md         # ::: nrcatalogtools.utils
```

Each page uses mkdocstrings directives like:

```markdown
## CatalogBase

::: nrcatalogtools.catalog.CatalogBase
    options:
      members:
        - load
        - get
        - get_metadata
        - get_parameters
        - simulations_list
        - simulations_dataframe
```

#### 4d: Update `mkdocs.yml` nav

```yaml
nav:
  - Home: index.md
  - Quick Start: quickstart.md
  - Catalog Reference: catalogs.md
  - API Reference:
      - Overview: api/index.md
      - catalog: api/catalog.md
      - rit: api/rit.md
      - sxs: api/sxs.md
      - maya: api/maya.md
      - waveform: api/waveform.md
      - metadata: api/metadata.md
      - registry: api/registry.md
      - lvc: api/lvc.md
      - utils: api/utils.md
  - Tutorials:
      - Loading a waveform: tutorials/load-waveform.md
      - Cross-catalog comparison: tutorials/cross-catalog.md
  - Architecture: architecture.md
  - Package Internals: package.md
  - Scientific Goal: goal.md
```

**Note**: The existing `docs/waveform.md` hand-written API reference should be converted to
prose-only (concepts, conventions, examples) and renamed `docs/waveform-guide.md`. The
auto-generated `docs/api/waveform.md` replaces the hand-maintained API table.

---

### REQ-DOC-5: Fix Stale References

The following stale references to the old `waveform.py` monolith must be updated:

| File | Location | Change |
|------|----------|--------|
| `docs/index.md` | Module tree (line ~59) | `waveform.py` → `waveform/` sub-package |
| `docs/package.md` | Module descriptions section | Update `waveform.py` references throughout |
| `docs/index.md` | Dependency table | Verify all dependency versions are current |

---

### REQ-DOC-6: Tutorials Section

Add `docs/tutorials/` with at least two complete, runnable end-to-end tutorials:

#### Tutorial 1: `load-waveform.md` — Loading and Plotting a Waveform

Cover:
1. Installing the package and dependencies
2. Loading the RIT catalog (`RITCatalog.load()`)
3. Picking a simulation and inspecting metadata
4. Loading the waveform (`catalog.get()`)
5. Getting the physical (2,2) mode (`wfm.get_mode(2, 2, ...)`)
6. Getting polarizations (`wfm.get_td_waveform(...)`)
7. Basic matplotlib plot of h+ vs time

#### Tutorial 2: `cross-catalog.md` — Cross-Catalog Comparison

Cover:
1. Loading the same physical system from RIT and SXS
2. Aligning time and phase at the (2,2) peak
3. Computing the noise-weighted mismatch
4. Interpreting the result in the context of NR accuracy

---

### REQ-DOC-7: GitHub-Readable Markdown

The documentation must render meaningfully on GitHub's markdown renderer, not just via MkDocs.

**Rules:**

1. **Math equations**: Every display equation in `.md` files must have a plain-text or
   AsciiMath fallback comment below it using HTML comment syntax:
   ```markdown
   $$\mathcal{M}(h_1, h_2) = 1 - \max_{t_c, \phi_c} \frac{\langle h_1 | h_2 \rangle}{\sqrt{\langle h_1|h_1\rangle\langle h_2|h_2\rangle}}$$
   <!-- Mismatch = 1 - max over (time shift, phase shift) of normalized inner product -->
   ```
   This ensures `goal.md` (heavy with LaTeX) still communicates on GitHub.

2. **Internal links**: Use relative links (e.g., `[Architecture](architecture.md)`) not
   absolute MkDocs site URLs. Relative links work both on GitHub and in MkDocs.

3. **Admonitions**: Replace all MkDocs `!!! note` admonitions with plain `> **Note:**`
   blockquotes **in files that are primarily reference material** (README, index.md). Keep
   admonitions in tutorial and guide pages where MkDocs rendering is expected.

4. **Code blocks**: Always specify the language identifier after triple backticks
   (e.g., ` ```python `, ` ```bash `, ` ```toml `). GitHub and MkDocs both support this.

5. **README.md**: Add a local preview section:
   ```markdown
   ## Building the docs locally
   pip install -e ".[docs]"
   mkdocs serve
   # Open http://127.0.0.1:8000
   ```

---

### REQ-DOC-8: CI — Docs Build Verification

The existing `docs.yml` workflow must be updated to:

1. Install `mkdocstrings[python]` and the docs extras group.
2. Run `mkdocs build --strict` (fails on warnings, including broken links).
3. Deploy to GitHub Pages on merges to `master`.

```yaml
# .github/workflows/docs.yml (key additions)
- name: Install docs dependencies
  run: pip install -e ".[docs]"

- name: Build docs (strict)
  run: mkdocs build --strict
```

---

## 4. Out-of-Scope Items (Deferred)

| Item | Reason |
|------|--------|
| Sphinx migration | Too large; MkDocs is already deployed and working |
| Doctest CI integration | Complex setup; can follow in a separate PR |
| Changelog / `CHANGELOG.md` | Separate effort |
| API versioning / deprecation notices | Separate effort |
| `scri` / `gwsurrogate` tutorial | Optional deps not installed in CI |

---

## 5. Implementation Plan

Work is organized into three phases. Each phase corresponds to one PR.

### Phase 1: Code Docstrings and Type Hints (PR #1)

**Goal**: Every public function and method has a complete docstring and type annotations.  
**Scope**: Pure Python changes; no docs site changes.

| Step | Files | Est. Changes |
|------|-------|-------------|
| 1.1 | `nrcatalogtools/waveform/units.py` | Module docstring, annotate constants |
| 1.2 | `nrcatalogtools/catalog.py` | Type hints + parameter docs for abstract methods |
| 1.3 | `nrcatalogtools/maya.py` | Module docstring; docstrings + type hints for all methods |
| 1.4 | `nrcatalogtools/rit.py` | `RITCatalogHelper` method docstrings; type hints on public API |
| 1.5 | `nrcatalogtools/lvc.py` | Docstrings for ~20 undocumented helpers; type hints on public fns |
| 1.6 | `nrcatalogtools/utils.py` | Type hints on all public functions |
| 1.7 | `nrcatalogtools/metadata.py` | Document module-level constants; type hints |
| 1.8 | `nrcatalogtools/registry.py` | Type hints (already well-documented) |
| 1.9 | `nrcatalogtools/__init__.py` | Module docstring listing public API |
| 1.10 | `nrcatalogtools/waveform/__init__.py` | Module docstring for sub-package |
| 1.11 | `nrcatalogtools/sxs.py` | Docstring for `_add_paths_to_metadata()`; type hints |

**Acceptance criteria:**
- `pydocstyle --convention=google nrcatalogtools/` passes with zero errors.
- `mypy nrcatalogtools/ --ignore-missing-imports` passes on public APIs (modules may have
  `# type: ignore` on third-party calls that are hard to annotate).

---

### Phase 2: MkDocs Auto-Generated API Reference (PR #2)

**Goal**: API docs auto-generated from docstrings; no hand-maintained API tables.  
**Prerequisite**: Phase 1 complete.

| Step | Files | Action |
|------|-------|--------|
| 2.1 | `pyproject.toml` | Add `[docs]` optional dependencies group |
| 2.2 | `mkdocs.yml` | Add `mkdocstrings` plugin config |
| 2.3 | `docs/api/*.md` | Create one page per module with `:::` directives |
| 2.4 | `mkdocs.yml` | Update nav with API Reference and Tutorials sections |
| 2.5 | `docs/waveform.md` | Convert to `waveform-guide.md` (prose only); remove API table |
| 2.6 | `.github/workflows/docs.yml` | Install docs extras; use `mkdocs build --strict` |

**Acceptance criteria:**
- `mkdocs build --strict` completes with zero warnings or errors.
- Every public class and function in `nrcatalogtools/` appears in at least one API page.
- No broken internal links reported by mkdocs.

---

### Phase 3: Content Fixes and Tutorials (PR #3)

**Goal**: Prose doc accuracy, GitHub rendering, tutorials.

| Step | Files | Action |
|------|-------|--------|
| 3.1 | `docs/index.md` | Fix module tree (`waveform/` sub-package); verify dependency table |
| 3.2 | `docs/package.md` | Update all `waveform.py` references to `waveform/` |
| 3.3 | `docs/goal.md` | Add plain-text fallbacks below each display equation |
| 3.4 | `README.md` | Add "Building the docs locally" section |
| 3.5 | `docs/tutorials/load-waveform.md` | Create Tutorial 1 (load and plot) |
| 3.6 | `docs/tutorials/cross-catalog.md` | Create Tutorial 2 (cross-catalog comparison) |
| 3.7 | All `.md` files | Audit: relative links, code block language tags, admonition style |

**Acceptance criteria:**
- All links in `docs/` resolve on GitHub (manual check or `markdown-link-check`).
- Tutorial code blocks are tested by copy-pasting in a conda environment with all deps.
- `mkdocs build --strict` still passes after Phase 3 changes.

---

## 6. Style Guide

All docstrings must follow Google style as used in `waveform/modes.py`. Key conventions:

```python
def example(param1: int, param2: str = "default") -> list[str]:
    """One-line summary ending in a period.

    Longer description if needed. Can span multiple paragraphs.

    Args:
        param1: Description of param1.
        param2: Description of param2. Defaults to ``"default"``.

    Returns:
        Description of return value.

    Raises:
        ValueError: If param1 is negative.
        FileNotFoundError: If the data file does not exist.

    Example:
        >>> result = example(3, "hello")
        >>> len(result)
        3
    """
```

**Additional rules:**

- Use `Args:` (not `Parameters:`) to match the existing convention in `modes.py`.
- Physical units always stated in parameter description (e.g., `total_mass (float): Total mass in solar masses.`).
- Dimensionless vs physical quantities always distinguished.
- Cross-references to related functions use `` `function_name` `` (backtick, no module prefix
  needed within the same module).

---

## 7. Success Metrics

| Metric | Target |
|--------|--------|
| Public functions with complete docstrings | 100% |
| Public functions with type annotations | 100% |
| `pydocstyle` errors | 0 |
| `mkdocs build --strict` warnings | 0 |
| Broken internal links | 0 |
| Tutorial pages | ≥ 2 complete, tested |
| `goal.md` equations with fallbacks | 100% |
