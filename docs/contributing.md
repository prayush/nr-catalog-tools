---
title: Contributing
nav_order: 11
---

# Contributing
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Development setup

```bash
git clone https://github.com/gwnrtools/nrcats.git
cd nrcats

# Option A: conda (recommended — ships lalsuite prebuilt)
conda env create -f environment.yml
conda activate nrcat

# Option B: pip into an existing environment
pip install -e .
pip install pytest coverage pre-commit
```

On Linux with pip, `lalsuite` needs system libraries first:
`sudo apt-get install libgsl-dev libfftw3-dev libhdf5-dev`.

## Code style

Formatting and linting are enforced by [pre-commit](https://pre-commit.com/)
(black + flake8 + merge-conflict checks):

```bash
pre-commit install                 # run automatically on every commit
pre-commit run --all-files         # run manually
```

CI runs the same hooks on every push and pull request.

## Running the tests

```bash
python -m pytest test/ -v -m "not cross_catalog"
```

- Tests marked `cross_catalog` require pre-downloaded waveform data from all three
  catalogs and run in a dedicated workflow (`test-cross-catalog.yml`); exclude them
  locally unless you have the data cached.
- Some tests skip gracefully when optional data or `waveformtools` is absent.

With coverage:

```bash
python -m coverage run --source=nrcats -m pytest test/ -m "not cross_catalog"
python -m coverage report
```

## Building the documentation

The docs are a [Jekyll](https://jekyllrb.com/) site using the
[just-the-docs](https://just-the-docs.com/) theme, with the API reference generated
from source docstrings by a griffe-based script.

```bash
# 1. Regenerate the API reference (only needs griffe, not the runtime deps)
pip install griffe
python bin/generate_api_docs.py

# 2. Build/serve the site
cd docs
bundle install                       # first time only; needs Ruby >= 3.0
bundle exec jekyll serve             # http://127.0.0.1:4000/nrcats/
```

Conventions:

- Pages live in `docs/` with just-the-docs front matter (`title`, `nav_order`,
  `parent`).
- `docs/api/*.md` (except `index.md`) are **generated — do not edit by hand**.
  Improve the source docstrings instead, then re-run
  `python bin/generate_api_docs.py`. CI fails if these pages are stale
  (`--check` mode).
- Write math as `$$...$$` (kramdown), both inline and display.
- Docstrings may be numpy- or Google-style; the generator auto-detects. New code
  should prefer **numpy-style**, which the codebase predominantly uses.

The site deploys automatically from `master` via `.github/workflows/docs.yml`.

## Adding a new catalog backend

New catalogs plug in via the [registry](api/registry.md) without touching core code:

1. Subclass [`CatalogBase`](api/catalog.md) and implement the `CatalogABC`
   contract (filename/filepath/URL/download methods for waveform, psi4, and
   metadata products).
2. Register it:

   ```python
   from nrcats.catalog import CatalogBase
   from nrcats.registry import register_catalog

   @register_catalog("MYCAT")
   class MyCatalog(CatalogBase):
       ...
   ```

3. If your catalog's metadata schema differs, add a key-mapping YAML under
   `nrcats/schemas/` (see `rit_keys.yaml` for the pattern) and wire it in
   [`metadata`](api/metadata.md).
4. Add the module to the `MODULES` list in `bin/generate_api_docs.py` so it appears
   in the API reference, and document catalog-specific conventions in
   [Catalog Reference](catalogs.md).

## Releases

Versioning is automatic via `setuptools_scm` from git tags. Publishing to PyPI is
triggered by creating a GitHub release (`python-publish.yml`). Update
`CHANGELOG.md` (Keep-a-Changelog format) as part of any user-facing change.
