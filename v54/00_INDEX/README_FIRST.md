# v53 bundle - start here

This bundle is append-only relative to sealed v52. All v52 root members remain at their historical paths and byte-exact. v53 adds a navigation layer; mirrored copies are conveniences, not replacements.

- `01_DOCS/` - orientation and current-document pointers.
- `02_VERSIONS/vNN/` - physical PDF/TEX mirrors by cumulative version.
- `03_SCRIPTS/vNN/` and `shared/` - executable/verifier/build scripts.
- `04_CERTIFICATS/vNN/` - certificates, reports, audits, manifests, tokens, ledgers.
- `05_SELLING/` - Selling-specific indexes and v53 typed Fig.2-6 functor.
- `06_CHIRES/` - CHIRES/Gamma3 indexes and canonical v53 obstruction data.
- `07_QA_SEAL/` - PDF QA, visual contact sheets, seal/postseal records.
- `08_FORMAL_CONTEXT/` - formal-context index.
- `09_ASSETS/` - heavy-asset index; heavy historical files remain single-copy at root.
- `00_INDEX/` - catalogs and bundle-tree audit.

Use `00_INDEX/ROOT_LEGACY_CATALOG.tsv` to locate original root members and `00_INDEX/NAVIGATION_MIRROR_CATALOG.tsv` for organized mirrors.
