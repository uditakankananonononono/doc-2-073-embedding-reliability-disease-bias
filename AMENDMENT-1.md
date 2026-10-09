# AMENDMENT-1 (2026-10-09 IST) - made before any outcome access

Reason: free Colab T4 refused ("Cannot connect to GPU backend: usage limits"). No labels, DMS scores or ESM-2 outputs have been read.
Change 1: primary scoring is ESM-2 wild-type-marginals (log p(mut) - log p(wt) from one unmasked forward pass per protein), not masked-marginals. This runs on free Colab CPU.
Change 2: model is ESM-2 650M if it fits CPU RAM, else ESM-2 150M; the choice is recorded in DATA_MANIFEST.tsv before labels are read. If a T4 becomes available, masked-marginal 650M is run as a secondary analysis only.
All gates, splits and the failure policy in PROTOCOL.md are unchanged.
