# AMENDMENT-2 (2026-10-09 IST) - pre-outcome; supersedes the lines named below

Applies the independent gate's three required fixes. Still no labels, DMS scores or ESM-2 outputs read.

## A. Lines of PROTOCOL.md that are SUPERSEDED (and by what)
1. "Model: ESM-2 650M, masked-marginal scoring, fp16 on Colab T4" is replaced by: ESM-2 wild-type-marginal scoring, one unmasked fp32 forward pass per window, on Colab CPU. fp16 and bf16 are not used.
2. "Outcome per assay: Reliability = Spearman rho between ESM-2 masked-marginal score and DMS fitness score" is replaced by: Reliability = Spearman rho between the wt-marginal score s(variant) = log p(mut aa | wt context) - log p(wt aa | wt context) and the DMS fitness score. The estimand is the reliability of wt-marginal ESM-2 scoring. It is related to, but not identical to, masked-marginal reliability, and every result will be worded that way.
3. "WT pseudo-perplexity" in G3 and in the mediator check is replaced by the feature WTLL, defined in B. True pseudo-perplexity (L masked passes) is NOT used anywhere.

## B. Feature WTLL (exact definition)
WTLL(protein) = mean over scored positions of log p(wt aa_i | wt sequence, unmasked), from the same single forward pass used for scoring. It is a different quantity from pseudo-perplexity and is inflated by the unmasked wt token being visible; it is used only as a covariate and mediator, in G3 and the mediator check. Same model and windows as scoring.

## C. Windowing for sequences > 1022 aa
Windows of 1022 residues, stride 511, one pass per window. For each position, use the log-probabilities from the window in which it is farthest from a window edge (ties: the earlier window). The wt-marginal score and WTLL use only these selected values. Sequences <= 1022 use one window.

## D. Single model choice (objective, same for ALL assays)
Rule: measure peak resident memory of one fp32 forward pass of ESM-2 650M on a random 1022-residue dummy sequence in the Colab CPU runtime used for the run. If peak RSS <= 9 GB, ESM-2 650M is used for every assay; otherwise ESM-2 150M is used for every assay. The measurement, the runtime RAM and the chosen model are written to DATA_MANIFEST.tsv before any label is read. No per-protein or per-assay model choice, ever. A later T4 masked-marginal run is a separate exploratory stratum, never swapped in.

## E. External cutoff check (done before the lock commit)
ESM-2 pretraining used UniRef built from the UniProt 2021_04 release (ESM maintainers, https://github.com/facebookresearch/esm/discussions/388). The 2022-01-01 external cutoff therefore postdates the training data release, so the external DMS assays were published after the model's training data. Caveat: protein sequences can predate the cutoff; the external test is of temporal DMS labels, not of unseen sequences.

## F. Lock commit
Before any label or score is read: scoring code, DATA_MANIFEST.tsv (versions, sha256, model choice, RAM measurement) and this amendment are committed together and tagged `lock-1`. Labels are fetched only after the tag exists.

All gates G1-G3, splits and the failure policy in PROTOCOL.md are otherwise unchanged.
