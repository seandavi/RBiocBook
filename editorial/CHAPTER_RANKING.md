# Chapter harmonization ranking

Worst-first priority list for the pedagogical harmonization effort. Scores are the
**mean of the 7-dimension 1–5 rubric** in [`STYLE_GUIDE.md`](STYLE_GUIDE.md) §4
(5 = matches the class-tested intro chapters). Full per-dimension scorecards for
the three pilot chapters are in [`PILOT_REVIEWS.md`](PILOT_REVIEWS.md).

*Scored 2026-06-01 by a panel of parallel reviewers, calibrated to the three
pilots (`dataframes_intro` 3.3, `t-stats-and-tests` 2.4, `bioc-summarizedexperiment`
2.0). A snapshot — re-score if chapters change substantially.*

## Done (merged to main)

| Chapter | Was | Now (self-est.) |
|---|:--:|:--:|
| `bioc-summarizedexperiment` | 2.0 | harmonized |
| `t-stats-and-tests` | 2.4 | harmonized |
| `norm` | (cold) | harmonized |
| `factors` | 1.7 | ~4.4 |
| `vectors` | 2.4 | ~4.5 |
| `matrices` | 2.4 | ~4.4 |
| `lists` | 3.1 | ~4.4 |
| `machine_learning/models` | 1.4 | ~4.7 |
| `visualization_guide` | 1.6 | ~4.7 |
| `ranges_and_signals` | 1.6 | ~4.3 |
| `intro` | 1.9 | ~4.7 |
| `kmeans` | 1.9 | ~4.5 |
| `eda_and_univariate_brfss` | 2.0 | ~4.5 |
| `atac-seq/atac-seq` | 1.6 | ~4.5 |
| `geoquery` | 1.7 | ~4.4 |
| `genomic_ranges_tutorial` | 2.1 | ~4.6 |

## Remaining — worst first

Centrality = teaching centrality (how core to the learning path). Priority ≈
`(5 − mean) × centrality`.

| Mean | Chapter | Centrality | Biggest problems |
|:--:|---|:--:|---|
| 2.1 | `machine_learning/mlr3verse_intro` | med | No frontmatter; Palmer Penguins (no biology); outputs uninterpreted |
| 2.3 | `reading_and_writing` | high | No objectives/exercises; non-bio examples; `df` variable name |
| 2.3 | `machine_learning/intro` | high | Biology motivation = 2 citations then generic spam/housing examples |
| 2.3 | `310_microbiome` | med | Installs a GitHub-HEAD package (fragile); ends on a link dump |
| 2.4 | `dplyr_intro_msleep` | high | msleep has no genomics motivation; duplicate objectives headers |
| 2.4 | `ranges_exercises` | med | Global `results='hide'`; solutions are code-only |
| 2.6 | `single_cell/setup` | high | Malformed callout `{.callout-tip)`; no objectives/exercises |
| 2.7 | `ggplot2/intro` | high | Zero exercises in an all-worked-example chapter; title as `#` |
| 2.9 | `machine_learning/machine_learning_mlr3` | high | Pervasive `=`; a truncated sentence; strong biology though |
| 3.3 | `dataframes_intro` | high | Pilot — not yet harmonized, but has a ready-to-paste draft in `PILOT_REVIEWS.md` (**quick win**) |

## Natural batches

- **ML part** — `mlr3verse_intro` (2.1), `ml/intro` (2.3),
  `machine_learning_mlr3` (2.9): `models` (1.4) is now harmonized; the rest remain.
  This is now the largest unharmonized cluster and the natural next batch.
- **Ranges / Bioc** — `ranges_exercises` (2.4) is the only one left;
  `ranges_and_signals`, `genomic_ranges_tutorial`, `atac-seq`, and `geoquery` are
  all harmonized.
- **High-leverage singletons** — done: `intro`, `kmeans`, `eda_and_univariate_brfss`.
  Remaining high-centrality targets: `reading_and_writing` (2.3),
  `dplyr_intro_msleep` (2.4), `ggplot2/intro` (2.7), `single_cell/setup` (2.6),
  and the quick-win pilot `dataframes_intro` (3.3, draft ready in `PILOT_REVIEWS.md`).

## Process reminders (from `STYLE_GUIDE.md` / `PILOT_REVIEWS.md`)

- One author per chapter (don't split a chapter's prose — it fractures voice);
  parallelize *across* chapters. Anchor voice to `r_basics.qmd`, structure to an
  already-harmonized chapter (e.g. `norm.qmd`).
- Before editing a chapter, run `Rscript tools/check-deps.R` and add any missing
  package to `DESCRIPTION` (editing a chapter forces CI re-execution that needs
  its deps).
- Audience is **not** very stats-savvy — bias toward accessible prose; prefer
  inline `r ...` over hardcoded numbers; seed random chunks for reproducible
  freeze.
