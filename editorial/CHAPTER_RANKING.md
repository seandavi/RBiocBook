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

## Remaining — worst first

Centrality = teaching centrality (how core to the learning path). Priority ≈
`(5 − mean) × centrality`.

| Mean | Chapter | Centrality | Biggest problems |
|:--:|---|:--:|---|
| 1.6 | `atac-seq/atac-seq` | med | Broken `\@ref(fig:)` crossrefs; inline `*Exercise*` stubs; separate `atac.bib` |
| 1.7 | `geoquery` | med | Opens with a pasted journal abstract; `=` everywhere; no frontmatter title |
| 1.9 | `intro` | high | Book's **first** chapter; dry marketing register; no objectives/exercises/summary — sets the wrong tone (high leverage) |
| 1.9 | `kmeans` | high | `=` + forbidden `{r #fig-}` labels; preprocessing wall-of-text; results never interpreted (last stats-part chapter) |
| 2.0 | `eda_and_univariate_brfss` | high | Splits into a polished half + an unscaffolded slide-outline half; hidden `path` var |
| 2.1 | `genomic_ranges_tutorial` | high | Blog-post register; no callouts; no exercises |
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
- **Ranges / Bioc** — `genomic_ranges_tutorial` (2.1), `ranges_exercises` (2.4),
  `atac-seq` (1.6); `ranges_and_signals` (1.6) is now harmonized.
- **High-leverage singletons** — `intro` (the face of the book, 1.9) and `kmeans`
  (1.9, finishes the stats part). With the top three (`models`,
  `visualization_guide`, `ranges_and_signals`) done, `intro` and `kmeans` are the
  highest-leverage next picks.

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
