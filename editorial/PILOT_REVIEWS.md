# Pilot Pedagogical Reviews (3 chapters)

*Calibration pass. Three chapters chosen to span the difficulty gradient you
described — one foundational, one statistical, one advanced Bioconductor — so you
can judge the output format before we commit to all ~30 instructional chapters.*

- `dataframes_intro.qmd` — foundational data structures (early-mid)
- `t-stats-and-tests.qmd` — statistics (mid-late, conceptually demanding)
- `bioc-summarizedexperiment.qmd` — Bioconductor (late, vignette-derived)

Scored against the 7-dimension rubric in `STYLE_GUIDE.md` (1–5; 5 = matches the
proven early chapters). Drafts below are **ready to paste** and written in the
house voice; line references are current as of 2026-05-31. Correctness bugs
already cataloged in `CHAPTER_REVIEW_FINDINGS.md` are noted but not re-explained.

---

## The gradient, quantified

| Chapter | Scaffold | Design | Voice | Adult-fit | Examples | Code↔prose | Consistency | **Mean** |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| `dataframes_intro.qmd` | 4 | 3 | 3 | 4 | 3 | 4 | 2 | **3.3** |
| `t-stats-and-tests.qmd` | 2 | 2 | 2 | 3 | 3 | 3 | 2 | **2.4** |
| `bioc-summarizedexperiment.qmd` | 2 | 1 | 2 | 2 | 2 | 3 | 2 | **2.0** |

The decline is exactly as you predicted: the further from the dice/intro
material, the more the chapter reverts to textbook/vignette register and drops
the scaffolding (objectives, motivation, callouts, exercises-with-solutions) that
makes the early chapters work. **Consistency is uniformly low** (the `=`-vs-`<-`
drift and missing callouts are book-wide), which is good news — it's a
mechanical, high-leverage fix.

---

## Chapter 1 — `dataframes_intro.qmd` (mean 3.3)

**One-line:** Solid bones and a great real dataset, but it reads flatter than the
chapters around it — no callouts, a duplicated objectives header, embedded
questions instead of real exercises, and an abrupt stop.

### Per-section

| Section / line | Issue | Action |
|---|---|---|
| Top (1–12) | Strong conceptual open (data.frame vs matrix) but no *motivating* hook — why should a biologist care today? | Add a 2-sentence biological hook before the definition (draft A). |
| `## Learning goals` (15) **and** `## Learning objectives` (23) | **Two** objective blocks back-to-back — redundant and confusing; one lists "goals," the other "objectives." | Merge into one `## What you'll learn` (draft B). |
| Whole chapter | **Zero callouts** in a book that otherwise uses them heavily — the `View()` note, the `NA` gotcha, and the `$`-vs-`[[ ]]` aside are all natural callouts. | Convert 3 asides to the callout taxonomy (draft C shows the `NA` one). |
| 179–196 "Some data exploration" | Bare imperative questions ("Make histograms of…", "What does `table()` do?") with no solutions — not framed as exercises, no payoff. | Promote to a real `## Exercises` block with collapsible solutions (draft D). |
| 250–257 | Simulation uses `=` for assignment (`smoker = …`, `x = rnorm(100)`) — conflicts with the taught `<-`. Plus typo "varibles". | Convert to `<-`; fix typo. |
| End (280) | Chapter just stops after `write.csv` — no recap, no reflection. | Add `## Summary` + `## Reflection` (draft E). |
| 18, 193, 257 | typos "examing" / "look a the" / "varibles" (see findings). | Fix inline. |

### Ready-to-paste drafts

**Draft A — biological hook (insert before line 5):**
```markdown
Almost every dataset you will meet in genomics arrives as a table: genes down the
side, samples across the top, measurements in the middle, and a few columns of
annotation bolted on. In R, the structure built for exactly this shape is the
`data.frame`. Get comfortable with it and most of the rest of the book — loading
expression data, filtering genes, summarizing by condition — becomes variations
on a theme you already know.
```

**Draft B — merged objectives (replace lines 15–28):**
```markdown
## What you'll learn

By the end of this chapter you will be able to:

- Explain how a `data.frame` differs from a matrix, and why that difference matters.
- Load tabular data from a file or URL with `read.csv()`.
- Inspect a data frame with `head()`, `dim()`, `summary()`, and friends.
- Pull out columns and rows by position, name, and logical condition.
- Summarize measurements by category with `aggregate()`.
- Build a `data.frame` from scratch and save it back to disk.

We'll do all of this on a real yeast gene-expression dataset.
```

**Draft C — the `NA` gotcha as a callout (wrap lines ~207–212):**
```markdown
::: {.callout-warning}
## Logical subsetting and missing values

When you filter with a condition like `ydat$symbol == 'LEU1'`, any row whose
`symbol` is `NA` returns `NA` for the test — not `FALSE` — so R keeps it and you
get phantom rows full of `NA`. This is one of the most common surprises in R.
The fix is to exclude the missing values explicitly with `!is.na()`:

\`\`\`{r}
head(ydat[ydat$symbol == 'LEU1' & !is.na(ydat$symbol), ])
\`\`\`
:::
```

**Draft D — real exercises (replace the loose questions at 179–196):**
```markdown
## Exercises

1. Which columns of `ydat` are numeric? Check `rate` and `expression`.
2. Make a histogram of the `expression` values, then of the `rate` values.
3. `rate` has only a handful of distinct values. Use `table()` to count how many
   rows fall at each rate. Which `rate` corresponds to the most nutrient-starved
   condition?

::: {.callout-note collapse="true"}
## Solutions
\`\`\`{r}
class(ydat$rate); class(ydat$expression)   # both numeric
hist(ydat$expression)
hist(ydat$rate)
table(ydat$rate)                            # the smallest rate = most starved
\`\`\`
Growth rate is the *limiting nutrient* dialed down, so the **lowest** `rate`
(0.05) is the most starved condition.
:::
```

**Draft E — close (append at end):**
```markdown
## Summary

You can now load a real dataset into a `data.frame`, inspect its shape and
contents, pull out the rows and columns you care about (by position, name, or a
logical test), summarize values by category, and write results back to disk.
These are the moves you'll repeat on every dataset for the rest of the book.

## Reflection

- Can you state one thing a `data.frame` can do that a matrix cannot?
- Given a column `ydat$rate`, can you write the subset that keeps only rows where
  `rate < 0.1`?
- When would you reach for `aggregate()` instead of subsetting?
```

---

## Chapter 2 — `t-stats-and-tests.qmd` (mean 2.4)

**One-line:** The simulation-driven approach to the t-distribution is genuinely
excellent pedagogy buried under textbook register — no objectives, no biological
motivation, a muddled mid-chapter "summary," `=` everywhere, the `df` name
collision, and at least one wall-of-text that needs to become steps.

### Per-section

| Section / line | Issue | Action |
|---|---|---|
| `## Background` (5–26) | Opens with a dense definition of the t-test before the reader knows *why they're here*; no objectives. | Add a motivating open + `## What you'll learn` (drafts A, B). |
| Throughout | Assignment with `=` (`zscore = …`, `df = data.frame(…)`, `samp = …`); the book teaches `<-`. | Convert all to `<-`. |
| 55, 146, 164 | Data frames named `df` — collides conceptually with the `df=` (degrees of freedom) argument used at line 200, in the same chapter. | Rename to `pdat` / `tdat`. |
| 116, 131 | Broken LaTeX (mismatched braces) in the SD formulas (see findings); also no plain-language gloss of what the formula *says*. | Fix braces + add a one-line "in words" gloss (draft C). |
| 104–170 "The t-distribution" | Strong content, but it's a lot of math with no callout signposting the key intuition (estimating σ → fatter tails). | Add a `callout-important` stating the one idea to keep (draft D). |
| 203–291 "Experiment" | Best part of the chapter (learn-by-simulation), but framed as exposition, not as an invitation to *do*. The payoff sentence ("Now the p-values are uniformly distributed…") deserves emphasis. | Light reframing as a guided exercise; keep the code. |
| 293–304 | **Wall of text** — the entire qqplot explanation is one 8-sentence paragraph mixing definition, usage, and interpretation. | Break into 3 short paragraphs / steps (draft E). |
| 311–324 "Summary of t-distribution vs normal" | A full summary appears *mid-chapter*, then the chapter continues into `t.test` and power — so the reader hits "Summary" and is surprised there's more. | Rename to `## Checkpoint: t vs normal` and add a real chapter-end `## Summary`. |
| End (491–493) | Ends on a bare "Resources" link; no recap tying back, no exercises. | Add `## Summary` + `## Exercises` with solutions. |
| 288, 459 | `hist(pt(t10k,5))` df should be 4; "is now easily" → "not" (see findings). | Fix inline. |

### Ready-to-paste drafts

**Draft A — motivating open (insert at top, before `## Background`):**
```markdown
Suppose you measure a gene's expression in five treated samples and five
controls, and the treated mean is higher. Is that a real effect, or the kind of
gap you'd see by chance from just five noisy measurements? With only a handful of
samples you can't appeal to the normal distribution as-is, because you had to
*estimate* the variability from the same tiny sample. The t-distribution is the
honest accounting for that extra uncertainty — and in this chapter we'll *build*
it by simulation so you can see exactly where it comes from, then use R's
`t.test()` to apply it.
```

**Draft B — objectives (insert after the open):**
```markdown
## What you'll learn

- Explain why estimating the standard deviation from a small sample forces us off
  the normal distribution and onto the t-distribution.
- Read a Z-score and a t-score and say what each measures.
- Use simulation to show that t-based p-values are calibrated (uniform under the
  null) when the normal-based ones are not.
- Run one-sample, two-sample, paired, and formula-based tests with `t.test()`.
- Estimate statistical power with `power.t.test()` and by simulation.
```

**Draft C — formula gloss (after the fixed `@eq-sample-sd`, ~line 131):**
```markdown
In words: the sample standard deviation $s$ is the typical distance of a data
point from the sample mean. The only change from the population formula is the
divisor — $N-1$ instead of $N$ — which corrects for the fact that we used the
data twice: once to estimate the mean, and again to measure spread around it.
```

**Draft D — the one idea, as a callout (after ~line 136):**
```markdown
::: {.callout-important}
## The whole point of the t-distribution

Because we *estimate* the standard deviation from a small sample, our test
statistic is more uncertain than a Z-score. The t-distribution captures that by
having **fatter tails** for small samples — so it takes a more extreme result to
reach significance. As the sample grows, the estimate sharpens and the
t-distribution slides back toward the normal.
:::
```

**Draft E — break the qqplot wall (replace lines 293–304):**
```markdown
A **Q-Q plot** plots the quantiles of two distributions against each other. If
the two distributions are identical, the points fall on a straight line; where
they differ, the points peel away from it.

Here we compare our simulated z-scores against our simulated t-scores. Because
the t-scores carry the extra uncertainty from estimating the standard deviation,
we expect them to disagree with the z-scores in the tails.

```{r}
qqplot(z10k, t10k)
abline(0, 1)
```

The points bow away from the line at both ends — the t-distribution's fatter
tails, made visible. What would happen if you raised the sample size from 5 to
50? Try it: the bow should shrink as the t-distribution approaches the normal.
```

---

## Chapter 3 — `bioc-summarizedexperiment.qmd` (mean 2.0)

**One-line:** Reads like the package vignette it was lifted from (there's even a
leftover HTML comment saying so) — accurate and well-sequenced, but with no
objectives, no motivation, no analogy for a genuinely abstract object, no
summary, and no exercises. This is the clearest example of the gradient you
described, and the biggest opportunity.

### Per-section

| Section / line | Issue | Action |
|---|---|---|
| 1–3 | Title is a `#` heading, not `title:` frontmatter (other late chapters do this too); line 3 is a leftover `<!-- taken from the vignette -->` comment. | Move title to frontmatter; delete the comment; reword the open so it doesn't read as a vignette. |
| 1–30 | Jumps straight into "stores rectangular matrices of experimental results" — assumes the reader already wants a `SummarizedExperiment` and knows `GRanges`, `DataFrame`, `ExpressionSet`. No hook, no objectives, no prerequisite bridge. | Add motivating open + `## What you'll learn` + a one-line bridge to the ranges chapters (drafts A, B). |
| 31–61 "Anatomy" | The three-part structure (assays / rowData / colData kept in sync) is *the* concept and is genuinely abstract — it needs an analogy, not just a diagram. | Add a `callout-note` analogy (draft C). |
| Throughout | Mixes `=` and `<-` (`sub_se = …`, `se <- airway`). | Convert to `<-`. |
| 63–147 | Good accessor tour, but several outputs printed without interpretation ("here's `colData(se)`" → then nothing about what to notice). | Add one interpretive sentence after the key prints. |
| 215–341 "Constructing" | Strong worked example (DeRisi), but it's the end of the chapter and then it just stops after a `hist()`. | Add `## Summary` + an exercise (draft D). |
| 47, 68, 208, 262 | typos / range bug: `"child""`, "experimental and assays", `100000:1100000`, "luorescence" (see findings). | Fix inline. |
| Whole chapter | No exercises anywhere. | Add a `## Exercises` block (draft D). |

### Ready-to-paste drafts

**Draft A — frontmatter + motivating open (replace lines 1–9):**
```markdown
---
title: "Storing experiments with SummarizedExperiment"
---

An RNA-seq or microarray experiment hands you several things at once: a matrix of
numbers (counts or intensities), a description of every feature down the side
(which gene is in each row), and a description of every sample across the top
(which condition, which patient). Keep these in three separate variables and they
*will* drift out of sync — you drop a bad sample from the count matrix, forget to
drop it from the sample table, and now your "treated vs control" labels point at
the wrong columns. Mismatches exactly like this have produced retracted papers.

`SummarizedExperiment` is Bioconductor's answer: one object that holds the
assay matrix, the feature data, and the sample data together, and keeps them
aligned for you whenever you subset.
```

**Draft B — objectives + bridge (insert after the open):**
```markdown
## What you'll learn

- Describe the three coordinated parts of a `SummarizedExperiment`: `assays()`,
  `rowData()`/`rowRanges()`, and `colData()`.
- Subset an experiment by sample or feature and trust that the metadata follows.
- Reach into assays, sample data, and experiment-wide metadata with the accessor
  functions.
- Build a `SummarizedExperiment` from raw matrices and data frames.

This chapter assumes you've met `GRanges` (see @sec-ranges) and `DataFrame`; if
"genomic ranges" isn't familiar yet, skim that chapter first.
```

**Draft C — the sync analogy (insert in "Anatomy", ~after line 61):**
```markdown
::: {.callout-note}
## A mental model: three linked sheets

Picture a spreadsheet workbook with three sheets that are *locked together*:

- the **assay** sheet — the matrix of numbers, genes × samples;
- the **row** sheet (`rowData`/`rowRanges`) — one row of annotation per gene,
  lined up with the assay's rows;
- the **column** sheet (`colData`) — one row of annotation per sample, lined up
  with the assay's columns.

The magic is the lock: when you delete a sample, all three sheets update together,
so the numbers and the labels can never disagree. Every accessor below is just a
way of looking at one of these sheets.
:::
```

**Draft D — summary + exercise (append at end):**
```markdown
## Summary

A `SummarizedExperiment` bundles an assay matrix with aligned feature data
(`rowData`/`rowRanges`) and sample data (`colData`) into one object that stays
internally consistent when you subset. You've toured the accessors, subset by
phenotype, and built one from scratch out of the DeRisi microarray data.

## Exercises

1. From the `airway` object, subset to just the untreated samples
   (`se$dex == "untrt"`) and confirm the assay and `colData` both shrink.
2. How many genes (rows) and samples (columns) does `airway` have? Which accessor
   tells you the sample treatment?

::: {.callout-note collapse="true"}
## Solutions
\`\`\`{r}
untrt <- se[, se$dex == "untrt"]
dim(untrt)          # columns drop; rows unchanged
dim(se)             # genes x samples
colData(se)$dex     # treatment per sample
\`\`\`
Subsetting the object subsets `colData` in lockstep — that's the whole point.
:::
```

---

## Cross-cutting observations (from these 3 + the proven chapters)

1. **Consistency is the cheapest win.** `=`→`<-`, the callout taxonomy, the
   single-objectives-block rule, and cross-ref syntax are mechanical and apply
   book-wide. Doing this pass first makes every later prose edit land in a
   consistent frame.
2. **The scaffolding, not the content, is what decays.** The later chapters are
   technically correct and well-sequenced; what's missing is the *wrapper* —
   hook, objectives, analogy, callouts, summary, exercises-with-solutions. That's
   a templatable transformation, which is what makes a multi-pass pipeline viable.
3. **Reuse the existing datasets for continuity.** yeast (`dataframes_intro`),
   airway + DeRisi (`bioc-summarizedexperiment`), BRFSS, msleep. Calling back to
   a dataset the reader already met is itself good scaffolding.
4. **Exercises with solutions are nearly absent past the intro.** Adding them is
   the single highest-value pedagogical upgrade and the most in-scope per your
   "adding examples/materials is in scope."

---

## Proposed implementation methodology (answering "multiple passes? subagents?")

**Recommendation: passes, not parallel authors — per chapter.** Prose voice
fractures if several agents co-write one chapter, so within a chapter a *single
author* makes sequential passes. Parallelism happens *across* chapters and in the
*verification* role, where independent perspectives help rather than hurt.

### Per-chapter pipeline (one author agent, sequential)

| Pass | Name | What it does | Why this order |
|---|---|---|---|
| 0 | **Correctness** | Apply the `CHAPTER_REVIEW_FINDINGS.md` fixes for this chapter (typos, broken code/refs). | Clean base before restructuring. |
| 1 | **Structure** | Insert objectives block, motivating open, section reorder, summary/reflection skeleton, exercise stubs. | Cheapest to do on a clean file; defines the frame. |
| 2 | **Explanation** | Fill scaffolding gaps, add analogies + callouts, interpret outputs. | Needs the structure from pass 1 to know where asides go. |
| 3 | **Voice** | Harmonize cadence to the house style: second person, break walls of text, settle `=`→`<-` and naming. | Polishing layer over settled content. |
| 4 | **Examples & exercises** | Write worked examples + exercises with collapsible solutions. | Last, so they target the now-final content. |
| 5 | **Render & verify** | `quarto render <chapter>.qmd`; fix any execution error; check cross-refs; commit the updated `_freeze/<chapter>/`. | Gate; nothing ships unrendered. |

### Orchestration across chapters

- **Fan out by chapter, isolated.** Each chapter runs its 0–5 pipeline in its own
  git worktree so parallel edits don't collide; `STYLE_GUIDE.md` is the shared
  input to every author.
- **Specialize the *critics*, not the authors.** After a chapter renders, run
  independent verification subagents, each with one lens — an **andragogy critic**
  (scores against the rubric), a **render/code verifier** (did it execute; are
  outputs sane), and, once several chapters are done, a **consistency critic**
  (terminology/voice/notation drift across chapters, redundant coverage, broken
  cross-refs between chapters). Their findings feed a touch-up pass.
- **Stay in the loop between phases.** Render output and rubric re-scores come
  back to you per chapter (or per batch) so you can approve voice/scope before the
  next batch — this is the calibration intent of the pilot, continued.

### Toolchain — confirmed working

`quarto 1.9.37`, `R 4.5.1`, and the committed `_freeze/` cache are all present, so
the render-and-verify gate (pass 5) is real, not aspirational. Per the repo's
freeze convention, each implemented chapter ships with its regenerated
`_freeze/<chapter>/` files.

---

## Suggested next step

Pick **one** of the three pilots and I'll run the full 0–5 pipeline on it
end-to-end — edits applied, chapter rendered, `_freeze/` updated, on a branch — so
you can review a *finished* chapter (not just drafts) and confirm the voice and
depth before we scale to the rest. My recommendation:
**`bioc-summarizedexperiment.qmd`**, because it's the lowest-scoring and most
representative of the work the back half of the book needs.
