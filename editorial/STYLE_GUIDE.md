# The RBioc Book — House Style Guide & Pedagogical Rubric

*Derived from the class-tested chapters (`r_basics.qmd`, `r_intro_mechanics.qmd`,
`packages_and_dice.qmd`, `intro_to_rstudio.qmd`) and adult-learning (andragogy)
principles. This is the target every chapter harmonizes toward. It is the shared
input for every implementation pass and the reference for the scorecards in
`PILOT_REVIEWS.md`.*

> **Not a rewrite mandate.** The early chapters already embody most of this. The
> guide exists to (a) make the implicit voice explicit so later chapters can be
> raised to it, and (b) settle the handful of inconsistencies that even the good
> chapters have (`=` vs `<-`, `qplot`, duplicate objective headers).

---

## 1. The voice (what makes the proven chapters work)

These are observed, recurring moves in the chapters that have taught well. Treat
them as the house voice.

| Move | Example from the book | Why it works for adults |
|---|---|---|
| **Open with a concrete goal, not a definition** | *"We're going to create a pair of virtual dice that can generate random numbers."* (`r_basics`) | Adults are problem-centered, not subject-centered. A task to accomplish beats a taxonomy to memorize. |
| **Direct second-person address** | *"You type R code into the bottom line… click Enter to run it."* | Creates a coaching relationship; reader is the doer. |
| **Name and defuse anxiety** | *"If you ever see an error message, don't panic."* | Adult learners carry fear of looking incompetent; naming it lowers the affective filter. |
| **Analogy before mechanism** | the jar/urn for sampling with/without replacement; the lightbulb for install-vs-`library` | Anchors the abstract to the familiar. |
| **Code → output → explanation, in small steps** | the `die`/`sample`/`replace=TRUE` build-up | One new idea per chunk; the reader is never parsing two unknowns at once. |
| **Callouts to separate the aside from the spine** | `callout-note` "What just happened?", `callout-warning` "Capitalization matters" | Lets a confident reader skim the spine and a struggling reader get the scaffolding. |
| **Earn a small win, then name it** | *"Congratulate yourself; you've just run your first simulation in R!"* | Spaced wins sustain motivation. |
| **Metacognitive close** | the `## Reflection` self-check list in `r_intro_mechanics` | Adults self-direct; a checklist lets them gauge readiness to move on. |

**Register & cadence targets**

- Reading level ≈ grade 9–11. Most sentences under ~25 words. Vary length;
  follow a long explanatory sentence with a short one.
- Second person ("you") as default. First person singular ("I") only for a brief
  instructor anecdote (e.g., the cluster story in `intro`). First person plural
  ("we") for shared journey ("let's…").
- Break any paragraph that runs past ~6 lines or bundles more than one idea. (The
  qqplot paragraph in `t-stats-and-tests` lines 293–304 is the canonical
  anti-pattern — one breath, eight sentences, three concepts.)
- Prefer active voice and concrete verbs.

---

## 2. The chapter skeleton (target structure)

Every instructional chapter should have these elements, in this order. Items
marked *(opt)* are encouraged but not mandatory for short chapters.

1. **Title** via YAML `title:` frontmatter (not a leading `#` heading — see §4).
2. **Motivating opening** — 1–2 short paragraphs: the concrete problem this
   chapter lets you solve, ideally with a biological framing. No definitions yet.
3. **`## What you'll learn`** — 3–6 objectives written as observable actions
   ("Load a CSV with `read.csv()`", "Subset a data frame by a logical condition").
   *One* objectives block — not the duplicated "Learning goals" + "Learning
   objectives" pair seen in `dataframes_intro`.
4. **Body sections** — each new concept: motivation → minimal code → its output →
   plain-language explanation. Use callouts (§3) for asides, gotchas, and deeper
   dives. Introduce an analogy for any genuinely abstract idea.
5. **Worked example on real data** *(opt for short chapters)* — ideally
   biological; the yeast (`dataframes_intro`), airway, and DeRisi datasets are
   good house examples to reuse for continuity.
6. **`## Exercises`** — 2–5 framed tasks, each with a **collapsible solution**
   (pattern in §3). Embedded rhetorical questions ("What does `table()` do?")
   become numbered exercises with answers.
7. **`## Summary`** — short recap that ties back to the objectives in §3
   ("You can now…").
8. **`## Reflection`** *(opt)* — self-check questions, as in `r_intro_mechanics`.
9. **`## Resources`** *(opt)* — further reading / links.

---

## 3. Conventions (settle the inconsistencies)

**Code**

- **Assignment is always `<-`.** Never `=` for assignment. The book *teaches*
  `<-` (`r_intro_mechanics` callout) then several later chapters use `=`
  (`t-stats-and-tests`, parts of `dataframes_intro`, `bioc-summarizedexperiment`).
  Harmonize all to `<-`.
- **Avoid reserved/overloaded names** the book warns against: don't name objects
  `df`, `c`, `T`, `data`, `mean`. (`t-stats-and-tests` uses `df` as a data-frame
  name *and* `df=` as the degrees-of-freedom argument in the same chapter —
  confusing. Rename the variable.)
- **`library()` calls are explicit** at first use in a chapter; never rely on a
  package being attached by a previous chapter (knitr sessions are per-chapter).
- Prefer current APIs: replace `qplot()` (`r_basics`, `packages_and_dice`) with
  `ggplot()`; `tidyr::gather()` → `tidyr::pivot_longer()`.
- Keep chunks small enough to explain in one or two sentences. Interpret the
  output — don't print and move on.

**Callout taxonomy** (use consistently):

- `callout-note` — an aside or "what just happened" deeper dive.
- `callout-tip` — a practical shortcut or pro move (e.g., googling errors).
- `callout-important` — a must-grasp concept (e.g., assignment vs expression).
- `callout-warning` — a common error or gotcha (e.g., case sensitivity, `=` vs `==`).

**Collapsible solution pattern** (house standard for exercises):

```markdown
::: {.callout-note collapse="true"}
## Solution
\`\`\`{r}
# answer code
\`\`\`
A sentence on *why*, not just the code.
:::
```

**Cross-references** — use `@fig-…`, `@tbl-…`, `@eq-…`, `@sec-…` (no stray braces
like `{@fig-…}` or `\@ref{…}`). Give knitr chunk labels without `#`
(`{r fig-foo}`, not `{r #fig-foo}`).

**Markup** — `title:` in frontmatter; sentence-case headings; define jargon on
first use and link prerequisites to their home chapter with `@sec-`.

---

## 4. The evaluation rubric (how every chapter is scored)

Each chapter is rated **1–5** on seven dimensions. 5 = matches the proven
chapters; 3 = serviceable but uneven; 1 = needs substantial work.

| # | Dimension | What a 5 looks like |
|---|---|---|
| 1 | **Scaffolding & level** | Each concept introduced before use; prerequisites bridged; no unexplained leaps; no time wasted over-explaining the trivial. |
| 2 | **Learning design** | Motivating open + single objectives block + logical flow + summary tied back to objectives. |
| 3 | **Prose cadence & voice** | Second-person, varied sentence length, no walls of text, the warm coaching register of the early chapters. |
| 4 | **Adult-learner fit** | Relevance/motivation surfaced (esp. biological); respects experience; problem-centered; anxiety defused. |
| 5 | **Examples & exercises** | Incremental worked examples on real data; framed exercises **with solutions**; output interpreted. |
| 6 | **Code↔prose integration** | Code narrated; outputs explained; chunks small; current APIs. |
| 7 | **Consistency & conventions** | `<-`, callout taxonomy, naming rules, cross-ref syntax, heading/frontmatter conventions all honored. |

A chapter's **priority for revision** ≈ (5 − mean score) × (teaching centrality).
Core early chapters that score low are the highest priority; optional appendices
that score low are the lowest.
