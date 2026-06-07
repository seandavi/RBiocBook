# JOSE submission, contributor credits, and the licensing cleanup

*Working notes captured 2026-06-07. **Revisit after the current course concludes —
target ~2026-06-28** (adjust to the real end date). Nothing here is urgent for the
course; it's the prep for a JOSE submission and a public CC0 claim.*

This file is the single place to pick this thread back up. It records (1) the plan
to submit the book to JOSE, (2) the full contributor picture, (3) the **licensing
conflict that must be resolved before submission**, and (4) title/landing-page
ideas.

---

## 1. JOSE (Journal of Open Source Education)

JOSE is the education sibling of JOSS: a short, peer-reviewed, citable paper for a
teaching resource, reviewed in the open on a GitHub thread against a checklist, with
a DOI'd paper + Zenodo archive on acceptance. Submission is a `paper.md` (YAML
front matter + ~1,000–3,000 words). Cost: zero. Turnaround: weeks to a few months.

**What the paper documents** (not a research result): what the material teaches and
to whom; learning objectives; the instructional/pedagogical design; evidence of use
in practice; how others adopt/contribute.

**Review criteria where submissions stall:**
- *Substantial scholarly effort* — a full multi-chapter book clears this easily.
- *Open license* on content and code (this is our problem — see §3).
- *Documentation for reuse* — install/setup, prerequisites (we have these).
- *Evidence of use* — even informal ("used with N trainees across M sessions").
- *Pedagogically designed vs. reference docs* — lean on the adult-learning framing
  (`index.qmd`, `editorial/STYLE_GUIDE.md`): objectives blocks, callout taxonomy,
  exercises-with-solutions, andragogy/autonomy language. This is our strong suit.

**Fit:** strong. This session already added much of the evidence — a Zenodo DOI
typed as a book, `CITATION.cff`, the biology primer (lowers the barrier), the
literate-programming chapter, and the harmonization pass (intentional design).

**TODO before submitting:**
- [ ] Pull the *current* JOSE author guide + `paper.md` template (the checklist
      drifts) — don't work from memory.
- [ ] Resolve the licensing conflict (§3). **This gates the submission.**
- [ ] Draft `paper.md`: what/for-whom, objectives, instructional design
      (adult-learning), evidence of use (CSHL / Bioconductor / AiM / med-school —
      get concrete "N trainees across M offerings"), adoption/contribution.
- [ ] Settle the author list (§2).

---

## 2. Contributors (full picture)

Now reflected in per-chapter bylines, `DESCRIPTION` `Authors@R`, and the
`contributors.qmd` page (added 2026-06-07).

**Authors (wrote/shaped material in the book):**
- Sean Davis — primary, most chapters.
- Lori Shepherd Kern — comparison_options, control_statements, factors,
  rproj_and_save_load; edits throughout (11 commits as `lshep`).
- Martin Morgan — comparison_options; much of eda_and_univariate_brfss (had been
  credited only in a hidden HTML comment — fixed).
- Stephen Turner — dplyr_intro_msleep; Data Frames data/material via his
  Bioconnector tutorials (byline added 2026-06-07).
- Garrett Grolemund — dice material + figures in r_basics / packages_and_dice /
  data_structures_overview, from *Hands-On Programming with R*.

**Adapted sources (credited, not co-authors) — see §3 for the license issue:**
- Garrett Grolemund, *Hands-On Programming with R* (CC BY-NC-ND) — dice chapters,
  `hopr_*` figures.
- Robert Kabacoff, *Modern Data Visualization with R* (CC BY-NC 4.0) — ggplot2
  chapter.
- Stephen Turner, Bioconnector tutorials — Data Frames / dplyr.
- Alboukadel Kassambara / STHDA — penalized-regression part of models.qmd.

**Other:** Janani Ravi — one `_quarto.yml` one-line commit (git contributor, not a
content author).

**Author-list decision for the paper (open):** likely Davis, Shepherd Kern, Morgan
as authors; Turner/Grolemund/Kabacoff/Kassambara as cited/adapted sources +
acknowledgments. Confirm whether "contributor" should include adapted-source
authors.

---

## 3. ⚠️ The licensing conflict (the load-bearing cleanup)

**The book is declared CC0** (license.qmd, About, CITATION.cff, the Zenodo
deposit), **but two strands of content are adapted from NonCommercial sources:**

- `ggplot2/intro.qmd` — adapted from Kabacoff, **CC BY-NC 4.0** (stated in-chapter).
- `r_basics` / `packages_and_dice` / `data_structures_overview` — adapted from
  Grolemund's *Hands-On Programming with R* and embed its figures
  (`hopr_0103`–`0107`, `hopr_0306`). HOPR's online edition is **CC BY-NC-ND**
  (confirm exact variant).

**Why it matters:**
- You cannot relicense NC/ND material as CC0 — those parts keep their terms.
- NC is **not "open"** by the Open Definition that JOSE/JOSS check → a reviewer
  will flag NC-derived content.
- The Zenodo deposit + CITATION.cff currently assert CC0 over the whole work, which
  is inaccurate for those chapters.

**Fair use is not the fix.** For a free course book, using these figures *with
attribution* is plausibly already within the CC BY-NC(-ND) licenses (NC = ✓ for a
free book; ND only forbids modifying the figures). But fair use / the NC license
does **not** make the work open or let us claim CC0 — so it doesn't satisfy JOSE.
Not legal advice; for a public CC0 claim, either get a real review or (simpler)
replace the assets.

**Cleanup options:**
- **(A) Make it genuinely open (recommended for JOSE).** Rewrite the ggplot2
  chapter and the dice sections to remove the derivation, and **replace the
  HOPR/Kabacoff figures with originals or openly-licensed equivalents**. Then the
  whole book is uniformly CC0/open. Much of this text may already be substantially
  rewritten — the figures are the concrete remaining dependency.
- (B) Carve out the adapted chapters under their NC licenses and change the claim to
  "CC0 except where noted" — honest, but those chapters then aren't open, which JOSE
  may reject.
- (C) Relicense the whole book — doesn't fix the NC restriction.

**Interim state (done 2026-06-07):** `contributors.qmd` now attributes the adapted
sources *with* their licenses and notes the book is "CC0 except for the adapted
material," so the current claim is at least honest while the cleanup waits.

**TODO:**
- [ ] Decide A/B/C (A recommended).
- [ ] If A: list every `hopr_*` and Kabacoff-derived figure; commission/redraw
      originals; confirm the ggplot2 and dice *text* is independent of the source.
- [ ] Re-credit the `r_basics` figure captions (currently uncredited inline; HOPR
      origin only in the byline + contributors page).
- [ ] Once uniformly open: update license.qmd / CITATION.cff / .zenodo.json to a
      clean single license, and re-cut a Zenodo release.

---

## 4. Title + landing page

`rbiocbook` is a handle, not a name. For visibility, give it a real
descriptive-with-subtitle title and a proper home page.

**Scope/audience (for the paper + positioning):** a path from *first line of R* →
data structures → dplyr/ggplot2 → statistics → machine learning → genomics/
Bioconductor, for **biologists and biomedical trainees who know the biology but are
new to computing** (plus a biology primer for the reverse case).

**Title candidates** (lean descriptive for discoverability; keep `RBiocBook` as the
URL handle):
1. **Data Science for Biologists: A Hands-On Introduction with R and Bioconductor**
   *(recommended)*
2. Computational Biology with R and Bioconductor: A Practical Introduction
3. R and Bioconductor for Biological Data Science
4. (named) From Bench to Bioinformatics: Learning Data Science with R and Bioconductor

**Positioning sentence:** "A hands-on, self-paced path from your first line of R to
real genomic data analysis — for biologists and biomedical trainees who know the
biology but are new to programming, statistics, and computation."

**Landing page (above the fold):** title + positioning → who-it's-for / what
you'll be able to do → the parts list → trust signals already in place (authors +
affiliations, Zenodo *How to cite* + DOI, license, "last updated", "built with
Quarto"). Gap: `index.qmd` opens with schema-theory philosophy — great as a
*preface*, weak as a *front door*. Restructure: hook + TOC + cite block above the
fold, philosophy below.

**TODO:**
- [ ] Pick a title; set it in `_quarto.yml` (and cover / `<title>` / citation).
- [ ] Restructure `index.qmd` into a landing page.
