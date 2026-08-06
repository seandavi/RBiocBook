<!--
F1000Research SOFTWARE TOOL ARTICLE draft — target: Bioconductor gateway
(the gateway explicitly accepts "teaching labs"). Redirected here 2026-07-01
because JOSE paused submissions (board deliberating eligibility changes).
The JOSE micro-paper is preserved untouched in paper/paper.md.
Resolved 2026-07-01: grant = Bioc U24 CA289073; author order senior-last
(Shepherd Kern, Davis); Martin Morgan moved to acknowledgements (author list
trimmed to two); competing interests = none disclosed. Abstract ~235 words (<300).

Revised 2026-08-06 (pre-submission review pass). All four 2026-07-01 TODOs were
already closed; these are the issues found on a reviewer-style read:
  1. Readership claim mixed scopes (site-wide views vs landing-page-only uniques
     in one clause). Restated as an explicit bound: 876 (landing page, lower
     bound) to 4,090 (sum of per-page uniques, upper bound because multi-chapter
     readers are counted once per chapter). This makes the claim rigorous AND
     removes the dependency on the blocked GA4 site-wide pull — ga4_pull.py is
     now a nice-to-have, not a submission blocker.
  2. Paper title ("Data Science for Biologists...") differs from the artifact's
     actual name everywhere else ("The RBioc Book" — site, Zenodo, CITATION.cff).
     Relationship now stated explicitly at first mention in the Introduction.
  3. "Feature-complete" (abstract) contradicted "continues to be developed in the
     open" (Use cases). Now "complete, classroom-tested".
  4. The representative-use-case paragraph was abstract ("a representative
     chapter..."). Now names the single-cell chapter concretely — which is also
     the highest-engagement page in the GA4 data, so evidence and pedagogy claim
     land on the same chapter.
  5. Added a Scope and limitations subsection — F1000 reviewers routinely ask,
     and the honest answer (no controlled learning-outcomes study) is better
     volunteered than extracted.
All figures re-verified against paper/data/ on 2026-08-06:
7,513 views / 60 pages / 876 landing-page uniques / 4,090 summed uniques /
single-cell 464.1s=7.7min, r-basics 262.7s=4.4min, simulation 361.9s=6.0min.
-->

---
title: "Data Science for Biologists: A Hands-On Introduction with R and Bioconductor"
article_type: Software Tool Article
gateway: Bioconductor
keywords:
  - R
  - Bioconductor
  - data science
  - genomics
  - bioinformatics
  - statistics
  - reproducible research
  - open educational resource
authors:
  - name: Lori Shepherd Kern
    orcid: 0000-0002-5910-4010
    affiliation: 2
  - name: Sean Davis
    orcid: 0000-0002-8991-6458
    affiliation: 1
    corresponding: true
affiliations:
  - name: University of Colorado Anschutz School of Medicine, Aurora, CO, USA
    index: 1
  - name: Roswell Park Comprehensive Cancer Center, Buffalo, NY, USA (Bioconductor Core Team)
    index: 2
date: 1 July 2026
bibliography: paper.bib
---

# Abstract

**Background.** Biomedical researchers are increasingly expected to analyze their
own data reproducibly, yet most have had no formal training in statistics or data
science, and generic "intro to R" materials teach syntax without the conceptual
grounding that biological data analysis requires.

**Methods.** *The RBioc Book* is an open-access, fully executable textbook, built as
a Quarto book on the `knitr` engine so that every chapter runs its own R code at
render time and prose, code, and output are always in sync. It is designed around
adult-learning principles, carries biological framing from the first chapters, and
builds a cumulative conceptual scaffold including R programming, data structures, and
statistical reasoning that the genomics and single-cell chapters later reuse. It is
openly licensed (CC-BY 4.0), versioned, and archived with a citable DOI.

**Results.** The book takes a reader from a first line of R to the analysis of real
genomic data with Bioconductor, teaching against public datasets.
It has been taught in four formal
course offerings (~100 trainees, 2025–2026) and, as an online resource, recorded
7,513 page views across 60 pages between May 2024 and June 2026, from between 876 and
4,090 unique readers, with per-chapter active-engagement times indicating readers
work through rather than skim the material. Instructors can adopt a chapter as a lab, fork it for a local course, or
assign it for self-paced study; every result is reproducible by re-running the code.

**Conclusions.** *The RBioc Book* is a complete, classroom-tested, openly-licensed
teaching resource that closes the concept-and-reproducibility gap for biologists endeavoring
to understand data science, statistics, and reproducible computational research. It is
also a template for large-scale executable, community-maintained teaching
materials in the Bioconductor ecosystem.

# Introduction

Biologists increasingly need to understand statistical and data science principles,
but many biologists have limited formal coursework in statistics or
data science [@biocteaching2025]. At the same time they are expected to drive
analyses of their own data — and to do so reproducibly. Many generic "intro to R"
courses show *how* to perform an analysis but do not motivate it with concept and
rationale, leaving learners able to copy a workflow but not to reason about it or
adapt it to their own problems.

*The RBioc Book* — the resource described here, published under that name at
<https://seandavi.github.io/RBiocBook> and archived under it at Zenodo — is an
open-access, executable textbook that takes a reader from
their first line of R [@rcore2024] to the analysis of real genomic data with
Bioconductor [@huber2015]. It is written for biomedical researchers who
know the biology and have been introduced to programming, statistics, and
computation but want a more formal and thorough treatment of the concepts and
skills; a biology primer makes it usable from the other direction as well, for
computational readers approaching biological data for the first time. Rather than
teaching programming in the abstract and leaving the biology for "later," the book
carries biological framing from the first chapters and builds the conceptual
scaffolding that the genomics chapters later depend on. Because the book is fully
executable and openly licensed, an instructor can adopt a chapter as a lab, fork it
for a local course, or assign it for self-paced study, and a learner can reproduce
every result by running the same code.

# Implementation

The book is built as a Quarto book [@quarto2024]
so that rendering *executes* the R code in each chapter: the prose, the code, and the
output the reader sees are produced from the same source and are always in sync. This
executable design is at the core of the book. It makes the entire book a reproducible
artifact rather than a static document, and it means that keeping the book correct is
the same act as keeping the code runnable.

The material is organized around adult-learning principles [@knowles2015adult]. Each
chapter opens with explicit learning objectives, motivates concepts before
formalizing them, and pairs worked examples with exercises and solutions. A
consistent callout taxonomy (objectives, notes, warnings, exercises) gives the reader
predictable signposts throughout. The narrative is intentionally cumulative: early
chapters earn the vocabulary and data structures that the statistics, machine-learning,
and genomics chapters reuse. The genomics material is taught against real, public
datasets and the Bioconductor classes practitioners actually use, so the skills
transfer directly to the reader's own analyses.

# Operation

**Minimal requirements.** Reading the book online requires only a web browser
(<https://seandavi.github.io/RBiocBook>). Building or extending it locally requires R
with the packages declared in the repository's `DESCRIPTION`, and
[Quarto](https://quarto.org/); the `knitr` engine executes each chapter's R at render
time.

**Workflow.** From a clone of the repository:

```sh
quarto render                       # all formats (html, pdf, epub) into _book/
quarto preview                      # live-reload server while editing
quarto render <file>.qmd --to html  # single chapter, fast iteration
```

**Adoption modes.** Because each chapter is a self-contained, executable `.qmd`
document, an instructor can (i) assign the online edition for self-paced study,
(ii) adopt an individual chapter as a lab, or (iii) fork the repository and adapt
chapters for a local course, re-rendering to regenerate every result. The book is
versioned and each release is archived with a citable DOI, so a course can pin to a
specific version.

# Use cases

*The RBioc Book* is designed to support both synchronous and asynchronous teaching:
learners ingest material at their own pace, freeing instructors to engage
participants directly rather than lecture syntax.

**Executable chapter (input → output).** The single-cell RNA-seq chapter is
representative. It loads a published mouse-brain dataset as a `SingleCellExperiment`,
log-normalizes the counts, selects highly variable genes, runs PCA and UMAP, clusters
the cells, and identifies and visualizes the marker genes that define each cluster —
producing every table and figure inline at render time, with no hidden state, so a
reader reproduces the chapter by rendering it. It also shows the cumulative design
concretely: `SingleCellExperiment` is introduced as an extension of the
`SummarizedExperiment` class taught in an earlier chapter, and the chapter tells the
reader to revisit that chapter first if the container's structure is not fresh. This
input (chapter source + declared R environment) → output (rendered prose, code, and
result) round-trip is the unit the book is built from.

**Formal course use.** The book has been taught in four formal offerings — the Cold
Spring Harbor Laboratory data-analysis course (2025 and 2026) and the BigCARE program
(2025 and 2026), each with roughly 25 participants per session (about 100 trainees
over two years).

**Open readership.** Beyond the classroom, between May 2024 and June 2026 the online
edition recorded 7,513 page views across 60 distinct pages. Unique readership over
that window is bounded rather than point-estimated: Google Analytics 4 recorded 876
unique readers on the landing page alone, which is a lower bound on the site-wide
figure, while the sum of per-page unique readers is 4,090, which is an upper bound
because a reader who works through several chapters is counted once per chapter. The
true site-wide count therefore lies between these values. Engagement is substantial
where the material is hardest: the single-cell chapter averaged 7.7 minutes of active
engagement per reader, and the R-basics and simulation chapters 4.4 and 6.0 minutes
respectively, indicating that readers work through the material rather than skim it.
The underlying per-page export is archived in the repository (`paper/data/`) so these
bounds can be recomputed. The book continues to be developed in the open and
contains social coding best practice guidance to accept contributions.

# Summary

*The RBioc Book* closes the gap between "how to run R" and "how to reason about, and
reproducibly perform, biological data analysis," in one continuous path from a first
line of R to real Bioconductor genomics. Its executable, openly-licensed,
version-archived design makes it a reusable teaching resource — adoptable whole or in
part — and a template for community-maintained, reproducible teaching materials in
the Bioconductor ecosystem.

**Scope and limitations.** The book is a teaching resource, not a methods reference:
it teaches the reasoning and the standard workflow for each topic rather than
surveying alternatives, and it deliberately stops short of specialized or
fast-moving areas — among them multi-omic integration, spatial transcriptomics,
and workflow-level pipeline orchestration — that would date quickly or require
prerequisites the book does not build. The evidence of use reported here is
adoption and engagement data; the book has not been evaluated against learning
outcomes in a controlled setting, and such a study would be the natural next step
for assessing pedagogical effect rather than reach. Being executable also carries an
ongoing maintenance obligation: chapters are re-rendered against current package
versions, and upstream Bioconductor changes require corresponding updates — a cost
the design converts from silent documentation rot into a visible build failure, but
does not eliminate.

# Data availability

The book teaches against publicly available datasets, cited and linked in place
within the relevant chapters. Website usage figures reported here were derived from
Google Analytics 4 for the book site over May 2024–June 2026; the underlying
per-page export used for the reported totals is included in the source repository
(`paper/data/`).

# Software and resource availability

- **Source code / book repository:** <https://github.com/seandavi/RBiocBook>
- **Live book:** <https://seandavi.github.io/RBiocBook>
- **Archived version (concept DOI):** <https://doi.org/10.5281/zenodo.20574829>
- **License:** Original content (prose, original figures, example code) is released
  under CC-BY 4.0; example code is additionally reusable as such. Third-party figures,
  datasets, and adapted passages retain their own open licenses, each documented on
  the book's *Contributors & acknowledgments* page.
- **Machine-readable citation:** `CITATION.cff` in the repository.

# Author contributions

Lori Shepherd Kern — authored *Base R versus the tidyverse*, *Control statements*,
*Factors*, and *Organizing, saving, and loading your work*, with edits throughout.
Sean Davis — conceived the book; primary author and editor; wrote most chapters. Both
authors reviewed and approved the manuscript.

# Competing interests

No competing interests were disclosed.

# Grant information

This work was supported in part by the US National Institutes of Health / National
Cancer Institute Informatics Technology for Cancer Research (ITCR) program through the
Bioconductor U24 award (U24 CA289073).

# Acknowledgements

We thank Martin Morgan (Bioconductor Core Team) for contributions to the *Base R
versus the tidyverse* chapter and the *BRFSS case study*; the Cold Spring Harbor
Laboratory and BigCARE course participants whose questions shaped the material; the
Bioconductor community; and the contributors credited on the book's *Contributors &
acknowledgments* page.

# References
