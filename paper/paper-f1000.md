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
     in one clause). SUPERSEDED same day — see the 2026-08-06 (later) note below.
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
Revised again 2026-08-06 (later), once GA4 credentials landed. The real
de-duplicated reader count is now in the paper, replacing the bound above.

CAUTION FOR ANYONE UPDATING THESE FIGURES — the trap that has now bitten this
paper three times is SCOPE. The GA4 property "Sean Davis — web" hosts much more
than the book (BigCARE course site, talks, campus-llm-kb, the agentic-AI
workshop site, the root). Property-wide activeUsers for this window is 3,233;
the BOOK's is 1,475. Quoting the property figure would overstate readership by
2.2x. ga4_pull.py now filters on pagePath BEGINS_WITH /RBiocBook/ and prints
both numbers under explicit labels so they cannot be swapped again.

Revised 2026-08-06 (final pre-submission pass), after Lori Shepherd Kern's
sign-off cleared the co-author gate. Checked against F1000Research's Software
Tool Article guidelines and fixed two conformance gaps:
  - The guidelines call it "good scholarly practice to mention previously
    developed tools that address similar needs, and why the current tool is
    needed." The draft had NO related-work discussion at all -- the single
    most likely reviewer objection ("how is this different from OSCA?").
    Added a "Relation to existing resources" subsection positioning the book
    against R4DS, Modern Statistics for Modern Biology, and OSCA, and stating
    the contribution as continuity across the span they leave open rather than
    depth in any one of them.
  - Guidelines specify Methods with Implementation and Operation as
    SUBsections; these were top-level. Restructured.
Also cited the three tools the book teaches but never referenced
(GenomicRanges, tidyverse, ggplot2). All 11 bib entries are now cited; no
missing keys, no orphans. The two new references were verified against
multiple independent sources before being added (Holmes & Huber 2019, CUP,
ISBN 978-1-108-70529-5; Wickham/Cetinkaya-Rundel/Grolemund 2023, O'Reilly 2nd
ed., ISBN 978-1-4920-9740-2) -- do not add a reference here without doing the
same, see the amc-ai-governance citation audit for why.
Concept DOI 10.5281/zenodo.20574829 verified resolving (-> record 20818328).
APC at time of writing: USD 1,595, payable on acceptance regardless of the
peer-review outcome.

Window extended to 2026-08-01 (was 06-22) since the pull was being redone.
All figures verified against paper/data/ga4_book_2024-05-26_to_2026-08-01.csv:
  9,098 views / 62 pages / 1,475 de-duplicated readers
  bounds check: 990 (landing page) <= 1,475 <= 4,933 (sum of per-page) OK
  single-cell 444.0s=7.4min, ggplot2 399.7s=6.7min,
  simulation 361.9s=6.0min, r-basics 262.7s=4.4min
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
9,098 page views across 62 pages from 1,475 unique readers between May 2024 and
August 2026, with per-chapter active-engagement times indicating readers
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

**Relation to existing resources.** Excellent open teaching resources already exist
at both ends of this path, and the gap *The RBioc Book* addresses is the span between
them. *R for Data Science* [@wickham2023r4ds] is the standard introduction to data
analysis in R, but it teaches against general-purpose datasets and does not carry the
reader toward biological data or the Bioconductor classes that biology requires.
*Modern Statistics for Modern Biology* [@holmes2019msmb] treats statistical reasoning
for biological data with far greater depth than we attempt, but assumes a reader
already fluent in R and comfortable with mathematical statistics. *Orchestrating
Single-Cell Analysis with Bioconductor* [@amezquita2020] is the reference work for one
domain and, by design, begins where a competent Bioconductor user begins. A learner
starting from no R at all must therefore assemble a path across three resources with
different notational conventions, different assumed backgrounds, and no shared
running example — and, in our teaching experience, most stall at the seams rather
than inside any one text.

*The RBioc Book* is deliberately not a deeper treatment of any of these. Its
contribution is continuity: one notation, one cumulative scaffold, and one executable
artifact spanning from a first line of R to real Bioconductor genomics, with the
biological motivation present from the beginning rather than deferred. It is intended
to precede rather than replace the resources above, and it points to them where a
reader should go deeper. This positioning also reflects a broader need identified
across the Bioconductor teaching community [@biocteaching2025].

# Methods

## Implementation

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
datasets and the Bioconductor classes practitioners actually use — including the
genomic-interval infrastructure of `GenomicRanges` [@lawrence2013] — so the skills
transfer directly to the reader's own analyses. Data manipulation and visualization
are taught with the `tidyverse` [@tidyverse2019] and `ggplot2` [@wickham2016ggplot2]
alongside base R, so that readers encounter the idioms actually used in practice
rather than a single house style.

## Operation

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

**Open readership.** Beyond the classroom, between May 2024 and August 2026 the online
edition recorded 9,098 page views across 62 distinct pages from 1,475 unique readers.
The reader count is Google Analytics 4's own de-duplicated figure for pages under the
book's path prefix, so a reader who works through several chapters is counted once; it
is not a sum of per-page counts, which would over-count the same reader to 4,933. For
scale, the landing page alone accounts for 990 of those readers. Engagement is
substantial where the material is hardest: the single-cell chapter averaged 7.4 minutes
of active engagement per reader, the `ggplot2` and simulation chapters 6.7 and 6.0
minutes, and the R-basics chapter 4.4 minutes, indicating that readers work through
the material rather than skim it. The underlying export is archived in the repository
(`paper/data/`), and the script that produced it is included (`paper/ga4_pull.py`), so
every figure here can be recomputed. The book continues to be developed in the open and
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
Google Analytics 4, restricted to pages under the book's path prefix, over
May 2024–August 2026. Both the export used for the reported totals and the script
that generated it are included in the source repository (`paper/data/` and
`paper/ga4_pull.py`), so the figures can be reproduced or recomputed for a later
window.

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
