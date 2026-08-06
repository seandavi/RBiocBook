---
title: "Data Science for Biologists: A Hands-On Introduction with R and Bioconductor"
tags:
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
  - name: Martin Morgan
    orcid: 0000-0002-5874-8148
    affiliation: 2
  - name: Sean Davis
    orcid: 0000-0002-8991-6458
    affiliation: 1
affiliations:
  - name: University of Colorado Anschutz School of Medicine, Aurora, CO, USA
    index: 1
  - name: Roswell Park Comprehensive Cancer Center, Buffalo, NY, USA (Bioconductor Core Team)
    index: 2
date: 23 June 2026
bibliography: paper.bib
---

# Summary

*The RBioc Book* is an open-access, executable textbook that takes a reader from
their first line of R to the analysis of real genomic data with Bioconductor
[@huber2015]. It is written for biologists and biomedical trainees who know the
biology and have been introduced to programming, statistics, and computation but 
want a more formal and thorough treatment of concepts and skills. The book also
includes a biology primer, for computational readers approaching biological data for
the first time. The book is built as a Quarto book [@quarto2024] in which every
chapter executes its own R [@rcore2024] code at render time, so the prose, the code,
and the output the reader sees are always in sync.

Unlike many other online materials for teaching R and statistics, the contents of the
book are meant to explore contents, provide the "why" for material, and to shore up
statistical and data science reasoning through clear descriptions, callouts, and 
explanatory figures. The goal is to allow learners to ingest material at their own 
paces and to allow instructors time to engage with course participants in the context
of a synchronous or asynchronous teaching environment. 
The book is published online, versioned,
and archived with a citable DOI.

# Statement of need

Biologists increasingly need to understand statistical and data science principles, 
but many generic "intro to R" courses show HOW to 
perform analysis but don't motivate with concept and rationale. Yet, most biomedical
researchers have never taken a formal course in statistics or data science. 
Yet, biomedical researchers are expected to drive analyses of their own data, and to 
do so in a reproducible manner. This book is meant for adult learners who want to
address well documented gaps [@biocteaching2025] in their training and to gain a more 
formal understanding of the principles of data science and statistics and to equip 
them with the skills to perform reproducible analyses of their own data.

*The RBioc Book* is designed to close that gap with one continuous path. Rather than
teaching programming in the abstract and leaving the biology for "later," it carries
biological framing from the first chapters and builds the conceptual scaffolding —
data structures, the tidyverse, statistical reasoning — that the genomics chapters
later depend on. Because the book is fully executable and openly licensed, an
instructor can adopt a chapter as a lab, fork it for a local course, or assign it for
self-paced study, and a learner can reproduce every result by running the same code.
It is feature-complete and has been used, in whole and in part, in formal courses and
by a broad online readership (see *Evidence of use*).

# Instructional design and content

The book is organized around adult-learning principles [@knowles2015adult]: each
chapter opens with explicit learning objectives, motivates concepts before formalizing
them, and pairs worked examples with exercises and solutions. A consistent callout
taxonomy (objectives, notes, warnings, exercises) gives the reader predictable signposts
throughout. The narrative is intentionally cumulative — early chapters earn the
vocabulary and data structures that the statistics, machine-learning, and genomics
chapters reuse — so that by the time a reader reaches `SummarizedExperiment` or
single-cell analysis, the container abstractions feel like an extension of skills they
already have rather than new machinery.

The genomics material is taught against real, public datasets and the Bioconductor
classes practitioners actually use, so the skills transfer directly to the reader's own
analyses. A biology primer lowers the barrier for readers without a biology background,
making the book usable from either direction. All materials provide substantial context
and rationale for the analyses, so that readers understand not just how to perform an
analysis but why. 

# Evidence of use

*The RBioc Book* has been taught in four formal offerings: the Cold Spring Harbor
Laboratory data-analysis course (2025 and 2026) and the BigCARE program (2025 and
2026), each with approximately 25 participants per session — roughly 100 trainees over
two years. Beyond the classroom, the online edition has reached a broad audience:
between May 2024 and August 2026, the book recorded 9,098 page views across 62
chapter pages from 1,475 unique readers (Google Analytics 4, de-duplicated and
restricted to the book's path prefix), with substantial per-chapter engagement (for
example, the single-cell chapter averaged over seven minutes of active
engagement per reader), indicating that readers work through the material rather than
skim it. Authorship and per-chapter contributions reflect both the primary author and
members of the Bioconductor Core Team, and the book continues to be developed in the
open with community contributions.

# Acknowledgements

We thank the Cold Spring Harbor Laboratory and BigCARE course participants whose
questions shaped the material, the Bioconductor community, and the contributors
credited on the book's *Contributors & acknowledgments* page. Original content is
released under CC-BY 4.0; third-party material retains its own (open) license as
documented in the book.

# References
