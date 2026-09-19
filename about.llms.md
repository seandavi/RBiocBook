# About this book

Code

Author

Affiliation

[Sean Davis](https://seandavi.github.io/) [](mailto:seandavi@gmail.com)

[University of Colorado  
Anschutz School of Medicine](https://medschool.cuanschutz.edu/)

Published

June 1, 2024

Modified

September 19, 2026

This page is a short, honest account of what *The RBioc Book* is, how it was made, and — because it would be strange not to say so in a book that encourages you to use AI — how AI was used in writing it.

## What this is

The book is a collection of chapters for learning introductory R, statistics, and Bioconductor, aimed at biologists and others coming to data science for the first time. It is built with [Quarto](https://quarto.org/), executed with `knitr` (so the code you read is the code that actually ran), and published openly to GitHub Pages. The source lives at [github.com/seandavi/RBiocBook](https://github.com/seandavi/RBiocBook). The book’s original content is released under [CC-BY 4.0](license.llms.md) — reuse it freely with attribution. Third-party material retains its own license, as noted on the [Contributors](contributors.llms.md) page.

## On the use of AI

I have used AI tools, variably, throughout the production of this book, and I want to be transparent about it rather than pretend otherwise.

The use has been *variable* by design — some chapters are long-standing class material lightly touched, others were drafted, restructured, or substantially edited with AI assistance. Where I have leaned into AI, I have leaned into it for three things in particular:

- **Editing.** Tightening prose, smoothing transitions, catching inconsistencies, and bringing chapters written over many years into a single, steadier voice.
- **Producing content.** Drafting explanations, worked examples, exercises, and figures — always reviewed, corrected, and signed off by a human author.
- **Supporting the content.** Checking that what a chapter claims is *well supported*, both in its **material** (is the statistics correct? does the code run? is the biology right?) and in its **presentation** — that each chapter reflects established **adult-learning best practices** rather than just dumping information.

That last point deserves a word. This book is built on a deliberate pedagogy ([Knowles et al. 2005](#ref-knowles_adult_2005); [Center 2016](#ref-author_lifelong_2016)): adults learn best when material is problem-centered, builds on what they already know, respects their autonomy, and gives them early, concrete wins. AI has been a genuinely useful collaborator for holding every chapter to that standard — asking, chapter by chapter, whether the motivation is concrete, whether jargon is defined before it is used, and whether the exercises actually let a reader check their understanding.

> **IMPORTANT:**
>
> None of this displaces human judgment, and that is the point. Large language models produce fluent text that can *feel* authoritative while being subtly wrong, and over-reliance on them can erode the very critical thinking this book is trying to build ([Messeri and Crockett 2024](#ref-doi:10.1038/s41586-024-07146-0); [Lee et al. 2025](#ref-doi:10.1145/3706598.3713778)). So every AI contribution here was reviewed against trusted sources, every code chunk was executed, and the responsibility for what the book says rests with its human authors. AI helped us get there faster and present it better; it did not get the final word.

If you would like to use these tools in your *own* learning — which I encourage — the [AI Tools appendix](ai_tools.llms.md) is a practical guide to doing so effectively and responsibly, and the [preface](index.llms.md) says more about the learning philosophy behind the book.

## Citing the book

If you use this material, a citation is appreciated. The book is archived on Zenodo with a DOI:

> Davis, S. *The RBioc Book*. <https://doi.org/10.5281/zenodo.20574829>

![](https://zenodo.org/badge/DOI/10.5281/zenodo.20574829.svg)

DOI

That DOI ([10.5281/zenodo.20574829](https://doi.org/10.5281/zenodo.20574829)) is the *concept* DOI — it always resolves to the latest version, and each release also gets its own version-specific DOI. A machine-readable citation lives in the repository’s [`CITATION.cff`](https://github.com/seandavi/RBiocBook/blob/main/CITATION.cff) (GitHub’s “Cite this repository” button exports APA and BibTeX). The book is a living document; chapters are added and revised over time.

Center, Pew Research. 2016. “Lifelong Learning and Technology.” In *Pew Research Center: Internet, Science & Tech*. <https://www.pewresearch.org/internet/2016/03/22/lifelong-learning-and-technology/>.

Knowles, Malcolm S., Elwood F. Holton, and Richard A. Swanson. 2005. *The Adult Learner: The Definitive Classic in Adult Education and Human Resource Development*. 6th ed. Elsevier.

Lee, Hao-Ping (Hank), Advait Sarkar, Lev Tankelevitch, et al. 2025. “The Impact of Generative AI on Critical Thinking: Self-Reported Reductions in Cognitive Effort and Confidence Effects From a Survey of Knowledge Workers.” *Proceedings of the 2025 CHI Conference on Human Factors in Computing Systems*, April 25, 1–22. <https://doi.org/10.1145/3706598.3713778>.

Messeri, Lisa, and M. J. Crockett. 2024. “Artificial Intelligence and Illusions of Understanding in Scientific Research.” *Nature* 627 (8002): 49–58. <https://doi.org/10.1038/s41586-024-07146-0>.
