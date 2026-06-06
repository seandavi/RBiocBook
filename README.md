# The RBioc Book

[**The RBioc Book**](https://seandavi.github.io/RBiocBook) — a [Quarto](https://quarto.org/)
book teaching introductory R, statistics, and Bioconductor for biologists and
others coming to data science for the first time.

📖 **Read it online:** <https://seandavi.github.io/RBiocBook>

## Building locally

The book is rendered with the `knitr` engine, so rendering executes the R in each
chapter — you need R with the packages declared in `DESCRIPTION`.

```sh
quarto render            # all formats (html, pdf, epub) into _book/
quarto preview           # live-reload server while editing
quarto render <file>.qmd --to html   # single chapter, fast iteration
```

## Citing

If you use this book, please cite it. A machine-readable citation is in
[`CITATION.cff`](CITATION.cff) — GitHub turns it into a **"Cite this repository"**
button (with APA and BibTeX export).

<!-- After the first Zenodo release, replace the line below with the DOI badge:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

> Davis, S. *The RBioc Book*. <https://seandavi.github.io/RBiocBook>

## License

Released under [CC0 1.0](license.qmd) — a public-domain dedication. Use it freely.
