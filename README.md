# The RBioc Book

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20574829.svg)](https://doi.org/10.5281/zenodo.20574829)

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

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20574829.svg)](https://doi.org/10.5281/zenodo.20574829)

The DOI above is the *concept* DOI — it always resolves to the latest version.

> Davis, S. *The RBioc Book*. <https://doi.org/10.5281/zenodo.20574829>

## License

Released under [CC0 1.0](license.qmd) — a public-domain dedication. Use it freely.
