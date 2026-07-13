# Rebuild the submission PDFs from their markdown sources.
# Both PDFs are committed, so re-run `make` in this directory after editing
# either .md or the PDFs go stale against the text.
#
#   paper.md        -> JOSE micro-paper (dormant; JOSE paused submissions)
#   paper-f1000.md  -> F1000Research Software Tool Article (live target)

PANDOC := pandoc
FLAGS  := --citeproc --bibliography=paper.bib

all: paper.pdf paper-f1000.pdf

%.pdf: %.md paper.bib
	$(PANDOC) $< -o $@ $(FLAGS)

clean:
	rm -f paper.pdf paper-f1000.pdf

.PHONY: all clean
