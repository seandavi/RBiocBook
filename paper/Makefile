# Rebuild the submission PDFs from their markdown sources.
# Both PDFs are committed, so re-run `make` in this directory after editing
# either .md or the PDFs go stale against the text.
#
#   paper.md        -> JOSE micro-paper (dormant; JOSE paused submissions)
#   paper-f1000.md  -> F1000Research Software Tool Article (live target)

PANDOC := pandoc
FLAGS  := --citeproc --bibliography=paper.bib

# paper-f1000.docx is the file F1000Research actually accepts: the submission
# portal takes DOC, DOCX or RTF only -- NOT PDF. The PDFs are for reading and
# for the repo record. Rebuild the docx whenever the markdown changes, or the
# submitted manuscript silently goes stale against the source.
all: paper.pdf paper-f1000.pdf paper-f1000.docx

%.pdf: %.md paper.bib
	$(PANDOC) $< -o $@ $(FLAGS)

%.docx: %.md paper.bib
	$(PANDOC) $< -o $@ $(FLAGS)

clean:
	rm -f paper.pdf paper-f1000.pdf paper-f1000.docx

.PHONY: all clean
