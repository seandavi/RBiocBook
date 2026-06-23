# ADR 0001: License the book's original content under CC-BY 4.0

- **Status:** Accepted
- **Date:** 2026-06-23
- **Deciders:** Sean Davis
- **Related:** `editorial/PUBLICATION_READINESS.md` (Phase 0), `editorial/JOSE_AND_LICENSING.md` §3

## Context

The book is currently declared **CC0-1.0** in four places: `license.qmd`,
`about.qmd`, `CITATION.cff`, and `.zenodo.json` (and the Zenodo deposit). A
provenance audit (2026-06-22/23) found that the book also contains third-party
material under more restrictive terms — CC BY-NC (Kabacoff), CC BY-NC-ND
(Grolemund/HOPR), all-rights-reserved (STHDA, publisher figures), and several
properly-attributable CC-BY / CC-BY-SA / BSD assets. A blanket CC0 claim is
therefore inaccurate: an author cannot waive rights to material they do not hold,
and CC0 is an unusual and needlessly strict choice for a textbook that quotes,
adapts, and wants to be *cited*.

We need a license that (a) is recognized as **open** by the Open Definition that
JOSE/JOSS apply, (b) is conventional and well-understood for an educational book,
(c) ensures the authors receive **attribution**, and (d) coexists cleanly with
third-party content that is openly licensed but carries its own terms.

## Decision

**Original content of the book is licensed under [Creative Commons Attribution 4.0
International (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/).**

Third-party material is **carved out**: it retains its own license, and every such
item is enumerated with its source and license in `contributors.qmd`. After the
publication-readiness cleanup (`editorial/PUBLICATION_READINESS.md`), the only
remaining carve-outs are themselves open (CC-BY family, BSD, public domain) or
defensible fair-use (e.g. product screenshots), so the work as a whole is uniformly
*open*.

This supersedes the prior CC0-1.0 declaration for original content.

## Consequences

**Positive**
- Conventional, instantly-recognized open license for a textbook; attribution to the
  authors is guaranteed.
- Satisfies the Open Definition, clearing the JOSE licensing blocker once NC/ND/ARR
  third-party content is removed (the rest of the cleanup plan).
- CC-BY is one-way compatible with reuse of CC-BY and (when share-alike is honored)
  coexists with the CC-BY-SA Wikimedia figures already used in the biology primer.

**Negative / trade-offs**
- Downstream reusers must now attribute the book (CC0 imposed no such requirement) —
  an acceptable and intended cost.
- Machine-readable metadata (`CITATION.cff`, `.zenodo.json`) currently encode
  `CC0-1.0` and must change to `CC-BY-4.0`. To keep the SPDX field truthful for the
  *aggregate*, flip those fields **only after** the NC/ND/all-rights-reserved
  third-party content is replaced (Phase 1 of the readiness plan); update the
  human-readable `license.qmd`/`about.qmd` to the accurate carve-out statement now.
- A fresh Zenodo release should be cut once the license metadata changes, producing a
  new version DOI.

## Alternatives considered

- **Keep CC0 for original content (status quo).** Rejected: needlessly strict for a
  textbook, gives authors no attribution, and still requires the same third-party
  carve-out work; offers no advantage over CC-BY here.
- **Dual/blanket relicense to a single license including the adapted parts.** Not
  possible — NC/ND terms cannot be relicensed by us; only the rights-holder can.
- **"CC0 except where noted" carve-out without removing NC/ND.** Honest, but leaves
  NC-derived chapters that JOSE may reject as not open; the cleanup removes them
  regardless, so CC-BY + open carve-outs is strictly better.

## Implementation

Tracked as Phase 0 of `editorial/PUBLICATION_READINESS.md`:
1. `license.qmd` — replace CC0 text with CC-BY 4.0 + carve-out statement; fix stale title.
2. `about.qmd` — match.
3. `CITATION.cff` / `.zenodo.json` — change `CC0-1.0` → `CC-BY-4.0` (after Phase 1).
4. `contributors.qmd` — canonical third-party registry; license note updated to match.
