# Publication-readiness: provenance & licensing cleanup

Goal: bring **The RBioc Book** to a **clean, ready-to-publish** state where every
piece of content is either (a) original, (b) third-party material that is openly
licensed (CC-BY family / BSD / public domain) **with correct attribution**, or
(c) defensible fair-use (e.g. product screenshots) **with a source caption** — and
the book's stated license is accurate and consistent everywhere.

This extends `editorial/JOSE_AND_LICENSING.md` §3 (which framed options A/B/C). The
decision below is a hybrid of A+B, enabled by the scope call on 2026-06-23: *hotlinks
and screenshots are fine if appropriately licensed or fair use.*

## The decision rule (applies to every item)

| Source license | Mitigation |
|---|---|
| **NC** (NonCommercial), **ND** (NoDerivatives), **all-rights-reserved**, or **unknown/closed** | **Replace** with an original/openly-licensed equivalent, or **rewrite** to remove the derivation |
| **CC-BY / CC-BY-SA / BSD / Artistic-2.0 / public domain** | **Keep** + ensure correct attribution (source + license in caption/credit) |
| **Verbatim copied text**, any license | **Reword** to original voice; add a credit if the source license requires attribution |
| **Product screenshot** (IDE/tool UI) | **Keep** as fair-use/educational + add a source caption; optionally replace with our own screenshot |

**Definition of done (ready to publish):**
- [ ] No NC / ND / all-rights-reserved / unknown-license third-party content remains.
- [ ] Every retained third-party asset carries source + license attribution.
- [ ] No uncredited verbatim copying of prose.
- [ ] `license.qmd`, `about.qmd`, `CITATION.cff`, `.zenodo.json` state the same, accurate license.
- [ ] `contributors.qmd` lists only *retained* third-party sources, all open-licensed.
- [ ] Book renders clean in html + pdf + epub; CI green on `main`.

---

## Phase 0 — License framing & decisions (do first; unblocks Phases 1–3)

- [x] **0.1 Licensing model — DECIDED: CC-BY 4.0 for original content.** Recorded in
  [`docs/adr/0001-license-original-content-cc-by-4-0.md`](../docs/adr/0001-license-original-content-cc-by-4-0.md)
  (2026-06-23). Original content → **CC-BY 4.0**; third-party content is **carved out**
  and retains its own (open) license, enumerated in `contributors.qmd`. After Phase 1 the
  only carve-outs left are CC-BY-family/BSD/PD, so the book is uniformly *open*.
- [ ] **0.2 Fix `license.qmd`.** Correct the stale title ("Statistical analysis of
  functional genomics dataa" → the real book title) and replace the CC0 text with:
  "Original content is licensed **CC-BY 4.0**; third-party material is listed in
  [Contributors] and retains its own license." Link to `contributors.qmd`. *(Safe to do
  now — accurate as an interim carve-out statement regardless of Phase 1 progress.)*
- [ ] **0.3 Align metadata.** Update `about.qmd` to match `license.qmd` now. For the
  machine-readable SPDX fields `CITATION.cff:22` and `.zenodo.json:14`, change
  `CC0-1.0` → `CC-BY-4.0` **after Phase 1** (so the aggregate SPDX field is truthful once
  NC/ND/ARR content is gone), then cut a fresh Zenodo release (Phase 5.3).
- [ ] **0.4 Make `contributors.qmd` the canonical third-party registry.** It already is
  (lines 27–49); ensure every retained third-party asset links back to it, and prune it
  in Phase 4 as items get replaced.

---

## Phase 1 — Must replace/rewrite: NC / ND / all-rights-reserved / closed

These cannot be made open by attribution; they are the load-bearing work.

- [ ] **1.1 `data_structures_overview.qmd:26` — `images/hopr_0306.png` (HOPR, CC BY-NC-ND).**
  Replace with an original drawio schematic of R's core data structures (vector / matrix
  / array / list / data.frame). Use the existing drawio→svg/png workflow (`drawings/`,
  `images/README.md`). Then remove HOPR from `contributors.qmd`.
- [ ] **1.2 `ggplot2/intro.qmd` — adapted from Kabacoff *Modern Data Visualization* (CC BY-NC 4.0).**
  Biggest single item. (a) Confirm/rewrite the prose to be independent of the source;
  (b) remove the carried-over captions `source: http://mosaic-web.org/` (L222) and
  `source: .../MLwR` (L253); (c) replace any Kabacoff-origin figures with our own
  ggplot output; (d) verify the *insurance* dataset (from Lantz, *Machine Learning with R*)
  is freely redistributable, else swap for an open dataset (e.g. `palmerpenguins`,
  `msleep`, or a genomics table already used in the book). Remove Kabacoff from
  `contributors.qmd` when done.
- [ ] **1.3 `models.qmd:221` — penalized-regression section "adapted from STHDA."**
  STHDA is all-rights-reserved. Rewrite the ridge/LASSO/elastic-net explanation in
  original voice (the equations are standard and not copyrightable), remove the
  "adapted from STHDA" line and the sthda.com link. Remove Kassambara/STHDA from
  `contributors.qmd`.
- [ ] **1.4 `single_cell/setup.qmd:85` — Nature Methods OSCA figure (Amezquita 2019), hotlinked, ARR.**
  Replace with an original SingleCellExperiment schematic (extend the existing
  "three-linked-sheets" SE figure / `single_cell/imgs/` drawio workflow). Remove the
  `media.springernature.com` hotlink.
- [ ] **1.5 `atac-seq/atac-seq.qmd` — Buenrostro 2013 *Nature Methods* figures (ARR).**
  - L78–83 `greenleaf_chip_dnase_atac.png` (Fig 4, reproduced locally) → replace with an
    original schematic or our own data-derived figure.
  - L259–263 `fig-fraglength-nucleosome` → **remove the NCBI CDN hotlink**
    (`cdn.ncbi.nlm.nih.gov/pmc/blobs/...`); redraw the fragment-length/nucleosome
    schematic ourselves or generate it from the chapter's own BAM data.
- [ ] **1.6 Stephen Turner / Bioconnector — `dataframes_intro.qmd` + `dplyr_intro_msleep.qmd`.**
  First **confirm the original license** of the Bioconnector tutorials (Phase 0 research).
  - If **CC-BY** → move to Phase 2 (keep + proper attribution + license note; the byline
    plus a license line suffices).
  - If **NC / unstated / closed** → rewrite both chapters to be independent of Turner's
    text (keep the public Brauer/`msleep` datasets, which are independently licensed), or
    obtain explicit permission. Update the `author:` byline and `contributors.qmd`
    accordingly.

---

## Phase 2 — Must attribute or reword: open-licensed but uncredited / verbatim

- [ ] **2.1 `genomic_ranges.qmd` — verbatim GenomicRanges vignette prose (Artistic-2.0).**
  Lines 85–87, 102–107, 385–387, 429–430, 516–521, 706–729 (incl. the strand-orientation
  table). Artistic-2.0 is open but attribution-required, and verbatim copying should be
  reduced regardless. **Reword** the copied sentences into the book's voice **and** add a
  credit ("This chapter draws on the Bioconductor *GenomicRanges* vignette by Aboyoun,
  Lawrence & Morgan"). Cite the vignette, not only the Lawrence 2013 paper.
- [ ] **2.2 `r_intro_mechanics.qmd:115–161` — Carpentries `weight_kg`/2.2/"animal's weight" block (CC-BY).**
  **Rewrite** as an original, genomics-flavored example consistent with the new intro arc
  (e.g. carry the DNA/sequence or a gene-measurement example), which also removes the
  off-theme animal-weight framing. (Rewrite preferred over adding Carpentries attribution,
  since it's small and self-contained — but attribution is an acceptable fallback.)
- [ ] **2.3 `norm.qmd:104–603` — "P/D/Q/R is for…" figures (blog, license unclear).**
  The figures are R-generated, so **re-author** them: change the titles/labels/visual
  treatment to the book's own, and **fix the "D is for Distribution" mislabel** (dnorm is
  *density*) — that error is the copying tell. Once genuinely re-authored, no attribution
  is owed; if any distinctive device is kept, credit the *Stats from Stardust* post.
- [ ] **2.4 `ranges_exercises.qmd` — exercise design from Kasper Hansen's Bioconductor course.**
  Confirm the course license. If reusable → add a credit line ("exercises adapted from
  Kasper Hansen's Bioconductor course"). If not → redesign with different AnnotationHub
  datasets/IDs.
- [ ] **2.5 `310_microbiome.qmd:116–121` — TSE figure + near-verbatim slot caption (OMA/mia docs).**
  Confirm the OMA/`mia` doc license. Reword the "five additional slots / rowTree:…"
  caption into original wording; if the figure itself is from OMA/`TreeSummarizedExperiment`,
  redraw it or attribute the specific figure with its license. (Chapter already credits OMA
  generally — make the figure credit specific.)

---

## Phase 3 — Verify & annotate: keep with proper credit (per the fair-use/licensed scope)

- [ ] **3.1 Product screenshots — keep as fair-use, add source captions.**
  `intro_to_rstudio.qmd:106` (VS Code; docs are CC-BY), `:110` (Positron/Posit),
  `git_and_github.qmd:34` (RStudio terminal, Posit docs). Add a "Source: …" caption to
  each; optionally replace with our own screenshots for full control. Confirm each
  vendor permits doc screenshots (VS Code docs = CC-BY ✓).
- [ ] **3.2 `git_and_github.qmd:107` — YouTube thumbnail.** Weakest fair-use case as an
  embedded image. **Replace** with a plain hyperlink to the video (drop the thumbnail).
- [ ] **3.3 `machine_learning/intro.qmd` — scikit-learn cheat-sheet figure (BSD).**
  Confirm it's the actual sklearn asset (BSD, reusable). Add an explicit "Source:
  scikit-learn, BSD-3" caption. Keep.
- [ ] **3.4 `atac-seq/atac-seq.qmd:69–74` — Tsompana & Buck 2014 figure.**
  Confirm the journal is CC-BY (*Epigenetics & Chromatin* = BMC, CC-BY 4.0). Add a
  "© authors, CC BY 4.0" caption with the DOI. Keep.
- [ ] **3.5 `biology_primer.qmd` — Wikimedia (CC BY-SA) + NIH BioArt (PD) figures.**
  Already attributed; verify each caption has source link + license, and that CC-BY-SA is
  acknowledged as a carve-out under the Phase 0 license statement (share-alike is open and
  fine). No replacement needed.
- [ ] **3.6 `visualization_guide.qmd:56` — Tufte book-cover JPG.** Book covers are
  copyrighted/trademarked. Replace the image with a styled quote/callout (no cover needed),
  or a clearly fair-use small thumbnail with citation. Recommend: drop the image.
- [ ] **3.7 `bioc-summarizedexperiment.qmd` — `se.png` + L129–136.** Confirm `se.png` is
  original (not the SummarizedExperiment vignette diagram); if from the vignette, redraw.
  Lightly reword the L129–136 class-description echo.

---

## Phase 4 — Cleanup & consistency sweep

- [ ] **4.1 Delete orphaned `images/hopr_0103.png … hopr_0107.png`** (no longer referenced
  after the dice rework).
- [ ] **4.2 Orphan drafts** (`simulation_basics.qmd`, `protein_simulation_basics.qmd` — both
  `author: Garrett Grolemund`; not in `_quarto.yml`): either delete, or if kept for future
  use, strip the Grolemund byline and rewrite. `talks/atacseq_slides.qmd` and `ml_practical.qmd`
  (out-of-book, many publisher hotlinks) — leave, or clean if ever promoted.
- [ ] **4.3 Re-run the mechanical sweep** to confirm zero regressions: no external image
  hotlinks except those approved in Phase 3; no `author:` byline outside
  {Sean Davis, Lori Shepherd Kern, Martin Morgan} unless that source is open + credited;
  no remaining `hopr_*` references.
- [ ] **4.4 Prune `contributors.qmd`** to reflect the final state (remove HOPR, Kabacoff,
  STHDA once replaced; keep only retained open sources). Update the "note on licensing"
  callout to match the Phase 0 statement.

---

## Phase 5 — Final verification & release

- [ ] **5.1** Full `quarto render` (html + pdf + epub) clean locally.
- [ ] **5.2** CI green on `main` after merge.
- [ ] **5.3** Re-confirm `license.qmd` / `about.qmd` / `CITATION.cff` / `.zenodo.json` are
  mutually consistent; cut a fresh Zenodo release if the license statement changed.
- [ ] **5.4** (Hands-off to `JOSE_AND_LICENSING.md`) With provenance clean, the JOSE
  blocker in §3 is cleared; title/landing-page work (§4) can proceed independently.

---

## Quick-win order (suggested execution sequence)

1. Phase 0 (license framing) — small, unblocks the honest-claim state immediately.
2. Phase 4.1 + 2.2 + 3.2 + 3.6 — fastest wins (delete orphans; rewrite the small
   Carpentries block; drop the YouTube thumbnail; drop the Tufte cover).
3. Phase 2.1 + 2.3 — reword genomic_ranges prose; re-author norm.qmd figures.
4. Phase 1 figure replacements (1.1, 1.4, 1.5) via drawio — batch the schematic work.
5. Phase 1.2 (ggplot2/Kabacoff) + 1.3 (STHDA) + 1.6 (Turner) — the larger rewrites.
6. Phase 3 verify-and-annotate; then Phase 4 sweep + Phase 5 release.

## License-research verdicts (resolved 2026-06-23)

- [x] **Bioconnector / Stephen Turner** — **CC-BY-NC 4.0** (bdsr) / CC-BY-NC-SA-4.0
  (workshops). NonCommercial → **must rewrite** `dataframes_intro.qmd` + `dplyr_intro_msleep.qmd` (1.6 = rewrite, not attribute).
- [x] **Kasper Hansen course** — CC-BY-NC-ND/-SA (conflicting), NonCommercial → **must
  redesign** `ranges_exercises.qmd` (2.4 = redesign, not attribute).
- [x] **OMA book prose/figures** — **CC-BY-NC** → **must replace** the TSE figure +
  reword caption (2.5). The `mia` (Artistic-2.0) / `TreeSummarizedExperiment` (GPL)
  *code* is fine to use with attribution.
- [x] ***Stats from Stardust* blog** — no license → **all-rights-reserved** → re-author
  `norm.qmd` figures independently (2.3).
- [x] **Lantz *insurance* dataset** — reusable only via Kaggle **ODbL-1.0/DbCL-1.0**;
  cite the Kaggle mirror, not the unlicensed GitHub copies (1.2).
- [x] **VS Code docs** — **CC-BY-3.0-US** (+ MS screenshot policy: OK for education, no
  splash/beta/partial screens). **Posit/RStudio** — no blanket grant; plain UI screenshots
  are defensible nominative fair use, **but keep logos/wordmarks out** (3.1).
- [x] **scikit-learn ml_map** — **BSD-3-Clause**, reusable with the BSD notice (3.3).
- [x] **GenomicRanges vignette** — **Artistic-2.0**: reproducing/adapting vignette prose
  **with attribution** is permitted → reword + credit (2.1, no replacement needed).
