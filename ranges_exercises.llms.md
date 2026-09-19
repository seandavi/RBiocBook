# 36  Ranges Exercises

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

The two ranges chapters — *Genomic ranges and features* and *Genomic ranges with real ChIP-seq peaks* — walked you through the `GenomicRanges` toolkit, first on small hand-built objects and then on a real CTCF peak set. This chapter is their hands-on companion: a set of short problems that ask you to *use* those tools rather than read about them.

The exercises come in two halves. The first half (Exercises 1–3) stays entirely on small `GRanges` you build by hand, so you can run them anywhere with nothing but `GenomicRanges` loaded — no downloads, no annotation packages. These drill the core mechanics: construction and accessors, the intra- and inter-range operations, and the `chr1`-vs-`1` naming trap that quietly breaks overlaps. The second half (Exercises 4–5) turns you loose on the **real** ENCODE CTCF peaks from the ChIP-seq chapter and a human gene model, to answer the kind of question a genomicist actually asks: which peaks sit in a promoter, and which gene is each peak closest to?

Each exercise is a short prompt. Try it yourself first; then expand the **Solution** box to reveal the code, its output, and a sentence or two on what to notice.

## 36.1 What you’ll learn

- Build a `GRanges` from scratch and pull it apart with the coordinate accessors ([`seqnames()`](https://rdrr.io/pkg/Seqinfo/man/seqinfo.html), [`start()`](https://rdrr.io/r/stats/start.html), [`width()`](https://rdrr.io/pkg/BiocGenerics/man/start.html), [`strand()`](https://rdrr.io/pkg/BiocGenerics/man/strand.html)) and [`mcols()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html).
- Reshape ranges with the intra-range operations ([`flank()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html), [`shift()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html), [`resize()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html)) and summarize a set with the inter-range operations ([`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html), [`disjoin()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html), [`gaps()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html), [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html)).
- Relate two range sets with [`findOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html), `%over%`, [`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html), [`countOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html), [`nearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html), and [`distanceToNearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html).
- Recognize the silent zero-overlap bug a `chr1`-vs-`1` `seqlevelsStyle` mismatch causes, and fix it.
- Apply the whole toolkit to a real ENCODE peak set and a `TxDb` gene model.

``` downlit
library(GenomicRanges)
```

## 36.2 Exercise 1: build a GRanges and take it apart

Before analysing anyone else’s data, it pays to be fluent in *constructing* a `GRanges` and reading every piece back out. In this exercise you’ll build a small catalogue of imaginary gene features by hand and practice the accessors.

**Build the object.** Construct a `GRanges` named `feat` with five ranges on `chr1`, with these starts and ends, strands, and a `gene` metadata column:

| name | start | end  | strand | gene |
|------|-------|------|--------|------|
| f1   | 100   | 250  | \+     | AAA  |
| f2   | 300   | 470  | \+     | BBB  |
| f3   | 600   | 760  | \-     | CCC  |
| f4   | 720   | 900  | \-     | DDD  |
| f5   | 1000  | 1180 | \+     | EEE  |

> **NOTE:**
>
> ``` downlit
> feat <- GRanges(
>   seqnames = "chr1",
>   ranges   = IRanges(start = c(100, 300, 600, 720, 1000),
>                      end   = c(250, 470, 760, 900, 1180)),
>   strand   = c("+", "+", "-", "-", "+"),
>   gene     = c("AAA", "BBB", "CCC", "DDD", "EEE"))
> names(feat) <- c("f1", "f2", "f3", "f4", "f5")
> feat
> ```
>
>     GRanges object with 5 ranges and 1 metadata column:
>          seqnames    ranges strand |        gene
>             <Rle> <IRanges>  <Rle> | <character>
>       f1     chr1   100-250      + |         AAA
>       f2     chr1   300-470      + |         BBB
>       f3     chr1   600-760      - |         CCC
>       f4     chr1   720-900      - |         DDD
>       f5     chr1 1000-1180      + |         EEE
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> The [`GRanges()`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html) constructor takes the coordinate pieces by name; any extra named argument (`gene` here) becomes a metadata column on the right of the `|`. A single `seqnames = "chr1"` is recycled to all five ranges.

**What’s on the left of the `|`, and what’s on the right?** Pull out the `seqnames`, the `strand`, the `width` of each range, and then the metadata.

> **NOTE:**
>
> ``` downlit
> seqnames(feat)
> ```
>
>     factor-Rle of length 5 with 1 run
>       Lengths:    5
>       Values : chr1
>     Levels(1): chr1
>
> ``` downlit
> strand(feat)
> ```
>
>     factor-Rle of length 5 with 3 runs
>       Lengths: 2 2 1
>       Values : + - +
>     Levels(3): + - *
>
> ``` downlit
> width(feat)        # end - start + 1, because GRanges are 1-based and inclusive
> ```
>
>     [1] 151 171 161 181 181
>
> ``` downlit
> mcols(feat)        # everything on the right of the |
> ```
>
>     DataFrame with 5 rows and 1 column
>               gene
>        <character>
>     f1         AAA
>     f2         BBB
>     f3         CCC
>     f4         DDD
>     f5         EEE
>
> ``` downlit
> mcols(feat)$gene   # one metadata column
> ```
>
>     [1] "AAA" "BBB" "CCC" "DDD" "EEE"
>
> The coordinate parts each have their own accessor; everything you attached lives in [`mcols()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html). Note `width(feat)[1]` is 151, not 150 — both endpoints count.

**Subset two ways.** Keep only the genes on the `+` strand, and — separately — keep only the ranges wider than 160 bp.

> **NOTE:**
>
> ``` downlit
> feat[strand(feat) == "+"]
> ```
>
>     GRanges object with 3 ranges and 1 metadata column:
>          seqnames    ranges strand |        gene
>             <Rle> <IRanges>  <Rle> | <character>
>       f1     chr1   100-250      + |         AAA
>       f2     chr1   300-470      + |         BBB
>       f5     chr1 1000-1180      + |         EEE
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> ``` downlit
> feat[width(feat) > 160]
> ```
>
>     GRanges object with 4 ranges and 1 metadata column:
>          seqnames    ranges strand |        gene
>             <Rle> <IRanges>  <Rle> | <character>
>       f2     chr1   300-470      + |         BBB
>       f3     chr1   600-760      - |         CCC
>       f4     chr1   720-900      - |         DDD
>       f5     chr1 1000-1180      + |         EEE
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> A `GRanges` subsets like a vector: any logical vector the same length as the object selects ranges. Here `strand(feat) == "+"` and `width(feat) > 160` are each such a vector.

## 36.3 Exercise 2: reshape and summarize one range set

This exercise stays with `feat` and works through the single-set operations: the *intra-range* ones that reshape each range on its own, and the *inter-range* ones that look at the whole set together.

**Promoters with [`flank()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html).** Grab the 100 bp immediately *upstream* of each feature. Then look carefully at `f3`/`f4` (the `-`-strand genes) and explain where their flanks landed.

> **NOTE:**
>
> ``` downlit
> proms <- flank(feat, 100)
> proms
> ```
>
>     GRanges object with 5 ranges and 1 metadata column:
>          seqnames    ranges strand |        gene
>             <Rle> <IRanges>  <Rle> | <character>
>       f1     chr1      0-99      + |         AAA
>       f2     chr1   200-299      + |         BBB
>       f3     chr1   761-860      - |         CCC
>       f4     chr1  901-1000      - |         DDD
>       f5     chr1   900-999      + |         EEE
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> [`flank()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html) is strand-aware. For the `+` genes the upstream flank sits at *lower* coordinates (just before the start); for the `-` genes `f3` and `f4`, “upstream” runs the other way, so their flanks sit at *higher* coordinates, just past their ends. Every flank is 100 bp wide — only the side changes.

**Standardize the widths with [`resize()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html).** Make every feature a 200 bp window centered on its midpoint, and confirm the widths.

> **NOTE:**
>
> ``` downlit
> win <- resize(feat, width = 200, fix = "center")
> width(win)
> ```
>
>     [1] 200 200 200 200 200
>
> ``` downlit
> all(width(win) == 200)
> ```
>
>     [1] TRUE
>
> `fix = "center"` holds each range’s midpoint still and grows or shrinks it to the requested width — the standard way to build comparable fixed-width windows.

**Collapse overlaps with [`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html).** Features `f3` (600–760) and `f4` (720–900) overlap. Merge any overlapping or touching ranges into a non-redundant set, then report how many ranges survive. (Try it once with the strands as-is, and once with `ignore.strand = TRUE`, and explain the difference.)

> **NOTE:**
>
> ``` downlit
> reduce(feat)
> ```
>
>     GRanges object with 4 ranges and 0 metadata columns:
>           seqnames    ranges strand
>              <Rle> <IRanges>  <Rle>
>       [1]     chr1   100-250      +
>       [2]     chr1   300-470      +
>       [3]     chr1 1000-1180      +
>       [4]     chr1   600-900      -
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> ``` downlit
> length(reduce(feat))
> ```
>
>     [1] 4
>
> ``` downlit
> reduce(feat, ignore.strand = TRUE)
> ```
>
>     GRanges object with 4 ranges and 0 metadata columns:
>           seqnames    ranges strand
>              <Rle> <IRanges>  <Rle>
>       [1]     chr1   100-250      *
>       [2]     chr1   300-470      *
>       [3]     chr1   600-900      *
>       [4]     chr1 1000-1180      *
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> [`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html) merges by chromosome *and strand* by default, so `f3` and `f4` — both on `-` — collapse into one range, leaving four. Because none of the ranges sit on opposite strands in the same place, `ignore.strand = TRUE` gives the same merge here, but in general it would also fuse ranges that overlap across strands.

**Per-base coverage.** Compute the [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html) of `feat` and find the single highest depth anywhere on `chr1`. Where does it occur, and why?

> **NOTE:**
>
> ``` downlit
> cov <- coverage(feat)
> cov                       # an RleList, one element per chromosome
> ```
>
>     RleList of length 1
>     $chr1
>     integer-Rle of length 1180 with 10 runs
>       Lengths:  99 151  49 171 129 120  41 140  99 181
>       Values :   0   1   0   1   0   1   2   1   0   1
>
> ``` downlit
> max(cov$chr1)             # the tallest pile-up
> ```
>
>     [1] 2
>
> [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html) counts, base by base, how many ranges cover each position. The maximum is 2, in the 720–760 stretch where `f3` and `f4` overlap; everywhere else at most one feature is present, so the depth is 1 (or 0 in the gaps).

## 36.4 Exercise 3: the `chr1`-vs-`1` overlap trap

Overlap operations match on the *exact* sequence name, so a range on `chr1` will never overlap a range on `1` — and R reports zero overlaps with no error or warning. This exercise reproduces that silent bug deliberately on two tiny hand-built objects, then fixes it, so you recognize it when it bites you for real.

**Build two overlapping sets in different naming styles.** Build `ucsc` with two ranges on `"chr1"`, and `ens` with two ranges on `"1"` that clearly overlap them by coordinate.

> **NOTE:**
>
> ``` downlit
> library(GenomeInfoDb)
> ucsc <- GRanges("chr1", IRanges(c(100, 500), c(300, 700)))
> ens  <- GRanges("1",    IRanges(c(250, 650), c(400, 800)))
> ucsc
> ```
>
>     GRanges object with 2 ranges and 0 metadata columns:
>           seqnames    ranges strand
>              <Rle> <IRanges>  <Rle>
>       [1]     chr1   100-300      *
>       [2]     chr1   500-700      *
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> ``` downlit
> ens
> ```
>
>     GRanges object with 2 ranges and 0 metadata columns:
>           seqnames    ranges strand
>              <Rle> <IRanges>  <Rle>
>       [1]        1   250-400      *
>       [2]        1   650-800      *
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> ``` downlit
> seqlevelsStyle(ucsc)
> ```
>
>     [1] "UCSC"
>
> ``` downlit
> seqlevelsStyle(ens)
> ```
>
>     [1] "NCBI"    "Ensembl" "MSU6"    "AGPvF"  
>
> By coordinate, `ucsc` range 1 (100–300) overlaps `ens` range 1 (250–400), and `ucsc` range 2 (500–700) overlaps `ens` range 2 (650–800). But the two objects use different naming conventions: `ucsc` is **UCSC** style (`chr1`) and `ens` is **Ensembl** style (`1`).

**Count the overlaps. What do you get, and is it right?**

> **NOTE:**
>
> ``` downlit
> sum(ucsc %over% ens)
> ```
>
>     Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
>       suppressWarnings() to suppress this warning.)
>
>     [1] 0
>
> ``` downlit
> countOverlaps(ucsc, ens)
> ```
>
>     Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
>       suppressWarnings() to suppress this warning.)
>
>     [1] 0 0
>
> Zero — even though we built the ranges to overlap. Nothing errors; the names `chr1` and `1` simply never match, so every overlap test fails silently. This is the trap: an overlap count of 0 that *should* be positive almost always means a `seqlevelsStyle` mismatch, not a biological result.

**Fix it and re-count.** Put both objects in the same naming style, then count again.

> **NOTE:**
>
> ``` downlit
> seqlevelsStyle(ens) <- "UCSC"   # rewrite "1" -> "chr1" in place
> ens                             # same ranges, new sequence names
> ```
>
>     GRanges object with 2 ranges and 0 metadata columns:
>           seqnames    ranges strand
>              <Rle> <IRanges>  <Rle>
>       [1]     chr1   250-400      *
>       [2]     chr1   650-800      *
>       -------
>       seqinfo: 1 sequence from an unspecified genome; no seqlengths
>
> ``` downlit
> sum(ucsc %over% ens)
> ```
>
>     [1] 2
>
> ``` downlit
> countOverlaps(ucsc, ens)
> ```
>
>     [1] 1 1
>
> Assigning a new `seqlevelsStyle` rewrites every sequence name without touching the coordinates. Once both objects are UCSC style the names match, and the overlap count is the 2 we expected all along.

> **WARNING:**
>
> Whenever two `GRanges` come from different sources — one from UCSC, one from Ensembl, one from a peak caller, one from an annotation package — confirm they share a `seqlevelsStyle` before you compare them. The mismatch produces no error, just wrong (usually zero) answers. `seqlevelsStyle(x) <- "UCSC"` (or `"Ensembl"`) harmonizes them.

## 36.5 Exercise 4: real CTCF peaks meet a gene model

Now we leave the hand-built toys behind. In the ChIP-seq chapter you imported a real set of ENCODE CTCF peaks; here you’ll load the same file and ask how those peaks relate to human genes. CTCF often binds near gene boundaries, so promoters are a natural place to look.

The peak file is the one the book keeps in its repository (ENCODE accession `ENCFF960ZGP`), read straight from its URL. The gene model is the hg38 UCSC knownGene `TxDb`, which is already declared as a dependency of the book.

``` downlit
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

bed_url <- "https://raw.githubusercontent.com/seandavi/RBiocBook/main/data/ctcf_peaks_ENCFF960ZGP.bed.gz"
ctcf <- import(bed_url, format = "narrowPeak")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
```

**Sanity-check the naming styles before you compare anything.** This is the lesson of Exercise 3, applied for real: confirm the peaks and the `TxDb` use the same `seqlevelsStyle`.

> **NOTE:**
>
> ``` downlit
> seqlevelsStyle(ctcf)
> ```
>
>     [1] "UCSC"
>
> ``` downlit
> seqlevelsStyle(genes(txdb))
> ```
>
>     [1] "UCSC"
>
> Both are **UCSC** style (`chr1`, `chr2`, …) — the peak file came from ENCODE in UCSC naming and the `TxDb` was built on the UCSC genome — so overlaps between them are safe. If either had come back `"Ensembl"`, you would harmonize first.

**Build a clean set of promoter regions.** Pull the transcripts from the `TxDb`, take a 2 kb promoter window around each transcription start site with [`promoters()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html), then collapse the overlapping windows into a non-redundant set with [`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html). How many promoter regions remain?

> **NOTE:**
>
> ``` downlit
> tx    <- transcripts(txdb)
> proms <- promoters(tx, upstream = 2000, downstream = 200)
> ```
>
>     Warning in valid.GenomicRanges.seqinfo(x, suggest.trim = TRUE): GRanges object contains 330 out-of-bound ranges located on sequences
>       chr1_GL383518v1_alt, chr1_KI270762v1_alt, chr2_GL383522v1_alt,
>       chr2_KI270774v1_alt, chr3_KI270777v1_alt, chr3_KI270781v1_alt,
>       chr4_GL000257v2_alt, chr4_KI270788v1_alt, chr5_GL339449v2_alt,
>       chr5_KI270795v1_alt, chr5_KI270898v1_alt, chr6_GL000250v2_alt,
>       chr6_GL000254v2_alt, chr6_KI270797v1_alt, chr6_KI270798v1_alt,
>       chr6_KI270801v1_alt, chr7_GL383534v2_alt, chr7_KI270803v1_alt,
>       chr7_KI270806v1_alt, chr7_KI270809v1_alt, chr8_KI270815v1_alt,
>       chr9_GL383540v1_alt, chr9_GL383541v1_alt, chr9_GL383542v1_alt,
>       chr9_KI270823v1_alt, chr10_GL383546v1_alt, chr11_KI270831v1_alt,
>       chr11_KI270902v1_alt, chr12_GL383551v1_alt, chr12_GL383553v2_alt,
>       chr12_KI270834v1_alt, chr13_KI270838v1_alt, chr14_KI270847v1_alt,
>       chr15_KI270848v1_alt, chr15_KI270850v1_alt, chr15_KI270851v1_alt,
>       chr15_KI270906v1_alt, chr16_GL383556v1_alt, chr16_GL383557v1_alt,
>       chr16_KI270854v1_alt, chr17_JH159146v1_alt, chr17_JH159147v1_alt,
>       chr17_KI270857v1_alt, chr17_KI270860v1_alt, chr17_KI270910v1_alt,
>       chr19_GL383575v2_alt, chr19_GL383576v1_alt, chr19_KI270866v1_alt,
>       chr19_KI270868v1_alt, chr19_KI270884v1_alt, chr19_KI270885v1_alt,
>       chr19_KI270889v1_alt, chr19_KI270890v1_alt, chr19_KI270891v1_alt,
>       chr19_KI270915v1_alt, chr19_KI270916v1_alt, chr19_KI270919v1_alt,
>       chr19_KI270922v1_alt, chr19_KI270923v1_alt, chr19_KI270929v1_alt,
>       chr19_KI270930v1_alt, chr19_KI270931v1_alt, chr19_KI270932v1_alt,
>       chr19_KI270933v1_alt, chr20_KI270869v1_alt, chr21_GL383581v2_alt,
>       chr21_KI270872v1_alt, chr22_KB663609v1_alt, chr22_KI270876v1_alt,
>       chr22_KI270879v1_alt, chrX_KI270880v1_alt, chr1_KI270706v1_random,
>       chr3_GL000221v1_random, chr4_GL000008v2_random, chr14_KI270726v1_random,
>       chr16_KI270728v1_random, chr22_KI270731v1_random, chrUn_KI270748v1,
>       chrUn_KI270750v1, chr4_ML143349v1_fix, chr5_KV575244v1_fix,
>       chr7_KZ208912v1_fix, chr10_MU273367v1_fix, chr11_KZ559108v1_fix,
>       chr11_MU273371v1_fix, chr15_ML143370v1_fix, chr16_ML143373v1_fix,
>       chr17_KV575245v1_fix, chr17_KV766196v1_fix, chr17_MU273381v1_fix,
>       chr17_MU273383v1_fix, chr19_MU273384v1_fix, chr20_MU273389v1_fix,
>       chrY_MU273398v1_fix, chr1_KQ458384v1_alt, chr2_KQ983256v1_alt,
>       chr4_KV766193v1_alt, chr17_KV766198v1_alt, chr19_KV575259v1_alt, and
>       chr19_MU273387v1_alt. Note that ranges located on a sequence whose length is
>       unknown (NA) or on a circular sequence are not considered out-of-bound (use
>       seqlengths() and isCircular() to get the lengths and circularity flags of the
>       underlying sequences). You can use trim() to trim these ranges. See
>       ?`trim,GenomicRanges-method` for more information.
>
> ``` downlit
> prom_regions <- reduce(proms)
> length(proms)          # one window per transcript (many overlap)
> ```
>
>     [1] 412044
>
> ``` downlit
> length(prom_regions)   # non-redundant promoter regions after reduce()
> ```
>
>     [1] 130668
>
> [`promoters()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html) builds a strand-aware window around each transcript’s start. Because many transcripts share a start neighbourhood, those windows overlap heavily; [`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html) merges each overlapping cluster into one region so nothing is double-counted downstream.

**How many CTCF peaks fall in a promoter region?** Use `%over%` for the count, and [`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html) to pull out the peaks themselves.

> **NOTE:**
>
> ``` downlit
> sum(ctcf %over% prom_regions)
> ```
>
>     [1] 12421
>
> ``` downlit
> prom_peaks <- subsetByOverlaps(ctcf, prom_regions)
> length(prom_peaks)
> ```
>
>     [1] 12421
>
> `%over%` returns a logical vector — `TRUE` for each peak that hits at least one promoter region — so [`sum()`](https://rdrr.io/r/base/sum.html) counts them. [`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html) returns those same peaks as a `GRanges`, ready for further work.

**Is that overlap more than chance?** If CTCF peaks were scattered over the genome *independently* of promoters, how many would you expect to land in a promoter by chance? Compare that expectation to what you observed.

> **NOTE:**
>
> ``` downlit
> # Fraction of the genome covered by promoter regions. Some contigs have unknown
> # (NA) lengths, so drop those when summing the genome size.
> genome_bp  <- sum(as.numeric(seqlengths(txdb)), na.rm = TRUE)
> prom_frac  <- sum(as.numeric(width(prom_regions))) / genome_bp
>
> expected <- prom_frac * length(ctcf)   # peaks expected in promoters by chance
> observed <- sum(ctcf %over% prom_regions)
> c(expected = expected, observed = observed, ratio = observed / expected)
> ```
>
>        expected    observed       ratio 
>      4424.94630 12421.00000     2.80704 
>
> Promoter regions are a small fraction of the genome, so the chance expectation is modest — yet the observed count is several-fold larger. That gap is the biology: CTCF binding is **not** randomly placed with respect to promoters. The overlap is meaningful precisely because it dwarfs the null expectation.

## 36.6 Exercise 5: which gene is each peak closest to?

Counting overlaps answers “is this peak *in* a feature?” But many peaks fall *between* genes, and there the useful question is “which gene is this peak nearest, and how far away?” That’s a job for [`nearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html) and [`distanceToNearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html).

**Annotate each peak with its nearest gene.** Get the per-gene ranges from the `TxDb`, then for every CTCF peak find the index of the nearest gene with [`nearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html). Attach the matched gene id to the peaks as a new metadata column.

> **NOTE:**
>
> ``` downlit
> gn <- genes(txdb)
> idx <- nearest(ctcf, gn)            # index into gn for each peak (NA if none)
> # Keep only peaks that found a nearest gene (e.g. drop peaks on contigs gn lacks).
> ok  <- !is.na(idx)
> ctcf_ann <- ctcf[ok]
> ctcf_ann$nearest_gene <- names(gn)[idx[ok]]
> head(ctcf_ann$nearest_gene)
> ```
>
>     [1] "8864"      "6615"      "100505534" "79648"     "105376325" "401335"   
>
> [`nearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html) returns, for each peak, the index of the single closest gene — overlapping, upstream, or downstream, whichever is nearest. Indexing `names(gn)` with that vector turns the indices into gene ids, which we store on the peaks.

**How far is each peak from its nearest gene?** Use [`distanceToNearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html) and summarize the distances. How many peaks sit *inside* a gene (distance 0)?

> **NOTE:**
>
> ``` downlit
> d2n  <- distanceToNearest(ctcf, gn)
> dist <- mcols(d2n)$distance
> summary(dist)
> ```
>
>        Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
>           0       0       0   12588    2838 2433030 
>
> ``` downlit
> sum(dist == 0)            # peaks overlapping a gene
> ```
>
>     [1] 29931
>
> ``` downlit
> mean(dist == 0)           # ...as a fraction of all matched peaks
> ```
>
>     [1] 0.6824213
>
> [`distanceToNearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html) returns a `Hits` object with a `distance` column: 0 means the peak overlaps its nearest gene, and any positive number is the base-pair gap to the closest gene edge. The summary shows the spread; the count and fraction of zeros tell you how much CTCF binding falls within genes versus in the intergenic space between them.

## 36.7 Summary

You’ve now exercised the whole range-analysis toolkit, from hand-built objects to real ENCODE data. You can:

- **Build a `GRanges`** from scratch and pull it apart with the coordinate accessors and [`mcols()`](https://rdrr.io/pkg/S4Vectors/man/Vector-class.html), and subset it like a vector.
- **Reshape** a set with [`flank()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html), [`resize()`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html), and [`reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html), and **summarize** it with [`coverage()`](https://rdrr.io/pkg/IRanges/man/coverage-methods.html) — reading strand-awareness and overlap depth correctly.
- **Recognize and fix** the silent zero-overlap bug a `chr1`-vs-`1` `seqlevelsStyle` mismatch causes.
- **Relate two real range sets** with `%over%`, [`subsetByOverlaps()`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html), [`nearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html), and [`distanceToNearest()`](https://rdrr.io/pkg/IRanges/man/nearest-methods.html), and judge an overlap against a chance expectation.

For the underlying concepts — how `GRanges` objects are built, the full catalog of set operations and overlap methods, and the coordinate-system details — return to the *Genomic ranges and features* chapter. For the real-data import workflow with `rtracklayer`, see *Genomic ranges with real ChIP-seq peaks*. This chapter was the place to put those tools to work.
