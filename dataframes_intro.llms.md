# 13  Data Frames

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

Almost every dataset you will meet in genomics arrives as a table: genes down the side, samples across the top, measurements in the middle, and a few columns of annotation bolted on. In R, the structure built for exactly this shape is the `data.frame`. Get comfortable with it and most of the rest of the book — loading expression data, filtering genes, summarizing by condition — becomes variations on a theme you already know.

A data.frame looks a bit like an R matrix in that it has two dimensions, rows and columns. The difference that matters is this: a matrix **must** hold a single data type, while a data.frame may hold a *different* data type in each column — numbers in one, gene names in another, a logical flag in a third. That is exactly what a real dataset looks like, which is why the data.frame, not the matrix, is the workhorse for tabular data in R.

## 13.1 What you’ll learn

By the end of this chapter you will be able to:

- Explain how a `data.frame` differs from a matrix, and why that difference matters.
- Build a `data.frame` from scratch with [`data.frame()`](https://rdrr.io/r/base/data.frame.html).
- Inspect a data frame with [`str()`](https://rdrr.io/r/utils/str.html), [`head()`](https://rdrr.io/r/utils/head.html), [`dim()`](https://rdrr.io/r/base/dim.html), [`nrow()`](https://rdrr.io/r/base/nrow.html), [`ncol()`](https://rdrr.io/r/base/nrow.html), and [`summary()`](https://rdrr.io/r/base/summary.html).
- Pull out columns with `$` and `[[ ]]`, and rows and columns with `[rows, columns]`.
- Add a new column and select rows by a logical condition.
- Save a data frame back to disk with [`write.csv()`](https://rdrr.io/r/utils/write.table.html).

We’ll do all of this on a small gene-expression table that we build ourselves, so there is nothing to download and every value is in plain sight.

## 13.2 Building a data.frame from scratch

The most direct way to understand a data.frame is to make one. Recall from the [vectors](vectors.llms.md) chapter that a column of a data frame is just a vector. So a data.frame is really a set of equal-length vectors stood up side by side, each one becoming a named column. The [`data.frame()`](https://rdrr.io/r/base/data.frame.html) function does exactly that assembly.

Imagine we have measured a handful of genes: each gene has a symbol, the chromosome it sits on, an average expression value, the log₂ fold change between a tumor and a normal sample, and a flag for whether it is a known cancer gene. Each of those is a vector, and they are all the same length — one entry per gene.

``` downlit
genes <- data.frame(
  symbol     = c("TP53", "BRCA1", "EGFR", "MYC", "PTEN", "KRAS", "ALK", "RB1"),
  chromosome = c("17", "17", "7", "8", "10", "12", "2", "13"),
  expression = c(8.2, 5.1, 11.4, 9.8, 6.7, 7.3, 4.2, 5.9),
  log2fc     = c(-1.4, -0.8, 2.1, 3.0, -2.2, 1.1, 0.3, -1.7),
  is_cancer  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
)
genes
```

      symbol chromosome expression log2fc is_cancer
    1   TP53         17        8.2   -1.4      TRUE
    2  BRCA1         17        5.1   -0.8      TRUE
    3   EGFR          7       11.4    2.1      TRUE
    4    MYC          8        9.8    3.0      TRUE
    5   PTEN         10        6.7   -2.2      TRUE
    6   KRAS         12        7.3    1.1      TRUE
    7    ALK          2        4.2    0.3     FALSE
    8    RB1         13        5.9   -1.7      TRUE

Notice that R lined up our five vectors into five columns and numbered the rows `1` through `8` down the left edge. Those row numbers are the default *row names*. Each column kept its own type: `symbol` and `chromosome` are character, `expression` and `log2fc` are numeric, and `is_cancer` is logical. That mix of types in one object is the whole point of a data.frame — and the one thing a matrix could never do.

We can confirm that R really did give us a `data.frame`:

``` downlit
class(genes)
```

    [1] "data.frame"

## 13.3 Inspecting a data.frame

Our table is small enough to print in full, but real datasets have thousands or tens of thousands of rows and you will not be able to look at them all at once. Fortunately, R gives you a toolkit for inspecting a data.frame a piece at a time.

- Structure at a glance
  - [`str()`](https://rdrr.io/r/utils/str.html) to see the columns, their types, and a preview of each
- Overviews of content
  - [`head()`](https://rdrr.io/r/utils/head.html) to show the first few rows
  - [`tail()`](https://rdrr.io/r/utils/head.html) to show the last few rows
- Size
  - [`dim()`](https://rdrr.io/r/base/dim.html) for dimensions (rows, columns)
  - [`nrow()`](https://rdrr.io/r/base/nrow.html) for the number of rows
  - [`ncol()`](https://rdrr.io/r/base/nrow.html) for the number of columns
- Names and summaries
  - [`colnames()`](https://rdrr.io/r/base/colnames.html) for the column names
  - [`rownames()`](https://rdrr.io/r/base/colnames.html) for the row “names” (often just the row numbers)
  - [`summary()`](https://rdrr.io/r/base/summary.html) for a per-column thumbnail of the data

The single most useful one is [`str()`](https://rdrr.io/r/utils/str.html), which you met in the [lists](lists.llms.md) chapter. It reports the structure: how many rows (“obs.”), how many columns (“variables”), and for each column its name, type, and first few values.

``` downlit
str(genes)
```

    'data.frame':   8 obs. of  5 variables:
     $ symbol    : chr  "TP53" "BRCA1" "EGFR" "MYC" ...
     $ chromosome: chr  "17" "17" "7" "8" ...
     $ expression: num  8.2 5.1 11.4 9.8 6.7 7.3 4.2 5.9
     $ log2fc    : num  -1.4 -0.8 2.1 3 -2.2 1.1 0.3 -1.7
     $ is_cancer : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...

[`head()`](https://rdrr.io/r/utils/head.html) and [`tail()`](https://rdrr.io/r/utils/head.html) show the top and bottom of the table — invaluable when the full thing is too long to read:

``` downlit
head(genes, 3)
```

      symbol chromosome expression log2fc is_cancer
    1   TP53         17        8.2   -1.4      TRUE
    2  BRCA1         17        5.1   -0.8      TRUE
    3   EGFR          7       11.4    2.1      TRUE

``` downlit
tail(genes, 3)
```

      symbol chromosome expression log2fc is_cancer
    6   KRAS         12        7.3    1.1      TRUE
    7    ALK          2        4.2    0.3     FALSE
    8    RB1         13        5.9   -1.7      TRUE

The size functions answer “how big is this?” at a glance. [`dim()`](https://rdrr.io/r/base/dim.html) returns the rows and columns together; [`nrow()`](https://rdrr.io/r/base/nrow.html) and [`ncol()`](https://rdrr.io/r/base/nrow.html) return them one at a time, which is handy inside other calculations:

``` downlit
dim(genes)
```

    [1] 8 5

``` downlit
nrow(genes)
```

    [1] 8

``` downlit
ncol(genes)
```

    [1] 5

The name functions tell you what the rows and columns are called. [`colnames()`](https://rdrr.io/r/base/colnames.html) is the one you’ll reach for constantly — it’s a menu of everything you can ask for:

``` downlit
colnames(genes)
```

    [1] "symbol"     "chromosome" "expression" "log2fc"     "is_cancer" 

``` downlit
rownames(genes)
```

    [1] "1" "2" "3" "4" "5" "6" "7" "8"

Finally, [`summary()`](https://rdrr.io/r/base/summary.html) walks across the columns and gives each an appropriate thumbnail: numeric ranges and quartiles for the measurement columns, simple counts for the others.

``` downlit
summary(genes)
```

           symbol      chromosome   expression         log2fc       is_cancer      
     Length   :8   Length   :8    Min.   : 4.200   Min.   :-2.200   Mode :logical  
     N.unique :8   N.unique :7    1st Qu.: 5.700   1st Qu.:-1.475   FALSE:1        
     N.blank  :0   N.blank  :0    Median : 7.000   Median :-0.250   TRUE :7        
     Min.nchar:3   Min.nchar:1    Mean   : 7.325   Mean   : 0.050                  
     Max.nchar:5   Max.nchar:2    3rd Qu.: 8.600   3rd Qu.: 1.350                  
                                  Max.   :11.400   Max.   : 3.000                  

> **NOTE:**
>
> In RStudio there is one more inspection tool, [`View()`](https://rdrr.io/r/utils/View.html) (note the capital “V”), that opens the data frame in a spreadsheet-style window you can scroll and sort. It is a great way to *get a feel* for a dataset before you start writing code against it. We don’t run it here because it only works in interactive RStudio, not when the book is rendered.
>
> ``` downlit
> View(genes)
> ```

## 13.4 Accessing columns

Because pulling out a single column is such a common operation, R gives data.frames two shortcuts for it. The first, the `$` operator, takes the data frame on the left and a column name on the right:

``` downlit
genes$symbol
```

    [1] "TP53"  "BRCA1" "EGFR"  "MYC"   "PTEN"  "KRAS"  "ALK"   "RB1"  

``` downlit
genes$expression
```

    [1]  8.2  5.1 11.4  9.8  6.7  7.3  4.2  5.9

Each of these hands you back a plain vector — exactly the vector you put in when you built the data frame. So everything you learned about vectors applies directly. You can, for example, take the mean of the expression column:

``` downlit
mean(genes$expression)
```

    [1] 7.325

The second shortcut comes from a fact worth remembering: in R, a data.frame is also a list, where each column is one element. So the list notation `[[ ]]` works too, either by position or by name:

``` downlit
genes[[3]]               # the third column, by position
```

    [1]  8.2  5.1 11.4  9.8  6.7  7.3  4.2  5.9

``` downlit
genes[["expression"]]    # the same column, by name
```

    [1]  8.2  5.1 11.4  9.8  6.7  7.3  4.2  5.9

> **NOTE:**
>
> Both pull out a single column as a plain vector, so why two notations? `$` is the quick, everyday form when you know the column name as you type it (`genes$expression`). `[[ ]]` is the form to reach for when the column is chosen by a *value* — a number (`genes[[3]]`) or a name stored in a variable (`col <- "expression"; genes[[col]]`), which `$` cannot do. Reach for `$` by default and `[[ ]]` when you’re programming with the column choice.

## 13.5 Subsetting with `[rows, columns]`

A data.frame has two dimensions, so to pull out a rectangular piece we use the square brackets with **two** values inside, separated by a comma: `[rows, columns]`. Leaving one side of the comma blank means “all of them.”

To get the first three rows and *all* columns, fill in the rows and leave the columns blank:

``` downlit
genes[1:3, ]
```

      symbol chromosome expression log2fc is_cancer
    1   TP53         17        8.2   -1.4      TRUE
    2  BRCA1         17        5.1   -0.8      TRUE
    3   EGFR          7       11.4    2.1      TRUE

To get specific rows *and* specific columns, fill in both:

``` downlit
genes[1:3, c("symbol", "expression")]
```

      symbol expression
    1   TP53        8.2
    2  BRCA1        5.1
    3   EGFR       11.4

And to get one whole column while leaving the rows blank, name the column on the right of the comma:

``` downlit
genes[, "chromosome"]
```

    [1] "17" "17" "7"  "8"  "10" "12" "2"  "13"

> **TIP:**
>
> With a one-dimensional vector you write `x[i]` — one value in the brackets. With a two-dimensional data.frame you write `df[i, j]` — two values, separated by a comma, meaning `[rows, columns]`. If you ever see an error like *“incorrect number of dimensions,”* you’ve almost always forgotten the comma.

## 13.6 Adding a column

Adding a column works just like assigning to a column that doesn’t exist yet: name it with `$` on the left of the assignment arrow, and give it a vector of the right length. Here we add an `abs_log2fc` column holding the *absolute* fold change, which is often what we sort on when we care about the size of a change regardless of its direction:

``` downlit
genes$abs_log2fc <- abs(genes$log2fc)
genes
```

      symbol chromosome expression log2fc is_cancer abs_log2fc
    1   TP53         17        8.2   -1.4      TRUE        1.4
    2  BRCA1         17        5.1   -0.8      TRUE        0.8
    3   EGFR          7       11.4    2.1      TRUE        2.1
    4    MYC          8        9.8    3.0      TRUE        3.0
    5   PTEN         10        6.7   -2.2      TRUE        2.2
    6   KRAS         12        7.3    1.1      TRUE        1.1
    7    ALK          2        4.2    0.3     FALSE        0.3
    8    RB1         13        5.9   -1.7      TRUE        1.7

The data frame now has six columns instead of five. Because [`abs()`](https://rdrr.io/r/base/MathFun.html) is vectorized (see the [vectors](vectors.llms.md) chapter), one short line computed the new value for every row at once.

## 13.7 Selecting rows by a logical condition

This is where data.frames earn their keep. Recall the logical-indexing pattern from the vectors chapter: ask a yes/no question of a column, get back a `TRUE`/`FALSE` vector, then use it to keep only the rows where the answer is `TRUE`. With a data.frame, that logical vector goes in the *row* position of `[rows, columns]`.

For example, keep only the genes with expression above 7:

``` downlit
genes[genes$expression > 7, ]
```

      symbol chromosome expression log2fc is_cancer abs_log2fc
    1   TP53         17        8.2   -1.4      TRUE        1.4
    3   EGFR          7       11.4    2.1      TRUE        2.1
    4    MYC          8        9.8    3.0      TRUE        3.0
    6   KRAS         12        7.3    1.1      TRUE        1.1

Read the inside first: `genes$expression > 7` produces a logical vector — one `TRUE`/`FALSE` per gene — and `genes[..., ]` keeps the rows where it is `TRUE`, with the blank after the comma keeping all the columns.

You can combine conditions with `&` (“and”) and `|` (“or”). Here we ask for cancer genes whose expression is also above 7:

``` downlit
genes[genes$is_cancer & genes$expression > 7, ]
```

      symbol chromosome expression log2fc is_cancer abs_log2fc
    1   TP53         17        8.2   -1.4      TRUE        1.4
    3   EGFR          7       11.4    2.1      TRUE        2.1
    4    MYC          8        9.8    3.0      TRUE        3.0
    6   KRAS         12        7.3    1.1      TRUE        1.1

And of course you can select rows and columns at the same time — keep the up-regulated genes (`log2fc > 0`) but show only their symbol and fold change:

``` downlit
genes[genes$log2fc > 0, c("symbol", "log2fc")]
```

      symbol log2fc
    3   EGFR    2.1
    4    MYC    3.0
    6   KRAS    1.1
    7    ALK    0.3

> **TIP:**
>
> Writing `genes$` in front of every column name gets repetitive. The [`subset()`](https://rdrr.io/r/base/subset.html) function lets you refer to the columns by their bare names, because it knows to look them up inside the data frame:
>
> ``` downlit
> subset(genes, is_cancer & expression > 7)
> ```
>
>       symbol chromosome expression log2fc is_cancer abs_log2fc
>     1   TP53         17        8.2   -1.4      TRUE        1.4
>     3   EGFR          7       11.4    2.1      TRUE        2.1
>     4    MYC          8        9.8    3.0      TRUE        3.0
>     6   KRAS         12        7.3    1.1      TRUE        1.1
>
> It reads almost like English: *subset the genes where `is_cancer` is true and `expression` exceeds 7.* Use whichever form is clearer to you; they do the same job.

## 13.8 Saving a data.frame

Once you have a data.frame worth keeping, you’ll want to save it. The most portable way is to write it to a `.csv` file (comma-separated values), which any other program — R, Python, even a spreadsheet — can read. The [`write.csv()`](https://rdrr.io/r/utils/write.table.html) function does this. We pass `row.names = FALSE` so the default row numbers aren’t written as an extra unnamed column:

``` downlit
write.csv(genes, "genes.csv", row.names = FALSE)
```

Reading it back later is just as direct, with [`read.csv()`](https://rdrr.io/r/utils/read.table.html):

``` downlit
genes <- read.csv("genes.csv")
```

The function names follow a pattern worth noticing: [`read.csv()`](https://rdrr.io/r/utils/read.table.html) and [`write.csv()`](https://rdrr.io/r/utils/write.table.html) are mirror images, and the `.csv` in each name tells you the file format. As always, [`help('read.csv')`](https://rdrr.io/r/utils/read.table.html) is worth a look for the many options these functions accept.

## 13.9 Exercises

These all use the `genes` data frame we built earlier in the chapter (including the `abs_log2fc` column we added).

1.  **How many genes, how many columns?** Report the number of rows and the number of columns in `genes` two ways: with a single call to [`dim()`](https://rdrr.io/r/base/dim.html), and with [`nrow()`](https://rdrr.io/r/base/nrow.html) and [`ncol()`](https://rdrr.io/r/base/nrow.html) separately.

    > **NOTE:**
    > ``` downlit
    > dim(genes)
    > ```
    >
    >     [1] 8 6
    >
    > ``` downlit
    > nrow(genes)
    > ```
    >
    >     [1] 8
    >
    > ``` downlit
    > ncol(genes)
    > ```
    >
    >     [1] 6
    >
    > [`dim()`](https://rdrr.io/r/base/dim.html) returns both numbers as a length-2 vector (rows first, then columns); [`nrow()`](https://rdrr.io/r/base/nrow.html) and [`ncol()`](https://rdrr.io/r/base/nrow.html) return each one on its own, which is handier inside other calculations.

2.  **Extract a column and summarize it.** Pull the `log2fc` column out as a plain vector and compute its mean. Confirm the extracted object is numeric.

    > **NOTE:**
    > ``` downlit
    > fc <- genes$log2fc
    > class(fc)
    > ```
    >
    >     [1] "numeric"
    >
    > ``` downlit
    > mean(fc)
    > ```
    >
    >     [1] 0.05
    >
    > `$` extracts the column as the plain numeric vector you put in, so ordinary vector functions like [`mean()`](https://rdrr.io/r/base/mean.html) work on it directly.

3.  **Subset by condition.** Keep only the rows where `is_cancer` is `TRUE` *and* `log2fc` is negative (down-regulated cancer genes). How many are there?

    > **NOTE:**
    > ``` downlit
    > down_cancer <- genes[genes$is_cancer & genes$log2fc < 0, ]
    > down_cancer
    > ```
    >
    >       symbol chromosome expression log2fc is_cancer abs_log2fc
    >     1   TP53         17        8.2   -1.4      TRUE        1.4
    >     2  BRCA1         17        5.1   -0.8      TRUE        0.8
    >     5   PTEN         10        6.7   -2.2      TRUE        2.2
    >     8    RB1         13        5.9   -1.7      TRUE        1.7
    >
    > ``` downlit
    > nrow(down_cancer)
    > ```
    >
    >     [1] 4
    >
    > The two conditions are joined with `&`; the combined logical vector goes in the *row* position of `[rows, columns]`, and the blank after the comma keeps every column. [`nrow()`](https://rdrr.io/r/base/nrow.html) then counts the surviving rows.

4.  **Add a column.** Add a logical column called `strong` that is `TRUE` when a gene’s absolute fold change (`abs_log2fc`) is at least 2. Then count how many genes are `strong`.

    > **NOTE:**
    > ``` downlit
    > genes$strong <- genes$abs_log2fc >= 2
    > sum(genes$strong)
    > ```
    >
    >     [1] 3
    >
    > Assigning to `genes$strong` creates the new column from a logical comparison. Because `TRUE` counts as 1, [`sum()`](https://rdrr.io/r/base/sum.html) of a logical column tallies the `TRUE`s — the same count-by-summing trick from the vectors chapter.

## 13.10 Summary

You can now build a `data.frame` from a set of equal-length vectors, inspect its shape and contents, pull out the rows and columns you care about, add new columns, and write the result back to disk. Specifically, you can:

- **Build** a data frame with `data.frame(name = vector, ...)`, mixing column types freely — the one thing a matrix cannot do.
- **Inspect** it with [`str()`](https://rdrr.io/r/utils/str.html) (structure), [`head()`](https://rdrr.io/r/utils/head.html)/[`tail()`](https://rdrr.io/r/utils/head.html) (content), [`dim()`](https://rdrr.io/r/base/dim.html) / [`nrow()`](https://rdrr.io/r/base/nrow.html) / [`ncol()`](https://rdrr.io/r/base/nrow.html) (size), and [`summary()`](https://rdrr.io/r/base/summary.html) (per-column thumbnail).
- **Extract a column** with `$` or `[[ ]]`, getting back a plain vector.
- **Subset** a rectangle with `[rows, columns]`, where a blank means “all.”
- **Filter rows** by a logical condition in the row position — the single most important move, and the same logical-indexing pattern you learned for vectors.
- **Add** a column by assigning to `df$newname`, and **save** with [`write.csv()`](https://rdrr.io/r/utils/write.table.html).

The logical-indexing pattern — ask a yes/no question of a column, then use the answer to keep rows — is the move you’ll repeat on every dataset for the rest of the book.

## 13.11 Reflection

- Can you state one thing a `data.frame` can do that a matrix cannot?
- Given the column `genes$expression`, can you write the subset that keeps only rows where `expression < 6`?
- What’s the difference between `genes[2]`, `genes[[2]]`, and `genes[, 2]`?
