# 8  Reading and writing data files

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

Almost no analysis starts inside R. Your data lives somewhere else first — an expression matrix exported from an instrument, a sample sheet a collaborator typed into Excel, a public dataset posted on a government website. Before you can do anything useful, you have to get that tabular data *into* R, and once you’ve cleaned or summarized it, you’ll often need to get it back *out* to share with someone who does not use R.

This chapter is about that round trip. You’ll write a small table of genes and their expression values to disk, read it back, pull a real dataset straight from a URL, and do the same with Excel files. None of it is hard once you know which function to reach for — and knowing that is most of the battle.

## 8.1 What you’ll learn

- Write a data frame to a CSV file with [`write_csv()`](https://readr.tidyverse.org/reference/write_delim.html) and read it back with [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).
- Read a CSV directly from a URL, without downloading it by hand.
- Read and write Excel files with [`read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html) and [`write_xlsx()`](https://docs.ropensci.org/writexl//reference/write_xlsx.html).
- Recognize when other formats (Google Sheets, JSON, databases) need a different package.

This chapter assumes you’ve already organized your analysis into an RStudio Project (see [sec-rstudio-projects](#sec-rstudio-projects)), so that the relative paths below resolve from your project folder.

## 8.2 CSV files

A CSV (“comma-separated values”) file is the lingua franca of tabular data: plain text, one row per line, columns separated by commas. Nearly every tool — Excel, Python, databases, sequencing pipelines — can read and write CSV, which is exactly why it’s so useful for moving data around.

Base R reads and writes CSV with [`read.csv()`](https://rdrr.io/r/utils/read.table.html) and [`write.csv()`](https://rdrr.io/r/utils/write.table.html). We’ll use the `readr` package instead, which provides [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) and [`write_csv()`](https://readr.tidyverse.org/reference/write_delim.html). They are faster, give cleaner output, and guess column types more sensibly. Install it once:

``` downlit
install.packages("readr")
```

Then load it at the start of your session:

``` downlit
library(readr)
```

### 8.2.1 Writing a CSV file

We need something to write, so let’s build a small biological data frame: a handful of genes, each with a symbol and a measured expression value.

``` downlit
expr <- data.frame(
  gene_id    = c("YAL001C", "YAL002W", "YAL003W", "YAL005C", "YAL007C"),
  symbol     = c("TFC3", "VPS8", "EFB1", "SSA1", "ERP2"),
  expression = c(8.2, 5.1, 11.7, 9.4, 6.8)
)
expr
```

      gene_id symbol expression
    1 YAL001C   TFC3        8.2
    2 YAL002W   VPS8        5.1
    3 YAL003W   EFB1       11.7
    4 YAL005C   SSA1        9.4
    5 YAL007C   ERP2        6.8

Each row is one gene; the columns mix text (`gene_id`, `symbol`) and numbers (`expression`), which is the everyday shape of real expression data. Now write it to disk with [`write_csv()`](https://readr.tidyverse.org/reference/write_delim.html):

``` downlit
write_csv(expr, "data.csv")
```

That call created a file named `data.csv` in your working directory. You can confirm it landed where you expect:

``` downlit
# Where am I working right now?
getwd()
```

    [1] "/Users/davsean/Documents/git/RBiocBook/.claude/worktrees/biology-primer"

``` downlit
# Did the file get created?
dir(pattern = "data.csv")
```

    [1] "data.csv"

If you wanted the file somewhere else, you’d pass a fuller path (for example `"outputs/tables/data.csv"`). Inside an RStudio Project, that path is relative to the project folder.

### 8.2.2 Reading a CSV file

Reading the file back is the mirror image — [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html):

``` downlit
expr_in <- read_csv("data.csv")
```

    Rows: 5 Columns: 3
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (2): gene_id, symbol
    dbl (1): expression

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Notice the message `readr` prints: it tells you how it guessed each column’s type (text vs. number). Let’s look at what came back:

``` downlit
expr_in
```

    # A tibble: 5 × 3
      gene_id symbol expression
      <chr>   <chr>       <dbl>
    1 YAL001C TFC3          8.2
    2 YAL002W VPS8          5.1
    3 YAL003W EFB1         11.7
    4 YAL005C SSA1          9.4
    5 YAL007C ERP2          6.8

The genes and their expression values round-tripped intact. The `gene_id` and `symbol` columns came back as character, `expression` as a number (`dbl`) — `readr` inferred all of that automatically.

### 8.2.3 Reading a CSV straight from the web

[`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) is not limited to files on your own disk. Hand it a URL and it will fetch and parse the file in one step. This matters because a great deal of useful data lives on public websites, and reading it directly — rather than downloading by hand — keeps your analysis reproducible: anyone who runs your script pulls the same data from the same place.

As an example, the `palmerpenguins` project publishes its penguin body measurements as a public CSV on GitHub:

``` downlit
penguins <- read_csv("https://raw.githubusercontent.com/allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv")
```

    Rows: 344 Columns: 8
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    chr (3): species, island, sex
    dbl (5): bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g, year

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` downlit
penguins
```

    # A tibble: 344 × 8
       species island    bill_length_mm bill_depth_mm flipper_length_mm body_mass_g
       <chr>   <chr>              <dbl>         <dbl>             <dbl>       <dbl>
     1 Adelie  Torgersen           39.1          18.7               181        3750
     2 Adelie  Torgersen           39.5          17.4               186        3800
     3 Adelie  Torgersen           40.3          18                 195        3250
     4 Adelie  Torgersen           NA            NA                  NA          NA
     5 Adelie  Torgersen           36.7          19.3               193        3450
     6 Adelie  Torgersen           39.3          20.6               190        3650
     7 Adelie  Torgersen           38.9          17.8               181        3625
     8 Adelie  Torgersen           39.2          19.6               195        4675
     9 Adelie  Torgersen           34.1          18.1               193        3475
    10 Adelie  Torgersen           42            20.2               190        4250
    # ℹ 334 more rows
    # ℹ 2 more variables: sex <chr>, year <dbl>

The dataset — bill and flipper measurements for three penguin species — is documented [here](https://allisonhorst.github.io/palmerpenguins/). The point is not the penguins themselves: it’s that [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) detected the file’s structure over the network and handed you a data frame, exactly as it did for your local file. Any public CSV, biological or otherwise, comes in the same way.

## 8.3 Excel files

Excel spreadsheets (`.xlsx`) are everywhere in biology — sample sheets, plate layouts, hand-curated gene lists. Unlike CSV, an Excel file can hold multiple sheets, formulas, and formatting, so it needs a dedicated package. Reading is handled by `readxl`; writing by `writexl`.

### 8.3.1 Reading an Excel file

Install and load `readxl`:

``` downlit
install.packages("readxl")
```

``` downlit
library(readxl)
```

We don’t have an Excel file lying around, and — unlike [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) — `readxl` can’t read straight from a URL. So we first *download* the file, then read it from disk. Here’s a small public example spreadsheet:[^1]

``` downlit
download.file(
  "https://raw.githubusercontent.com/seandavi/RBiocBook/main/data/SaleData.xlsx",
  "SaleData.xlsx",
  mode = "wb"
)
```

The `mode = "wb"` (“write binary”) matters for non-text files like `.xlsx`; without it, the download can be corrupted on some platforms. Confirm it arrived:

``` downlit
dir(pattern = "SaleData.xlsx")
```

    [1] "SaleData.xlsx"

Now read it with [`read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html):

``` downlit
sales <- read_excel("SaleData.xlsx")
sales
```

    # A tibble: 45 × 8
       OrderDate           Region  Manager SalesMan  Item  Units Unit_price Sale_amt
       <dttm>              <chr>   <chr>   <chr>     <chr> <dbl>      <dbl>    <dbl>
     1 2018-01-06 00:00:00 East    Martha  Alexander Tele…    95       1198   113810
     2 2018-01-23 00:00:00 Central Hermann Shelli    Home…    50        500    25000
     3 2018-02-09 00:00:00 Central Hermann Luis      Tele…    36       1198    43128
     4 2018-02-26 00:00:00 Central Timothy David     Cell…    27        225     6075
     5 2018-03-15 00:00:00 West    Timothy Stephen   Tele…    56       1198    67088
     6 2018-04-01 00:00:00 East    Martha  Alexander Home…    60        500    30000
     7 2018-04-18 00:00:00 Central Martha  Steven    Tele…    75       1198    89850
     8 2018-05-05 00:00:00 Central Hermann Luis      Tele…    90       1198   107820
     9 2018-05-22 00:00:00 West    Douglas Michael   Tele…    32       1198    38336
    10 2018-06-08 00:00:00 East    Martha  Alexander Home…    60        500    30000
    # ℹ 35 more rows

[`read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html) has options for choosing a specific sheet (`sheet =`), a cell range (`range =`), or how column types are guessed — handy when a spreadsheet has more than one tab or some messy header rows.

### 8.3.2 Writing an Excel file

To go the other direction, use [`write_xlsx()`](https://docs.ropensci.org/writexl//reference/write_xlsx.html) from `writexl`. Install and load it:

``` downlit
install.packages("writexl")
```

``` downlit
library(writexl)
```

Then write our gene expression data frame to an `.xlsx` file:

``` downlit
write_xlsx(expr, "data.xlsx")
```

And confirm it was created:

``` downlit
getwd()
```

    [1] "/Users/davsean/Documents/git/RBiocBook/.claude/worktrees/biology-primer"

``` downlit
dir(pattern = "data.xlsx")
```

    [1] "data.xlsx"

If you open `data.xlsx` in Excel, you’ll see the same five genes and their expression values you started with.

## 8.4 Additional formats

CSV and Excel cover most of what you’ll meet, but R reads and writes many other formats too. When you hit one, reach for the matching package:

- **Google Sheets** — `googlesheets4` reads and writes sheets directly from your Google account.
- **JSON** — `jsonlite` converts between R objects and JSON, common for data pulled from web APIs.
- **Databases** — `DBI` together with a backend like `RSQLite` lets you query a database and pull results straight into a data frame.

The pattern is always the same: find the package for the format, then look for its `read_*` and `write_*` functions.

## 8.5 Exercises

Use the `expr` data frame from this chapter (five yeast genes with a `symbol` and an `expression` value) unless a problem says otherwise.

1.  **Round-trip a CSV.** Write `expr` to a file called `genes.csv`, then read it into a new object called `genes_in`. Confirm the file exists with [`dir()`](https://rdrr.io/r/base/list.files.html).

    > **NOTE:**
    > ``` downlit
    > write_csv(expr, "genes.csv")
    > genes_in <- read_csv("genes.csv")
    > ```
    >
    >     Rows: 5 Columns: 3
    >     ── Column specification ────────────────────────────────────────────────────────
    >     Delimiter: ","
    >     chr (2): gene_id, symbol
    >     dbl (1): expression
    >
    >     ℹ Use `spec()` to retrieve the full column specification for this data.
    >     ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    >
    > ``` downlit
    > dir(pattern = "genes.csv")
    > ```
    >
    >     [1] "genes.csv"
    >
    > [`write_csv()`](https://readr.tidyverse.org/reference/write_delim.html) saves it; [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) brings it back. The file name is just a string, so you can call it anything you like.

2.  **Read a column back.** From `genes_in`, print just the `symbol` column.

    > **NOTE:**
    > ``` downlit
    > genes_in$symbol
    > ```
    >
    >     [1] "TFC3" "VPS8" "EFB1" "SSA1" "ERP2"
    >
    > The `$` operator pulls one column out of a data frame as a vector — the same trick you used in the data frames chapter.

3.  **Save to Excel.** Write `expr` to `genes.xlsx` and confirm it was created.

    > **NOTE:**
    > ``` downlit
    > write_xlsx(expr, "genes.xlsx")
    > dir(pattern = "genes.xlsx")
    > ```
    >
    >     [1] "genes.xlsx"
    >
    > [`write_xlsx()`](https://docs.ropensci.org/writexl//reference/write_xlsx.html) from `writexl` writes a `.xlsx` file with no need for Excel to be installed.

4.  **Why a relative path?** A colleague sends you a script whose first line is `read.csv("/Users/dana/project/data/counts.csv")`. It fails on your machine. Why, and what would you change it to?

    > **NOTE:**
    > That absolute path points to a folder on Dana’s computer that doesn’t exist on yours. Inside an RStudio Project, a relative path like `read.csv("data/counts.csv")` works for anyone who opens the project, because the working directory is the project folder.

## 8.6 Summary

You can now move tabular data in and out of R in both directions:

- **Organize first.** As covered in [sec-rstudio-projects](#sec-rstudio-projects), an RStudio Project gives you a stable working directory so relative paths (`data/file.csv`) stay portable; absolute paths do not.
- **CSV** — [`write_csv()`](https://readr.tidyverse.org/reference/write_delim.html) saves a data frame, [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) reads it back, and [`read_csv()`](https://readr.tidyverse.org/reference/read_delim.html) will even read straight from a URL for reproducible access to public data.
- **Excel** — [`read_excel()`](https://readxl.tidyverse.org/reference/read_excel.html) reads `.xlsx` (download first; it can’t read a URL), and [`write_xlsx()`](https://docs.ropensci.org/writexl//reference/write_xlsx.html) writes one.
- **Other formats** — Google Sheets, JSON, and databases each have a package (`googlesheets4`, `jsonlite`, `DBI`/`RSQLite`) that follows the same `read_*`/`write_*` pattern.

With data flowing freely in and out, you’re ready to spend your time on the part that matters: the analysis.

[^1]: This spreadsheet originally comes from a [w3resource pandas exercise](https://www.w3resource.com/python-exercises/pandas/excel/SaleData.xlsx). We keep a copy in this book’s repository and download it from there, so the example doesn’t depend on a third-party host staying up.
