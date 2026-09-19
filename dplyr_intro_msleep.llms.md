# 16  Reshaping tables with dplyr

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

Open almost any biological analysis and a table is sitting at the center of it. A sample sheet has one row per sample. A count matrix has one row per gene. A differential-expression result has one row per transcript, with columns for the fold change, the p-value, and the adjusted p-value. The data changes from study to study, but the *moves* you make on those tables hardly change at all: keep the rows that matter, throw away the columns that don’t, sort by some number, compute a derived quantity, and boil many rows down to a handful of summary numbers.

The **dplyr** package gives each of those moves a name and a function. The functions are small, they each do one job, and they snap together. Once you can do these moves on a table of mammal sleep times, you can do them — unchanged — on a table of gene counts. That is the whole bet of this chapter: learn the moves on something easy to picture, then carry them straight over to data you actually care about.

## 16.1 What you’ll learn

- Keep and drop columns with [`select()`](https://dplyr.tidyverse.org/reference/select.html).
- Keep rows that pass a test with [`filter()`](https://dplyr.tidyverse.org/reference/filter.html).
- Sort rows with [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html), including largest-first with [`desc()`](https://dplyr.tidyverse.org/reference/desc.html).
- Build new columns from old ones with [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).
- Reduce a column to one number with [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html).
- Compute those summaries one-per-group with [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) — the *split-apply-combine* idea.
- Chain all of these together with the pipe, `|>`, so an analysis reads top to bottom.

> **TIP:**
>
> We’ll practice on `msleep`, a small table of how long different mammals sleep. It’s a stand-in, picked because it’s tiny enough to see all at once and the columns need no explaining. But hold the genomics version in your head the whole time: filtering to mammals that sleep more than 16 hours is the *same operation* as filtering to genes with an adjusted p-value below 0.05; pulling out two columns is the same as pulling `sample_id` and `condition` from a sample sheet; and averaging sleep within each taxonomic order is the same as averaging expression within each treatment group. Same verbs, bigger stakes.

## 16.2 A small table to practice on

The `msleep` dataset records sleep times, brain weights, and body weights for 83 mammals, drawn from a study by Savage and West[^1]. It ships inside **ggplot2**, so there is nothing to download — install ggplot2 once if you don’t already have it, then load it to make the data available.

``` downlit
install.packages("ggplot2")
```

``` downlit
library(ggplot2)
data(msleep)
```

A few of the columns will do most of the work for us: `name` (the common name), `order` and `vore` (taxonomic order and feeding type — carnivore, herbivore, omnivore, or insectivore), `sleep_total` and `sleep_rem` (hours of total and REM sleep), and `bodywt` and `brainwt` (weights in kilograms). The full set of 11 columns is described on the help page:

``` downlit
?msleep
```

Glancing at the first few rows is a good habit before you touch any table:

``` downlit
head(msleep)
```

    # A tibble: 6 × 11
      name    genus vore  order conservation sleep_total sleep_rem sleep_cycle awake
      <chr>   <chr> <chr> <chr> <chr>              <dbl>     <dbl>       <dbl> <dbl>
    1 Cheetah Acin… carni Carn… lc                  12.1      NA        NA      11.9
    2 Owl mo… Aotus omni  Prim… <NA>                17         1.8      NA       7  
    3 Mounta… Aplo… herbi Rode… nt                  14.4       2.4      NA       9.6
    4 Greate… Blar… omni  Sori… lc                  14.9       2.3       0.133   9.1
    5 Cow     Bos   herbi Arti… domesticated         4         0.7       0.667  20  
    6 Three-… Brad… herbi Pilo… <NA>                14.4       2.2       0.767   9.6
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

The verbs we’re about to meet live in dplyr. Install it once, then load it:

``` downlit
install.packages("dplyr")
```

``` downlit
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

## 16.3 `select()`: choose columns

[`select()`](https://dplyr.tidyverse.org/reference/select.html) returns a table with only the columns you ask for, in the order you ask for them. Notice that you write the column names bare — no quotes — as if they were already variables in scope:

``` downlit
select(msleep, name, sleep_total)
```

    # A tibble: 83 × 2
       name                       sleep_total
       <chr>                            <dbl>
     1 Cheetah                           12.1
     2 Owl monkey                        17  
     3 Mountain beaver                   14.4
     4 Greater short-tailed shrew        14.9
     5 Cow                                4  
     6 Three-toed sloth                  14.4
     7 Northern fur seal                  8.7
     8 Vesper mouse                       7  
     9 Dog                               10.1
    10 Roe deer                           3  
    # ℹ 73 more rows

A leading minus sign means “everything *except* this.” To drop the `name` column and keep the rest:

``` downlit
head(select(msleep, -name))
```

    # A tibble: 6 × 10
      genus      vore  order    conservation sleep_total sleep_rem sleep_cycle awake
      <chr>      <chr> <chr>    <chr>              <dbl>     <dbl>       <dbl> <dbl>
    1 Acinonyx   carni Carnivo… lc                  12.1      NA        NA      11.9
    2 Aotus      omni  Primates <NA>                17         1.8      NA       7  
    3 Aplodontia herbi Rodentia nt                  14.4       2.4      NA       9.6
    4 Blarina    omni  Soricom… lc                  14.9       2.3       0.133   9.1
    5 Bos        herbi Artioda… domesticated         4         0.7       0.667  20  
    6 Bradypus   herbi Pilosa   <NA>                14.4       2.2       0.767   9.6
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

A colon picks a contiguous run of columns by name, the same way `1:5` picks a run of numbers:

``` downlit
head(select(msleep, name:order))
```

    # A tibble: 6 × 4
      name                       genus      vore  order       
      <chr>                      <chr>      <chr> <chr>       
    1 Cheetah                    Acinonyx   carni Carnivora   
    2 Owl monkey                 Aotus      omni  Primates    
    3 Mountain beaver            Aplodontia herbi Rodentia    
    4 Greater short-tailed shrew Blarina    omni  Soricomorpha
    5 Cow                        Bos        herbi Artiodactyla
    6 Three-toed sloth           Bradypus   herbi Pilosa      

When column names follow a pattern, dplyr’s selection helpers save typing. To grab every column whose name begins with `"sleep"`:

``` downlit
head(select(msleep, starts_with("sleep")))
```

    # A tibble: 6 × 3
      sleep_total sleep_rem sleep_cycle
            <dbl>     <dbl>       <dbl>
    1        12.1      NA        NA    
    2        17         1.8      NA    
    3        14.4       2.4      NA    
    4        14.9       2.3       0.133
    5         4         0.7       0.667
    6        14.4       2.2       0.767

There are companions to [`starts_with()`](https://tidyselect.r-lib.org/reference/starts_with.html) — [`ends_with()`](https://tidyselect.r-lib.org/reference/starts_with.html), [`contains()`](https://tidyselect.r-lib.org/reference/starts_with.html) (a literal substring), and [`matches()`](https://tidyselect.r-lib.org/reference/starts_with.html) (a regular expression) — that all select columns by their names. They matter more than they look: a real count matrix might have hundreds of columns, and naming them one at a time is not an option.

## 16.4 `filter()`: choose rows

Where [`select()`](https://dplyr.tidyverse.org/reference/select.html) works on columns, [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) works on rows. You give it a logical test, and it keeps the rows where the test is `TRUE`. Here are the mammals that sleep at least 16 hours a day:

``` downlit
filter(msleep, sleep_total >= 16)
```

    # A tibble: 8 × 11
      name    genus vore  order conservation sleep_total sleep_rem sleep_cycle awake
      <chr>   <chr> <chr> <chr> <chr>              <dbl>     <dbl>       <dbl> <dbl>
    1 Owl mo… Aotus omni  Prim… <NA>                17         1.8      NA       7  
    2 Long-n… Dasy… carni Cing… lc                  17.4       3.1       0.383   6.6
    3 North … Dide… omni  Dide… lc                  18         4.9       0.333   6  
    4 Big br… Epte… inse… Chir… lc                  19.7       3.9       0.117   4.3
    5 Thick-… Lutr… carni Dide… lc                  19.4       6.6      NA       4.6
    6 Little… Myot… inse… Chir… <NA>                19.9       2         0.2     4.1
    7 Giant … Prio… inse… Cing… en                  18.1       6.1      NA       5.9
    8 Arctic… Sper… herbi Rode… lc                  16.6      NA        NA       7.4
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

The test is an ordinary comparison of the kind you met in the Vectors chapter — `sleep_total >= 16` is `TRUE` or `FALSE` for each row, and [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) keeps the `TRUE` ones. Any comparison works: `>`, `<`, `>=`, `<=`, `==`, `!=`.

List several conditions, separated by commas, and a row must satisfy *all* of them — the comma acts as “and.” Mammals that both sleep a lot *and* are not tiny:

``` downlit
filter(msleep, sleep_total >= 16, bodywt >= 1)
```

    # A tibble: 3 × 11
      name    genus vore  order conservation sleep_total sleep_rem sleep_cycle awake
      <chr>   <chr> <chr> <chr> <chr>              <dbl>     <dbl>       <dbl> <dbl>
    1 Long-n… Dasy… carni Cing… lc                  17.4       3.1       0.383   6.6
    2 North … Dide… omni  Dide… lc                  18         4.9       0.333   6  
    3 Giant … Prio… inse… Cing… en                  18.1       6.1      NA       5.9
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

To test membership in a set of allowed values, reach for `%in%`. It asks, for each row, “is this value one of these?” Here we keep two taxonomic orders at once:

``` downlit
filter(msleep, order %in% c("Perissodactyla", "Primates"))
```

    # A tibble: 15 × 11
       name   genus vore  order conservation sleep_total sleep_rem sleep_cycle awake
       <chr>  <chr> <chr> <chr> <chr>              <dbl>     <dbl>       <dbl> <dbl>
     1 Owl m… Aotus omni  Prim… <NA>                17         1.8      NA       7  
     2 Grivet Cerc… omni  Prim… lc                  10         0.7      NA      14  
     3 Horse  Equus herbi Peri… domesticated         2.9       0.6       1      21.1
     4 Donkey Equus herbi Peri… domesticated         3.1       0.4      NA      20.9
     5 Patas… Eryt… omni  Prim… lc                  10.9       1.1      NA      13.1
     6 Galago Gala… omni  Prim… <NA>                 9.8       1.1       0.55   14.2
     7 Human  Homo  omni  Prim… <NA>                 8         1.9       1.5    16  
     8 Mongo… Lemur herbi Prim… vu                   9.5       0.9      NA      14.5
     9 Macaq… Maca… omni  Prim… <NA>                10.1       1.2       0.75   13.9
    10 Slow … Nyct… carni Prim… <NA>                11        NA        NA      13  
    11 Chimp… Pan   omni  Prim… <NA>                 9.7       1.4       1.42   14.3
    12 Baboon Papio omni  Prim… <NA>                 9.4       1         0.667  14.6
    13 Potto  Pero… omni  Prim… lc                  11        NA        NA      13  
    14 Squir… Saim… omni  Prim… <NA>                 9.6       1.4      NA      14.4
    15 Brazi… Tapi… herbi Peri… vu                   4.4       1         0.9    19.6
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

This is exactly the operation behind “keep only the genes that passed significance” or “keep only the samples from these three batches.”

## 16.5 The pipe: reading an analysis top to bottom

So far each call has stood alone. Real work strings these verbs together, and nesting the calls — `head(select(filter(msleep, ...), ...))` — gets unreadable fast, because you have to read it inside-out. The **pipe** fixes this. Written `|>`, it takes whatever is on its left and feeds it as the *first argument* to the function on its right.

> **TIP:**
>
> `x |> f() |> g()` means: start with `x`, **and then** apply `f`, **and then** apply `g`. The steps appear in the order they actually happen, top to bottom, instead of turned inside-out by nesting. Older code uses `%>%` from the magrittr package for the same idea; `|>` is built into R itself, and the two behave the same for everything we do here.

Compare the nested form with the piped form. Nested, read inside-out:

``` downlit
head(select(msleep, name, sleep_total))
```

    # A tibble: 6 × 2
      name                       sleep_total
      <chr>                            <dbl>
    1 Cheetah                           12.1
    2 Owl monkey                        17  
    3 Mountain beaver                   14.4
    4 Greater short-tailed shrew        14.9
    5 Cow                                4  
    6 Three-toed sloth                  14.4

Piped, read top to bottom — take `msleep`, **and then** select two columns, **and then** show the head:

``` downlit
msleep |>
    select(name, sleep_total) |>
    head()
```

    # A tibble: 6 × 2
      name                       sleep_total
      <chr>                            <dbl>
    1 Cheetah                           12.1
    2 Owl monkey                        17  
    3 Mountain beaver                   14.4
    4 Greater short-tailed shrew        14.9
    5 Cow                                4  
    6 Three-toed sloth                  14.4

Both produce the same answer. The piped version wins the moment a chain grows past two steps, so we’ll use it from here on.

## 16.6 `arrange()`: sort rows

[`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html) reorders the rows by one or more columns, smallest first by default. Sort by taxonomic order:

``` downlit
msleep |> arrange(order) |> head()
```

    # A tibble: 6 × 11
      name    genus vore  order conservation sleep_total sleep_rem sleep_cycle awake
      <chr>   <chr> <chr> <chr> <chr>              <dbl>     <dbl>       <dbl> <dbl>
    1 Tenrec  Tenr… omni  Afro… <NA>                15.6       2.3      NA       8.4
    2 Cow     Bos   herbi Arti… domesticated         4         0.7       0.667  20  
    3 Roe de… Capr… herbi Arti… lc                   3        NA        NA      21  
    4 Goat    Capri herbi Arti… lc                   5.3       0.6      NA      18.7
    5 Giraffe Gira… herbi Arti… cd                   1.9       0.4      NA      22.1
    6 Sheep   Ovis  herbi Arti… domesticated         3.8       0.6      NA      20.2
    # ℹ 2 more variables: brainwt <dbl>, bodywt <dbl>

Name more than one column and [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html) sorts by the first, breaking ties with the second. Here we select three columns, then sort by order and, within each order, by total sleep:

``` downlit
msleep |>
    select(name, order, sleep_total) |>
    arrange(order, sleep_total) |>
    head()
```

    # A tibble: 6 × 3
      name     order        sleep_total
      <chr>    <chr>              <dbl>
    1 Tenrec   Afrosoricida        15.6
    2 Giraffe  Artiodactyla         1.9
    3 Roe deer Artiodactyla         3  
    4 Sheep    Artiodactyla         3.8
    5 Cow      Artiodactyla         4  
    6 Goat     Artiodactyla         5.3

Wrap a column in [`desc()`](https://dplyr.tidyverse.org/reference/desc.html) to sort it from largest to smallest. The sleepiest mammals within each order, heaviest sleepers first:

``` downlit
msleep |>
    select(name, order, sleep_total) |>
    arrange(order, desc(sleep_total)) |>
    head()
```

    # A tibble: 6 × 3
      name     order        sleep_total
      <chr>    <chr>              <dbl>
    1 Tenrec   Afrosoricida        15.6
    2 Pig      Artiodactyla         9.1
    3 Goat     Artiodactyla         5.3
    4 Cow      Artiodactyla         4  
    5 Sheep    Artiodactyla         3.8
    6 Roe deer Artiodactyla         3  

Sorting a results table by p-value, or a count table by total expression, is this same verb.

## 16.7 `mutate()`: build new columns

[`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) adds columns computed from the ones already there, leaving the original columns in place. REM sleep as a fraction of total sleep is a natural quantity to derive:

``` downlit
msleep |>
    mutate(rem_fraction = sleep_rem / sleep_total) |>
    select(name, sleep_total, sleep_rem, rem_fraction) |>
    head()
```

    # A tibble: 6 × 4
      name                       sleep_total sleep_rem rem_fraction
      <chr>                            <dbl>     <dbl>        <dbl>
    1 Cheetah                           12.1      NA         NA    
    2 Owl monkey                        17         1.8        0.106
    3 Mountain beaver                   14.4       2.4        0.167
    4 Greater short-tailed shrew        14.9       2.3        0.154
    5 Cow                                4         0.7        0.175
    6 Three-toed sloth                  14.4       2.2        0.153

You can compute several new columns in one call, separated by commas, and a later one may even use a column the same [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) just created. Here we convert body weight to grams and flag the heavyweights:

``` downlit
msleep |>
    mutate(bodywt_g = bodywt * 1000,
           is_large = bodywt_g > 10000) |>
    select(name, bodywt, bodywt_g, is_large) |>
    head()
```

    # A tibble: 6 × 4
      name                        bodywt bodywt_g is_large
      <chr>                        <dbl>    <dbl> <lgl>   
    1 Cheetah                     50        50000 TRUE    
    2 Owl monkey                   0.48       480 FALSE   
    3 Mountain beaver              1.35      1350 FALSE   
    4 Greater short-tailed shrew   0.019       19 FALSE   
    5 Cow                        600       600000 TRUE    
    6 Three-toed sloth             3.85      3850 FALSE   

Turning raw counts into counts-per-million, or a fold change into a log fold change, is [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) doing exactly this job on a bigger table.

## 16.8 `summarise()`: reduce a column to one number

[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html) collapses a whole column down to a single summary value. The mean amount of sleep across all 83 mammals:

``` downlit
msleep |>
    summarise(avg_sleep = mean(sleep_total))
```

    # A tibble: 1 × 1
      avg_sleep
          <dbl>
    1      10.4

Ask for several summaries at once and you get a one-row table with one column per answer. Alongside [`mean()`](https://rdrr.io/r/base/mean.html), the usual reducers — [`min()`](https://rdrr.io/r/base/Extremes.html), [`max()`](https://rdrr.io/r/base/Extremes.html), [`median()`](https://rdrr.io/r/stats/median.html), [`sd()`](https://rdrr.io/r/stats/sd.html), [`sum()`](https://rdrr.io/r/base/sum.html), and [`n()`](https://dplyr.tidyverse.org/reference/context.html) (which simply counts the rows) — are all fair game:

``` downlit
msleep |>
    summarise(avg_sleep = mean(sleep_total),
              min_sleep = min(sleep_total),
              max_sleep = max(sleep_total),
              n_mammals = n())
```

    # A tibble: 1 × 4
      avg_sleep min_sleep max_sleep n_mammals
          <dbl>     <dbl>     <dbl>     <int>
    1      10.4       1.9      19.9        83

That single row reports the average, the extremes, and the sample size together.

## 16.9 `group_by()`: one summary per group

[`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html) on its own gives one number for the whole table. Pair it with [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) and you get one number *per group* instead — and almost every interesting summary is a per-group summary.

> **NOTE:**
>
> [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) doesn’t change the data; it tags the rows with which group each belongs to. Picture the table sorted into piles, one pile per taxonomic order. The next [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html) runs *separately* on each pile (**apply**) and stacks the per-pile answers back into one tidy table (**combine**). Splitting the data, applying a calculation to each part, and combining the results is such a common shape that it has a name: **split-apply-combine**. Base R does it with [`aggregate()`](https://rdrr.io/r/stats/aggregate.html) and [`tapply()`](https://rdrr.io/r/base/tapply.html); dplyr just makes it read like a sentence.

Group by taxonomic order, then ask for the same summaries as before. The result has one row per order:

``` downlit
msleep |>
    group_by(order) |>
    summarise(avg_sleep = mean(sleep_total),
              min_sleep = min(sleep_total),
              max_sleep = max(sleep_total),
              n_mammals = n())
```

    # A tibble: 19 × 5
       order           avg_sleep min_sleep max_sleep n_mammals
       <chr>               <dbl>     <dbl>     <dbl>     <int>
     1 Afrosoricida        15.6       15.6      15.6         1
     2 Artiodactyla         4.52       1.9       9.1         6
     3 Carnivora           10.1        3.5      15.8        12
     4 Cetacea              4.5        2.7       5.6         3
     5 Chiroptera          19.8       19.7      19.9         2
     6 Cingulata           17.8       17.4      18.1         2
     7 Didelphimorphia     18.7       18        19.4         2
     8 Diprotodontia       12.4       11.1      13.7         2
     9 Erinaceomorpha      10.2       10.1      10.3         2
    10 Hyracoidea           5.67       5.3       6.3         3
    11 Lagomorpha           8.4        8.4       8.4         1
    12 Monotremata          8.6        8.6       8.6         1
    13 Perissodactyla       3.47       2.9       4.4         3
    14 Pilosa              14.4       14.4      14.4         1
    15 Primates            10.5        8        17          12
    16 Proboscidea          3.6        3.3       3.9         2
    17 Rodentia            12.5        7        16.6        22
    18 Scandentia           8.9        8.9       8.9         1
    19 Soricomorpha        11.1        8.4      14.9         5

Swap “order” for “treatment condition” and “sleep” for “expression” and you have a standard table from a genomics paper, produced with the same three lines.

## 16.10 Exercises

For each exercise, with ggplot2 and dplyr loaded, try writing the pipe yourself before opening the solution.

1.  **The right two verbs.** You want only the mammals that sleep *fewer than* 5 hours, showing just `name` and `sleep_total`. Which two verbs, and in what order?

    > **NOTE:**
    > ``` downlit
    > msleep |>
    >     filter(sleep_total < 5) |>
    >     select(name, sleep_total)
    > ```
    >
    >     # A tibble: 11 × 2
    >        name             sleep_total
    >        <chr>                  <dbl>
    >      1 Cow                      4  
    >      2 Roe deer                 3  
    >      3 Asian elephant           3.9
    >      4 Horse                    2.9
    >      5 Donkey                   3.1
    >      6 Giraffe                  1.9
    >      7 Pilot whale              2.7
    >      8 African elephant         3.3
    >      9 Sheep                    3.8
    >     10 Caspian seal             3.5
    >     11 Brazilian tapir          4.4
    >
    > [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) chooses the rows; [`select()`](https://dplyr.tidyverse.org/reference/select.html) chooses the columns. Either order returns the same answer, but filtering first means [`select()`](https://dplyr.tidyverse.org/reference/select.html) has fewer rows to carry — the habit that pays off on large tables.

2.  **Sleepiest carnivores.** Keep only the carnivores (`vore == "carni"`), then sort them from most to least total sleep, showing `name`, `vore`, and `sleep_total`.

    > **NOTE:**
    > ``` downlit
    > msleep |>
    >     filter(vore == "carni") |>
    >     select(name, vore, sleep_total) |>
    >     arrange(desc(sleep_total))
    > ```
    >
    >     # A tibble: 19 × 3
    >        name                       vore  sleep_total
    >        <chr>                      <chr>       <dbl>
    >      1 Thick-tailed opposum       carni        19.4
    >      2 Long-nosed armadillo       carni        17.4
    >      3 Tiger                      carni        15.8
    >      4 Northern grasshopper mouse carni        14.5
    >      5 Lion                       carni        13.5
    >      6 Domestic cat               carni        12.5
    >      7 Arctic fox                 carni        12.5
    >      8 Cheetah                    carni        12.1
    >      9 Slow loris                 carni        11  
    >     10 Jaguar                     carni        10.4
    >     11 Dog                        carni        10.1
    >     12 Red fox                    carni         9.8
    >     13 Northern fur seal          carni         8.7
    >     14 Genet                      carni         6.3
    >     15 Gray seal                  carni         6.2
    >     16 Common porpoise            carni         5.6
    >     17 Bottle-nosed dolphin       carni         5.2
    >     18 Caspian seal               carni         3.5
    >     19 Pilot whale                carni         2.7
    >
    > [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) selects the rows, [`select()`](https://dplyr.tidyverse.org/reference/select.html) trims the columns, and `arrange(desc(...))` sorts largest-first. Reading the pipe top to bottom tells the whole story.

3.  **A derived column.** Add a `rem_fraction` column (REM sleep divided by total sleep), then sort so the mammals with the *highest* REM fraction come first, and look at their `bodywt`. Do the heaviest animals sit at the top?

    > **NOTE:**
    > ``` downlit
    > msleep |>
    >     mutate(rem_fraction = sleep_rem / sleep_total) |>
    >     select(name, bodywt, rem_fraction) |>
    >     arrange(desc(rem_fraction)) |>
    >     head()
    > ```
    >
    >     # A tibble: 6 × 3
    >       name                   bodywt rem_fraction
    >       <chr>                   <dbl>        <dbl>
    >     1 European hedgehog       0.77         0.347
    >     2 Thick-tailed opposum    0.37         0.340
    >     3 Giant armadillo        60            0.337
    >     4 Tree shrew              0.104        0.292
    >     5 Dog                    14            0.287
    >     6 North American Opossum  1.7          0.272
    >
    > [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) builds the new column and `arrange(desc(...))` ranks by it. The heaviest mammals are *not* at the top, hinting that large-bodied animals spend a smaller fraction of their sleep in REM — a pattern worth a proper plot to confirm.

4.  **Summarize by group.** For each feeding type (`vore`), compute the average total sleep and the number of species. Which feeding type sleeps the most on average?

    > **NOTE:**
    > ``` downlit
    > msleep |>
    >     group_by(vore) |>
    >     summarise(avg_sleep = mean(sleep_total),
    >               n = n())
    > ```
    >
    >     # A tibble: 5 × 3
    >       vore    avg_sleep     n
    >       <chr>       <dbl> <int>
    >     1 carni       10.4     19
    >     2 herbi        9.51    32
    >     3 insecti     14.9      5
    >     4 omni        10.9     20
    >     5 <NA>        10.2      7
    >
    > Split-apply-combine: split by `vore`, apply [`mean()`](https://rdrr.io/r/base/mean.html) and [`n()`](https://dplyr.tidyverse.org/reference/context.html) to each group, combine into one row per group. The `NA` row gathers mammals whose feeding type is unrecorded — a standing reminder to check for missing values in real data.

## 16.11 Summary

You now have the core of dplyr, and every verb maps onto something you’ll do to a biological table:

- [`select()`](https://dplyr.tidyverse.org/reference/select.html) keeps the columns you want (e.g. `sample_id` and `condition` from a sample sheet).
- [`filter()`](https://dplyr.tidyverse.org/reference/filter.html) keeps the rows that pass a test (e.g. genes below a p-value cutoff).
- [`arrange()`](https://dplyr.tidyverse.org/reference/arrange.html) sorts rows, with [`desc()`](https://dplyr.tidyverse.org/reference/desc.html) for largest-first.
- [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) adds derived columns (e.g. log fold changes).
- [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html) reduces a column to a summary statistic.
- [`group_by()`](https://dplyr.tidyverse.org/reference/group_by.html) turns [`summarise()`](https://dplyr.tidyverse.org/reference/summarise.html) into one-summary-per-group — the split-apply-combine pattern at the heart of countless analyses.
- The pipe, `|>`, threads these verbs so an analysis reads top to bottom as a sequence of “and then” steps.

Master these six verbs on `msleep` and they carry over, unchanged, to count matrices, sample sheets, and results tables.

[^1]: A quantitative, theoretical framework for understanding mammalian sleep. Van M. Savage, Geoffrey B. West. Proceedings of the National Academy of Sciences Jan 2007, 104 (3) 1051-1056; DOI: [10.1073/pnas.0610080104](https://doi.org/10.1073/pnas.0610080104)
