# 14  Factors

Code

Author

Sean Davis & Lori Kern

Published

June 1, 2024

Modified

September 19, 2026

Open almost any biological dataset and you will meet columns that are not numbers but *labels*: a sample is **treated** or **control**, a mouse is **wild-type**, **heterozygous**, or **homozygous**, a tumor is stage **I**, **II**, **III**, or **IV**. These categorical variables are the backbone of experiments — they say which group each sample belongs to. R has a special tool built exactly for them: the **factor**.

You could store these labels as ordinary text, and R would let you. But the moment you want to count how many samples are in each group, run an analysis of variance, or fit a model that compares treatment to control, R needs to know that “treated” and “control” are *categories* — a fixed set of possible values — not free-form strings. A factor is how you tell it that. This chapter shows you what factors are, how they behave, and the one quirk that trips up nearly every beginner.

## 14.1 What you’ll learn

- Create a factor from a character vector with [`factor()`](https://rdrr.io/r/base/factor.html).
- Explain how a factor stores its data as integer codes plus a lookup table of **levels**.
- Use [`table()`](https://rdrr.io/r/base/table.html) to count how many observations fall in each category.
- Control the set of levels, and understand how missing values arise.
- Avoid the classic [`as.numeric()`](https://rdrr.io/r/base/numeric.html) trap that returns codes instead of labels.

## 14.2 A vector with a built-in dictionary

Imagine you are recording the country of origin for nine samples in a study. The raw data is just text:

``` downlit
# country of origin for nine samples
citizen <- c("uk", "us", "no", "au", "uk", "us", "us", "no", "au")
citizen
```

    [1] "uk" "us" "no" "au" "uk" "us" "us" "no" "au"

That is a plain character vector — nine strings, with the quotation marks R always uses to show text. Now turn it into a factor:

``` downlit
citizen_f <- factor(citizen)
citizen_f
```

    [1] uk us no au uk us us no au
    Levels: au no uk us

Two things changed. The values printed *without* quotation marks, and a new line appeared: `Levels: au no uk us`. That `Levels:` line is the heart of the whole idea.

> **IMPORTANT:**
>
> A factor *looks* like text but is stored as **integers** behind a **lookup table**. Think of a coat check at a theater: you hand over your coat and get a numbered ticket. The tickets (1, 2, 3, …) are small and easy to handle; the coats themselves stay on a rack, each hanging at the number on its ticket.
>
> A factor works the same way. The **levels** are the rack — the unique categories, stored once, in alphabetical order by default (`au`, `no`, `uk`, `us`). Each element of your data is just a *ticket number* pointing at one of those levels. So `"uk"` is not stored nine times; it is stored as the integer `3`, and R looks up level 3 to show you `uk`.

We can pull back the curtain and see the tickets directly. Every R object can carry hidden **attributes**, and a factor’s class and levels live there:

``` downlit
attributes(citizen_f)
```

    $levels
    [1] "au" "no" "uk" "us"

    $class
    [1] "factor"

``` downlit
class(citizen_f)
```

    [1] "factor"

The [`attributes()`](https://rdrr.io/r/base/attributes.html) call shows the `levels` (the rack) and confirms the object’s `class` is `"factor"`. If we strip that class away with [`unclass()`](https://rdrr.io/r/base/class.html), the disguise drops and the raw ticket numbers appear:

``` downlit
unclass(citizen_f)
```

    [1] 3 4 2 1 3 4 4 2 1
    attr(,"levels")
    [1] "au" "no" "uk" "us"

Read that output against the levels `au no uk us`. The first sample was `"uk"`, which is the **3rd** level, so its code is `3`. The third sample was `"no"`, the **2nd** level, so its code is `2`. Each number is simply a pointer into the lookup table.

## 14.3 Counting categories with `table()`

The single most useful thing you will do with a factor is count how many samples fall in each category. That is exactly what [`table()`](https://rdrr.io/r/base/table.html) does:

``` downlit
table(citizen_f)
```

    citizen_f
    au no uk us 
     2  2  2  3 

At a glance you can see the shape of your sample set: two samples each from Australia and Norway, two from the UK, and three from the US. For a real experiment this is how you check that your treatment and control groups are balanced before you analyze anything.

## 14.4 Converting a factor back

You will often need to turn a factor back into something else — and here is where people get burned. There are two ways to convert, and they give very different answers.

To recover the original labels, use [`as.character()`](https://rdrr.io/r/base/character.html):

``` downlit
as.character(citizen_f)
```

    [1] "uk" "us" "no" "au" "uk" "us" "us" "no" "au"

That gives back the text you started with — `"uk"`, `"us"`, and so on. But watch what happens if you reach for [`as.numeric()`](https://rdrr.io/r/base/numeric.html) instead:

``` downlit
as.numeric(citizen_f)
```

    [1] 3 4 2 1 3 4 4 2 1

> **WARNING:**
>
> [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a factor returns the **ticket numbers**, not the values. You get `3 4 2 1 3 4 4 2 1` — the integer codes — *not* the countries, and certainly not any number the labels might have looked like.
>
> This bites hardest when a factor’s labels are themselves numbers. Suppose a dose column reads `"10"`, `"20"`, `"50"` but was accidentally stored as a factor. [`as.numeric()`](https://rdrr.io/r/base/numeric.html) would return `1, 2, 3` (the codes), not `10, 20, 50` (the doses) — a silent, dangerous error. The safe recipe when labels are numeric is to go through character first: `as.numeric(as.character(dose_f))`.

## 14.5 Controlling the levels

By default, the levels are every unique value in your data, sorted alphabetically. But you are allowed to specify the levels yourself — and this is genuinely useful. Maybe you only care about the US and UK samples for one analysis:

``` downlit
citizen_f2 <- factor(citizen, levels = c("us", "uk"))
citizen_f2
```

    [1] uk   us   <NA> <NA> uk   us   us   <NA> <NA>
    Levels: us uk

Look closely: the values that were not in your chosen levels — the `no` and `au` samples — turned into `<NA>`, R’s marker for a missing value. By naming only `us` and `uk` as legal levels, you told R that every other value has no category to belong to, so it became missing.

``` downlit
table(citizen_f2)
```

    citizen_f2
    us uk 
     3  2 

By default, [`table()`](https://rdrr.io/r/base/table.html) quietly drops those missing values, so the count reflects only the four samples that matched. If you *want* to see the missing ones, add them as an explicit level with [`addNA()`](https://rdrr.io/r/base/factor.html):

``` downlit
table(addNA(citizen_f2))
```


      us   uk <NA> 
       3    2    4 

Now the `<NA>` category appears with its own count of 5 — the five samples that fell outside your chosen levels.

> **NOTE:**
>
> Setting the levels yourself does more than filter. It also fixes their **order**, which controls how categories appear in tables and plots and which group a model treats as the baseline. For a treatment study you might write `factor(group, levels = c("control", "treated"))` so that “control” is the reference everything else is compared against. The default settings are not always what your analysis needs — knowing what they are, and overriding them on purpose, is part of doing the analysis correctly.

## 14.6 Exercises

1.  **Genotypes.** A small genotyping study records these calls for eight mice: `"wt"`, `"het"`, `"hom"`, `"wt"`, `"wt"`, `"het"`, `"hom"`, `"wt"`. Make a factor from this vector and use [`table()`](https://rdrr.io/r/base/table.html) to count how many mice carry each genotype.

    > **NOTE:**
    > ``` downlit
    > genotype <- c("wt", "het", "hom", "wt", "wt", "het", "hom", "wt")
    > genotype_f <- factor(genotype)
    > table(genotype_f)
    > ```
    >
    >     genotype_f
    >     het hom  wt 
    >       2   2   4 
    >
    > [`table()`](https://rdrr.io/r/base/table.html) counts each category for you: four `wt`, two `het`, and two `hom`. The factor turned a list of labels into a tidy summary of the sample set.

2.  **Read the codes.** Using the `genotype_f` factor from Exercise 1, predict what `unclass(genotype_f)` will print *before* you run it. Then run it and check.

    > **NOTE:**
    > ``` downlit
    > unclass(genotype_f)
    > ```
    >
    >     [1] 3 1 2 3 3 1 2 3
    >     attr(,"levels")
    >     [1] "het" "hom" "wt" 
    >
    > The levels are alphabetical — `het`, `hom`, `wt` — so `het` is code `1`, `hom` is `2`, and `wt` is `3`. The first mouse (`wt`) prints as `3`. Each integer is a ticket pointing into the levels table, not a measurement.

3.  **The numeric trap.** A colleague stored a dose column as a factor with labels `"5"`, `"10"`, `"20"` and ran [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on it to get the doses back. Why are their numbers wrong, and what should they do instead?

    > **NOTE:**
    > ``` downlit
    > dose_f <- factor(c("5", "10", "20", "5", "20"))
    > as.numeric(dose_f)               # WRONG: returns the integer codes
    > ```
    >
    >     [1] 3 1 2 3 2
    >
    > ``` downlit
    > as.numeric(as.character(dose_f)) # RIGHT: recovers the real doses
    > ```
    >
    >     [1]  5 10 20  5 20
    >
    > [`as.numeric()`](https://rdrr.io/r/base/numeric.html) returns the level codes, and because levels are sorted, `"10"` sorts before `"5"` — so the codes do not even match the order you expect. Converting to character first turns each label back into its text, which [`as.numeric()`](https://rdrr.io/r/base/numeric.html) can then read as the true number.

4.  **Set the reference group.** You have a `group` vector of `"treated"` and `"control"` labels. Build a factor whose levels are ordered so that `control` comes first (the baseline a model would compare against).

    > **NOTE:**
    > ``` downlit
    > group <- c("treated", "control", "treated", "control", "treated")
    > group_f <- factor(group, levels = c("control", "treated"))
    > group_f
    > ```
    >
    >     [1] treated control treated control treated
    >     Levels: control treated
    >
    > By naming `control` first in `levels`, you make it level 1 — the reference category. Most modeling functions in R compare every other group against the first level, so setting it deliberately keeps your results interpretable.

## 14.7 Summary

You can now recognize and build the categorical variables at the center of nearly every biological experiment:

- **Create** a factor from a character vector with [`factor()`](https://rdrr.io/r/base/factor.html), and recognize it by the `Levels:` line and the missing quotation marks.
- **Picture** it as codes plus a lookup table — integer tickets pointing at a rack of levels — which explains every behavior in this chapter.
- **Count** categories with [`table()`](https://rdrr.io/r/base/table.html) to check group sizes and spot imbalance.
- **Control** the levels to filter categories, set a reference group, and decide how missing values are handled.
- **Convert** safely: [`as.character()`](https://rdrr.io/r/base/character.html) gives the labels back, while [`as.numeric()`](https://rdrr.io/r/base/numeric.html) gives the codes — so reach for `as.numeric(as.character(x))` whenever a factor’s labels are really numbers.

Keep the coat-check picture in mind and factors stop being mysterious: a factor is just a vector of tickets with a table of what each ticket stands for.
