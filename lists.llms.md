# 12  Lists: a catch-all container

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

So far in our journey through R’s data structures, we’ve dealt with vectors and matrices. These are fantastic tools, but they have one strict rule: all their elements must be of the *same data type*. You can have a vector of numbers or a matrix of characters, but you can’t mix and match.

But what about real-world biological data? A single experiment can generate a dizzying variety of information. Imagine you’re studying a particular gene. You might have:

- The gene’s name (text).
- Its expression level across several samples (a set of numbers).
- A record of whether it’s a known cancer-related gene (a simple TRUE/FALSE).
- The raw fluorescence values from your qPCR machine (a matrix of numbers).
- Some personal notes about the experiment (a paragraph of text).

How could you possibly store all of this related, yet different, information together? You could create many separate variables, but that would be clunky and hard to manage. This is exactly the problem that **lists** are designed to solve.

A list in R is like a flexible, multi-compartment container. It’s a single object that can hold a collection of *other* R objects, and those objects can be of any type, length, or dimension. You can put vectors, matrices, logical values, and even other lists inside a single list. This makes them one of the most fundamental and powerful data structures for bioinformatics analysis. In fact, almost every rich result R hands you back — a fitted model, the output of a statistical test, a whole data frame — is a list under the hood. Learn lists and you’ve learned how to take those results apart.

The key features of lists are:

- **Flexibility**: They can contain a mix of any data type.
- **Organization**: You can and should *name* the elements of a list, making your data self-describing.
- **Hierarchy**: Because lists can contain other lists, you can create complex, nested data structures to represent sophisticated relationships in your data.

## 12.1 What you’ll learn

- Create a list with [`list()`](https://rdrr.io/r/base/list.html) and give its elements descriptive names.
- Inspect a list’s contents with [`str()`](https://rdrr.io/r/utils/str.html), [`length()`](https://rdrr.io/r/base/length.html), [`names()`](https://rdrr.io/r/base/names.html), and [`class()`](https://rdrr.io/r/base/class.html).
- Pull a single element *out* of a list with `[[ ]]` or `$`.
- Take a smaller *sub-list* out with single `[ ]`, and explain why that is different.
- Add, update, and remove list elements.
- Use a list to build a self-contained record for a gene.

## 12.2 Creating a list

You create a list with the [`list()`](https://rdrr.io/r/base/list.html) function. The best practice is to name the elements as you create them. This makes your code infinitely more readable and your data easier to work with.

Let’s create a list to store the information for our hypothetical gene study.

``` downlit
# An experiment tracking list for the gene TP53
experiment_data <- list(
  experiment_id = "EXP042",
  gene_name = "TP53",
  read_counts = c(120, 155, 98, 210),
  is_control = FALSE,
  sample_matrix = matrix(1:4, nrow = 2, dimnames = list(c("Treated", "Untreated"), c("Replicate1", "Replicate2")))
)

print(experiment_data)
```

    $experiment_id
    [1] "EXP042"

    $gene_name
    [1] "TP53"

    $read_counts
    [1] 120 155  98 210

    $is_control
    [1] FALSE

    $sample_matrix
              Replicate1 Replicate2
    Treated            1          3
    Untreated          2          4

Notice how the printout walks through the list one named compartment at a time: `$experiment_id` holds a single string, `$read_counts` holds a four-number vector, and `$sample_matrix` holds an entire 2×2 matrix — all living happily inside one object.

> **NOTE:**
>
> The [`print()`](https://rdrr.io/r/base/print.html) function displays the contents of an R object in the console. For a list, it shows each element and its contents. It’s the default action when you just type the variable’s name and hit Enter, so `print(experiment_data)` and simply typing `experiment_data` do the same thing.

## 12.3 Inspecting your list: what’s inside?

When someone hands you a tube in the lab, the first thing you do is look at the label. When R gives you a complex object like a list, you need to do the same. R provides several “introspection” functions to help you understand the contents and structure of your lists.

### 12.3.1 `str()`: the structure function

This is arguably the most useful function for inspecting any R object, especially lists.

``` downlit
str(experiment_data)
```

    List of 5
     $ experiment_id: chr "EXP042"
     $ gene_name    : chr "TP53"
     $ read_counts  : num [1:4] 120 155 98 210
     $ is_control   : logi FALSE
     $ sample_matrix: int [1:2, 1:2] 1 2 3 4
      ..- attr(*, "dimnames")=List of 2
      .. ..$ : chr [1:2] "Treated" "Untreated"
      .. ..$ : chr [1:2] "Replicate1" "Replicate2"

The output of [`str()`](https://rdrr.io/r/utils/str.html) tells us everything we need to know: it’s a “List of 5”, and for each of the 5 elements, it shows the name (e.g., `experiment_id`), the data type (e.g., `chr` for character, `num` for numeric), and a preview of the content.

> **TIP:**
>
> [`str()`](https://rdrr.io/r/utils/str.html) provides a compact, human-readable summary of an object’s internal **str**ucture. Whenever a function hands you back something you don’t recognize — a list, a model, a complicated result — run [`str()`](https://rdrr.io/r/utils/str.html) on it first. It’s the fastest way to see what’s inside without printing pages of output.

### 12.3.2 `length()`, `names()`, and `class()`

These functions give you more specific information about the list itself.

``` downlit
# How many top-level elements does the list contain?
length(experiment_data)
```

    [1] 5

``` downlit
# What are the elements called?
names(experiment_data)
```

    [1] "experiment_id" "gene_name"     "read_counts"   "is_control"   
    [5] "sample_matrix"

``` downlit
# What kind of object is this?
class(experiment_data)
```

    [1] "list"

Read together, these tell a story: the list has 5 elements ([`length()`](https://rdrr.io/r/base/length.html)), they are named `experiment_id`, `gene_name`, and so on ([`names()`](https://rdrr.io/r/base/names.html)), and the object as a whole is a `list` ([`class()`](https://rdrr.io/r/base/class.html)). The [`names()`](https://rdrr.io/r/base/names.html) output is especially handy — it’s a menu of everything you’re allowed to ask for.

## 12.4 Accessing list elements: getting things out

Okay, you’ve packed your experimental data into a list. Now, how do you get specific items out? This is a critical concept, and R has a few ways to do it, each with a distinct purpose.

### 12.4.1 The mighty `[[ ]]` and `$` for single items

To pull out a *single element* from a list in its original form, you use either double square brackets `[[ ]]` or the dollar sign `$` (for named lists). Think of this as carefully reaching into a specific compartment of your container and taking out the item itself.

Let’s use our `experiment_data` list.

``` downlit
# Get the gene name using [[ ]]
gene <- experiment_data[["gene_name"]]
print(gene)
```

    [1] "TP53"

``` downlit
class(gene) # It's a character vector, just as it was when we put it in.
```

    [1] "character"

``` downlit
# Get the read counts using the $ shortcut. This is often easier to read.
reads <- experiment_data$read_counts
print(reads)
```

    [1] 120 155  98 210

``` downlit
class(reads) # It's a numeric vector.
```

    [1] "numeric"

``` downlit
# The [[ ]] has a neat trick: you can use a variable to specify the name.
element_to_get <- "read_counts"
experiment_data[[element_to_get]]
```

    [1] 120 155  98 210

The key takeaway is that `[[ ]]` and `$` **extract the element**. The result is the object that was stored inside the list — a plain character vector, a plain numeric vector — ready to be used in calculations. Notice that `$` needs the name typed out literally, while `[[ ]]` can take a name stored in a variable, which is invaluable when you write loops over a list’s elements.

### 12.4.2 The subsetting `[ ]` for new lists

The single square bracket `[ ]` behaves differently. It always returns a *new, smaller list* that is a subset of the original list. It’s like taking a whole compartment, label and all, out of your larger container — the item stays wrapped inside its list.

``` downlit
# Get the gene name using [ ]
gene_sublist <- experiment_data["gene_name"]

print(gene_sublist)
```

    $gene_name
    [1] "TP53"

``` downlit
# Note the class!
class(gene_sublist)
```

    [1] "list"

See the difference? `experiment_data[["gene_name"]]` gave us a `"character"` vector, but `experiment_data["gene_name"]` gives us back a `"list"` that still has `gene_name` inside it.

> **IMPORTANT:**
>
> This single-versus-double-bracket distinction is the most common source of confusion for R beginners — so it’s worth burning in:
>
> - **`[[ ]]` and `$` reach in and pull the item *out*.** You get the vector, matrix, or value itself.
> - **`[ ]` returns a smaller *list*.** The item is still wrapped inside a list.
>
> This matters the moment you try to compute. To take the [`mean()`](https://rdrr.io/r/base/mean.html) of the read counts you must *extract* them:
>
> ``` downlit
> mean(experiment_data[["read_counts"]])   # works: operating on a numeric vector
> ```
>
>     [1] 145.75
>
> But `mean(experiment_data["read_counts"])` would fail, because you can’t take the mean of a *list*. When in doubt, ask yourself: do I want the thing itself (`[[ ]]`), or a smaller list (`[ ]`)?

## 12.5 Modifying lists

Your data is rarely static. You can easily add, remove, or update elements in a list after you’ve created it.

### 12.5.1 Adding and updating elements

You can add a new element or change an existing one by using the `$` or `[[ ]]` assignment syntax. If the name already exists, R updates it; if it’s new, R adds it.

``` downlit
# Add the date of the experiment
experiment_data$date <- "2024-06-05"

# Add some notes using the [[ ]] syntax
experiment_data[["notes"]] <- "Initial pilot experiment. High variance in read counts."

# Let's update the control status
experiment_data$is_control <- TRUE

# Let's look at the structure now
str(experiment_data)
```

    List of 7
     $ experiment_id: chr "EXP042"
     $ gene_name    : chr "TP53"
     $ read_counts  : num [1:4] 120 155 98 210
     $ is_control   : logi TRUE
     $ sample_matrix: int [1:2, 1:2] 1 2 3 4
      ..- attr(*, "dimnames")=List of 2
      .. ..$ : chr [1:2] "Treated" "Untreated"
      .. ..$ : chr [1:2] "Replicate1" "Replicate2"
     $ date         : chr "2024-06-05"
     $ notes        : chr "Initial pilot experiment. High variance in read counts."

The list has grown from 5 elements to 7 — `date` and `notes` are new — and `is_control` has flipped from `FALSE` to `TRUE`. Same syntax, two different jobs: a name that already exists gets overwritten, a name that doesn’t gets created.

### 12.5.2 Removing elements

To remove an element from a list, you simply assign `NULL` to it. `NULL` is R’s special object representing nothingness.

``` downlit
# We've decided the matrix isn't needed for this summary object.
experiment_data$sample_matrix <- NULL

# See the final structure of our list
str(experiment_data)
```

    List of 6
     $ experiment_id: chr "EXP042"
     $ gene_name    : chr "TP53"
     $ read_counts  : num [1:4] 120 155 98 210
     $ is_control   : logi TRUE
     $ date         : chr "2024-06-05"
     $ notes        : chr "Initial pilot experiment. High variance in read counts."

The `sample_matrix` element is gone, and the list is back down to 6 elements. Assigning `NULL` is R’s idiom for “delete this compartment entirely.”

## 12.6 A biological example: a self-contained gene record

Let’s put this all together. Lists are perfect for creating self-contained records that you can easily pass to functions or combine into larger lists.

``` downlit
brca1_gene <- list(
  gene_symbol = "BRCA1",
  full_name = "BRCA1 DNA repair associated",
  chromosome = "17",
  expression_log2 = log2(c(45, 50, 30, 88, 120)),
  related_diseases = c("Breast Cancer", "Ovarian Cancer")
)

# Now we can easily work with this structured information
cat("Analyzing gene:", brca1_gene$gene_symbol, "\n")
```

    Analyzing gene: BRCA1 

``` downlit
cat("Located on chromosome:", brca1_gene$chromosome, "\n")
```

    Located on chromosome: 17 

``` downlit
# Calculate the average log2 expression
avg_expression <- mean(brca1_gene$expression_log2)
cat("Average log2 expression:", avg_expression, "\n")
```

    Average log2 expression: 5.881784 

Each [`cat()`](https://rdrr.io/r/base/cat.html) line prints a label followed by a value pulled straight out of the list with `$`, and the average log2 expression comes out around 5.9.

> **NOTE:**
>
> - **[`cat()`](https://rdrr.io/r/base/cat.html)** concatenates and prints its arguments to the console. Unlike [`print()`](https://rdrr.io/r/base/print.html), it lets you seamlessly join text and variables on one line, and the `"\n"` character adds a newline (a line break).
> - **[`log2()`](https://rdrr.io/r/base/Log.html)** calculates the base-2 logarithm. It’s very common in gene-expression analysis to transform skewed count data with [`log2()`](https://rdrr.io/r/base/Log.html) to make it more symmetric and easier to model.

This simple `brca1_gene` list is now a complete, portable record. You could imagine creating a list of these gene records, building a powerful, hierarchical database for your entire project.

## 12.7 Exercises

For these exercises, keep the `brca1_gene` list from above handy. Each one practices a move you’ll use constantly when pulling apart analysis results.

1.  **Extract by name.** Pull the `related_diseases` element *out* of `brca1_gene` so that you end up with a plain character vector (not a list). Confirm its class.

    > **NOTE:**
    > ``` downlit
    > diseases <- brca1_gene[["related_diseases"]]
    > diseases
    > ```
    >
    >     [1] "Breast Cancer"  "Ovarian Cancer"
    >
    > ``` downlit
    > class(diseases)
    > ```
    >
    >     [1] "character"
    >
    > Using `[[ ]]` (or `brca1_gene$related_diseases`) extracts the element itself, so [`class()`](https://rdrr.io/r/base/class.html) reports `"character"`. Had you used single `[ ]`, you’d have gotten a one-element *list* back instead.

2.  **Add an element.** A collaborator tells you BRCA1 sits on the *minus* strand. Add a new element called `strand` with the value `"-"` to `brca1_gene`, then confirm it’s there.

    > **NOTE:**
    > ``` downlit
    > brca1_gene$strand <- "-"
    > names(brca1_gene)
    > ```
    >
    >     [1] "gene_symbol"      "full_name"        "chromosome"       "expression_log2" 
    >     [5] "related_diseases" "strand"          
    >
    > Assigning to a name that doesn’t yet exist *adds* a new compartment. Checking [`names()`](https://rdrr.io/r/base/names.html) is a quick way to confirm the addition landed.

3.  **Nest a vector inside a list.** Add an element called `samples` holding the character vector `c("Tumor", "Normal", "Tumor")`. Then count how many samples are tumors.

    > **NOTE:**
    > ``` downlit
    > brca1_gene$samples <- c("Tumor", "Normal", "Tumor")
    > sum(brca1_gene$samples == "Tumor")
    > ```
    >
    >     [1] 2
    >
    > A list element can be a whole vector. Once you *extract* it with `$`, it behaves like any ordinary vector — here we compare it to `"Tumor"` and [`sum()`](https://rdrr.io/r/base/sum.html) the `TRUE`s to get 2.

4.  **Spot the bracket bug.** A labmate runs `mean(brca1_gene["expression_log2"])` and gets an error. Explain why, then fix it.

    > **NOTE:**
    > ``` downlit
    > mean(brca1_gene[["expression_log2"]])
    > ```
    >
    >     [1] 5.881784
    >
    > Single `[ ]` returns a *list* containing the vector, and [`mean()`](https://rdrr.io/r/base/mean.html) can’t average a list. Switching to `[[ ]]` (or `$`) extracts the numeric vector itself, so [`mean()`](https://rdrr.io/r/base/mean.html) works.

## 12.8 Summary

You can now reach for a list whenever you need to keep mixed, named things together in one object — which is exactly what a real analysis result looks like. Specifically, you can:

- **Build** a named list with `list(name = value, ...)`.
- **Inspect** it with [`str()`](https://rdrr.io/r/utils/str.html) (the structure), [`length()`](https://rdrr.io/r/base/length.html), [`names()`](https://rdrr.io/r/base/names.html), and [`class()`](https://rdrr.io/r/base/class.html).
- **Extract** a single element with `[[ ]]` or `$` (you get the item itself), and take a **sub-list** with single `[ ]` (you get a smaller list).
- **Modify** a list by assigning to `$name` / `[["name"]]` to add or update, and by assigning `NULL` to remove.
- **Assemble** a self-contained gene record and pull values back out for calculation.

The one idea to carry forward: `[[ ]]` reaches in and hands you the thing; `[ ]` hands you a smaller list. Get that straight, and lists — and all the model and test results built on them — stop being mysterious.
