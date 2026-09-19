# 15  Control statements

Code

Author

Lori Kern

Published

June 1, 2024

Modified

September 19, 2026

So far your code has run straight down the page: every line executes, once, in order. Real analyses rarely stay that simple. You want R to *decide* — drop a sample only if its read count falls below a threshold — and to *repeat* — run the same quality check across every one of a hundred samples without writing the line a hundred times. **Control statements** are how you express those decisions and repetitions.

This chapter introduces the two workhorses, conditionals (`if` / `else`) and loops (`for`, `while`, `repeat`), along with the keywords that steer them ([`break`](https://rdrr.io/r/base/Control.html), [`next`](https://rdrr.io/r/base/Control.html), `return`) and a first look at error handling with [`tryCatch()`](https://rdrr.io/r/base/conditions.html). We close with a note on when *not* to write a loop at all — because in R, a vectorized function is usually the better tool.

## 15.1 What you’ll learn

- Branch your code with `if`, `if`/`else`, and the vectorized [`ifelse()`](https://rdrr.io/r/base/ifelse.html).
- Repeat work with `for`, `while`, and `repeat` loops, and choose the right one.
- Steer a loop with [`break`](https://rdrr.io/r/base/Control.html), [`next`](https://rdrr.io/r/base/Control.html), and `return`.
- Handle errors and warnings gracefully with [`tryCatch()`](https://rdrr.io/r/base/conditions.html).
- Recognize when a vectorized function or an [`apply()`](https://rdrr.io/r/base/apply.html)-family call beats a loop.

> **IMPORTANT:**
>
> Control statements lean heavily on punctuation. Every open parenthesis `(` needs a matching close `)`, and every open brace `{` needs a matching `}`. A single missing bracket is the most common reason a control statement refuses to run, so keep them balanced.

> **TIP:**
>
> Once a control statement opens a [`{ }`](https://rdrr.io/r/base/Paren.html) block, indent everything inside it. The indentation has no effect on R, but it makes it obvious at a glance which lines belong to the block — which matters enormously once you start nesting one statement inside another.

## 15.2 Conditional statements

### 15.2.1 `if`

An `if` statement evaluates an expression and, **only if** the result is `TRUE`, runs the block of code in its braces.

**Syntax:**

``` downlit
if (expression) {
  # code that runs only when expression is TRUE
  ...
}
```

The expression in parentheses must boil down to a single logical value — `TRUE` or `FALSE` — and that value decides whether the braced code runs.

**Example:**

``` downlit
x <- 12
if (x > 0) {
  message(x, " is greater than 0")
  x <- 0
}
```

    12 is greater than 0

``` downlit
x
```

    [1] 0

> **NOTE:**
>
> When the condition is `TRUE`, *all* the lines in the braces execute — including the assignment that reset `x` to `0`. That’s why printing `x` afterward shows `0`, not `12`.

### 15.2.2 `if`/`else`

An `if`/`else` statement adds a second block that runs when the expression is `FALSE`. Exactly one of the two blocks always runs.

**Syntax:**

``` downlit
if (expression) {
  # code that runs when expression is TRUE
  ...
} else {
  # code that runs when expression is FALSE
  ...
}
```

Read it aloud: *if* the expression is true, run this code; *else* (otherwise) run that other code.

**Example:**

``` downlit
if (x > 0) {
  message(x, " is greater than 0")
  x <- 0
} else {
  message(x, " is not greater than 0")
  x <- x + 2
}
```

    0 is not greater than 0

``` downlit
x
```

    [1] 2

Because `x` was `0` going in, the condition `x > 0` is `FALSE`, so the `else` branch runs and `x` becomes `2`.

### 15.2.3 `ifelse()`

[`ifelse()`](https://rdrr.io/r/base/ifelse.html) is a *vectorized* cousin of `if`/`else`. Instead of testing one value, it tests an entire vector at once and returns a vector of results — one “yes” or “no” value per element.

**Syntax:**

``` downlit
ifelse(test_expression, yes_value, no_value)
```

**Example:**

``` downlit
num_vec <- -3:3
ifelse(num_vec >= 0, "positive", "negative")
```

    [1] "negative" "negative" "negative" "positive" "positive" "positive" "positive"

> **WARNING:**
>
> A plain `if` expects a *single* `TRUE`/`FALSE`. Hand it a whole vector and, in recent versions of R, it errors. When you want to test many values at once — a column of a data frame, say — reach for [`ifelse()`](https://rdrr.io/r/base/ifelse.html).

## 15.3 Loops

Loops repeat a block of code: a fixed number of times, once for each element of a collection, or until some condition is met.

### 15.3.1 `for`

A `for` loop runs its body once for each element of a vector or list, assigning the current element to a variable you can use inside the loop.

**Syntax:**

``` downlit
for (value in vector) {
  # code that runs once per element; `value` holds the current one
  ...
}
```

**Examples:**

Here we loop over a vector of names, printing how many characters each one has:

``` downlit
names <- c("Donna", "John", "Bradley", "Kara")
for (nm in names) {
  print(paste(nm, "has", nchar(nm), "letters"))
}
```

    [1] "Donna has 5 letters"
    [1] "John has 4 letters"
    [1] "Bradley has 7 letters"
    [1] "Kara has 4 letters"

Here, for each value `1` through `5`, we add it to a running total. Notice how `x` is updated on every pass:

``` downlit
x <- 0
for (i in 1:5) {
  print(paste("add", i, "to", x))
  x <- x + i
}
```

    [1] "add 1 to 0"
    [1] "add 2 to 1"
    [1] "add 3 to 3"
    [1] "add 4 to 6"
    [1] "add 5 to 10"

``` downlit
x
```

    [1] 15

We can also loop over the elements of a list. [`seq_along()`](https://rdrr.io/r/base/seq.html) gives us the positions `1, 2, 3, ...`, which we use both to pull each element out and to look up its name:

``` downlit
my_list <- list(
  people  = names,
  ages    = c(54, 78, 40, 5, 25),
  animals = c("dog", "fish")
)
for (i in seq_along(my_list)) {
  print(paste("List element", names(my_list)[i], "contains",
              length(my_list[[i]]), "values"))
}
```

    [1] "List element people contains 4 values"
    [1] "List element ages contains 5 values"
    [1] "List element animals contains 2 values"

### 15.3.2 `while`

A `while` loop keeps running *as long as* its condition stays `TRUE`. It checks the condition before each pass, so you must change something inside the loop that eventually makes the condition `FALSE`.

**Syntax:**

``` downlit
while (expression) {
  # code that runs while expression is TRUE
  # remember to update the variable the expression tests
  ...
}
```

**Example:**

Starting at `1`, we print the value and add `1` each pass, stopping once the value passes `5`:

``` downlit
value <- 1
while (value <= 5) {
  print(value)
  value <- value + 1
}
```

    [1] 1
    [1] 2
    [1] 3
    [1] 4
    [1] 5

> **WARNING:**
>
> A `while` loop only ends when its condition becomes `FALSE`. If nothing inside the loop ever makes that happen, the loop runs forever — an **infinite loop**. Always make sure the body updates the variable the condition tests. (If you do get stuck in one, press `Esc` or the stop button to interrupt R.)

### 15.3.3 `repeat`

A `repeat` loop has no condition of its own — it runs forever until a [`break`](https://rdrr.io/r/base/Control.html) statement stops it. That makes [`break`](https://rdrr.io/r/base/Control.html) mandatory.

**Syntax:**

``` downlit
repeat {
  # code to evaluate
  if (condition) {
    break
  }
}
```

**Example:**

Here we repeat ourselves until we’ve done it the set number of times. The body increments `i`, and the `if` checks whether we’ve hit the limit:

``` downlit
i <- 0
times <- 3
repeat {
  print("I am repeating myself")
  i <- i + 1
  if (i == times) {
    break
  }
}
```

    [1] "I am repeating myself"
    [1] "I am repeating myself"
    [1] "I am repeating myself"

## 15.4 Steering a loop: `break`, `next`, and `return`

These keywords let you change a loop’s flow from the inside. They’re especially handy in longer or nested statements.

### 15.4.1 `break`

We’ve already seen [`break`](https://rdrr.io/r/base/Control.html) with `repeat`. It stops the loop immediately and exits the moment it’s reached — no further iterations run.

### 15.4.2 `next`

[`next`](https://rdrr.io/r/base/Control.html) skips the rest of the *current* iteration and jumps straight to the next one, without ending the loop. Here we print only the even numbers from 1 to 10, using [`next`](https://rdrr.io/r/base/Control.html) to skip the odd ones:

``` downlit
for (i in 1:10) {
  if (i %% 2 != 0) {
    next
  }
  print(i)
}
```

    [1] 2
    [1] 4
    [1] 6
    [1] 8
    [1] 10

### 15.4.3 `return`

`return` exits the surrounding **function** immediately and hands back a value. It’s used inside functions rather than bare loops.

**Syntax:**

``` downlit
return(expression)
```

This function takes an argument `x`. If `x` is `0` it returns `"zero"`. Otherwise it adds `4`; if the result is `0` or less it returns that value, and if not it returns the value doubled:

``` downlit
classify <- function(x) {
  if (x == 0) {
    return("zero")
  }
  x <- x + 4
  if (x <= 0) {
    return(x)
  } else {
    return(x * 2)
  }
}

classify(0)
```

    [1] "zero"

``` downlit
classify(-8)
```

    [1] -4

``` downlit
classify(6)
```

    [1] 20

## 15.5 Other useful patterns

### 15.5.1 Nesting

You can place one control statement inside another — we already did, when the `repeat` loop held an `if` block. Any statement can nest, to any depth. A classic example is a pair of `for` loops walking the rows and columns of a matrix.

Here we build a numeric matrix with 5 rows and 3 columns, filled with `1` through `15` by column, then loop over it row by row and print each cell:

``` downlit
mat <- matrix(1:15, ncol = 3)
mat
```

         [,1] [,2] [,3]
    [1,]    1    6   11
    [2,]    2    7   12
    [3,]    3    8   13
    [4,]    4    9   14
    [5,]    5   10   15

``` downlit
for (i in seq(nrow(mat))) {
  for (j in seq(ncol(mat))) {
    print(mat[i, j])
  }
}
```

    [1] 1
    [1] 6
    [1] 11
    [1] 2
    [1] 7
    [1] 12
    [1] 3
    [1] 8
    [1] 13
    [1] 4
    [1] 9
    [1] 14
    [1] 5
    [1] 10
    [1] 15

> **TIP:**
>
> It’s tempting to use `r` and `c` for “row” and “column”, but `c` is also R’s function for combining values (`c(1, 2, 3)`). Naming a variable `c` shadows that function and invites confusing bugs. Using `i` and `j` keeps the built-in [`c()`](https://rdrr.io/r/base/c.html) available and is the conventional choice for loop counters.

### 15.5.2 `try()` / `tryCatch()`

Not strictly a control statement, but it fits the theme: [`tryCatch()`](https://rdrr.io/r/base/conditions.html) lets you *handle* code that might throw an error or a warning, instead of letting it halt your whole script. This matters when you’re looping over many files or samples and don’t want one bad case to stop the rest.

**Syntax:**

``` downlit
tryCatch(
  expr,
  error = function(e) {
    # what to do if expr raises an error
  },
  warning = function(w) {
    # what to do if expr raises a warning
  },
  finally = {
    # code that runs no matter what
  }
)
```

You can supply any subset of `error`, `warning`, and `finally`. `expr` is the code you’re attempting. The `error` handler runs if `expr` fails; its argument `e` is an error object carrying the details. The `warning` handler works the same way for warnings, via `w`. The `finally` block runs regardless of the outcome — which makes it the right place to close a file or database connection that `expr` may have opened.

### 15.5.3 Vectorization and the `apply()` family

Here’s the twist: a lot of the work people reach for loops to do, R can do *without* a loop — and faster. Many functions are already **vectorized**, meaning they operate on a whole vector at once. To add `1` to every element of a vector, you don’t loop; you just write `x + 1`.

``` downlit
x <- 1:5
x + 1          # vectorized: no loop needed
```

    [1] 2 3 4 5 6

``` downlit
mean(x)        # already summarizes the whole vector
```

    [1] 3

When there isn’t a ready-made vectorized function, the [`apply()`](https://rdrr.io/r/base/apply.html) family applies a function across a structure for you: [`apply()`](https://rdrr.io/r/base/apply.html) over the rows or columns of a matrix, and [`lapply()`](https://rdrr.io/r/base/lapply.html), [`sapply()`](https://rdrr.io/r/base/lapply.html), [`mapply()`](https://rdrr.io/r/base/mapply.html), and [`tapply()`](https://rdrr.io/r/base/tapply.html) over vectors and lists. Before writing a loop, it’s worth asking whether a vectorized function or an [`apply()`](https://rdrr.io/r/base/apply.html) call would be clearer and quicker.

## 15.6 Exercises

1.  **Even or odd.** Write an `if`/`else` statement that prints `"even"` if a number `n` is even and `"odd"` otherwise. (Hint: `n %% 2` gives the remainder after dividing by 2.)

    > **NOTE:**
    > ``` downlit
    > n <- 7
    > if (n %% 2 == 0) {
    >   print("even")
    > } else {
    >   print("odd")
    > }
    > ```
    >
    >     [1] "odd"
    >
    > `n %% 2` is `0` for even numbers and `1` for odd ones, so testing it against `0` tells the two apart.

2.  **Sum a vector with a loop.** Using a `for` loop, add up the numbers in `c(4, 8, 15, 16, 23, 42)`. Then check your answer with the built-in [`sum()`](https://rdrr.io/r/base/sum.html).

    > **NOTE:**
    > ``` downlit
    > nums <- c(4, 8, 15, 16, 23, 42)
    > total <- 0
    > for (n in nums) {
    >   total <- total + n
    > }
    > total
    > ```
    >
    >     [1] 108
    >
    > ``` downlit
    > sum(nums)   # the vectorized way — same answer, far less typing
    > ```
    >
    >     [1] 108
    >
    > The loop accumulates the running total one element at a time; [`sum()`](https://rdrr.io/r/base/sum.html) does the same thing in one vectorized call. This is exactly the kind of task where a vectorized function beats a hand-written loop.

3.  **Skip the missing values.** Loop over `c(3, NA, 7, NA, 2)` and print only the values that are *not* `NA`, using [`next`](https://rdrr.io/r/base/Control.html) to skip the missing ones. (Hint: `is.na(x)` is `TRUE` when `x` is missing.)

    > **NOTE:**
    > ``` downlit
    > vals <- c(3, NA, 7, NA, 2)
    > for (v in vals) {
    >   if (is.na(v)) {
    >     next
    >   }
    >   print(v)
    > }
    > ```
    >
    >     [1] 3
    >     [1] 7
    >     [1] 2
    >
    > When `v` is `NA`, [`next`](https://rdrr.io/r/base/Control.html) jumps to the following iteration before the [`print()`](https://rdrr.io/r/base/print.html) line is reached, so only the real numbers are printed.

## 15.7 Summary

You can now control the flow of your code rather than running straight down the page:

- **Branch** with `if` and `if`/`else` for single values, and [`ifelse()`](https://rdrr.io/r/base/ifelse.html) to test a whole vector at once.
- **Repeat** with `for` (once per element), `while` (until a condition flips), and `repeat` (until you [`break`](https://rdrr.io/r/base/Control.html)).
- **Steer** loops with [`break`](https://rdrr.io/r/base/Control.html) (stop), [`next`](https://rdrr.io/r/base/Control.html) (skip), and `return` (exit a function with a value).
- **Guard** risky code with [`tryCatch()`](https://rdrr.io/r/base/conditions.html) so one error doesn’t sink an entire run.
- **Reach for vectorization first.** Many tasks that look like they need a loop are better done with a vectorized function or an [`apply()`](https://rdrr.io/r/base/apply.html)-family call.
