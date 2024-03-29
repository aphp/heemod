## set test sets
tempdf <- expand.grid(arg1 = c("A", "B", "C"), arg2 = 1:4, arg3 = 1:5)
tempdf$value <- 1:60


just_numeric <- 1:10
names(just_numeric) <- 10 * 0:9

oned_df <- data.frame(age = 10 * 0:9, value = 1:10)
oned_df_othername <- oned_df
names(oned_df_othername)[2] <- "decade"

multiple_outputs_df <- oned_df_othername
multiple_outputs_df$rev_decade <- 11 - multiple_outputs_df$decade

set.seed(4)
tempdf_rand <- tempdf[sample(1:60), ]
oned_df_rand <- oned_df[sample(1:10), ]

test_that(
  "Warnings on non-matched values", {
    
    expect_identical(
      suppressWarnings(
        look_up(
          tempdf, arg1 = c("A"),
          arg2 = c(3), arg3 = c(0.5),
          bin  = c("arg2", "arg3")
        )
      ),
      NA_integer_
    )
    
    expect_warning(
      look_up(
        tempdf, arg1 = c("A"),
        arg2 = c(3), arg3 = c(0.5),
        bin  = c("arg2", "arg3")
      ),
      "Some values were not found, returning missing data"
    )
    
    tempdf <- expand.grid(arg1 = c("A", "B", "C"), arg2 = 1:4, arg3 = 1:5)
    tempdf$value <- 1:60
    
    expect_warning(
      look_up(
        data = tempdf,
        value = "value",
        arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(-1, 1, 1, 2, 3)
      ),
      "arguments to look_up:
arg1 : A, C, A
arg2 : 1, 3.2, 5
arg3 : -1, 1, 3"
    )
    
  }
)

test_that(
  "Errors on improper input.", {
    
    expect_error(
      look_up(
        tempdf, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 1, 1, 1),
        numeric_cols  = c("arg2", "arg3")
      )
    )
    
    expect_error(
      look_up(
        tempdf, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "xxx")
      ),
      "Names in 'bin' not found in source data: xxx.",
      fixed = TRUE
    )
    expect_error(
      look_up(
        tempdf, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg1")
      ),
      "Some variables in 'bin' are not numeric in the selection data: arg1.",
      fixed = TRUE
    )
    tempdf2 <- tempdf
    tempdf2$arg2 <- as.character(tempdf2$arg2)
    expect_error(
      look_up(
        tempdf2, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      "Some variables in 'bin' are not numeric in the source data: arg2.",
      fixed = TRUE
    )
    tempdf3 <- rbind(
      tempdf, data.frame(
        arg1 = "A", arg2 = 1, arg3 = 1, value = pi))
    expect_error(
      look_up(
        tempdf3, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      "Some rows in 'data' are duplicates: 61",
      fixed = TRUE
    )
    expect_error(
      look_up(
        tempdf, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      "All variables passed to 'look_up()' must be named.",
      fixed = TRUE
    )
    tempdf4 <- tempdf
    tempdf4$arg2[60] <- Inf
    expect_error(
      look_up(
        tempdf4, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      "infinite values in look_up table element: arg2",
      fixed = TRUE
    )
    
    expect_error(
      look_up(
        just_numeric, age = c(0, 30, 40)
      ),
      "'data' must be a data.frame",
      fixed = TRUE
    )
    
    expect_error(
      look_up(oned_array, age = c(0, 30, 40)),
      "not found"
    )
    
    ## unspecified output_col with multiple possibilities
    expect_error(
      look_up(multiple_outputs_df, age = c(0, 30, 40)),
      "Names passed to 'look_up()' not found in 'data': value.",
      fixed = TRUE
    )
    expect_error(
      look_up(
        multiple_outputs_df, age = c(0, 30, 40), 
        value = "random"
      ),
      "Names passed to 'look_up()' not found in 'data': random.",
      fixed = TRUE
    )
  }
)

test_that(
  "Correct values with correct input.", {
    
    expect_equal(
      look_up(
        tempdf, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      c(1, 2, 9, 20, 34)
    )
    
    expect_equal(
      look_up(
        oned_df, age = c(0, 10, 20)
      ),
      c(1,2,3)
    )
    
    expect_equal(
      look_up(
        oned_df_othername, age = c(0.5, 23, 42.7),
        bin = "age", value = "decade"
      ),
      c(1, 3, 5)
    )
    
    expect_equal(
      look_up(
        oned_df_othername, age = c(0, 10, 20), value = "decade"
      ),
      c(1, 2, 3)
    )
    
    expect_equal(
      look_up(
        oned_df, age = c(0.5, 23, 42.7),
        bin = "age"
      ),
      c(1, 3, 5)
    )
    
    expect_equal(
      look_up(
        multiple_outputs_df,
        age = c(0, 10, 20),
        value = "rev_decade"
      ),
      c(10, 9, 8)
    )
    
    expect_equal(
      look_up(
        multiple_outputs_df,
        age = c(0.5, 23, 42.7),
        bin = "age", value = "decade"
      ),
      c(1, 3, 5)
    )
    
    
    expect_equal(
      look_up(
        tempdf_rand, arg1 = c("A", "B", "C", "B", "A"),
        arg2 = c(1, 1, 3.2, 3.0, 5), 
        arg3 = c(1, 1, 1, 2, 3),
        bin  = c("arg2", "arg3")
      ),
      c(1, 2, 9, 20, 34)
    )
    
    expect_equal(
      look_up(
        oned_df_rand, age = c(0, 10, 20)
      ),
      c(1, 2, 3)
    )
    
    expect_equal(
      look_up(
        oned_df_rand, age = c(0.5, 23, 42.7),
        bin = "age"
      ),
      c(1, 3, 5)
    )
  }
)
