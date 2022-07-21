# This function requires the following packages:
#   - crayon
#   - glue
#   - stringr

generate_octopus <- function() {

    # Initialize art, split by newlines ---------------------------------------
    octopus_ascii <- "                        ___\n                     .-'   `'.\n                    /         \\\n                    |         ;\n                    |         |           ___.--,\n           _.._     |0) ~ (0) |    _.---'`__.-( (_.\n    __.--'`_.. '.__.\\    '--. \\_.-' ,.--'`     `\"\"`\n   ( ,.--'`   ',__ /./;   ;, '.__.'`    __\n   _`) )  .---.__.' / |   |\\   \\__..--\"\"  \"\"\"--.,_\n  `---' .'.''-._.-'`_./  /\\ '.  \\ _.-~~~````~~~-._`-.__.'\n        | |  .' _.-' |  |  \\  \\  '.               `~---`\n         \\ \\/ .'     \\  \\   '. '-._)\n          \\/ /        \\  \\    `=.__`~-.\n          / /\\         `) )    / / `\"\".`\\\n    , _.-'.'\\ \\        / /    ( (     / /\n     `--~`   ) )    .-'.'      '.'.  | (\n            (/`    ( (`          ) )  '-;\n             `      '-;         (-'"
    octopus_ascii <- stringr::str_split(octopus_ascii, "\\n")[[1]]


    # Adjustments to ASCII art ------------------------------------------------

    # wink
    if (stats::runif(1) < 1/5) {
        octopus_ascii <-
            stringr::str_replace(octopus_ascii, "(?<=\\()0", "-") |>
            stringr::str_replace("~", " ")

        # blink
        if (stats::runif(1) < 1/3) {
            octopus_ascii <- stringr::str_replace(octopus_ascii, "0", "-")

            # sleep
            if (stats::runif(1) < 1/2) {
                octopus_ascii[2] <- stringr::str_replace(octopus_ascii[2], "^                  ", "{crayon::white('              Z  Z')}")
                octopus_ascii[3] <- stringr::str_replace(octopus_ascii[3], "^                  ", "{crayon::white('             Z   Z')}")
                octopus_ascii[4] <- stringr::str_replace(octopus_ascii[4], "^                  ", "{crayon::white('            Z Z   ')}")
                octopus_ascii[5] <- stringr::str_replace(octopus_ascii[5], "^                  ", "{crayon::white('             Z  Z ')}")
            }
        }

    } else {

        # angry
        if (stats::runif(1) < 1/5) {
            octopus_ascii <-
                stringr::str_replace_all(octopus_ascii, "0", '{crayon::bold(crayon::red("0"))}') |>
                stringr::str_replace("~", "v")
        }

    }

    # glue::glue() for various expressions
    octopus_ascii <- vapply(octopus_ascii, glue::glue, FUN.VALUE = character(1))


    # Assign colors -----------------------------------------------------------

    color_options <- c(
        "solid",
        "stripe",
        "rainbow"
    )

    color_funs <- list(
        crayon::cyan,
        crayon::blue,
        crayon::magenta,
        crayon::green
    )

    color_funs <-switch(
        sample(color_options, 1),
        "solid" = color_funs[[sample(1:length(color_funs), 1)]],
        "stripe" = color_funs[sample(1:length(color_funs), 2)],
        "rainbow" = color_funs
    )

    color_ascii <- function(string) {

        if (length(color_funs) == 1) {
            color_fun <- color_funs
        } else {
            color_fun <- color_funs[[sample(1:length(color_funs), 1)]]
        }

        color_fun(string)
    }

    # Collapse back to vector of length 1, w/ newlines
    octopus_ascii <- lapply(octopus_ascii, color_ascii) |>
        paste(collapse = "\n")



    # cat, returning string invisibly -----------------------------------------
    cat(octopus_ascii)
    invisible(octopus_ascii)
}



.First <- function(){
    generate_octopus()
    library(devtools, quietly = TRUE)
    library(testthat, quietly = TRUE, warn.conflicts = FALSE)
    library(checkmate, quietly = TRUE)
    library(covr, quietly = TRUE)
    library(beepr, quietly = TRUE)
    library(fortunes, quietly = TRUE)
    print(fortunes::fortune(12))
    print(fortunes::fortune())
    beep(6)
    # rm(generate_octopus)
}

.Last <- function(){
    cat("\nGoodbye at ", date(), "\n")
}
