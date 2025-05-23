`%||%` <- function(x, y) {
    if (length(x) <= 0L) {
        y
    } else {
        x
    }
}

`%|||%` <- function(x, y) {
    if (is.null(x)) {
        y
    } else {
        x
    }
}

`%||NA%` <- function(x, y) {
    if (anyNA(x)) {
        y
    } else {
        x
    }
}

`%||nf%` <- function(x, y) {
    if (length(x) <= 0L || anyNA(x)) {
        y
    } else {
        x
    }
}

if_any <- function(condition, x, y) {
    if (any(condition)) {
        x
    } else {
        y
    }
}
