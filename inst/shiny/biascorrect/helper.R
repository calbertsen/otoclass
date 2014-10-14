summary_wna <- function (object, ..., digits = max(3L, getOption("digits") - 
    3L)) 
{
        nas <- is.na(object)
        object <- object[!nas]
        qq <- stats::quantile(object)
        qq <- signif(c(qq[1L:3L], mean(object), qq[4L:5L]), digits)
        names(qq) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", 
            "Max.")
        if (any(nas)) 
            c(qq, `NA's` = sum(nas))
        else c(qq, `NA's` = 0 )
}

