##' @export
groupProportion <- function(x,...)
    UseMethod("groupProportion")

##' @method groupProportion mlld
##' @export
groupProportion.mlld <- function(x, tabulate = FALSE,...){
    ssdr <- summary(x$sdr)
    tIn <- ssdr[rownames(ssdr)=="theta_logit",]
    reshp <- function(d) do.call("structure", c(list(.Data = unname(d),
                                                     dimnames = list(x$muNames[[3]],
                                                                     levels(x$proportionGroup))),
                                                     attributes(x$rp$theta_logit)))
    Estimate <- reshp( 1 / (1.0 + exp(-tIn[,1])) )
    ciLow <- reshp( 1 / (1.0 + exp(-(tIn[,1] - 2.0 * tIn[,2]))) )
    ciHigh <- reshp( 1 / (1.0 + exp(-(tIn[,1] + 2.0 * tIn[,2]))) )
    ## rownames(Estimate) <- rownames(ciLow) <- rownames(ciHigh) <- x$muNames[[3]]
    ## colnames(Estimate) <- colnames(ciLow) <- colnames(ciHigh) <- levels(x$proportionGroup)
    if(!tabulate){
        return(list(Estimate = Estimate,
                    ciLow = ciLow,
                    ciHigh = ciHigh))
    }
    toTab <- function(xx, indx){
        sp <- do.call("rbind",strsplit(colnames(xx),":"))
        tab <- xtabs(Val ~ ., data = cbind(data.frame(Val=xx[indx,]),sp))
        names(attr(tab,"dimnames")) <- NULL
        tabCount <- xtabs(Val ~ ., data = cbind(data.frame(Val=rep(1,nrow(sp))),sp))
        tab[tabCount == 0] <- NA
        tab
    }
    r <- lapply(1:nrow(Estimate), function(ii){
        list(Estimate = toTab(Estimate, ii),
             ciLow = toTab(ciLow, ii),
             ciHigh = toTab(ciHigh, ii))
    })
    names(r) <- x$muNames[[3]]
    r
}

##' @export
groupMean <- function(x,...)
    UseMethod("groupMean")

##' @method groupMean mlld
##' @export
groupMean.mlld <- function(x,...){
    ssdr <- summary(x$sdr)
    mIn <- ssdr[rownames(ssdr) == "muUse",]
    atr <- attributes(x$rp$muUse)
    nm <- x$muNames
    est <- do.call("structure",c(list(.Data = unname(mIn[,1]), dimnames = nm), atr))
    ciLow <- do.call("structure",c(list(.Data = unname(mIn[,1] - 2 * mIn[,2]), dimnames = nm), atr))
    ciHigh <- do.call("structure",c(list(.Data = unname(mIn[,1] + 2 * mIn[,2]), dimnames = nm), atr))
    list(Estimate = est,
         ciLow = ciLow,
         ciHigh = ciHigh)
}


##' @export
commonMean <- function(x,...)
    UseMethod("commonMean")

##' @method groupProportion mlld
##' @export
commonMean.mlld <- function(x,...){
    ssdr <- summary(x$sdr)
    mIn <- ssdr[rownames(ssdr) == "commonMu",]
    atr <- attributes(x$rp$commonMu)
    nm <- x$commonMuNames
    est <- do.call("structure",c(list(.Data = unname(mIn[,1]), dimnames = nm), atr))
    ciLow <- do.call("structure",c(list(.Data = unname(mIn[,1] - 2 * mIn[,2]), dimnames = nm), atr))
    ciHigh <- do.call("structure",c(list(.Data = unname(mIn[,1] + 2 * mIn[,2]), dimnames = nm), atr))
    list(Estimate = est,
         ciLow = ciLow,
         ciHigh = ciHigh)
}
