##' @importFrom stats nlminb optimHess
##' @importFrom utils head
##' @export
biascorrect <- function(assigned, cv, method = c("mle","mi")){
    if(!is.factor(assigned))
        stop("assigned must be a factor of assigned groups.")
    if(!inherits(cv,"crossval_confusion"))
        stop("cv must be the output of the crossval function")
    method <- match.arg(method)
    M <- t(cv$confusion_matrix)
    if(method == "mi")
        return(solve(M, table(assigned) / length(assigned)))

    X <- as.numeric(table(assigned))

    if(length(X) != nrow(M) | length(X) != ncol(M))
        stop("Assignments and cross validation does not match")

    f <- function(logitP){
        .Call("biascorrect_logistic",
              x=logitP,
              PACKAGE = "otoclass")
    }
    grf <- function(logitP){
        .Call("biascorrect_logistic_gradient",
              x=logitP,
              PACKAGE = "otoclass")
    }

    objective <- function(logitP){
        .Call("biascorrect_objective",
              x=logitP,
              X = X,
              M = M,
              PACKAGE = "otoclass")
    }

    gradient <- function(logitP){
      .Call("biascorrect_gradient",
              x=logitP,
              X = X,
              M = M,
              PACKAGE = "otoclass")
    }


    opt <- stats::nlminb(utils::head(X*0.0,-1), objective, gradient)
    he <- stats::optimHess(opt$par,objective, gradient)
    gr <- grf(opt$par)
    covar <- gr %*% solve(he) %*% t(gr)
    tab <- cbind(f(opt$par), sqrt(diag(covar)))
    rownames(tab) <- levels(assigned)
    rownames(covar) <- colnames(covar) <- levels(assigned)
    colnames(tab) <- c("Estimate","Std.Err")
    list(estimate_table = tab,
         full_covariance = covar,
         mle = opt$par,
         hessian = he)
}

