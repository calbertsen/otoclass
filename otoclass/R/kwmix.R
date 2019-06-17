
dkwmix <- Vectorize(function(x, p, df, a, b, log = FALSE){
    logres <- .Call("logdkwmix", x=x,p=qlogis(p),df=df,a=a,b=b)
    if(log)
        return(logres)
    return(exp(logres))
})

pkwmix <- Vectorize(function(x, p, df, a, b, log = FALSE){
    logres <- .Call("logpkwmix", x=x,p=qlogis(p),df=df,a=a,b=b)
    if(log)
        return(logres)
    return(exp(logres))
})

dkwmixE0 <- Vectorize(function(x, p, df, a, b, log = FALSE){
    logres <- .Call("logdkwmixE0", x=x,p=qlogis(p),df=df,a=a,b=b)
    if(log)
        return(logres)
    return(exp(logres))
})

pkwmixE0 <- Vectorize(function(x, p, df, a, b, log = FALSE){
    logres <- .Call("logpkwmixE0", x=x,p=qlogis(p),df=df,a=a,b=b)
    if(log)
        return(logres)
    return(exp(logres))
})
