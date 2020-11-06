
load("~/work_data/basic_data_sgd/gnms.Rdata")
gnms <- cbind(rep(names(gnms), sapply(gnms, length)),
                            unlist(gnms))
norm.nam <- function(x, nf)
    {
        x <- gsub("SGD:", "", x)
        aux <- nf[match(x, nf[,2]),]
        r <- ifelse(x %in% nf[,1],
                    x, aux)
        r
    }

dot.product <- function(x,y=x) drop(crossprod(x,y))
eucl.norm <- function(x) sqrt(dot.product(x))

extract.sigs <- function(x, n)
{    
    sgs <- matrix(0, nrow=nrow(x[[2]]), ncol=n)
    for (i in 1:n) {
        aux <- as.matrix(x$u[,i]) %*% t(x$d[i] * x$v[,i])
        sgs[,i] <- aux[,which.max(apply(aux, 2, eucl.norm))]
    }
    sgs
}

## .. get percentage of variation for a mode 
get.pov <- function(x) x[[1]]^2/sum(x[[1]]^2)
