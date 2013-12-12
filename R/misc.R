mclapplyPb <- function(X, FUN, mc.cores = 2, ...) {
   env <- environment()
   pb_Total <- length(X)/mc.cores
   counter <- 0
   pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
   wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
   }
   res <- mclapply(X, wrapper, mc.cores = mc.cores, ...)
   #res <- lapply(X, wrapper,  ...)
   close(pb)
   res
}

# round2.max <- function(x) {
#   max.x  <- ceiling(x)
#   diff.x <- max.x - x
#     if (diff.x == 0.5) {
#        ceiling(x)
#     } else{
#        round(x)
#     }
# }
