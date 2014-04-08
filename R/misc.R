## mclapply with a progress bar

mclapplyPb <- function(X, FUN, mc.cores = getOption("mc.cores", 2L), ...) {
   env <- environment()
   pb_Total <- length(X) / mc.cores
   counter <- 0
   pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)
   wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1 ,envir = env)
      setTxtProgressBar(get("pb", envir = env), curVal + 1)
      FUN(...)
   }
   res <- mclapply(X, wrapper, mc.cores = mc.cores, ...)
   close(pb)
   res
}
