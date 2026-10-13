## Number of workers from `cores`: a number or a BiocParallel param object.
.nWorkers <- function(cores = 1) {
  if (is.numeric(cores) && length(cores) == 1 && !is.na(cores) && cores >= 1) {
    return(as.integer(cores))
  }
  if (is(cores, "BiocParallelParam")) {
    return(max(1L, as.integer(BiocParallel::bpnworkers(cores))))
  }
  stop("`cores` must be a positive number or a BiocParallel param object.",
    call. = FALSE
  )
}

## lapply() over `cores` forked workers (parallel::mclapply), forked once per
## call; sequential with one core or on Windows. BLAS runs one thread per
## worker: a multithreaded OpenBLAS would otherwise start its full thread pool
## in every forked worker.
#' @importFrom parallel mclapply
.parLapply <- function(X, FUN, ..., cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows" && length(X) > 1L) {
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      blasThreads <- RhpcBLASctl::blas_get_num_procs()
      RhpcBLASctl::blas_set_num_threads(1)
      on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
    }
    out <- parallel::mclapply(X, FUN, ...,
      mc.cores = min(cores, length(X)),
      mc.preschedule = FALSE
    )
    bad <- vapply(out, inherits, logical(1), what = "try-error")
    if (any(bad)) stop(attr(out[[which(bad)[1]]], "condition"))
    out
  } else {
    lapply(X, FUN, ...)
  }
}

## Random number streams, one per image: the streams BiocParallel gave each
## task with RNGseed = seed, so results do not depend on the number of cores.
.rngStreams <- function(n, seed) {
  oldSeed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  oldKind <- RNGkind()
  on.exit({
    RNGkind(oldKind[1], oldKind[2], oldKind[3])
    if (is.null(oldSeed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", oldSeed, envir = globalenv())
    }
  })
  RNGkind("default", "default", "default")
  set.seed(seed)
  RNGkind("L'Ecuyer-CMRG")
  s <- get(".Random.seed", envir = globalenv())
  streams <- vector("list", n)
  for (i in seq_len(n)) {
    streams[[i]] <- s
    s <- parallel::nextRNGSubStream(s)
  }
  streams
}

## Evaluate `expr` with the random number stream `stream`, then restore the
## caller's generator.
.withStream <- function(stream, expr) {
  oldSeed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  oldKind <- RNGkind()
  on.exit({
    RNGkind(oldKind[1], oldKind[2], oldKind[3])
    if (is.null(oldSeed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", oldSeed, envir = globalenv())
    }
  })
  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", stream, envir = globalenv())
  expr
}
