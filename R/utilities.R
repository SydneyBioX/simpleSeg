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

## A BiocParallel param for `cores` workers: forked workers (MulticoreParam),
## or socket workers (SnowParam) on Windows. A param object is returned
## unchanged.
.bpparam <- function(cores = 1L) {
  if (is(cores, "BiocParallelParam")) return(cores)
  cores <- .nWorkers(cores)
  if (cores <= 1L) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(cores)
  } else {
    BiocParallel::MulticoreParam(cores)
  }
}

## The task function and elements of the current .parLapply() call, found by
## forked workers in their copy of this environment.
.parallelEnv <- new.env(parent = emptyenv())
.runTask <- function(i) .parallelEnv$task(.parallelEnv$X[[i]])

## lapply(X, FUN) over the workers of `cores` (a number or a BiocParallel
## param).
## - Forked workers start once, after the task exists, so they share the
##   session's data without being sent it, and only task numbers and results
##   pass between processes.
## - BLAS runs one thread in each worker during the call. OpenBLAS otherwise
##   starts one thread per core in every worker. The setting is changed in this
##   session before the workers start, because changing it inside a worker
##   first restarts every thread.
## - Tasks are handed out in up to eight chunks per worker, so that quick and
##   slow tasks balance out.
#' @importFrom BiocParallel bpnworkers bpisup bplapply bptasks bptasks<- bpstart bpstop
.parLapply <- function(X, FUN, cores = 1L) {
  BPPARAM <- .bpparam(cores)
  workers <- BiocParallel::bpnworkers(BPPARAM)
  if (workers <= 1L || length(X) <= 1L) return(lapply(X, FUN))

  blasThreads <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1L)
  on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
  if (BiocParallel::bpisup(BPPARAM) || !is(BPPARAM, "MulticoreParam")) {
    # Workers that are already running, or are new R sessions, are sent the
    # task with each chunk.
    threadVars <- Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"), unset = NA)
    Sys.setenv(OPENBLAS_NUM_THREADS = 1L, OMP_NUM_THREADS = 1L)
    on.exit(.restoreEnv(threadVars), add = TRUE)
    return(BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM))
  }

  assign("task", FUN, envir = .parallelEnv)
  assign("X", X, envir = .parallelEnv)
  on.exit(rm(list = c("task", "X"), envir = .parallelEnv), add = TRUE)
  tasks <- BiocParallel::bptasks(BPPARAM)
  if (tasks == 0L) {
    BiocParallel::bptasks(BPPARAM) <- min(length(X), 8L * workers)
    on.exit(BiocParallel::bptasks(BPPARAM) <- tasks, add = TRUE, after = FALSE)
  }
  BiocParallel::bpstart(BPPARAM)
  on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE, after = FALSE)
  out <- BiocParallel::bplapply(seq_along(X), .runTask, BPPARAM = BPPARAM)
  names(out) <- names(X)
  out
}

## Put environment variables back as they were (NA: unset).
.restoreEnv <- function(values) {
  for (name in names(values)) {
    if (is.na(values[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(values[name]))
    }
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

## Gaussian blur with the kernel of EBImage::gblur() and wrap-around edges, as
## two one-dimensional passes in C++.
.gblur <- function(x, sigma) {
  x <- as.array(x)
  storage.mode(x) <- "double"
  .gaussianBlurCpp(x, sigma, TRUE)
}
