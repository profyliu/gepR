#' Load gepR.dll
#' @param dll_file the path and name of the dll file (or .dylib file) to be loaded.
gep_load_dll <- function(dll_file = "gepR.dll"){
  dyn.load(dll_file)
}

#' Unload gepR.dll
#' @param dll_file the path and name of the dll file (or .dylib file) to be unloaded.
gep_unload_dll <- function(dll_file = "gepR.dll"){
  dyn.unload(dll_file)
}


#' Train a GEP model
#'
#' Fit a composite linear-nonlinear regression model using Gene Expression Programming. The gepR.dll must be present in the work directory and loaded using the \link{gep_load_dll} function.
#' @param y numeric vector containing the response variable
#' @param x numeric dataframe containing the independent variables
#' @param px1 one-point crossover rate, between 0 and 1
#' @param px2 two-point crossover rate, between 0 and 1
#' @param pm mutation rate, between 0 and 1
#' @param maxiter maximum number of iteration (generations) per round
#' @param headlen head length of a gene
#' @param popsize population size (number of individuals per generation), default 100. This is where the OpenMP "parallel for" takes effect. So it is preferred to be an integer multiple of nthreads.
#' @param eliterate elite rate, the proportion in the population to survive to the next generation, between 0 and 1
#' @param goal the targeted R-square value, between 0 and 1, default 0.95. The program stops when R-square achieves this value.
#' @param rseed random seed to be used in the C program, integer value
#' @param nthreads number of parallel threads to be used in computation, ideally set to equal to the number of cores available
#' @param verbose verbose leve, 0 to 2, 0 terse, 2 verbose
#' @param fit_method 0: regression, 1: classification, default 0. Only 0 is implemented in this version.
#' @param sol_file character string of the path and name of the file where the GEP model is to be saved.
#' @return sol_file
#'
#' @examples
#' gep_load_dll()
#' gepmod <- gep_train(elecdemand[,1], elecdemand[,2:4], nthreads=20)
#' predval <- gep_score(elecdemand[,2:4], gepmod)
#' gep_unload_dll()
#'
gep_train <- function(y,x,px1=0.4,px2=0.1,pm=0.3,maxiter=1000,headlen=5,popsize=100,eliterate=0.1,goal=0.95,
                      rseed=8888,nthreads=4,verbose=1,fit_method=0,maxpass=3,sol_file='gep_sol.dat'){
  
  # Check if DLL is loaded
  if(!("gepR" %in% names(getLoadedDLLs()) || "gepR.dylib" %in% names(getLoadedDLLs())))
    stop("The gepR.dll (or gepR.dylib) is not loaded. Use the gep_load_dll() function to load it.")
  
  nrows = length(y)
  dim_x = dim(x)
  if(nrows != dim_x[1])
    stop("Lengths of x and y do not match")
  else
    nvars = dim_x[2]
  
  # A lot of input validity checks are needed here to prevent the DLL from crashing.
  if (any(is.na(y)))
    stop("y contains missing value")
  
  x <- as.double(unlist(x))
  if (any(is.na(x)))
    stop("x contains missing value")
  
  if (!is.numeric(x))
    stop("x must be a numeric vector")
  if (px1 < 0 | px1 > 1)
    stop("px1 value is invalid.")
  if (px2 < 0 | px2 > 1)
    stop("px2 value is invalid.")
  if (pm < 0 | pm > 1)
    stop("pm value is invalid.")
  if (maxiter <= 1)
    stop("maxiter value is invalid")
  if (headlen <= 0 | headlen > 100)
    stop("headlen value is invalid")
  if (popsize <= 0)
    stop("popsize value is invalid")
  if (eliterate < 0 | eliterate > 1)
    stop("eliterate value is invalid.")
  if (goal < 0 | goal > 1)
    stop("goal value is invalid.")
  if (rseed <=0)
    stop("rseed value is invalid")
  if (nthreads < 1)
    stop("nthreads value is invalid")
  if (verbose < 0 | verbose > 2)
    stop("verbose value is invalid")
  if (fit_method != 0 & fit_method != 1)
    stop("fit_method value is invalid")
  
  gepmod <- .C(getNativeSymbolInfo("gep_run"),
               nrows=as.integer(nrows),
               nvars=as.integer(nvars),
               y=as.double(unlist(y)),
               x=as.double(unlist(x)),
               px1=as.double(px1),
               px2=as.double(px2),
               pm=as.double(pm),
               maxiter=as.integer(maxiter),
               headlen=as.integer(headlen),
               popsize=as.integer(popsize),
               eliterate=as.double(eliterate),
               goal=as.double(goal),
               rseed=as.integer(rseed),
               nthreads=as.integer(nthreads),
               verbose=as.integer(verbose),
               fit_method=as.integer(fit_method),
               maxpass=as.integer(maxpass),
               sol_file
  )
  cat(paste("GEP model saved to file", sol_file, "\n"))
  return(sol_file)
}

#' Predict using a GEP model
#' @param x data frame containing the input data. Must have the same number of columns as in the training data
#' @param sol_file character string of the file name where the trained model is stored. Default = "gep_sol.dat".
#' @return numeric vector containing the predicted values
#'
gep_score <- function(x, sol_file = 'gep_sol.dat'){
  
  # Check if DLL is loaded
  if(!("gepR" %in% names(getLoadedDLLs()) || "gepR.dylib" %in% names(getLoadedDLLs())))
    stop("The gepR.dll (or gepR.dylib) is not loaded. Use the gep_load_dll() function to load it.")
  
  y <- rep(0, nrow(x))
  if (any(is.na(x)))
    stop("x contains missing value")
  score = .C(getNativeSymbolInfo("gep_score"),
             sol_file,
             x=as.double(unlist(x)),
             y=as.double(y),
             nrows=as.integer(nrow(x)),
             nvars=as.integer(ncol(x))
  )
  return(score$y)
}
