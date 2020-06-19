# f has to be a decreasing continuous function of lambda, and current_lambda between 0 and 1
#'@export
search_lambda <- function(current_lambda, f, objective, maxsteps = 1000, tolerance = 1e-2){
  if ((f(current_lambda) < objective)|| f(1) > objective){
    print("problem! there's no solution to the binary search")
  }
  attempt <- 1
  current_size <- (1 - current_lambda)/2
  fattempt <- f(attempt)
  istep <- 0
  while (!(fattempt >= objective && fattempt < objective+tolerance) && (istep < maxsteps)){
    istep <- istep + 1
    if (fattempt > objective){
      attempt <- attempt + current_size
      fattempt <- f(attempt)
      current_size <- current_size / 2
    } else {
      attempt <- attempt - current_size
      fattempt <- f(attempt)
      current_size <- current_size / 2
    }
  }
  return(list(x = attempt, f = fattempt))
}
