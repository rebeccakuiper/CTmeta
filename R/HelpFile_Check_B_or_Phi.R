
Check_B_or_Phi <- function(B){

  if(is.null(dim(B))){
    if(!is.null(length(B))){
      stop(paste0("The argument Drift or Phi is not a matrix. It should be a matrix of size q x q."))
    }else{
      stop(paste0("The argument Drift or Phi is not found: The lagged effects matrix Drift or Phi is unknown, but should be part of the input."))
    }
  }else if(length(dim(B)) < 2){
    stop(paste0("The lagged effects matrix Drift or Phi should be a q x q matrix."))
  }else if(length(dim(B)) > 2){
    stop(paste0("The lagged effects matrix Drift or Phi should be a q x q matrix."))
  }else if(dim(B)[1] == 0 | dim(B)[2] == 0){
    stop(paste0("The lagged effects matrix Drift or Phi should be a square matrix of size q x q, with q > 0. Currently, it is of size ", dim(B)[1], " x ", dim(B)[2]))
  }else if(dim(B)[1] != dim(B)[2]){
    stop(paste0("The lagged effects matrix Drift or Phi should be a square matrix of size q x q. Currently, it is of size ", dim(B)[1], " x ", dim(B)[2]))
  }

}
