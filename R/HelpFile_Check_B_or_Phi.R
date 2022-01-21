
Check_B_or_Phi <- function(B){

  if(is.null(dim(B))){
    if(!is.null(length(B))){
      print(paste0("The argument Drift or Phi is not a matrix. It should be a matrix of size q x q."))
      stop()
    }else{
      print(paste0("The argument Drift or Phi is not found: The lagged effects matrix Drift or Phi is unknown, but should be part of the input."))
      stop()
    }
  }else if(length(dim(B)) < 2){
    print(paste0("The lagged effects matrix Drift or Phi should be a q x q matrix."))
    stop()
  }else if(length(dim(B)) > 2){
    print(paste0("The lagged effects matrix Drift or Phi should be a q x q matrix."))
    stop()
  }else if(dim(B)[1] == 0 | dim(B)[2] == 0){
    print(paste0("The lagged effects matrix Drift or Phi should be a square matrix of size q x q, with q > 0. Currently, it is of size ", dim(B)[1], " x ", dim(B)[2]))
    stop()
  }else if(dim(B)[1] != dim(B)[2]){
    print(paste0("The lagged effects matrix Drift or Phi should be a square matrix of size q x q. Currently, it is of size ", dim(B)[1], " x ", dim(B)[2]))
    stop()
  }

}
