
Check_Phi_or_B <- function(Phi){

  if(is.null(dim(Phi))){
    if(!is.null(length(Phi))){
      print(paste0("The argument Phi or Drift is not a matrix. It should be a matrix of size q times q."))
      stop()
    }else{
      print(paste0("The argument Phi or Drift is not found: The lagged effects matrix Phi or Drift is unknown, but should be part of the input."))
      stop()
    }
  }else if(length(dim(Phi)) < 2){
    print(paste0("The lagged effects matrix Phi or Drift should be an q times q matrix."))
    stop()
  }else if(length(dim(Phi)) > 2){
    print(paste0("The lagged effects matrix Phi or Drift should be an q times q matrix. Currently, it is of size ", dim(Phi)))
    stop()
  }else if(dim(Phi)[1] == 0 | dim(Phi)[2] == 0){
    print(paste0("The lagged effects matrix Phi or Drift should be a square matrix of size q times q, with q > 0. Currently, it is of size ", dim(Phi)[1], " times ", dim(Phi)[2]))
    stop()
  }else if(dim(Phi)[1] != dim(Phi)[2]){
    print(paste0("The lagged effects matrix Phi or Drift should be a square matrix of size q times q. Currently, it is of size ", dim(Phi)[1], " times ", dim(Phi)[2]))
    stop()
  }

}
