
Check_WhichElts <- function(WhichElements, q){

  if(q == 0){
    print(paste0("The argument q = 0, it should be larger than 0. Notably, WhichElements is a q times q matrix."))
    stop()
  }

  # Checks on WhichElements
  if(length(WhichElements) != 1){
    if(length(dim(WhichElements)) < 2){
      print(paste0("The argument WhichElements should be a square matrix. It should be of size q times q, with q = ", q))
      stop()
    }else if(length(dim(WhichElements)) > 2){
      print(paste0("The argument WhichElements should be an q times q matrix, with q = ", q, ". Currently, it is of size ", dim(WhichElements)))
      stop()
    }else if(dim(WhichElements)[1] != dim(WhichElements)[2]){
      print(paste0("The argument WhichElements should be a square matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(WhichElements)[1], " times ", dim(WhichElements)[2]))
      stop()
    }else if(dim(WhichElements)[1] != q){
      print(paste0("The argument WhichElements should be a matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(WhichElements)))
      stop()
    }
  }else if(q != 1){
    print(paste0("The argument WhichElements is a scalar, but it should be of size q times q, with q = ", q))
    stop()
  }
  #
  if(any(WhichElements != 0 & WhichElements != 1)){
    print(paste0("The argument WhichElements should consist of solely 1s and 0s."))
    stop()
  }

}
