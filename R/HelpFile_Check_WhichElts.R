
Check_WhichElts <- function(WhichElements, q){

  if(q == 0){
    stop("The argument q = 0. q should be larger than 0 as, WhichElements is a q x q matrix.")
  }

  # Checks on WhichElements
  if(length(WhichElements) != 1){
    if(length(dim(WhichElements)) < 2){
      stop(paste0("The argument WhichElements should be a square matrix. It should be of size q x q, with q = ", q))
    }else if(length(dim(WhichElements)) > 2){
      stop(paste0("The argument WhichElements should be an q x q matrix, with q = ", q, "."))
    }else if(dim(WhichElements)[1] != dim(WhichElements)[2]){
      stop(paste0("The argument WhichElements should be a square matrix of size q x q, with q = ", q, ". Currently, it is of size ", dim(WhichElements)[1], " x ", dim(WhichElements)[2]))
    }else if(dim(WhichElements)[1] != q){
      stop(paste0("The argument WhichElements should be a matrix of size q x q, with q = ", q, ". Currently, it is of size ", dim(WhichElements)[1], " x ", dim(WhichElements)[2]))
    }
  }else if(q != 1){
    stop(paste0("The argument WhichElements is a scalar, but it should be of size q x q, with q = ", q))
  }
  #
  if(any(WhichElements != 0 & WhichElements != 1)){
    stop("The argument WhichElements should consist of solely 1s and 0s.")
  }

}
