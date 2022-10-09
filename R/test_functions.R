#' Test Me Out
#' @description A test function made while working out package and github setup.
#'
#' @param x The input object
#'
#' @return A character string.
#' @export
#'
#' @examples
#' test_me_out(1)
#' test_me_out("Yes")
test_me_out <- function(x){
  if(x == "Yes"){
    return("Yayyyyyyy")
  } else {
    return("Nooooooooooo")
  }
}