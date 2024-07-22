# cheapr_cast <- function(x, y){
#   if (identical(class(x), class(y))){
#     x
#   } else {
#     x[0] <- y[0]
#     `mostattributes<-`(x, attributes(y))
#   }
# }
