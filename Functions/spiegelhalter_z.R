##### SPIEGELHALTER_Z Function ####
# Returns z-score and calculates XXXXXXXXXXX

spiegelhalter_z <- function(y, prob){
  alpha = 0.05
  z_score = sum((y-prob)*(1-2*prob))/sqrt(sum(((1-2*prob)^2)*prob*(1-prob)))
  p_val = 2*(1-pnorm(abs(z_score)))
  return(c(z_score, p_val))
}