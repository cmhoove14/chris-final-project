snail_sei <- function(s_0, e_0, i_0, M, params, time){

#initialize vectors to fill    
  s_fill <- vector("numeric", length = time+1)
    s_fill[1] <- s_0
  e_fill <- vector("numeric", length = time+1)
    e_fill[1] <- e_0
  i_fill <- vector("numeric", length = time+1)
    i_fill[1] <- i_0
  n_fill <- vector("numeric", length = time+1)
    n_fill[1] <- s_0+e_0+i_0
    
#Run discrete time SEI model with parameter set    
  with(as.list(params),{
    for(i in 1:time){
      s_fill[i+1] <- s_fill[i]*(1+f_N*(s_fill[i]=z*e_fill[i])*(1-n_fill[i]/K)-mu_N-beta_0*(1-exp(-M/n_fill[i])))
      e_fill[i+1] <- e_fill[i]+s_fill[i]*(1-exp(-M/n_fill[i]))-e_fill[i]*mu_N-e_fill[i]*sigma
      i_fill[i+1] <- i_fill[i]+e_fill[i]*sigma-(mu_N+mu_I)*i_fill[i]
    }
  })
  return(cbind(s_fill, e_fill, i_fill, n_fill))
}

snail_sei()