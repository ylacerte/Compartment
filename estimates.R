est_model <- function (t, var, params) {
  X <- var[1]
  Y <- var[2]
  Z <- var[3]
  N <- X + Y + Z
  with(as.list(params),{
    dX <- -a*X*Y/N
    dY <- a*X*Y/N - b*Y
    dZ <- b*Y
    res <- c(dX,dY,dZ)
    list(res)
  })
}

est_predict <- function(parameters, population, times) {
  # solve an ODE system via the ode function from deSolve package
  out <- ode(
    func  = est_model,
    y     = population,
    times = times,
    parms = parameters
  )
  return(data.frame(out))
}

estimate_parameters <- function() {
  
  # generate a random sample based on a model
  # estimate parameters
  
  times <- seq(0, 30, by=1/10)
  parameters <- c(a=0.5, b=0.1)
  population <- c(X=9, Y=1, Z=0)
  N <- sum(population)
  
  out <- est_predict(parameters=parameters, 
                     population=population,
                     times=times)
  
  # sample
  s <- seq(1,nrow(out),20) 
  sample <- out[s,c(1,3)]
  sample$Y <- jitter(sample$Y, 70)
  sample <- sample[order(sample$time),]
  
  
  SS <- function(init, data=sample) {
    ### least squares function
    # compute sum of squares errors between data and model
    
    times <- data$time
    parameters <- init
    population <- c(X=9, Y=1, Z=0)
    pred <- est_predict(parameters=parameters, 
                        population=population, 
                        times=times)
    
    df <- cbind(pred, data=data$Y)
    df$SS <- ( pred$Y - data$Y )^2
    return(sum(df$SS))
  }
  #init <- c(a=0.004, b=0.5) ; SS(init)
  
  
  f_min <- function(init) { SS(init) }
  # init <- c(a=0.004, b=0.5)  ;  f_min(init)
  
  
  loglik <- function (params, data=sample) {
    # MLE ... maximum likelihood of the data 
    # given the model and its parameters
    
    times <- data$time
    parameters <- params[1:2]
    sigma <- 1
    population <- c(X=9, Y=1, Z=0)
    pred <- est_predict(parameters, population, times)
    LL <- dnorm(x=data$Y,mean=pred$Y,sd=sigma,log=TRUE) 
    return(-sum(LL))
  }
  # params <- c(a=.004, b=.5)  ;  loglik(params=params)
  
  
  # Estimate the parameter values with the function optim() 
  # Input:
  ## init = Initial values for the parameters to be optimized over.
  ## f_min = A function to be minimized (or maximized).
  ###        This would be sum of squares or MLE
  params <- c(a=0.004, b=0.5)
  Optim1 <- optim(params, f_min)
  Optim2 <- optim(params, loglik)
  opt <- data.frame(rbind(least_squares=Optim1$par, 
                          mle=Optim2$par))
  
  return(list(sample=sample, opt=opt))
  
}

est_test <- function() {
  est <- estimate_parameters()
  est$opt
  
  times <- est$sample$time
  parameters <- c(a=est$opt$a[1], b=est$opt$b[1])
  population <- c(X=9, Y=1, Z=0)
  out <- est_predict(parameters, population, times)
  
  plot_ly(showlegend=TRUE) %>%
    add_trace(type='scatter', mode='lines',  data=out, x=~time,y=~X, name="X") %>%
    add_trace(type='scatter', mode='lines',  data=out, x=~time,y=~Y, name="Y") %>%
    add_trace(type='scatter', mode='lines',  data=out, x=~time,y=~Z, name="Z") %>%
    add_trace(type='scatter', mode='markers', data=est$sample, x=~time,y=~Y, name="Y'") %>%
    layout(yaxis=list(title="aircraft"))
}
est_test()

