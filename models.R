library(deSolve)
library(plotly)

# compartment models with solutions using the function ode() 
# from the deSolve package


example.model <- function (t, x, params) {
# 3 compartments 
  # with exogenous flows, 
  # no population constraints
  X <- x[1]
  Y <- x[2]
  Z <- x[3]
  with(as.list(params),{
    dX <- x_in - x_out - alpha*X
    dY <- y_in - y_out + alpha*X - beta*Y
    dZ <- beta*Y
    res <- c(dX,dY,dZ)
    list(res)
  })
}



basic.model <- function (t, x, params) {
# 3 compartments 
  # no exogenous flows, 
  # population constraints
  S <- x[1]
  I <- x[2]
  R <- x[3]
  N <- S + I + R
  with(as.list(params),{
    dS <- -beta*S*I/N
    dI <- beta*S*I/N - delta*I
    dR <- delta*I
    res <- c(dS,dI,dR)
    list(res)
  })
}


birthDeath.model <- function (t, x, params) {
# 3 compartments 
  # exogenous flows, mu = birth rate = death rate
  # population constraints
  # Notice damped oscillations toward an equilibrium.  
  S <- x[1]
  I <- x[2]
  R <- x[3]
  N <- S + I + R
  with(as.list(params),{
    dS <- mu*(N-S) - beta*S*I/N
    dI <- beta*S*I/N - delta*I - mu*I
    dR <- delta*I - mu*R
    res <- c(dS,dI,dR)
    list(res)
  })
}


seasonal.model <- function (t, x, params) {
  # 3 compartments 
  # no exogenous flows, 
  # population constraints
  # beta varies over time
  with(
    as.list(c(x,params)),{
      N <- S + I + R
      beta <- beta0 + beta1*cos(2*pi*t)
      dS <- mu*(N-S) - beta*S*I/N
      dI <- beta*S*I/N - delta*I - mu*I
      dR <- delta*I - mu*R
      res <- c(dS,dI,dR)
      list(res)
    }
  )
}




# solve each SIR model

predict.example <- function (parameters, population, times) {
  # solve the ODE system via the ode function from deSolve package
  # parameters: c(alpha=0.004, beta=0.5, x_in=1, x_out=1, y_in=1, y_out=1)
  # population: c(X=999, Y=1, Z=0)
  # times: seq(0, 100, by=1/120)
  out <- ode(
    func=example.model,
    y=population,
    times=times,
    parms=parameters
  )
  return(data.frame(out))
}


predict.basic <- function (parameters, population, times) {
# solve the ODE system via the ode function from deSolve package
# parameters: c(beta=0.004, delta=0.5)
# population: c(S=999, I=1, R=0)
# times: seq(0, 100, by=1/120)
  out <- ode(
    func=basic.model,
    y=population,
    times=times,
    parms=parameters
  )
  return(data.frame(out))
}

predict.birthDeath <- function (parameters, population, times) {
  # solve the ODE system via the ode function from deSolve package
  # parameters: c(mu=.02, beta=1000, delta=30)
  # population: c(S=999, I=1, R=0)
  # times: seq(0, 100, by=1/120)
  out <- ode(
    func=birthDeath.model,
    y=population,
    times=times,
    parms=parameters
  )
  return(data.frame(out))
}


predict.seasonal <- function (parameters, population, times) {
  # solve the ODE system via the ode function from deSolve package
  # parameters: c(mu=.02, beta=1000, delta=30)
  # population: c(S=999, I=1, R=0)
  # times: seq(0, 100, by=1/120)
  out <- ode(
    func=seasonal.model,
    y=population,
    times=times,
    parms=parameters
  )
  return(data.frame(out))
}


test <- function(what) {
  if ( what == "example" ) {
    times <- seq(0, 20, by=.1)
    parameters <- c(alpha=0.5, beta=0.1, x_in=.1, x_out=1, y_in=1,  y_out=1)
    population <- c(X=9, Y=5, Z=10)
    predict <- predict.example(parameters=parameters, 
                               population=population, 
                               times=times)
    plot_ly(data=predict, type='scatter', mode='lines', showlegend=TRUE) %>%
      add_trace(x=~time,y=~X, name="X") %>%
      add_trace(x=~time,y=~Y, name="Y") %>%
      add_trace(x=~time,y=~Z, name="Z") %>%
      layout(yaxis=list(title="population"))
    
  } else if ( what == "basic" ) {
    times <- seq(0, 40, by=.1)
    parameters <- c(beta=0.5, delta=0.1)
    population <- c(S=9, I=1, R=0)
    predict <- predict.basic(parameters=parameters, 
                             population=population, 
                             times=times)
    plot_ly(data=predict, type='scatter', mode='lines', showlegend=TRUE) %>%
      add_trace(x=~time,y=~S, name="S") %>%
      add_trace(x=~time,y=~I, name="I") %>%
      add_trace(x=~time,y=~R, name="R") %>%
      layout(yaxis=list(title="population"))
  } else if ( what == "birth death" ) {
    times <- seq(0, 15, by=.001)
    parameters <- c(mu=0.02, beta=1000, delta=30)
    population <- c(S=60, I=1, R=0)
    predict <- predict.birthDeath(parameters=parameters, 
                                  population=population, 
                                  times=times)
    plot_ly(data=predict, type='scatter', mode='lines', showlegend=TRUE) %>%
      add_trace(x=~time,y=~S, name="S") %>%
      add_trace(x=~time,y=~I, name="I") %>%
      add_trace(x=~time,y=~R, name="R") %>%
      layout(yaxis=list(title="population (log scale)", type='log'))
  } else if ( what == "seasonal" ) {
    times <- seq(0, 100, by=1/100)
    parameters <- c(mu=0.02, delta=30, beta0=1000, beta1=400)
    population <- c(S=60, I=1, R=0)
    predict <- predict.seasonal(parameters=parameters, 
                                population=population, 
                                times=times)
    plot_ly(data=predict, type='scatter', mode='lines', showlegend=TRUE) %>%
      add_trace(x=~time,y=~S, name="S") %>%
      add_trace(x=~time,y=~I, name="I") %>%
      add_trace(x=~time,y=~R, name="R") %>%
      layout(yaxis=list(title="population (log scale)", type='log'))
  } else {
    print("ERROR")
  }
}

#test("seasonal")
#test(what="birth death")
#test(what="basic")
#test(what="example")



# simulate stochastic differential equation model
library(yuima)
sde.model <- function(parameters, population, times) {
  
  # population ... c(S=8, I=2, R=0)
  # parameters ... list(beta = 0.5, delta = 0.1)
  # times ... sequence
  
  drift <- c("-beta *S*I / (S+I+R)",
             "( beta *S*I / (S+I+R) ) - delta *I",
             "delta *I")
  
  diffusion <- t(matrix(c("-sqrt(beta *S*I/(S+I+R))", "0",
                          "sqrt(beta *S*I/(S+I+R))", "-sqrt(delta *I)",
                          "0", "sqrt(delta *I)"),
                        2, 3))
  sol <- c("S", "I", "R")
  model <- setModel(drift = drift,
                    diffusion = diffusion,
                    solve.variable = sol,
                    state.variable = sol,
                    xinit=population)
  
  ysamp <- setSampling(Terminal=times[length(times)], n=2000)
  yuima <- setYuima(model=model, sampling=ysamp)
  
  GOsim <- simulate(yuima, xinit=population, 
                    true.parameter=parameters)
  
  z <- GOsim@data@zoo.data$`Series 1`
  p.df <- data.frame(t=index(z), S=z)
  
  z <- GOsim@data@zoo.data$`Series 2`
  p.df <- cbind(p.df, data.frame(I=z))
  
  z <- GOsim@data@zoo.data$`Series 3`
  p.df <- cbind(p.df, data.frame(R=z))
  traj <- p.df[complete.cases(p.df),]

  if ( nrow(traj) == nrow(p.df) ) {
    est <- lse(GOsim, start=list(beta=.5, delta=.1))
#    qmle <- qmle(GOsim, start=list(beta=.5, delta=.1))
    return(list(traj=traj, est=est))  
  } else { 
    return(NULL) 
  }
  
}

test <- function() {
  times <- seq(0, 40, by=.1)
  parameters <- c(beta=0.5, delta=0.1)
  population <- c(S=8, I=2, R=0)
  sde <- sde.model(parameters, population, times)
  if ( ! is.null(sde$est) ) print(sde$est)
  
  predict <- predict.basic(parameters=parameters, 
                           population=population, 
                           times=times)
  
  plot_ly( type='scatter', mode='lines', showlegend=TRUE) %>%
    add_trace(data=predict, x=~time,y=~S, name="S") %>%
    add_trace(data=predict, x=~time,y=~I, name="I") %>%
    add_trace(data=predict, x=~time,y=~R, name="R") %>%
    layout(yaxis=list(title="Population")) %>%
    add_trace(data=sde$traj, x=~t, y=~I, name="I (sim)") 
}
#test()



