###############################
run=3

## Libraries
###############################
has_ggplot2 <- require(ggplot2)
has_coda <- require(coda)
generate_original_results <- FALSE
library(nimble)
library(readr)
###############################

# Read in dataset
data <- read_csv("completedata.csv", col_types = cols(X1 = col_skip(), time = col_date(format = "%Y%m%d")))

# Period defines test and training sets
Period <- cbind(c("2015-10-13","2015-07-14","2015-04-14","2015-01-13","2014-10-14","2014-07-15","2014-04-15","2014-01-14","2013-10-15","2013-07-16"),
                c("2020-12-29","2020-09-29","2020-06-30","2020-03-31","2019-12-31","2019-10-01","2019-07-02","2019-04-02","2019-01-01","2018-10-02"))

# List of covariates
vars = c("air", "soilm")

# Model specifier
models = c("Full", "NoSpace", "Neither")

time.start = proc.time()[3]

# which period are we working with
thisPeriod <- 1

# which model are we working with
thisModel <- 1

# Data preparation
data$drought <- as.ordered(data$drought)
data$drought <- as.integer(data$drought) - 1

# Subset data temporally to the period we are working with on this run
data <- data[which(data$time >= Period[thisPeriod, 1] & data$time <=Period[thisPeriod, 2] ),]

# Subset data spatially
# data <- data[which(data$grid %in% unique(data$grid[data$lon == -106.75])),] ## One column, I=34
## data <- data[which(data$grid %in% unique(data$grid[data$lon %in% c(-106.75, -106.25)])),] ## One column, I=34
## ## data <- data[which(data$cal == TRUE),] ## CA, I=169
## data <- data[which(data$grid %in% unique(data$grid[data$lon < -117.5])),] ## Less West, I=343
data <- data[which(data$grid %in% unique(data$grid[data$lon < -114.5])),] ## More West, I=540

# Standardize covariates for faster runs
data[,vars] <- scale(data[,vars])

dim(data)
grids <- unique(data$grid)
weeks <- unique(data$time)
I <- length(grids)
p <- 13
J <- length(weeks) - p
C <- length(vars)
Training <- data[which(data$time %in% weeks[1:J]),]
Test <-  data[which(data$time %in% weeks[-(1:J)]),]
Y <- matrix(Training$drought,I, J)
A <- matrix(Training$air ,I, J)
S <- matrix(Training$soilm ,I, J)
X <- array(data = c(A,S), dim = c(I,J,C))
w <- rbind(data$x_sin[which(data$grid == grids[1])],data$x_cos[which(data$grid == grids[1])])
## X.s2 <- sin(2*asin(w[1,]))
## X.c2 <- cos(2*acos(w[2,]))
## w <- rbind(w, X.s2, X.c2)
## w <- unname(w)
W <- nrow(w)
G <- matrix(NA, I, C)
for(i in 1:I){
  G[i,] <- c(mean(A[i,]), mean(S[i,]))
}

# Initialize X.p
X.p = array(dim = c(I,p,C))
for(c in 1:(C)){
  for(x in 1:p){
    X.p[,x,c] = X[,J,c]
  }
}

# Adjacency matrix
coords = cbind(data$lon[1:I],data$lat[1:I])
D = as.matrix(dist(coords))
D = 1*((D<0.51) & (D>0)) # 1 for rook contiguity

# Get needed vectors for ICAR model
D = as.carAdjacency(D)

# Training time
TT = length(unique(data$time)) - 13

DEcode <- nimbleCode({
  #Drought Variable
  for(i in 1:I){
    for(t in 1:J){
      Y[i,t] ~ dinterval(Z[i,t], cut[])
    }
  }
  #Latent Gaussian Variable - includes ICAR U
  for(i in 1:I){
    mu[i,1] <- intercept + inprod(X[i,1,],beta[]) + U[i]
    Z[i,1] ~ dnorm(mu[i,1], tau = tau.z)
    for(t in 2:J){
      mu[i,t] <- intercept + inprod(X[i,t,],beta[]) + U[i] +
        rho.z[i]*(Z[i,(t-1)] - (intercept + inprod(X[i,t-1,],beta[]) + U[i]))
      Z[i,t] ~ dnorm(mu[i,t], tau = tau.z)
    }
  }

  # ICAR prior for the spatial random effects
  U[1:I] ~ dcar_normal(adj[], weights[], num[], tau.U, zero_mean=1)

  #Covariate Models
  for(i in 1:I){
    # U[i] <- 0
    for(c in 1:C) {
      X[i,1,c] ~ dnorm(G[c,i] + inprod(delta[c,,i],w[,1]), tau = tau.x[c])
    }
    for(t in 2:J){
      for(c in 1:C) {
        X[i,t,c] ~ dnorm(G[c,i] + inprod(delta[c,,i],w[,t]) + rho.x[i,c]*(X[i,(t-1),c] - (G[c,i] + inprod(delta[c,,i],w[,(t-1)]))), tau = tau.x[c])
      }
    }
  }
  #Random Effects
  for(i in 1:I){
    for(c in 1:C){G[c,i] ~ dnorm(0, tau = tau.G[c])}
  }
  for(i in 1:I){rho.z[i] ~ T(dnorm(rho.mean, tau = tau.rho),-1,1)}

  #Priors
  for(c in 1:C){
    tau.G[c] ~ dgamma(0.01, 0.01)
    tau.x[c] ~ dgamma(0.01, 0.01)
  }
  tau.rho ~ dgamma(0.01, 0.01)
  tau.z ~ dgamma(0.01, 0.01)
  tau.U ~ dgamma(0.01, 0.01)

  for (i in 1:I) {
    for (c in 1:C){
      rho.x[i,c] ~ dunif(-1,1)
    }
  }
  rho.mean ~ dunif(-1,1)
  for(i in 1:I) {
    for(c in 1:C){
      for (x in 1:W){
        delta[c,x,i] ~ dnorm(0, tau = 0.01)
      }
    }
  }
  intercept ~ dnorm(0, tau = 0.01)
  for(c in 1:C){beta[c] ~ dnorm(0, tau = 0.01)}
  for(i in 1:4){d[i] ~ dnorm(0, tau = 1)}

  cut[1:5] <- c(0, rep(exp(d[1]), 4)) + c(0, 0, rep(exp(d[2]), 3)) + c(0, 0, 0, rep(exp(d[3]), 2)) + c(0,0,0, 0,exp(d[4]))

  #Predict Covariates, Generate Z's, Predict Drought
  for(i in 1:I){
    for(c in 1:C){
      X.p[i,1,c] ~ dnorm(G[c,i] + inprod(delta[c,,i],w[,(J+1)]) + rho.x[i,c]*(X[i,J,c] - (G[c,i] + inprod(delta[c,,i],w[,J]))), tau = tau.x[c])
    }
    Z.p[i,1] ~ dnorm(intercept + inprod(X.p[i,1,],beta[]) + U[i] + rho.z[i]*(Z[i,J] - (intercept + inprod(X[i,J,],beta[]) + U[i])), tau = tau.z)
    Y.p[i,1] <- 0 + 1*(Z.p[i,1]>cut[1]) + 1*(Z.p[i,1]>cut[2]) + 1*(Z.p[i,1]>cut[3]) + 1*(Z.p[i,1]>cut[4]) + 1*(Z.p[i,1]>cut[5])

    for(t in 2:p){
      for(c in 1:C){
        X.p[i,t,c] ~ dnorm(G[c,i] + inprod(delta[c,,i],w[,(J+t)]) + rho.x[i,c]*(X.p[i,(t-1),c] - (G[c,i] + inprod(delta[c,,i],w[,(J+t-1)]))), tau = tau.x[c])
      }
      Z.p[i,t] ~ dnorm(intercept + inprod(X.p[i,t,],beta[]) + U[i] + rho.z[i]*(Z.p[i,(t-1)] - (intercept + inprod(X.p[i,(t-1),],beta[]) + U[i])), tau = tau.z)
      Y.p[i,t] <- 0 + 1*(Z.p[i,t]>cut[1]) + 1*(Z.p[i,t]>cut[2]) + 1*(Z.p[i,t]>cut[3]) + 1*(Z.p[i,t]>cut[4]) + 1*(Z.p[i,t]>cut[5])
    }
  }
})
## ends nimble code

nimbleOptions(MCMCusePredictiveDependenciesInCalculations = TRUE)

DEconstants <- list(I=I, J=J, p=p, adj=D$adj, num=D$num, weights=D$weights, w=w)

DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants, dimensions = list(cut = 5, X = c(I,J,C), beta = C, X.p = c(I,p,C),
                                                                  delta = c(C,W,I), w = c(W,(J+p))))

# This sets the values and *flags the nodes as data*
DEmodel$setData(list(Y = Y, X = X))

DEinits <- list(d = rep(log(3), 4),
                intercept = 0,
                beta = rep(0,C),
                tau.z = 1,
                tau.U = 1,
                Z = (Y - 0.5)*3,
                tau.x = rep(1,C),
                rho.x = array(data = 0.9, dim = c(I, C)),
                delta = array(data = 0, dim = c(C, W, I)),
                tau.G = rep(1,C),
                G = t(G),
                rho.z = rep(0.9, I),
                rho.mean = 0.9,
                tau.rho = 1,
                Z.p = matrix((Y[,J]-0.5)*3, I, p),
                X.p = X.p,
                U=rep(0,I)
)


DEmodel$setInits(DEinits)

conf <- configureMCMC(DEmodel) # What does this do?
#"Creates a default MCMC configuration for a given model

conf$addMonitors("mu", "rho.z", "U", "Y.p", "Z.p", "Z", "cut", "G", "delta", "rho.x", "X.p", "tau.z", "tau.U", "tau.x", "tau.G", "rho.mean", "tau.rho")

DEmcmc <- buildMCMC(conf) # Returns an uncompiled executable MCMC function
cDEmodel <- compileNimble(DEmodel, resetFunctions=TRUE, showCompilerOutput=TRUE)
cDEmcmc <- compileNimble(DEmcmc, project = DEmodel)
# cDEmcmc$run(10000)
# samples1 <- as.matrix(cDEmcmc$mvSamples)

## save.image("LaptopTest.Rdata")
set.seed(run*100+run*10+run)

samples2 <- runMCMC(cDEmcmc,
                    niter = 25000,
                    nburnin = 10000,
                    thin = 5,
                    nchains = 1,
                    samplesAsCodaMCMC = TRUE)

output=samples2


mu5 <- paste("mu[",seq(1:I),", 5]", sep="")
mu50 <- paste("mu[",seq(1:I),", 50]", sep="")
mu150 <- paste("mu[",seq(1:I),", 150]", sep="")

Z5 <- paste("Z[",seq(1:I),", 5]", sep="")
Z50 <- paste("Z[",seq(1:I),", 50]", sep="")
Z150 <- paste("Z[",seq(1:I),", 150]", sep="")

Yp1 <- paste("Y.p[",seq(1:I),", 1]", sep="")
Yp6 <- paste("Y.p[",seq(1:I),", 6]", sep="")
Yp13 <- paste("Y.p[",seq(1:I),", 13]", sep="")

## Save Ys here
TraceVars <- c(Z5, Z50, Z150, mu5, mu50, mu150, Yp1, Yp6, Yp13, "intercept", "beta[1]", "beta[2]", "d[1]", "d[2]", "d[3]", "d[4]","tau.z", "tau.U", "tau.x[1]", "tau.x[2]", "tau.G[1]", "tau.G[2]", "rho.mean", "tau.rho")
MCMCSave <- output[,which(names(output[1,]) %in% TraceVars)]

# Used to check if mixing
## save(MCMCSave, file = paste("RerunInformativePriors",run,".Rdata", sep=""))
loc <- as.character(seq(1:I))
tim <- as.character(seq(1:p))
var <- as.character(seq(1:C))
wav <- as.character(seq(1:W))
mat <- matrix(NA, I, p)
for(i in 1:I){
  for(t in 1:p){
    mat[i,t] <- paste("[", loc[i],", ", tim[t], "]", sep = "")
  }
}

Y.l <- matrix(NA, I, p)
Y.m <- matrix(NA, I, p)
Y.u <- matrix(NA, I, p)
Z.l <- matrix(NA, I, p)
Z.m <- matrix(NA, I, p)
Z.u <- matrix(NA, I, p)

Y.1 <- matrix(NA, I, p)
Y.2 <- matrix(NA, I, p)
Y.3 <- matrix(NA, I, p)
Y.4 <- matrix(NA, I, p)
Y.5 <- matrix(NA, I, p)
Y.6 <- matrix(NA, I, p)

for(i in 1:I){
  for(t in 1:p){
    Y.l[i,t] <- quantile(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))], probs = 0.025, type = 1)
    Y.m[i,t] <- quantile(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))], probs = 0.5, type = 1)
    Y.u[i,t] <- quantile(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))], probs = 0.975, type = 1)
    Z.l[i,t] <- quantile(output[,which(names(output[1,])==paste("Z.p", mat[i,t], sep = ""))], probs = 0.025)
    Z.m[i,t] <- quantile(output[,which(names(output[1,])==paste("Z.p", mat[i,t], sep = ""))], probs = 0.5)
    Z.u[i,t] <- quantile(output[,which(names(output[1,])==paste("Z.p", mat[i,t], sep = ""))], probs = 0.975)
    Y.1[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==0)
    Y.2[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==1)
    Y.3[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==2)
    Y.4[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==3)
    Y.5[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==4)
    Y.6[i,t] <- mean(output[,which(names(output[1,])==paste("Y.p", mat[i,t], sep = ""))]==5)
  }
}
# Used to compute RPS
# save(Y.l, Y.m, Y.u, Z.l, Z.m, Z.u, Y.1, Y.2, Y.3, Y.4, Y.5, Y.6,
#     file = paste("PredictInformativePriors",run,".Rdata", sep=""))

rho.l = rep(NA, I)
rho.m = rep(NA, I)
rho.u = rep(NA, I)
rho.x.l = array(data = NA, dim = c(I,C))
rho.x.m = array(data = NA, dim = c(I,C))
rho.x.u = array(data = NA, dim = c(I,C))
G.l = array(data = NA, dim = c(I,C))
G.m = array(data = NA, dim = c(I,C))
G.u = array(data = NA, dim = c(I,C))
delta.l = array(data = NA, dim = c(I,C,W))
delta.m = array(data = NA, dim = c(I,C,W))
delta.u = array(data = NA, dim = c(I,C,W))
U.m = rep(NA,I)
U.l = rep(NA,I)
U.u = rep(NA,I)


for(i in 1:I){
  rho.l[i] <- quantile(output[,which(names(output[1,])==paste("rho.z[", loc[i], "]", sep = ""))], probs = 0.025)
  rho.m[i] <- quantile(output[,which(names(output[1,])==paste("rho.z[", loc[i], "]", sep = ""))], probs = 0.5)
  rho.u[i] <- quantile(output[,which(names(output[1,])==paste("rho.z[", loc[i], "]", sep = ""))], probs = 0.975)
  U.l[i] <- quantile(output[,which(names(output[1,])==paste("U[", loc[i], "]", sep = ""))], probs = 0.025)
  U.m[i] <- quantile(output[,which(names(output[1,])==paste("U[", loc[i], "]", sep = ""))], probs = 0.5)
  U.u[i] <- quantile(output[,which(names(output[1,])==paste("U[", loc[i], "]", sep = ""))], probs = 0.975)
  for(c in 1:C) {
    G.l[i,c] <- quantile(output[,which(names(output[1,])==paste("G[", var[c], ", " ,  loc[i], "]", sep = ""))], probs = 0.025)
    G.m[i,c] <- quantile(output[,which(names(output[1,])==paste("G[", var[c], ", " , loc[i], "]", sep = ""))], probs = 0.5)
    G.u[i,c] <- quantile(output[,which(names(output[1,])==paste("G[", var[c], ", " , loc[i], "]", sep = ""))], probs = 0.975)
    rho.x.l[i,c] <- quantile(output[,which(names(output[1,])==paste("rho.x[", loc[i], ", " , var[c], "]", sep = ""))], probs = 0.025)
    rho.x.m[i,c] <- quantile(output[,which(names(output[1,])==paste("rho.x[", loc[i], ", " , var[c], "]", sep = ""))], probs = 0.5)
    rho.x.u[i,c] <- quantile(output[,which(names(output[1,])==paste("rho.x[", loc[i], ", " , var[c], "]", sep = ""))], probs = 0.975)
    for (w in 1:W) {
      delta.l[i,c,w] <- quantile(output[,which(names(output[1,])==paste("delta[", var[c], ", " , wav[w], ", " , loc[i], "]", sep = ""))], probs = 0.025)
      delta.m[i,c,w] <- quantile(output[,which(names(output[1,])==paste("delta[", var[c], ", " , wav[w], ", " , loc[i], "]", sep = ""))], probs = 0.5)
      delta.u[i,c,w] <- quantile(output[,which(names(output[1,])==paste("delta[", var[c], ", " , wav[w], ", " , loc[i], "]", sep = ""))], probs = 0.975)
    }
  }
}
# Used for plotting
save(rho.l, rho.m, rho.u, U.l, U.m, U.u, G.l, G.m, G.u, rho.x.l, rho.x.m, rho.x.u,
     delta.l, delta.m, delta.u, MCMCSave, Y.l, Y.m, Y.u, Z.l, Z.m, Z.u, Y.1, Y.2, Y.3, Y.4, Y.5, Y.6,
     file = paste("InformativePriors",run,".Rdata", sep=""))

this.runtime = proc.time()[3] - time.start
print(this.runtime)



