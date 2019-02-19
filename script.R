N = 1000;
M_empty = 42600
M_pload = M_empty + 19900
Ra = 3e6
g = 9.81

# 1)
## N-Ã©chantillon de la vitesse V
V = runif(N, 226, 234);
V_moy = mean(V);
V_var = mean(V**2) - mean(V)**2;
hist(V,freq=FALSE,nclass=100)
curve(dunif(x, 226, 234),add=TRUE, col="blue")

## N-Ã©chantillon de la consomation SFC
SFC = rexp(N,3.45) + 17.23;
SFC_moy = mean(SFC);
SFC_var = mean(SFC**2) - mean(SFC)**2;
hist(SFC,freq=FALSE,nclass=100)
curve(dexp(x-17.23, 3.45),add=TRUE, col="blue")

## N-Ã©chantillon de la finesse f
f = (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7;
f_moy = mean(f);
f_var = mean(f**2) - mean(f)**2;
hist(f,freq=FALSE,nclass=50)
curve(dbeta((x - 18.7)/(19.05 - 18.7), 7, 2)/(19.05 - 18.7),add=TRUE, col="blue")

# 2)
## illustration de la LFGN
V1 = cumsum(V)/(1:N);
SFC1 = cumsum(SFC)/(1:N);
f1 = cumsum(f)/(1:N);
f2 = cumsum(c(f,(19.05 - 18.7)*rbeta(4000, 7, 2) + 18.7))/(1:5000);
plot(V1, type = "l", xlab = "n", ylab= "sum(V)/n",col = "blue", main = "Illustration de la LFGN pour V")
lines(0*1:1000 +230, type = "l")
plot(SFC1, type = "l", xlab = "n", ylab= "sum(SFC)/n", col = "blue", main = "Illustration de la LFGN pour SFC")
lines(0*1:1000 +1/3.45 + 17.23, type = "l")

layout(matrix(1:2, 1,2))
plot(f1, type = "l", xlab = "n=1000", ylab= "sum(f)/n", col = "blue")
lines(0*1:1000 +(19.05 - 18.7)*7/9 + 18.7, type = "l")
plot(f2, type = "l", xlab = "n=5000", ylab= "", col = "blue")
lines(0*1:5000 +(19.05 - 18.7)*7/9 + 18.7, type = "l")
layout(matrix(1))
title("Illustration de la LFGN pour F")

# 3)
X_V = c()
X_SFC = c()
X_f = c()
for (k in 1:1000) {
  Vk = runif(N, 226, 234)
  SFCk = rexp(N,3.45) + 17.23
  fk = (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7
  X_V = c(X_V,sum(Vk)/N)
  X_SFC = c(X_SFC,sum(SFCk)/N)
  X_f = c(X_f,sum(fk)/N)
};

hist(X_V,freq=FALSE, xlab = "v",nclass=50, main = "Illustration du TCL pour V (n=1000)")
curve(dnorm(x, 230, sqrt((234 - 226)**2/12/N)), add=TRUE, col="blue")
legend(229.78, 7, "N(230, 5.333e-3)", col = "blue", lty = 1)

hist(X_SFC,freq = FALSE, xlab = "sfc",nclass=50, main = "Illustration du TCL pour SFC (n=1000)")
curve(dnorm(x, 1/3.45 + 17.23, sqrt(1/3.45**2/N)) , add=TRUE, col="blue")
legend(17.491, 55, "N(17.52, 8.4e-5)", col = "blue", lty = 1)

hist(X_f,freq=FALSE, xlab = "f",nclass=50, main = "Illustration du TCL pour F (n=1000)")
curve(dnorm(x,(19.05 - 18.7)*7/9 + 18.7, sqrt(0.00212/N)), add=TRUE, col="blue")
legend(18.9664, 275, "N(18.972, 2.12e-6)", col = "blue", lty = 1)

# 4)
X_sigma = function(X,sigma) {
  return (X + sigma*rnorm(1000,0,1))
};

layout(matrix(1:4,2,2))
for (sigma in c(0.1, seq(1,3,1))) {
  M_fuel = (M_empty + M_pload)*(exp(X_sigma(SFC, sigma)*g*Ra*1e-6/(X_sigma(V, sigma)*(X_sigma(f, sigma)))) - 1)
  hist(M_fuel,freq=FALSE,nclass = 30, xlim = c(5000, 27000),main = 
         paste("Histogram of M_fuel \n for sigma =",sigma,", écart type =", round(sqrt(var(M_fuel)), 0), "\n moyenne =", round(mean(M_fuel),0)))
  
}
layout(matrix(1))

#1)
# Etude sans bruit
S_f = c()
for (n in 1:N) {
  fn = (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7
  V1n = runif(n, 226, 234)
  SFC1n = rexp(n,3.45) + 17.23
  
  V2n = runif(n, 226, 234)
  SFC2n = rexp(n,3.45) + 17.23
  
  y1n = (M_empty + M_pload)*(exp(SFC1n*g*Ra*1e-6/(V1n*fn)) - 1)
  y2n = (M_empty + M_pload)*(exp(SFC2n*g*Ra*1e-6/(V2n*fn)) - 1)
  
  S_f = c(S_f, (mean(y1n*y2n) - mean(y1n)*mean(y2n))/(mean(y1n**2) - mean(y1n)**2))
}

S_sfc = c()
for (n in 1:N) {
  SFCn = rexp(n,3.45) + 17.23
  
  V1n = runif(n, 226, 234)
  f1n = (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7
  
  V2n = runif(n, 226, 234)
  f2n = (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7
  
  y1n = (M_empty + M_pload)*(exp(SFCn*g*Ra*1e-6/(V1n*f1n)) - 1)
  y2n = (M_empty + M_pload)*(exp(SFCn*g*Ra*1e-6/(V2n*f2n)) - 1)
  
  S_sfc = c(S_sfc, (mean(y1n*y2n) - mean(y1n)*mean(y2n))/(mean(y1n**2) - mean(y1n)**2))
}

plot(S_f, type = 'l',xlab = 'n' ,main ="Indice de Sobol de M_fuel par rapport à F\n en fonction de n" ,col = "blue")
lines(0*1:1000 +mean(S_f[500:1000]), type = "l", lwd = 2, col = "red")
plot(S_sfc, type = 'l',xlab = 'n' ,main ="Indice de Sobol de M_fuel par rapport à SFC\n en fonction de n" ,col = "blue")
lines(0*1:1000 +mean(S_sfc[500:1000]), type = "l", lwd = 2, col ="red")

#2)
R_liste = c()
for (j in 1:10){
  f_index = c()
  for (val in (19.05 - 18.7)*rbeta(1000, 7, 2) + 18.7) {
    f_index = rbind(f_index, c(val, "f"))
  }
  
  SFC_index = c()
  for (val in rexp(N,3.45) + 17.23) {
    SFC_index = rbind(SFC_index, c(val, "sfc"))
  }
  
  f_SFC = rbind(f_index, SFC_index)
  f_SFC = f_SFC[order(f_SFC[, 1]), ]
  
  R = 0
  for (i in 1:1000) {
    if (f_SFC[i, 2]=="sfc"){
      R = R + i
    }
  }
  R_liste = c(R_liste, R)
}

R_liste_bruit = c()
for (sigma in seq(0.1,1,0.1)){
  f_index = c()
  for (val in X_sigma(f, sigma)) {
    f_index = rbind(f_index, c(val, "f"))
  }
  
  SFC_index = c()
  for (val in X_sigma(SFC, sigma)) {
    SFC_index = rbind(SFC_index, c(val, "sfc"))
  }
  
  f_SFC = rbind(f_index, SFC_index)
  f_SFC = f_SFC[order(f_SFC[, 1]), ]
  
  R_bruit = 0
  for (i in 1:1000) {
    if (f_SFC[i, 2]=="sfc"){
      R_bruit = R_bruit + i
    }
  }
  R_liste_bruit = c(R_liste_bruit, R_bruit)
}

#1)
#2)
sigma = 0.1
V_sigma = X_sigma(runif(N, 226, 234), sigma)
SFC_sigma = X_sigma(rexp(N,3.45) + 17.23, sigma)
f_sigma = X_sigma((19.05 - 18.7)*rbeta(N, 7, 2) + 18.7, sigma)

M_fuel = (M_empty + M_pload)*(exp(SFC_sigma*g*Ra*1e-6/(V_sigma*f_sigma)) - 1)
hist(M_fuel,freq=FALSE,nclass = 50,main = 
       paste("Histogram of M_fuel \n for sigma =",sigma,", ", mean(M_fuel)))

X1 = c()
for (i in 1:500) {
  X1 = rbind(X1,c(1, V_sigma[i], f_sigma[i], SFC_sigma[i]))
}
b1 = solve(aperm(X1)%*%X1)%*%(aperm(X1)%*%M_fuel[1:500])

X2 = c()
for (i in 501:1000) {
  X2 = rbind(X2,c(1, V_sigma[i], f_sigma[i], SFC_sigma[i]))
}
b2 = solve(aperm(X2)%*%X2)%*%(aperm(X2)%*%M_fuel[501:1000])

delta_b = abs((b2 - b1)/(b2 + b1))

y = M_fuel[1:500]
SST = y%*%y -500*mean(y)**2
SSR = aperm(b1)%*%aperm(X1)%*%y - 500*mean(y)**2

R_carrée = SSR/SST
a0 = b1[1]
a1 = b1[2]
a2 = b1[3]
a3 = b1[4]


#3)
S_fp = c()
for (n in 1:N) {
  fn = (19.05 - 18.7)*rbeta(n, 7, 2) + 18.7
  V1n = runif(n, 226, 234)
  SFC1n = rexp(n,3.45) + 17.23
  
  V2n = runif(n, 226, 234)
  SFC2n = rexp(n,3.45) + 17.23
  
  y1n = a0 + a1*V1n + a2*fn + a3*SFC1n
  y2n = a0 + a1*V1n + a2*fn + a3*SFC2n
  
  S_fp = c(S_fp, (mean(y1n*y2n) - mean(y1n)*mean(y2n))/(mean(y1n**2) - mean(y1n)**2))
}

S_sfcp = c()
for (n in 1:N) {
  SFCn = rexp(n,3.45) + 17.23
  
  V1n = runif(n, 226, 234)
  f1n = (19.05 - 18.7)*rbeta(n, 7, 2) + 18.7
  
  V2n = runif(n, 226, 234)
  f2n = (19.05 - 18.7)*rbeta(n, 7, 2) + 18.7
  
  y1n = a0 + a1*V1n + a2*f1n + a3*SFCn
  y2n = a0 + a1*V2n + a2*f2n + a3*SFCn
  
  S_sfcp = c(S_sfcp, (mean(y1n*y2n) - mean(y1n)*mean(y2n))/(mean(y1n**2) - mean(y1n)**2))
}

plot(S_fp, type = 'l',xlab = 'n' ,main ="Indice de Sobol de M_fuel par rapport à F\n en fonction de n" ,col = "blue")
lines(0*1:1000 +mean(S_fp[500:1000]), type = "l", lwd = 2, col = "red")
plot(S_sfcp, type = 'l',xlab = 'n' ,main ="Indice de Sobol de M_fuel par rapport à SFC\n en fonction de n" ,col = "blue")
lines(0*1:1000 +mean(S_sfcp[500:1000]), type = "l", lwd = 2, col ="red")
