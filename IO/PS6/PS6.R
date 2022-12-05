fpa <- read.table("~/Downloads/fpa.dat", quote="\"", comment.char="")

fpa$M1 <- apply(cbind(fpa$V2, fpa$V3, fpa$V4),1,max,na.rm=TRUE)
fpa$M2 <- apply(cbind(fpa$V1, fpa$V3, fpa$V4),1,max,na.rm=TRUE)
fpa$M3 <- apply(cbind(fpa$V1, fpa$V2, fpa$V4),1,max,na.rm=TRUE)
fpa$M4 <- apply(cbind(fpa$V1, fpa$V2, fpa$V3),1,max,na.rm=TRUE)

fpa$max <- apply(cbind(fpa$V1, fpa$V2, fpa$V3, fpa$V4),1,max,na.rm=TRUE)

h <- 20
grid <- seq(0, 160, 0.1)
kd_est <- data.frame(grid)
kd_est$G <- 0
kd_est$g <- 0
for (ind in 1:length(grid)){
  b <- grid[ind]
  kd_est$G1[ind] <- (1/(h))*mean((dnorm((b - fpa$V1)/h)*(fpa$M1 < b) ))
  kd_est$G2[ind] <- (1/(h))*mean((dnorm((b - fpa$V2)/h)*(fpa$M2< b) ))
  kd_est$G3[ind] <- (1/(h))*mean((dnorm((b - fpa$V3)/h)*(fpa$M3< b) ))
  kd_est$G4[ind] <- (1/(h))*mean((dnorm((b - fpa$V4)/h)*(fpa$M4< b) ))
  kd_est$g1[ind] <- (1/(h^2))*mean(dnorm((b - fpa$V1)/h)*dnorm((b - fpa$M1)/h))
  kd_est$g2[ind] <- (1/(h^2))*mean(dnorm((b - fpa$V2)/h)*dnorm((b - fpa$M2)/h))
  kd_est$g3[ind] <- (1/(h^2))*mean(dnorm((b - fpa$V3)/h)*dnorm((b - fpa$M3)/h))
  kd_est$g4[ind] <- (1/(h^2))*mean(dnorm((b - fpa$V4)/h)*dnorm((b - fpa$M4)/h))
}

plot(grid,kd_est$G1/kd_est$g1)

kd_est$markdown1 <- kd_est$G1/kd_est$g1
kd_est$markdown2 <- kd_est$G2/kd_est$g2
kd_est$markdown3 <- kd_est$G3/kd_est$g3
kd_est$markdown4 <- kd_est$G4/kd_est$g4

markdown1 <- approxfun(grid, y=kd_est$markdown1)
markdown2 <- approxfun(grid, y=kd_est$markdown2)
markdown3 <- approxfun(grid, y=kd_est$markdown3)
markdown4 <- approxfun(grid, y=kd_est$markdown4)

fpa$u1 <- fpa$V1 + markdown1(fpa$V1)
fpa$u2 <- fpa$V2 + markdown2$(fpa$V2)
fpa$u3 <- fpa$V3 + markdown3(fpa$V3)
fpa$u4 <- fpa$V4 + markdown4(fpa$V4)

plot(fpa$u1, fpa$V1, xlim=c(0, max(fpa$u1)), ylim=c(0, max(fpa$u1)), main="Bidder 1's bid as function of (implied) valuation", xlab="Valuation", ylab="Bid")
abline(0,1)
plot(fpa$u2, fpa$V2, xlim=c(0, max(fpa$u2)), ylim=c(0, max(fpa$u2)), main="Bidder 2's bid as function of (implied) valuation", xlab="Valuation", ylab="Bid")
abline(0,1)
plot(fpa$u3, fpa$V3, xlim=c(0, max(fpa$u3)), ylim=c(0, max(fpa$u3)), main="Bidder 3's bid as function of (implied) valuation", xlab="Valuation", ylab="Bid")
abline(0,1)
plot(fpa$u4, fpa$V4, xlim=c(0, max(fpa$u4)), ylim=c(0, max(fpa$u4)), main="Bidder 4's bid as function of (implied) valuation", xlab="Valuation", ylab="Bid")
abline(0,1)

#construct f_U
u_grid <- seq(0, 325, 0.2)

u1_L <- quantile(fpa$u1, probs=0.25)
u1_H <- quantile(fpa$u1, probs=0.75)
u2_L <- quantile(fpa$u2, probs=0.25)
u2_H <- quantile(fpa$u2, probs=0.75)
u3_L <- quantile(fpa$u3, probs=0.25)
u3_H <- quantile(fpa$u3, probs=0.75)
u4_L <- quantile(fpa$u4, probs=0.25)
u4_H <- quantile(fpa$u4, probs=0.75)

u1_vals <- c(u1_L, u1_H)
u2_vals <- c(u2_L, u2_H)
u3_vals <- c(u3_L, u3_H)
u4_vals <- c(u4_L, u4_H)
eval_pts <- expand.grid(u1_vals, u2_vals, u3_vals, u4_vals)
colnames(eval_pts) <- c("u1", "u2", "u3", "u4")
eval_pts$F_U <- 0



  
for (i in 1:16){
  u1 <- eval_pts$u1[i]
  u2 <- eval_pts$u2[i]
  u3 <- eval_pts$u3[i]
  u4 <- eval_pts$u4[i]
  eval_pts$F_U[i ] <- sum(fpa$u1 <= u1 & fpa$u2 <= u2 & fpa$u3 <= u3 & fpa$u4 <= u4)/length(fpa$u1)
}
  
#selected using rule of thumb - probably not optimal, but what can you do
h1 <- 17.27
h2 <- 18.45
h3 <- 17.57
h4 <- 18.13

marg_dist <- data.frame(matrix(nrow = length(u_grid), ncol=5))
colnames(marg_dist) <- c("u", "dens1", "dens2", "dens3", "dens4")
for (i in 1:length(u_grid)){
  u <- u_grid[i]
  dens1 <- (1/h1)*mean(dnorm((u - fpa$u1)/h_f))
  dens2 <- (1/h2)*mean(dnorm((u - fpa$u2)/h_f))
  dens3 <- (1/h3)*mean(dnorm((u - fpa$u3)/h_f))
  dens4 <- (1/h4)*mean(dnorm((u - fpa$u4)/h_f))
  marg_dist[i,] <- c(u, dens1, dens2, dens3, dens4)
}
plot(marg_dist$u, marg_dist$dens1, col="black", type='l', xlab="Valuation", ylab="Density", main="Estimated marginal densities of valuation by bidder")
lines(marg_dist$u,marg_dist$dens2, col="red")
lines(marg_dist$u, marg_dist$dens3, col="blue")
lines(marg_dist$u, marg_dist$dens4, col="green")
legend(220, 0.005,legend=c("Bidder 1", "Bidder 2", "Bidder 3", "Bidder 4"), col=c("black", "red", "blue", "green"), lty=1)

#using the ready-made density to check whether I'm correct and to also find the right bandwidth to use
est1 <- density(fpa$u1, from=0, to=max(u_grid))
est2 <- density(fpa$u2, from=0, to=max(u_grid))
est3 <- density(fpa$u3, from=0, to=max(u_grid))
est4 <- density(fpa$u4, from=0, to=max(u_grid))
plot(est1)
plot(est2)
plot(est3)
plot(est4)
lines(est1$x, est1$y, type="l")
lines(est2$x, est2$y)
lines(est3$x, est3$y)
lines(est4$x, est4$y)


#test for correlation
cor(fpa$u1, fpa$u2)
cor(fpa$u1, fpa$u3)
cor(fpa$u1, fpa$u4)
cor(fpa$u2, fpa$u3)
cor(fpa$u2, fpa$u4)
cor(fpa$u3, fpa$u4)


#plot of marginal distribution of valuations assuming independence and symmetry
f_indsym <- density(unlist(stack(fpa, select=c("u1", "u2", "u3", "u4"))[1]), from=0, to=max(u_grid))

plot(f_indsym, lwd=3, type='l', xlab="Valuation", ylab="Density", main="Marginal density assuming independence and symmetry \n (plotted with individual marginals for comparison)")
lines(marg_dist$u, marg_dist$dens1, col="orange")
lines(marg_dist$u,marg_dist$dens2, col="red")
lines(marg_dist$u, marg_dist$dens3, col="blue")
lines(marg_dist$u, marg_dist$dens4, col="green")
legend(200, 0.005,legend=c("Marginal distribution \n (under ind + sym) ", "Bidder 1 Marginal", "Bidder 2 Marginal", "Bidder 3 Marginal", "Bidder 4 Marginal"), col=c("black", "red", "orange", "blue", "green"), lty=1, lwd=c(3, 1,1,1,1))


f_dist <- approxfun(f_indsym$x, f_indsym$y, yleft=0, yright=0)
f_cdf <- apply(cbind(-Inf, f_indsym$x), 1, function(x) {integrate(f_dist, lower=x[1], upper=x[2])$value}) #ecdf(unlist(stack(fpa, select=c("u1", "u2", "u3", "u4"))[1]))

f_1 <- approxfun(est1$x, est1$y, yleft=0, yright=0)
cdf_1 <- apply(cbind(0, est1$x), 1, function(x) {integrate(f_1, lower=x[1], upper=x[2])$value})
f_2 <- approxfun(est2$x, est2$y, yleft=0, yright=0)
cdf_2 <- apply(cbind(0, est2$x), 1, function(x) {integrate(f_2, lower=x[1], upper=x[2])$value})
f_3 <- approxfun(est3$x, est3$y, yleft=0, yright=0)
cdf_3 <- apply(cbind(0, est3$x), 1, function(x) {integrate(f_3, lower=x[1], upper=x[2])$value})
f_4 <- approxfun(est4$x, est4$y, yleft=0, yright=0)
cdf_4 <- apply(cbind(0, est4$x), 1, function(x) {integrate(f_4, lower=x[1], upper=x[2])$value})

plot(f_indsym$x,f_cdf, type='l',lwd=5,  main ="CDF of valuations assuming independence and symmetry \n (plotted with individual marginals for comparison)")
lines(est1$x, cdf_1, col="orange")
lines(est2$x, cdf_2, col="red")
lines(est3$x, cdf_3, col="blue")
lines(est4$x, cdf_4, col="green")
legend(210, 0.4,legend=c("CDF of valuations \n (under ind + sym) ", "Bidder 1", "Bidder 2", "Bidder 3", "Bidder 4 "), col=c("black", "orange", "red", "blue", "green"), lty=1, lwd=c(5, 1,1,1,1))



#symmetric independent private values
h_g <- 20
#second_highest <- function(x) order(x, decreasing=TRUE)[2]
#fpa$second <- apply(cbind(fpa$V1, fpa$V2, fpa$V3, fpa$V4),1, function(x) x[second_highest(x)])
#g <- density(fpa$second, from=0, to=max(u_grid))
#g_int <- approxfun(g$x, g$y, yleft=0, yright=0)
#G_vals <- apply(cbind(-Inf, g$x), 1, function(x) {integrate(g_int, lower=x[1], upper=x[2])$value})
#G <- approxfun(g$x, G_vals, yleft=0, yright=1)

g <- f_indsym
G_vals <- f_cdf
G <- approxfun(g$x, G_vals, yleft=0, yright=1)

fpa$uhat_1 <- fpa$V1 + G(fpa$V1)/(3*g_int(fpa$V1))
fpa$uhat_2 <- fpa$V2 + G(fpa$V2)/(3*g_int(fpa$V2))
fpa$uhat_3 <- fpa$V3 + G(fpa$V3)/(3*g_int(fpa$V3))
fpa$uhat_4 <- fpa$V4 + G(fpa$V4)/(3*g_int(fpa$V4))

f1 <- density(fpa$uhat_1, from=0, to=max(u_grid))
f2 <- density(fpa$uhat_2, from=0, to=max(u_grid))
f3 <- density(fpa$uhat_3, from=0, to=max(u_grid))
f4 <- density(fpa$uhat_4, from=0, to=max(u_grid))

f_SIPV <- density(unlist(stack(fpa, select=c("uhat_1", "uhat_2", "uhat_3", "uhat_4"))[1]), from=0, to=max(u_grid))
f_SIPV_int <- approxfun(f_SIPV$x, f_SIPV$y, yleft=0, yright=0)
F_SIPV <- apply(cbind(-Inf, u_grid), 1, function(x) {integrate(f_SIPV_int, lower=x[1], upper=x[2])$value})
F_SIPV_fun <- approxfun(u_grid, F_SIPV, yleft=0, yright=1)


plot(f_SIPV, lwd=3, main ="Marginal distribution of valuations under initial SIPV assumption \n (plotted with previous marginal for comparison)")
lines(f_indsym$x,f_indsym$y, type="l", col="red", lwd=3)
legend(150, 0.008, legend=c("Marginal distribution (under init. SIPV assn.) ", "Previous marginal (not under initial SIPV)"), col=c("black", "red" ), lty=1, lwd=c(3, 3))

plot(F_SIPV_fun, xlim=c(0,300), type='l', col="black", lwd=3, main="CDF of valuations under initial SIPV assumption \n (plotted with previous CDF for comparison")
lines(f_indsym$x, f_cdf, lwd=3, col="red")
legend(120, 0.15, legend=c("CDF of valuations (under init. SIPV assn.) ", "Previous CDF (not under initial SIPV)"), col=c("black", "red" ), lty=1, lwd=c(3, 3))

