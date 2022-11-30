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
  m <- grid[ind]
  kd_est$G1[ind] <- (1/(h))*mean((dnorm((m - fpa$V1)/h)*(fpa$M1 < m) ))#+ dnorm((m - fpa$V2)/h)*(fpa$M2 < m) + dnorm((m - fpa$V3)/h)*(fpa$M3 < m)+ dnorm((m - fpa$V4)/h)*(fpa$M4 < m) ))
  kd_est$G2[ind] <- (1/(h))*mean((dnorm((m - fpa$V2)/h)*(fpa$M2< m) ))
  kd_est$G3[ind] <- (1/(h))*mean((dnorm((m - fpa$V3)/h)*(fpa$M3< m) ))
  kd_est$G4[ind] <- (1/(h))*mean((dnorm((m - fpa$V4)/h)*(fpa$M4< m) ))
  kd_est$g1[ind] <- (1/(h^2))*mean(dnorm((m - fpa$V1)/h)*dnorm((m - fpa$M1)/h))# + dnorm((m - fpa$V2)/h)*dnorm((m - fpa$M2)/h) + dnorm((m - fpa$V3)/h)*dnorm((m - fpa$M3)/h) + dnorm((m - fpa$V4)/h)*dnorm((m - fpa$M4)/h) )
  kd_est$g2[ind] <- (1/(h^2))*mean(dnorm((m - fpa$V2)/h)*dnorm((m - fpa$M2)/h))
  kd_est$g3[ind] <- (1/(h^2))*mean(dnorm((m - fpa$V3)/h)*dnorm((m - fpa$M3)/h))
  kd_est$g4[ind] <- (1/(h^2))*mean(dnorm((m - fpa$V4)/h)*dnorm((m - fpa$M4)/h))
}

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


#construct f_U
u_grid <- seq(0, 280, 0.2)

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
lines(marg_dist$u, marg_dist$dens1, col="black", type='l')
lines(marg_dist$u,marg_dist$dens2, col="red")
lines(marg_dist$u, marg_dist$dens3, col="blue")
lines(marg_dist$u, marg_dist$dens4, col="green")


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


#idea: if they follow the same distribution, just combine them


#symmetric independent private values
h_g <- 20
second_highest <- function(x) order(x, decreasing=TRUE)[2]
fpa$second <- apply(cbind(fpa$V1, fpa$V2, fpa$V3, fpa$V4),1, function(x) x[second_highest(x)])
G <- ecdf(fpa$second)
g <- density(fpa$second)
g_int <- approxfun(x=g$x, y=g$y, rule = 2)


Gmax <- ecdf(fpa$V2)
gmax <- density(fpa$V2)
gmax_int <- approxfun(x = gmax$x, y=g$y, rule=2)
fpa$max_md <- fpa$V2 + Gmax(fpa$V2)/(3*gmax_int(fpa$V2))
max_d <- density(fpa$max_md)

fpa$uhat_1 <- fpa$V1 + G(fpa$V1)/(3*g_int(fpa$V1))
fpa$uhat_2 <- fpa$V2 + G(fpa$V2)/(3*g_int(fpa$V2))
fpa$uhat_3 <- fpa$V3 + G(fpa$V3)/(3*g_int(fpa$V3))
fpa$uhat_4 <- fpa$V4 + G(fpa$V4)/(3*g_int(fpa$V4))

f1 <- density(fpa$uhat_1)
f2 <- density(fpa$uhat_2)
f3 <- density(fpa$uhat_3)
f4 <- density(fpa$uhat_4)

f <- density(unlist(stack(fpa, select=c("uhat_1", "uhat_2", "uhat_3", "uhat_4"))[1]))
plot(f)
lines(f1$x,f1$y, type="l")
lines(f2$x, f2$y, col="red")
lines(f3$x, f3$y, col="blue")
lines(f4$x, f4$y, col="green")     

ggplot(fpa) + stat_density(aes(uhat_1), color="red", fill="red", alpha=0.2) + stat_density(aes(uhat_2), color="blue", fill="blue", alpha=0.2) + stat_density(aes(uhat_3),  color="green",fill="green", alpha=0.2) + stat_density(aes(uhat_2), fill="purple", alpha=0.2)

