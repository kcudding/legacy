
rm(list=ls())
library(rootSolve)
library(deSolve)

#a = 2; c = 10; r = 1.3; p0 = 17; k = 20; m = .010; j=0.01; h=0.02
a = 30.5; b=0.8; c =5.3; f = .098; m = .9
j=0.3; h=0.05; mu=0.0005



Pars <- c(a = a, b=b, c=c, m=m, f=f, j=j, mu=mu)
init.state<-c(.163,.012, 0.01, j)


#numerical solution

RMenginedecay<-function(t,y,p){
  with(as.list(c(y,p)), {
    dxdt= y[1]*(1-y[1]/(1+y[3])-a*y[2]/(b+y[1]))
  #  dydt= y[2]*(c*y[1]/(b+y[1]) - m*(1-f*y[3]))
    dydt= y[2]*(c*y[1]/(b+y[1]) - m + f*y[3])
    dzdt= -y[4]*y[3]+ y[2]
    dmudt= mu
    list(c( dxdt,dydt, dzdt, dmudt))
  })
}

RMengine<-function(t,y,p){
  with(as.list(c(y,p)), {
    dxdt= y[1]*(1-y[1]/(1+y[3])-a*y[2]/(b+y[1]))
   # dydt= y[2]*(c*y[1]/(b+y[1]) - m*(1-f*y[3]))
    dydt= y[2]*(c*y[1]/(b+y[1]) - m + f*y[3])
    dzdt= -j*y[3]+ y[2]
    
    list(c( dxdt,dydt, dzdt))
  })
}


out.time<-0:5000

out.state2<-ode(y=init.state,times=out.time,func=RMenginedecay,parms=Pars,rtol=1e-4)
#par(mfrow=c(1,2))
plot(out.state2[,1],out.state2[,5])

Pars[7]=0.00005

out.state3<-ode(y=init.state,times=out.time,func=RMenginedecay,
                parms=Pars,rtol=1e-4, method="rk4")
#par(mfrow=c(1,2))
plot(out.state3[,1],out.state3[,5])

init.state2<-init.state[1:3]
out.state<-ode(y=init.state2,times=out.time,func=RMengine,parms=Pars,rtol=1e-4)

matplot(out.state[,1],out.state[,2:4],
        #   ylim=c(0,.051),
        xlim=c(0,1000),
        lwd=2,lty=1,
        type="l",ylab="density",xlab="time",main="Rosenzweig-MacArthur")
legend("topright", c("prey", "pred", "env"), col=c(1,2,3), 
       lwd=2,lty=1, bty="n")
lines(out.state2[,1],out.state2[,2], col="purple")
lines(out.state3[,1],out.state3[,2], col="orange")

#set up jacobian
dx<-expression(x*(1-(x/(1+z))) - (a*x*y)/(b+x))
dy<-expression((c*x*y)/(b+x) - m*(1-f*z)*y)
dz<-expression(-j*z+y)
dxx=D(dx, "x")
dxy=D(dx, "y")
dxz=D(dx, "z")
dyx=D(dy, "x")
dyy=D(dy, "y")
dyz=D(dy, "z")
dzx=D(dz, "x")
dzy=D(dz, "y")
dzz=D(dz, "z")

J<-expression(matrix(c(eval(dxx), eval(dxy), eval(dxz), eval(dyx), 
                       eval(dyy), eval(dyz), eval(dzx), eval(dzy), eval(dzz)),
                     nrow=3, byrow=TRUE))



#find eigenvalues for interior eq'm (2D)

mu=0.4
alpha=1.6

##approximate eq'm for starting conditions in rootsolver

RMengsolv <- function (x,Pars) {
  with(as.list(c(Pars)), {
    dxdt= x[1]*(1-x[1]/(1+x[3])-a*x[2]/(b+x[1]))
 #   dydt= x[2]*(c*x[1]/(b+x[1]) - m*(1-f*x[3]))
    dydt= x[2]*(c*x[1]/(b+x[1]) - m + f*x[3])
    dzdt= -j*x[3]+ x[2]
    return((c(dxdt, dydt,dzdt)))
  })
}



##### bifurcate decay
jseq=seq(from=0.01, to=0.6, by=0.005)
eqval=data.frame(eint=double(),ek=double(), yint=double(),yk=double(),
                 zint=double(),zk=double(),
                 s1=integer(), s2=integer(), 
                 mxx=double(), mnx=double())

for (i in seq_along(jseq)){
  Pars["j"] <- jseq[i]
  ss <- multiroot(f = RMengsolv, parms=Pars,jacfunc=J,
                 # start = c(eval(xapr), eval(yapr), eval(eapr)))
                 start = init.state2)
  ssk <- multiroot(f = RMengsolv, parms=Pars,start = c(1, 0,0))
  xeq=ss$root[1]
  yeq=ss$root[2]
  zeq=ss$root[3]
  x1=ssk$root[1]
  y2=ssk$root[2]
  e2=ssk$root[3]
  
  x=xeq
  y=yeq
  z=zeq
  
  
  
  out.time<-0:15000

  init.state<-c(.03,.012,0.18)
  init.stateQ<-c(xeq+0.0001,yeq+0.001,zeq+0.001)
  out.state<-ode(y=init.stateQ,times=out.time,func=RMengine,
                 parms=Pars,rtol=1e-7, method="ode23")
  indlen=c(nrow(out.state)-2000,nrow(out.state)) 
  mxx=max(out.state[indlen[1]:indlen[2],2])
  mnx=min(out.state[indlen[1]:indlen[2],2])
  eqval[i,]=c(ss$root[1], ssk$root[1], ss$root[2], ssk$root[2],
              ss$root[3], ssk$root[3],NA, NA, mxx, mnx)
  eqval[i,"s1"]=ifelse(all(Re(eigen(eval(J))$values)<0), 1,2)
  eqval[i,"s2"]=ifelse(all(Re(eigen(eval(J))$values)<0), 1,2)
  
  
}
stor=eqval
stor[,1]=Re(stor[,1])
stor[,2]=Re(stor[,2])



out.time<-0:10000
inx=which(round(jseq,4)==0.02)

init.stateQ2=c(stor$eint[inx]+0.001,stor$yint[inx]+0.001,
             stor$zint[inx]+0.001,jseq[inx] )
Pars["mu"]=0.00005
out.state2<-ode(y=init.stateQ2,times=out.time,func=RMenginedecay,parms=Pars,
                rtol=1e-7, method="ode23")
Pars["mu"]=0.0005
out.state3<-ode(y=init.stateQ2,times=out.time,func=RMenginedecay,parms=Pars,
                rtol=1e-7, method="ode23")
#matplot(out.state2[,1], cbind(out.state2[,2], out.state3[,2]), type="l")
Pars["mu"]=0.005
out.state4<-ode(y=init.stateQ2,times=out.time,func=RMenginedecay,parms=Pars,
                rtol=1e-7, method="ode23")
Pars["mu"]=0.05
out.state5<-ode(y=init.stateQ2,times=out.time,func=RMenginedecay,parms=Pars,
                rtol=1e-7, method="ode23")

idd=which(stor$s1==1)
bind=idd[1]
jpeg('baseR_figure0.jpeg', width =180, 
     height = 150,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

plot(x=jseq,stor$eint, lty=2, 
     type="l", col=1, lwd=2, las=1,
     ylim=c(0,0.7),
     xlim=c(0, 0.5),
     ylab="prey density (x)",
     xlab="environmental decay rate (j)",
     cex.lab=1.5, bty="L"
)
lines(jseq[1:bind+3], stor$mnx[1:bind+3],col=1, lty=2, type="l", lwd=2) 
lines(jseq[1:bind+3], stor$mxx[1:bind+3],col=1, lty=2,  type="l",  lwd=2) 

lines(x=jseq,stor$eint, lty=2, type="l", col=1, lwd=2)
lines(jseq[bind:length(jseq)], stor$eint[bind:length(jseq)],lty=1, lwd=2.1)
legend("topright", legend=c("unstable", "stable"), 
       lty=c(2,1), lwd=2, bty="n")

dev.off()

#########

jpeg('baseR_figure.jpeg', width =180, 
     height = 150,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)

plot(x=jseq,stor$eint, lty=2, 
     type="l", col=1, lwd=2, las=1,
     ylim=c(0,0.3),
     xlim=c(0, 0.5),
     ylab="prey density (x)",
     xlab="environmental decay rate (j)",
     cex.lab=1.5, bty="L"
)


lines(out.state2[,5],out.state2[,2], col="pale goldenrod", typ="l")
lines(out.state3[,5],out.state3[,2], 
      col="purple", typ="l")

lines(out.state4[,5],out.state4[,2], col="green", typ="l", lwd=2)
lines(out.state5[,5],out.state5[,2], col="red", typ="l", lwd=2)

points(jseq[1:bind+3], stor$mnx[1:bind+3],col=1,  type="l", lwd=2) 
points(jseq[1:bind+3], stor$mxx[1:bind+3],col=1,  type="l",  lwd=2) 

lines(x=jseq,stor$eint, lty=2, type="l", col=1, lwd=2)
lines(jseq[bind:length(jseq)], stor$eint[bind:length(jseq)],lty=1, lwd=2.1)


legend("bottomright", legend=c(0.00005,0.0005,0.005,0.05), lty=1, lwd=2,
       col=c("pale goldenrod", "purple", "green", "red"), bty="n", title="rate of change")

dev.off()



par(mfrow=c(1,3))

#########
inx=which(round(jseq,4)==0.02)
init.stateQ1=c(stor$eint[inx]+0.001,stor$yint[inx]+0.001,
               stor$zint[inx]+0.001)
Pars["j"]=jseq[inx] 
out.state<-ode(y=init.stateQ1,times=out.time,func=RMengine,parms=Pars,
                rtol=1e-7, method="ode23")
matplot(out.state[9800:10000,1], 
        out.state[9800:10000,2:4], pch=16, cex=0.25, 
        bty="l", ylim=c(0,1.8),
        type="l", lwd=2,cex.lab=1.5,
        xlab="time", ylab="density")
#legend("top", legend=c("prey (x)", "predator(y)", "environment (z)"),
 #      lwd=2,ncol=3,
  #     col=c(1,2,3), lty=c(1,2,3), bty="n")
######

#########
inx=35
init.stateQ1=c(stor$eint[inx]+0.001,stor$yint[inx]+0.001,
               stor$zint[inx]+0.001)
Pars["j"]=jseq[inx] 
out.state<-ode(y=init.stateQ1,times=out.time,func=RMengine,parms=Pars,
               rtol=1e-7, method="ode23")
matplot(out.state[9800:10000,1], 
        out.state[9800:10000,2:4], pch=16, cex=0.25, 
        bty="l", ylim=c(0,1.8),
        type="l", lwd=2,cex.lab=1.5,
        xlab="time", ylab="density")
#legend("top", legend=c("prey (x)", "predator(y)", "environment (z)"),
 #      lwd=2,ncol=3,
#       col=c(1,2,3), lty=c(1,2,3), bty="n")
######

#########
inx=98
init.stateQ1=c(stor$eint[inx]+0.001,stor$yint[inx]+0.001,
               stor$zint[inx]+0.001)
Pars["j"]=jseq[inx] 
out.state<-ode(y=init.stateQ1,times=out.time,func=RMengine,parms=Pars,
               rtol=1e-7, method="ode23")
matplot(out.state[9800:10000,1], 
        out.state[9800:10000,2:4], pch=16, cex=0.25, 
        bty="l", ylim=c(0,1.8),
        type="l", lwd=2,cex.lab=1.5,
        xlab="time", ylab="density")
legend("topright", legend=c("prey (x)", "predator(y)", "environment (z)"),
       lwd=2,
      col=c(1,2,3), lty=c(1,2,3), bty="n")
######

jpeg('baseR_figure3.jpeg', width =180, 
     height = 150,
     units = 'mm', pointsize = 12,
     quality = 75,
     res = 300)
plot(x=out.state2[,1], y=out.state2[,5], typ="l", col="pale goldenrod",
     xlab="time", ylab="environmental decay rate (j)",
     ylim=c(0,10), lwd=2, bty="L", cex.lab=1.5)
lines(x=out.state3[,1], y=out.state3[,5],col="purple",lwd=2)
lines(x=out.state4[,1], y=out.state4[,5],col="green",lwd=2)
lines(x=out.state[,1], y=out.state5[,5],col="red",lwd=2)
legend("topright", legend=c(0.00005,0.0005,0.005,0.05), lty=1, lwd=2,
       col=c("pale goldenrod", "purple", "green", "red"), bty="n", title="rate of change")
dev.off()
########
xstar=eval(xapr)
ystar=eval(yapr)
estar=eval(eapr)

init.state<-c(xeq,yeq, eeq)
out.state<-ode(y=init.state,times=out.time,func=RMengine,parms=Pars,rtol=1e-5)

plot(out.state[,1], out.state[,4])
matplot(out.state[,1], out.state[,2:4], typ="l", lwd=2)
init.state<-c(xeq,yeq, eeq, 0.0001)
out.state<-ode(y=init.state,times=out.time,func=RMenginedecay,parms=Pars,rtol=1e-5)

plot(out.state[,1], out.state[,5])
matplot(out.state[,1], out.state[,2:4], typ="l", lwd=2)

#########
x=x1
y=y2
e=e2
egk=eigen(eval(J))$values
ss0<- multiroot(f = RMengsolv, parms=Pars,start = c(0, 0,0))
x=0;y=0;e=0
eg0=eigen(eval(J))$values
#abline(h=xstar, col="red", lty=2, lwd=1.5)
#abline(h=ystar, col="blue", lty=2, lwd=1.5)
#abline(h=s3, col="dark green", lty=3, lwd=1.5)
#abline(h=s4, col="purple", lty=3, lwd=1.5)



#text(3000,1.8,paste("eig=",round(eg1[1],2),", ",round(eg1[2],2),round(eg1[2],3)))
#text(3000,1.7,paste("eig=",round(egk[1],2),", ",round(egk[2],2),round(egk[2],3)))
#legend("topright", legend=c("mu/(1-mu)", "mu/(alpha-mu)", "(1-(mu/(1-mu)))(1-(mu/(alpha-mu)))",  "(1-(mu/(alpha-mu)))(1-(mu/(alpha-mu)))"), col=c("blue","dark green", "red", "purple"), lty=c(2,3,2,3), bty="n" )

yiso=eval(expression(r/a*(1-x/(k+e)*(p0+x))))

###3D


library("plot3D")

lines3D(out.state[,2], out.state[,3],out.state[,4],colvar = out.state[,2], col = NULL, add = FALSE,theta = 5, phi = 25)
lines3D(out.state4[4500:5000,2], out.state4[4500:5000,3], out.state4[4500:5000,4],pch=16, cex=0.1)
lines3D(out.state[4500:5000,2], out.state[4500:5000,3], out.state[4500:5000,4],pch=16, cex=0.1)
