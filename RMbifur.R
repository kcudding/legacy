#RMbifur
rm(list=ls())
library(rootSolve)
library(deSolve)

#a = 2; c = 10; r = 1.3; p0 = 17; k = 20; m = .010
a = 1.; c = 2; r = .5; p0 = .16; k = 20; m = .2; mu=0.0005
Pars <- c(a = a, c = c, r = r, p0 = p0, k = k, m =m, mu=mu)
state <- c(x = 0.05, y = 0.1)
RMpredsim <- function (t,state,Pars) {
  with(as.list(c(state,Pars)), {
    dx = r*x[1]*(1-x[1]/k)-a*x[1]*x[2]/(p0+x[1])
    dy = a/c*x[1]*x[2]/(p0+x[1])-m*x[2]
    return(list(c(dx, dy)))
  })
}

rosenarthur<-function(t,y,p){
  with(as.list(c(y,p)), {
    dx= r*y[1]*(1-(y[1]/k)) - (y[1]*y[2]/(p0+y[1]))
    dy= y[2]*((a/c)*y[1]/(p0+y[1]) - m)
    
    list(c(dx,dy))
  })
}

#bifurcation(RMpredsim, state, Pars)

RMpredsolv <- function (x,Pars) {
  with(as.list(c(Pars)), {
    dx = r*x[1]*(1-x[1]/k)-a*x[1]*x[2]/(p0+x[1])
    dy = a/c*x[1]*x[2]/(p0+x[1])-m*x[2]
    return((c(dx, dy)))
  })
}

#set up jacobian

dx=expression(r*x*(1-x/k)-a*x*y/(p0+x))
dy = expression(a/c*x*y/(p0+x)-m*y)

dxx=D(dx, "x")
dxy=D(dx, "y")
dyx=D(dy, "x")
dyy=D(dy, "y")

J<-expression(matrix(c(eval(dxx), eval(dxy), eval(dyx), eval(dyy)), nrow=2, byrow=TRUE))


kl=seq(from=0.01, to=.8, by=0.01)
stor=setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
              c("xeq", "yeq", "eg1", "eg2", "k","x1","y2", "egk1", "egk2",
                "mxx", "mnx"))
for (i in seq_along(kl)) {
k=kl[i]
Pars <- c(a = a, c = c, r = r, p0 = p0, k = k, m =m)
xapr=expression(m*p0/(a*c-m))
xstar=eval(xapr)
yapr=expression(r/a*(1-xstar/k)*(p0+xstar))
                                           
ss <- multiroot(f = RMpredsolv, parms=Pars,start = c(eval(xapr), eval(yapr)))
ssk <- multiroot(f = RMpredsolv, parms=Pars,start = c(k, 0))
xeq=ss$root[1]
yeq=ss$root[2]
x1=ssk$root[1]
y2=ssk$root[2]

x=xeq
y=yeq
eg1=eigen(eval(J))$values
x=x1
y=y2
egk=eigen(eval(J))$values




out.time<-0:1000
init.state<-c(.01,.01)
out.state<-ode(y=init.state,times=out.time,func=rosenarthur,parms=Pars,rtol=1e-5)
indlen=c(nrow(out.state)-50,nrow(out.state)) 
mxx=max(out.state[indlen[1]:indlen[2],2])
mnx=min(out.state[indlen[1]:indlen[2],2])

stor[i,]=c(xeq,yeq,eg1, k,x1,y2, egk, mxx, mnx)
}

stor[,1]=Re(stor[,1])
stor[,2]=Re(stor[,2])
stor[,5]=Re(stor[,5])
stor$col=1
stor$col=ifelse(Re(stor[,3])<0 & Re(stor[,4])<0, stor$col,2)
stor$col2=1
stor$col2=ifelse(Re(stor$egk1)<0 & Re(stor$egk2)<0, stor$col2,2)
which(stor$col==2)
plot(x=kl,stor$x1, col=stor$col2, type="l", pch=16, cex=.5, lty=2)
kdex=which(stor$col2==2)
lines
points(kl, stor$xeq,col=stor$col, type="l", pch=16, cex=.5) 
points(kl, stor$mxx,col=stor$col, type="l", pch=16, cex=.5) 
points(kl, stor$mnx,col=stor$col, type="l", pch=16, cex=.5) 
#######

rosendecay<-function(t,y,p){
  with(as.list(c(y,p)), {

    dx= r*y[1]*(1-(y[1]/y[3])) - (y[1]*y[2]/(p0+y[1]))
    dy= y[2]*((a/c)*y[1]/(p0+y[1]) - m)
    dk=mu
    list(c(dx,dy, dk))
  })
}

state <- c(x = 0.05, y = 0.1)
out.time<-0:1000
init.state<-c(.1,.1, .01)
Pars[6]=0.2
out.state<-ode(y=init.state,times=out.time,func=rosendecay,
               parms=Pars,rtol=1e-5)
lines(out.state[,4],out.state[,2], col="green")
init.state<-c(.1,.1, .3)
Pars[6]=0.2
out.state<-ode(y=init.state,times=out.time,
      func=rosendecay,parms=Pars,rtol=1e-5, method="ode23")
lines(out.state[,4],out.state[,2], col="blue")
##library("nleqslv")
#https://staff.fnwi.uva.nl/a.m.deroos/projects/BifurcationTheory/43-HopfPoint-Rosenzweig.html