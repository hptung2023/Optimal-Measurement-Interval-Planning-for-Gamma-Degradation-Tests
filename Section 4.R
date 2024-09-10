library(dplyr)
library(MASS)
library(pracma)

#read data
data=read.csv("LIGHT INTENSITY DEGRADATION DATA OF 12 LEDS.csv",header=T)
#parameter estimation
est<-fitdistr(na.omit(data$incre), "gamma")
alpha<-as.numeric(est$estimate[1])
beta<-as.numeric(est$estimate[2])
ks.test(data$incre, "pgamma", alpha,beta ,exact=FALSE)

#cost setting
n<-12
m<-5
Ti<-250
tmin<-5
deltat<-50
ci<-64
cm<-0.9
co<-0.1
costor<-ci*n+cm*m*n+co*Ti

ci<-64/costor
cm<-0.9/costor
co<-0.1/costor
ct<-1

np<-1 #scale objective function, 1 is no scale
npd=1 #scale D objective function , 1 is no scale
t1<-50
alpha<-as.numeric(est$estimate[1])/t1
beta<-as.numeric(est$estimate[2])
ga<-log(alpha/(beta))
w0<-50

#A-optimality####
#grid search
A_GS<-function(ct,ci,cm,co){
  op_GS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci-co*tmin)/cm)
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      Ti<-(ct-ci*n-cm*m*n)/co
      if (Ti<=(m*tmin)) {
        break
      }
      
      sol_GS<-(1/((m-1)*tmin^2*trigamma(alpha*tmin)+(Ti-(m-1)*tmin)^2*trigamma(alpha*(Ti-(m-1)*tmin))-Ti/alpha)+1/(alpha*Ti))/(n/np)
      if (sol_GS<op_GS) {
        op_GS<-sol_GS
        op_GS_n<-n
        op_GS_m<-m
        op_GS_T<-Ti
      }
    }
  }
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 1
A_BS<-function(ct,ci,cm,co){
  op_BS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  tr_BS <- function(w) { #objective function
    Ti<-(ct-ci*n-cm*n+cm*w*n/tmin)/(cm*n/tmin+co)
    re<-1/((Ti-w)*tmin*trigamma(alpha*tmin)+w^2*trigamma(alpha*w)-Ti/alpha)+1/(alpha*Ti)
    re<-re/(n/np)
    return(re)
  }
  for (n in 1:nmax) {
    Timax<-(ct-ci*n-cm*n)/co
    op<-optimize(tr_BS,c(tmin,Timax))
    sol_BS<-op$objective
    if (sol_BS<op_BS) { #rounding m by Ti
      Ti<-(ct-ci*n-cm*n+cm*op$minimum*n/tmin)/(cm*n/tmin+co)
      m<-floor((Ti-op$minimum)/tmin+1)
      Ti_mrd<-(ct-ci*n-cm*m*n)/co
      Ti_mru<-(ct-ci*n-cm*(m+1)*n)/co
      if (Ti_mrd<(m*tmin)) {
        obj_mrd<-1e+15
      }else{
        w_mrd<-Ti_mrd-(m-1)*tmin
        obj_mrd<-tr_BS(w_mrd)
      }
      if (Ti_mru<((m+1)*tmin)) {
        obj_mru<-1e+15
      }else{
        w_mru<-Ti_mru-m*tmin
        obj_mru<-tr_BS(w_mru)
      }
      op_BS_tem<-min(obj_mrd,obj_mru)
      if (op_BS_tem<op_BS) {
        op_BS<-op_BS_tem
        op_BS_n<-n
        if (op_BS_tem==obj_mrd) {
          op_BS_m<-m
          op_BS_T<-Ti_mrd
        }else{
          op_BS_m<-m+1
          op_BS_T<-Ti_mru
        }
        
      }
    }
  }
  op<-list("objective"=op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_T)
  return(op)    
  
}

#equal interval
#grid search
A_EI_GS<-function(ct,ci,cm,co){
  op_GS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci)/(cm+co*tmin))
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      t<-(ct-ci*n-cm*m*n)/(co*m)
      if (t<tmin) {
        break
      }
      
      sol_GS<-1/((n/np)*m*t)*(alpha/(alpha*t*trigamma(alpha*t)-1)+1/alpha)
      if (sol_GS<op_GS) {
        op_GS_n<-n
        op_GS_m<-m
        op_GS<-sol_GS
      }
    }
  }
  op_GS_T<-(ct-ci*op_GS_n-cm*op_GS_m*op_GS_n)/(co)
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 2
A_EI_BS<-function(ct,ci,cm,co){
  op_BS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  for (n in 1:nmax) {
    mmax<-floor((ct-ci*n)/(cm*n+co*tmin))
    tr_BS <- function(m) { #objective function
      # browser()
      tn<-(ct-ci*n-cm*n*m)/(co*m)
      re<-1/((n/np)*m*tn)*(alpha/(alpha*tn*trigamma(alpha*tn)-1)+1/alpha)
      return(re)
    }
    if(n==nmax){
      sol_BS <- tr_BS(1)
      if(sol_BS<op_BS){
        op_BS_n <- n
        op_BS_m <- 1
        op_BS_t <- (ct-ci*n-cm*n)/(co)
        op_BS <- sol_BS
      }
    }else{
      op <- optimize(tr_BS,c(1,mmax))
      sol_BS <- op$objective
      mrd <- floor(op$minimum)
      mru <- mrd+1
      sol_BSd <- tr_BS(mrd)
      sol_BSu <- tr_BS(mru)
      ord <- order(c(sol_BSd,sol_BSu,op_BS))
      if(ord[1]==1){
        op_BS_n <- n
        op_BS_m <- mrd
        op_BS_t <- (ct-ci*n-cm*n*mrd)/(co*mrd)
        op_BS <- sol_BSd
      }else if(ord[1]==2){
        op_BS_n <- n
        op_BS_m <- mru
        op_BS_t <- (ct-ci*n-cm*n*mru)/(co*mru)
        op_BS <- sol_BSu
      } 
    }
  }
  op<-list("objective"=op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_t*op_BS_m)
  return(op)
}

#D-optimality####
#grid search
D_GS<-function(ct,ci,cm,co){
  op_GS<-0
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci-co*tmin)/cm)
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      Ti<-(ct-ci*n-cm*m*n)/co
      if (Ti<(m*tmin)) {
        break
      }
      sol_GS <-((n/npd)^2*(alpha*(Ti/npd)*((m-1)*tmin^2*trigamma(alpha*tmin)/npd+((Ti-(m-1)*tmin)/sqrt(npd))^2*trigamma(alpha*(Ti-(m-1)*tmin)))-(Ti/npd)^2))
      if (sol_GS>op_GS) {
        op_GS_n<-n
        op_GS_m<-m
        op_GS_T<-Ti
        op_GS<-sol_GS
      }
      
    }
  }
  
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 1
D_BS<-function(ct,ci,cm,co){
  op_BS<-0
  nmax<-floor((ct-co*tmin)/(ci+cm))
  # browser()
  detI_BS <- function(w) {
    Ti<-(ct-ci*n-cm*n+cm*w*n/tmin)/(cm*n/tmin+co)
    re<-alpha*Ti*((Ti-w)*tmin*trigamma(alpha*tmin)+w^2*trigamma(alpha*w))-Ti^2
    re<-re*n^2/np
    return(re)
  }
  
  for (n in 1:nmax) {
    Timax<-(ct-ci*n-cm*n)/co
    op<-optimize(detI_BS,c(tmin,Timax),maximum=TRUE)
    # print(op$maximum)
    sol_BS<-op$objective
    if (sol_BS>op_BS) {
      T<-(ct-ci*n-cm*n+cm*op$maximum*n/tmin)/(cm*n/tmin+co)
      m<-floor((T-op$maximum)/tmin+1)
      Ti_mrd<-(ct-ci*n-cm*m*n)/co
      Ti_mru<-(ct-ci*n-cm*(m+1)*n)/co
      if (Ti_mrd<(m*tmin)) {
        obj_mrd<-0
      }else{
        w_mrd<-Ti_mrd-(m-1)*tmin
        obj_mrd<-detI_BS(w_mrd)
      }
      if (Ti_mru<((m+1)*tmin)) {
        obj_mru<-0
      }else{
        w_mru<-Ti_mru-m*tmin
        obj_mru<-detI_BS(w_mru)
      }
      op_BS_tem<-max(obj_mrd,obj_mru)
      if (op_BS_tem>op_BS) {
        op_BS<-op_BS_tem
        op_BS_n<-n
        if (op_BS_tem==obj_mrd) {
          op_BS_m<-m
          op_BS_T<-Ti_mrd
        }else{
          op_BS_m<-m+1
          op_BS_T<-Ti_mru
        }
        
      }
    }
  }
  
  op<-list("objective"=op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_T)
  return(op)    
}

#equal interval
#grid search
D_EI_GS<-function(ct,ci,cm,co){
  op_GS<-0
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci)/(cm+co*tmin))
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      t<-(ct-ci*n-cm*m*n)/(co*m)
      if (t<tmin) {
        break
      }
      sol_GS<-((n^2)*(m^2)*(alpha*t^3*trigamma(alpha*t)-t^2))/np
      if (sol_GS>op_GS) {
        op_GS_n<-n
        op_GS_m<-m
        op_GS<-sol_GS
      }
    }
  }
  op_GS_t<-(ct-ci*op_GS_n-cm*op_GS_m*op_GS_n)/(co*op_GS_m)
  op_GS_T<-op_GS_m*op_GS_t
  op_GS<-(op_GS_n^2*(op_GS_m^2)*(alpha*op_GS_t^3*trigamma(alpha*op_GS_t)-op_GS_t^2))/np
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 2
D_EI_BS<-function(ct,ci,cm,co){
  op_BS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  for (n in 1:nmax) {
    mmax<-floor((ct-ci*n)/(cm*n+co*tmin))
    det_BS <- function(m) { #objective function
      # browser()
      tn<-(ct-ci*n-cm*n*m)/(co*m)
      re<-((n^2)*(m^2)*(alpha*tn^3*trigamma(alpha*tn)-tn^2))/np
      return(-re)
    }
    if(n==nmax){
      sol_BS <- det_BS(1)
      if(sol_BS<op_BS){
        op_BS_n <- n
        op_BS_m <- 1
        op_BS_t <- (ct-ci*n-cm*n)/(co)
        op_BS <- sol_BS
      }
    }else{
      # browser()
      op <- optimize(det_BS,c(1,mmax))
      sol_BS <- op$objective
      mrd <- floor(op$minimum)
      mru <- mrd+1
      sol_BSd <- det_BS(mrd)
      sol_BSu <- det_BS(mru)
      ord <- order(c(sol_BSd,sol_BSu,op_BS))
      if(ord[1]==1){
        op_BS_n <- n
        op_BS_m <- mrd
        op_BS_t <- (ct-ci*n-cm*n*mrd)/(co*mrd)
        op_BS <- sol_BSd
      }else if(ord[1]==2){
        op_BS_n <- n
        op_BS_m <- mru
        op_BS_t <- (ct-ci*n-cm*n*mru)/(co*mru)
        op_BS <- sol_BSu
      } 
    }
  }
  op<-list("objective"=-op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_t*op_BS_m)
  return(op)
}

#V-optimality####
#partial derivetive of xi_p
f<-function(xi_p){incgam(beta*w0, alpha*xi_p)/gamma(alpha*xi_p)-0.05}
q<-uniroot(f,c(1,1000))
xi_p<-q$root
f=expression(alpha^(alpha*xi_p)*x^(alpha*xi_p-1)*exp(-alpha*x*exp(-ga)-alpha*ga*xi_p)/gamma(alpha*xi_p))
f_xi_p<-function(x) {eval(D(f,'xi_p'))}
f_alpha<-function(x){eval(D(f,'alpha'))}
f_ga<-function(x){eval(D(f,'ga'))}
a<-integrate(f_xi_p,lower=w0,upper=Inf)
b<-integrate(f_alpha,w0,Inf)
c<-integrate(f_ga,w0,Inf)
h1<-b$value/a$value
h2<-c$value/a$value

#grid search
V_GS<-function(ct,ci,cm,co){
  op_GS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci-co*tmin)/cm)
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      Ti<-(ct-ci*n-cm*m*n)/co
      if (Ti<=(m*tmin)) {
        break
      }
      sol_GS<-(h1^2/((m-1)*tmin^2*trigamma(alpha*tmin)+(Ti-(m-1)*tmin)^2*trigamma(alpha*(Ti-(m-1)*tmin))-Ti/alpha)+h2^2/(alpha*Ti))/(n/np)
      if (sol_GS<op_GS) {
        op_GS<-sol_GS
        op_GS_n<-n
        op_GS_m<-m
        
      }
    }
  }
  op_GS_T<-(ct-ci*op_GS_n-cm*op_GS_m*op_GS_n)/co
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 1
V_BS<-function(ct,ci,cm,co){
  op_BS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  avar_BS <- function(w) {
    Ti<-(ct-ci*n-cm*n+cm*w*n/tmin)/(cm*n/tmin+co)
    re<-h1^2/((Ti-w)*tmin*trigamma(alpha*tmin)+w^2*trigamma(alpha*w)-Ti/alpha)+h2^2/(alpha*Ti)
    re<-re/(n/np)
    return(re)
  }
  for (n in 1:nmax) {
    Timax<-(ct-ci*n-cm*n)/co
    op<-optimize(avar_BS,c(tmin,Timax))
    sol_BS<-op$objective
    if (sol_BS<op_BS) {
      T<-(ct-ci*n-cm*n+cm*op$minimum*n/tmin)/(cm*n/tmin+co)
      m<-floor((T-op$minimum)/tmin+1)
      Ti_mrd<-(ct-ci*n-cm*m*n)/co
      Ti_mru<-(ct-ci*n-cm*(m+1)*n)/co
      if (Ti_mrd<(m*tmin)) {
        obj_mrd<-1e+15
      }else{
        w_mrd<-Ti_mrd-(m-1)*tmin
        obj_mrd<-avar_BS(w_mrd)
      }
      if (Ti_mru<((m+1)*tmin)) {
        obj_mru<-1e+15
      }else{
        w_mru<-Ti_mru-m*tmin
        obj_mru<-avar_BS(w_mru)
      }
      op_BS_tem<-min(obj_mrd,obj_mru)
      if (op_BS_tem<op_BS) {
        op_BS<-op_BS_tem
        op_BS_n<-n
        if (op_BS_tem==obj_mrd) {
          op_BS_m<-m
          op_BS_T<-Ti_mrd
        }else{
          op_BS_m<-m+1
          op_BS_T<-Ti_mru
        }
        
      }
    }
  }
  op<-list("objective"=op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_T)
  return(op)    
  
}

#equal interval
#grid search
V_EI_GS<-function(ct,ci,cm,co){
  op_GS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  mmax<-floor((ct-ci)/(cm+co*tmin))
  for (n in 1:nmax) {
    for (m in 1:mmax) {
      t<-(ct-ci*n-cm*m*n)/(co*m)
      if (t<tmin) {
        break
      }
      sol_GS<-1/((n/np)*m*t)*(h1^2*alpha/(alpha*t*trigamma(alpha*t)-1)+h2^2/alpha)
      if (sol_GS<op_GS) {
        op_GS_n<-n
        op_GS_m<-m
        op_GS<-sol_GS
      }
    }
  }
  op_GS_t<-(ct-ci*op_GS_n-cm*op_GS_m*op_GS_n)/(co*op_GS_m)
  op_GS_T<-op_GS_m*op_GS_t
  op<-list("objective"=op_GS,"n"=op_GS_n,"m"=op_GS_m,"T"=op_GS_T)
  return(op)
}

#algorithm 2
V_EI_BS<-function(ct,ci,cm,co){
  op_BS<-1e+15
  nmax<-floor((ct-co*tmin)/(ci+cm))
  for (n in 1:nmax) {
    mmax<-floor((ct-ci*n)/(cm*n+co*tmin))
    var_BS <- function(m) { #objective function
      # browser()
      tn<-(ct-ci*n-cm*n*m)/(co*m)
      re<-1/((n/np)*m*tn)*(h1^2*alpha/(alpha*tn*trigamma(alpha*tn)-1)+h2^2/alpha)
      return(re)
    }
    if(n==nmax){
      sol_BS <- var_BS(1)
      if(sol_BS<op_BS){
        op_BS_n <- n
        op_BS_m <- 1
        op_BS_t <- (ct-ci*n-cm*n)/(co)
        op_BS <- sol_BS
      }
    }else{
      op <- optimize(var_BS,c(1,mmax))
      sol_BS <- op$objective
      mrd <- floor(op$minimum)
      mru <- mrd+1
      sol_BSd <- var_BS(mrd)
      sol_BSu <- var_BS(mru)
      ord <- order(c(sol_BSd,sol_BSu,op_BS))
      if(ord[1]==1){
        op_BS_n <- n
        op_BS_m <- mrd
        op_BS_t <- (ct-ci*n-cm*n*mrd)/(co*mrd)
        op_BS <- sol_BSd
      }else if(ord[1]==2){
        op_BS_n <- n
        op_BS_m <- mru
        op_BS_t <- (ct-ci*n-cm*n*mru)/(co*mru)
        op_BS <- sol_BSu
      } 
    }
  }
  op<-list("objective"=op_BS,"n"=op_BS_n,"m"=op_BS_m,"T"=op_BS_t*op_BS_m)
  return(op)
}

cost<-function(op,ci,cm,co){
  ci*op$n+cm*op$m*op$n+co*op$T
}


#original experiment configuration
A_ob<-1/(n*m*deltat^2*trigamma(alpha*deltat)-n*Ti/alpha)+1/(n*alpha*Ti)
D_ob<-1/(n^2*alpha*Ti*m*deltat^2*trigamma(alpha*deltat)-(n*Ti)^2)
V_ob<-h1^2/(n*m*deltat^2*trigamma(alpha*deltat)-n*Ti/alpha)+h2^2/(n*alpha*Ti)




#optimal designs####
#A-optimality
A_op_GS<-A_GS(ct,ci,cm,co)
cost_GS<-cost(A_op_GS,ci,cm,co)

A_op_BS<-A_BS(ct,ci,cm,co)
cost_BS<-cost(A_op_BS,ci,cm,co)

A_op_EI_GS<-A_EI_GS(ct,ci,cm,co)
cost_EI_GS<-cost(A_op_EI_GS,ci,cm,co)

A_op_EI_BS<-A_EI_BS(ct,ci,cm,co)
cost_EI_BS<-cost(A_op_EI_BS,ci,cm,co)

opsum_A<-matrix(c(unlist(A_op_GS),cost_GS,
                  unlist(A_op_BS),cost_BS),
                ncol=5,byrow=TRUE,
                dimnames = list(c("A-grid search", "A-algorithm 1"),
                                c("obj", "n", "m","T","cost")))
opsum_A_EI<-matrix(c(unlist(A_op_EI_GS),cost_EI_GS,
                     unlist(A_op_EI_BS),cost_EI_BS),
                ncol=5,byrow=TRUE,
                dimnames = list(c("A-equal interval grid search","A-algorithm 2"),
                                c("obj", "n", "m","T","cost")))
#D-optimality
D_op_GS<-D_GS(ct,ci,cm,co)
cost_GS<-cost(D_op_GS,ci,cm,co)

D_op_BS<-D_BS(ct,ci,cm,co)
cost_BS<-cost(D_op_BS,ci,cm,co)

D_op_EI_GS<-D_EI_GS(ct,ci,cm,co)
cost_EI_GS<-cost(D_op_EI_GS,ci,cm,co)

D_op_EI_BS<-D_EI_BS(ct,ci,cm,co)
cost_EI_BS<-cost(D_op_EI_BS,ci,cm,co)

opsum_D<-matrix(c(unlist(D_op_GS),cost_GS,
                  unlist(D_op_BS),cost_BS),
                ncol=5,byrow=TRUE,
                dimnames = list(c("D-grid search", "D-algorithm 1"),
                                c("obj", "n", "m","T","cost")))
opsum_D_EI<-matrix(c(unlist(D_op_EI_GS),cost_EI_GS,
                     unlist(D_op_EI_BS),cost_EI_BS),
                   ncol=5,byrow=TRUE,
                   dimnames = list(c("D-equal interval grid search","D-algorithm 2"),
                                   c("obj", "n", "m","T","cost")))
opsum_D[,1]<-1/opsum_D[,1]
opsum_D_EI[,1]<-1/opsum_D_EI[,1]

#V-optimality
V_op_GS<-V_GS(ct,ci,cm,co)
cost_GS<-cost(D_op_GS,ci,cm,co)

V_op_BS<-V_BS(ct,ci,cm,co)
cost_BS<-cost(D_op_BS,ci,cm,co)

V_op_EI_GS<-V_EI_GS(ct,ci,cm,co)
cost_EI_GS<-cost(D_op_EI_GS,ci,cm,co)

V_op_EI_BS<-V_EI_BS(ct,ci,cm,co)
cost_EI_BS<-cost(V_op_EI_BS,ci,cm,co)

opsum_V<-matrix(c(unlist(V_op_GS),cost_GS,
                  unlist(V_op_BS),cost_BS),
                ncol=5,byrow=TRUE,
                dimnames = list(c("V-grid search", "V-algorithm 2"),
                                c("obj", "n", "m","T","cost")))
opsum_V_EI<-matrix(c(unlist(V_op_EI_GS),cost_EI_GS,
                     unlist(V_op_EI_BS),cost_EI_BS),
          ncol=5,byrow=TRUE,
          dimnames = list(c("V-equal interval grid search","V-algorithm 2"),
                c("obj", "n", "m","T","cost")))

opsum<-rbind(opsum_A,opsum_D,opsum_V)
opsum_EI<-rbind(opsum_A_EI,opsum_D_EI,opsum_V_EI)

#relative efficiency
A_eff<-opsum[2,1]/opsum_EI[2,1]*100
D_eff<-opsum[4,1]/opsum_EI[4,1]*100
V_eff<-opsum[6,1]/opsum_EI[6,1]*100

A_eff2<-opsum[2,1]/A_ob*100
D_eff2<-opsum[4,1]/D_ob*100
V_eff2<-opsum[6,1]/V_ob*100

effsum<-c(D_eff,A_eff,V_eff)
effsum2<-c(D_eff2,A_eff2,V_eff2)

print(opsum)
print(opsum_EI)
print(c(D_ob,A_ob,V_ob))
print(c(effsum,effsum2))
#run above code for table 1 ####




#sensitivity analysis####
FImatrix=function(par){
  n=par[1]
  m=par[2]
  ti=par[3]
  re=n*matrix(c((ti^2*trigamma(alpha*ti)*m-m*ti/alpha),0,0,alpha*m*ti),2,2)
  return(re)
}
estsd=sqrt(solve(FImatrix(c(12,5,50))))
is<-c(-2,-1,0,1,2)
xi_ps<-data.frame()
h1s<-data.frame()
h2s<-data.frame()
for (i in 1:5) { #calculate h1 and h2
  alpha<-as.numeric(est$estimate[1])/deltat
  alpha<-alpha+estsd[1,1]*is[i]
  for (j in 1:5) {
    ga<--log(as.numeric(est$estimate[2])/(as.numeric(est$estimate[1])/deltat))
    ga<-ga+estsd[2,2]*is[j]
    beta<-exp(-ga)*alpha
    
    f<-function(xi_p){incgam(beta*w0, alpha*xi_p)/gamma(alpha*xi_p)-0.05}
    q<-uniroot(f,c(1,1000))
    xi_ps[j,i]<-q$root
    xi_p<-xi_ps[j,i]
    f=expression(alpha^(alpha*xi_p)*x^(alpha*xi_p-1)*exp(-alpha*x*exp(-ga)-alpha*ga*xi_p)/gamma(alpha*xi_p))
    f_xi_p<-function(x) {eval(D(f,'xi_p'))}
    f_alpha<-function(x){eval(D(f,'alpha'))}
    f_ga<-function(x){eval(D(f,'ga'))}
    a<-integrate(f_xi_p,lower=w0,upper=Inf)
    b<-integrate(f_alpha,w0,Inf)
    c<-integrate(f_ga,w0,Inf)
    h1s[j,i]<-b$value/a$value
    h2s[j,i]<-c$value/a$value
  }
}

ob_A=function(par){
  n=par[1]
  m=par[2]
  Ti=par[3]
  re=(1/((m-1)*tmin^2*trigamma(alpha*tmin)+(Ti-(m-1)*tmin)^2*trigamma(alpha*(Ti-(m-1)*tmin))-Ti/alpha)+1/(alpha*Ti))/(n/np)
  return(re)
}

ob_D=function(par){
  n=par[1]
  m=par[2]
  Ti=par[3]
  re=((n/npd)^2*(alpha*(Ti/npd)*((m-1)*tmin^2*trigamma(alpha*tmin)/npd+((Ti-(m-1)*tmin)/sqrt(npd))^2*trigamma(alpha*(Ti-(m-1)*tmin)))-(Ti/npd)^2))
  return(re)
}

ob_V=function(par){
  n=par[1]
  m=par[2]
  Ti=par[3]
  re=(h1^2/((m-1)*tmin^2*trigamma(alpha*tmin)+(Ti-(m-1)*tmin)^2*trigamma(alpha*(Ti-(m-1)*tmin))-Ti/alpha)+h2^2/(alpha*Ti))/(n/np)
  return(re)
}

concl<-matrix(0, ncol = 5)
for (i in 1:5) {
  alpha<-as.numeric(est$estimate[1])/deltat
  alpha<-alpha+estsd[1,1]*is[i]
  for (j in 1:5) {
    ga<--log(as.numeric(est$estimate[2])/(as.numeric(est$estimate[1])/deltat))
    ga<-ga+estsd[2,2]*is[j]
    beta<-exp(-ga)*alpha
    xi_p<-xi_ps[j,i]
    h1<-h1s[j,i]
    h2<-h2s[j,i]
    
    A_op_BS<-A_BS(ct,ci,cm,co)
    A_op_BS_e<-ob_A(opsum_A[2,2:4])
    D_op_BS<-D_BS(ct,ci,cm,co)
    D_op_BS_e<-ob_D(opsum_D[2,2:4])
    V_op_BS<-V_BS(ct,ci,cm,co)
    V_op_BS_e<-ob_V(opsum_V[2,2:4])
    opsum<-matrix(c(is[i],is[j],
                    A_op_BS$objective/A_op_BS_e,
                    (1/D_op_BS$objective)/(1/D_op_BS_e),
                    V_op_BS$objective/V_op_BS_e),ncol=5,byrow=TRUE)
    concl<-rbind(concl,opsum)
  }
}
D_sensitive=round((concl[seq(2,26,by=5),4]),3)
A_sensitive=round((concl[seq(2,26,by=5),3]),3)
V_sensitive=round(matrix(c(concl[-1,5]),5,5),3)
print(D_sensitive)
print(A_sensitive)
print(V_sensitive)
#run above code for table 2 ####
