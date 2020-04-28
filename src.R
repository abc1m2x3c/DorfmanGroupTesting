
####################################################################
############functions###############################################
####################################################################

#given true_status (can be vector),Se,Sp, generate a diagnosis result
test<-function(true_status,Se,Sp){
  p = ifelse(true_status==1,Se,1-Sp)
  return (rbinom(length(p),1,p))
}

#########function to take symbolic first derivative of function "func" with respects to "vars=p,Se,Sp"
get.gradient<-function(func,vars){ #func is the expression of the formula, vars are the variables that are going to be taken derivative
  funcD <- sapply(vars, function(v) D(func, v))
  return (funcD)
}

#########function to take symbolic second derivative of function "func" with respects to "vars=p,Se,Sp"
get.hessian<-function(func,vars){ #func is the expression of the formula, vars are the variables that are going to be taken derivative
  funcD <- sapply(vars, function(v) D(func, v))
  funcDD <- matrix(list(), 3,3)
  for (i in 1:length(vars)) 
    funcDD[,i] <- sapply(vars, function(v) D(funcD[[i]], v))
  return (funcDD)
}

#generate Dorfman testing data
Dorfman.generate_data<-function(theta_true,grp_num,grp_size){
  p_true=theta_true[1]
  Se_true=theta_true[2]
  Sp_true=theta_true[3]
  c=grp_size
  n=grp_num
  N=n*c
  Y_ind_true=rbinom(N,1,p_true) #individual true binary status 
  #group testing stage
  Y_grp_true=ifelse(colSums(matrix(Y_ind_true,ncol=n,nrow=c))>0,1,0) #group binary true status
  Y_grp=test(Y_grp_true,Se_true,Sp_true)
  #second stage
  Z=rep(0,sum(Y_grp)) #Z records the number of positive individual diagnosis among those positive groups
  count=0
  for (j in 1:n){
    if (Y_grp[j] ==1){ #if the jth pool's diagnosis is positive
      count=count+1
      Z[count]=sum(test(Y_ind_true[((j-1)*c+1):(j*c)],Se_true,Sp_true)) #Z indicates the number of postive diagnosis after retesting each individual
    }
    #for each group, could be zero
  }
  #y = sum(Y_grp)
  return (Z)
}


#the log likelihood function for Dorfman tsting dataset with pool size of 'grpsize'
#theta is a vector representing (p, Se, Sp)
#grpsize is the pool size
#data is Dorfman testing data which is a list of 
#1. total number of pools
#2. number of postive pools
#3. a vector of all possible values of x, where x denotes the positive individuals in a pool in the
#   individual testing stage
#4. the frequency of pools corrsponding to each x
Dorfman.llf.each_grpsize <- function(theta,grpsize,data){
  n=data[[1]]
  y=data[[2]]
  Z_category=data[[3]]
  Z_freq=data[[4]]
  
  c=grpsize
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  if (y==0){
    return ((n-y)*log(1-(Se-(Se+Sp-1)*(1-p)^c)))
  }else{
    return ((n-y)*log(1-(Se-(Se+Sp-1)*(1-p)^c)) + sum(Z_freq*(log((1-Sp-Se)*(1-p)^c*(1-Sp)^Z_category*Sp^(c-Z_category)+Se*(p*Se+(1-p)*(1-Sp))^Z_category*((1-p)*Sp+p*(1-Se))^(c-Z_category)))))
  }
}


#score function
Dorfman.Score.each_grpsize <-function(theta,grpsize,data,Dorfman_loglikelihood_D_part1,Dorfman_loglikelihood_D_part2){
  n=data[[1]]
  y=data[[2]]
  Z_category=data[[3]]
  Z_freq=data[[4]]
  c=grpsize
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  gradient_part1=(n-y)*unlist(lapply(Dorfman_loglikelihood_D_part1, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)))
  gradient_part2=
    colSums(
      matrix(Z_freq,nrow=length(Z_category),ncol=length(theta))* 
        matrix(unlist(lapply(Dorfman_loglikelihood_D_part2, eval, env=list(Se=Se,Sp=Sp,p=p,c=c,Z_category=Z_category))),nrow=length(Z_category),ncol=length(theta))
    )
  if (y==0){
    return (gradient_part1)
  }else{
    return (gradient_part1+gradient_part2)
  }
}


#dorfman information matrix
#function to calculate the information matrix of Dorfman, does not depend on simulated data.
Dorfman.IM.each_grpsize <-function(theta,grp_size,data,Dorfman_loglikelihood_DD_part1,Dorfman_loglikelihood_DD_part2){ 
  n=data[[1]]
  y=data[[2]]
  Z_category=data[[3]]
  Z_freq=data[[4]]
  
  
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]  
  c=grp_size
  
  hessian1=matrix(sapply(Dorfman_loglikelihood_DD_part1, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)), length(theta))

  
  hessian2<-function(z){
    return (matrix(sapply(Dorfman_loglikelihood_DD_part2, eval, env=list(Z_category=z,Se=Se,Sp=Sp,p=p,c=c)), length(theta)))
  }
  
  #P(G=1)=(1-p)^c*(1-Sp)+(1-(1-p)^c)*Se
  P_G1=(1-p)^c*(1-Sp)+(1-(1-p)^c)*Se #Pr(G=1)
  
  #P(Z=z,G=1)=(1-Sp-Se)*(1-p)^c*(1-Sp)^Z*Sp^(c-Z)+Se*(p*Se+(1-p)*(1-Sp))^Z*((1-p)*Sp+p*(1-Se))^(c-Z)
  P_ZzG1_rmchoose<-function(z){ #Pr(G=1,Z=z) remove choose(c,seq_z)
    return ((1-Sp-Se)*(1-p)^c*(1-Sp)^z*Sp^(c-z)+Se*(p*Se+(1-p)*(1-Sp))^z*((1-p)*Sp+p*(1-Se))^(c-z))
  }
  
  hessian2_PZzG1_sum<-function(seq_z){
    res=0
    for (i in 1:length(seq_z)){
      res=res+hessian2(seq_z[i])*choose(c,seq_z[i])*P_ZzG1_rmchoose(seq_z[i])
    }
    return (res)
  }
  
  #hessian1*E(n-y)+E(y)*E(hessian2(Z)|Y=1)
  #hessian matrix expectation
  if (y==0){
    hessian_dorfman = hessian1*n*(1-P_G1)
  }else{
    hessian_dorfman = hessian1*n*(1-P_G1)+n*hessian2_PZzG1_sum(seq(0,c))
  }
  return (-hessian_dorfman)
}


##############################below, the input arguments of functions is not complete
##############################functions depend on outer variables such as res_mle, etc.
#score confidence region
Dorfman.S.CI<-function(theta){
  IM = Dorfman.IM.all(theta)
  Score = Dorfman.Score.all(theta)
  return (t(Score)%*%solve(IM)%*%Score-qchisq(0.95,3))
}


#theta=c(p,Se,Sp)
Dorfman.LR.CI<-function(theta){ #likelihood ratio confidence interval
  return (-2*(Dorfman.llf.all(theta)-res_mle$maximum)-qchisq(0.95,3))
}



#wald confidence interval of (p,Se,Sp)
Dorfman.W.CI<-function(theta){ #likelihood ratio confidence interval
  return (t(theta-res_mle$estimate)%*%Dorfman_IM_mle%*%(theta-res_mle$estimate)-qchisq(0.95,3))
}




#Wald confidence interval 
Dorfman.profileWald.CI<-function(theta){ #likelihood ratio confidence interval
  #return ((p-res_mle$estimate[1])^2*Dorfman_IM_mle[1,1]-qchisq(0.95,1))
  return (abs(theta-res_mle$estimate)/sqrt(diag(-solve(res_mle$hessian)))-qnorm(0.975))
}



#theta=c(p,Se,Sp)
Dorfman.profileLR.CI<-function(theta){ #likelihood ratio confidence interval
  llk.maximum=rep(0,3)
  for (i in 1:3){
    param=theta[i]
    max.func<-function(other_param){
      theta=append(other_param,param,i-1)
      return (Dorfman.llf.all(theta))
    }
    A=matrix(nrow=6,ncol=3,c( 1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1))
    B=c(0,1,-0.5,1,-0.5,1)
    A=A[-((i*2-1):(i*2)),-i]
    B=B[-((i*2-1):(i*2))]
    res=maxNM(max.func,start=res_mle$estimate[-i],constraints=list(ineqA=A, ineqB=B),control=list(printLevel=0))
    llk.maximum[i]=res$maximum
  }
  return (-2*(llk.maximum-res_mle$maximum)-qchisq(0.95,1))
}


#bisection function designed for profile CI of p, Se, and Sp
bisect.forProfile<-function(func_name,x0,x1,tol){ 
  fun=match.fun(func_name)
  count<-1
  y0=fun(x0)
  y1=fun(x1)
  if(sign(y0)==sign(y1)){
    #print("Starting vaules are not suitable")
    return(NA)
  }else
  {
    while(abs(x1-x0)>tol){
      x2<-(x0+x1)/2
      
      y2<-fun(x2)
      
      if(sign(y1)==sign(y2)){x1<-x2}else{x0<-x2}
      y0<-fun(x0)
      y1<-fun(x1)
      count<-count+1
      #cat('The ',count,'th iteration:',x2,"\n") #Print the value obtained in each iteration next line 
      
    }
    return (x2)
  }
  
}



Dorfman.profileLR.CI.EP<-function(){ #likelihood ratio confidence interval #EP=end points
  EP=matrix(nrow=3,ncol=2)
  for (i in 1:3){
    profle.func<-function(param){
      max.func<-function(other_param){
        theta=append(other_param,param,i-1)
        return (Dorfman.llf.all(theta))
      }
      A=matrix(nrow=6,ncol=3,c( 1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1))
      B=c(0,1,-0.5,1,-0.5,1)
      A=A[-((i*2-1):(i*2)),-i]
      B=B[-((i*2-1):(i*2))]
      res=maxNM(max.func,start=res_mle$estimate[-i],constraints=list(ineqA=A, ineqB=B),control=list(printLevel=0))
      return(2*(res_mle$maximum-res$maximum)-qchisq(0.95,1))
    }
    EP[i,1]=bisect.forProfile(profle.func,x0=res_mle$estimate[i],x1=0,tol=1e-6)
    if (is.na(EP[i,1])){
      EP[i,1]=0
    }
    EP[i,2]=bisect.forProfile(profle.func,x0=res_mle$estimate[i],x1=1,tol=1e-6)
    if (is.na(EP[i,2])){
      EP[i,2]=1
    }    
  }
  return (EP)
}
##################################################4/18/2020 above
#calculate the characteristic corresponding to expr (e.g. ET, EC, PPV, NPV) at mle of p,Se and Sp
get_characteristic_at_mle <- function(expr){
  p  = theta_mle[1]
  Se = theta_mle[2]
  Sp = theta_mle[3]
  gradient = get.gradient(expr,vars)
  gradient=unlist(lapply(gradient, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)))
  
  mu = unlist(lapply(expr, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)))
  sigma = sqrt(t(gradient)%*%solve(Dorfman_IM_mle)%*%gradient)
  res = c(mu-1.96*sigma,mu+1.96*sigma)
  return (res)
}

#calculate the true characteristic corresponding to expr (e.g. ET, EC, PPV, NPV)
get_characteristic_true <- function(expr){
  res = unlist(lapply(expr, eval, env=list(Se=Se_true,Sp=Sp_true,p=p_true,c=c)))
  return (res)
}  


#calculate the characteristic corresponding to expr (e.g. ET, EC, PPV, NPV) at mle of p,Se and Sp
get_characteristic_at_mle_realdata <- function(expr,c){
  p  = theta_mle[1]
  Se = theta_mle[2]
  Sp = theta_mle[3]
  gradient = get.gradient(expr,vars)
  gradient=unlist(lapply(gradient, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)))
  
  mu = unlist(lapply(expr, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)))
  sigma = sqrt(t(gradient)%*%solve(Dorfman_IM_mle)%*%gradient)
  res = c(mu-1.96*sigma,mu+1.96*sigma)
  return (res)
}






#func_name can be W_CI, LR_CI, S_CI
#start, end =c(p,Se,sp), is the x,y,z coordinate of the point
bisect<-function(func_name,start_point,end_point,origin,tol){ 
  FUN=match.fun(func_name)
  start_value=FUN(start_point)
  end_value=FUN(end_point)
  #exception throw
  # if (sign(start_value) == sign(end_value)){
  #   stop ('Function of two ends should have different sign')
  # }
  # count=1
  #this if prevents the start_point is already out of the sphere
  
  
  
  
  while (sqrt(sum((start_point-end_point)^2))>tol){
    mid_point=(start_point+end_point)/2
    mid_value=FUN(mid_point)
    if ((start_value > 0) | (is.na(start_value))){
      end_value = start_value
      end_point = start_point
      start_point = origin
      start_value = FUN(origin)
    }else if ((sign(start_value)<=0) & ( (sign(mid_value)>0) | (is.na(mid_value)) )){ #if start_value<0 and mid_value>0
      end_point = mid_point
      end_value = mid_value
    }else if((sign(mid_value)<=0) & ( (sign(end_value)>0) | (is.na(end_value)) ) ){
      start_point = mid_point
      start_value = mid_value
    } else{ #(sign(end_value)<=0)
      start_value = end_value
      start_point = end_point
      end_point = origin + 4*(end_point-origin)
      end_value = FUN(end_point)
    }
  }
  return (mid_point)
  
}


##3D plot
my.plot3D<-function(func_name,UL_limit,len,main_name=''){
  FUN=Dorfman.W.CI
  FUN=match.fun(func_name)
  record_array=array(0,c(len,len,len))
  L_p=UL_limit[1]
  U_p=UL_limit[2]
  L_Se=UL_limit[3]
  U_Se=UL_limit[4]
  L_Sp=UL_limit[5]
  U_Sp=UL_limit[6]
  
  seq_Se = seq(L_Se,U_Se,length=len)
  seq_Sp = seq(L_Sp,U_Sp,length=len)
  seq_p  = seq(L_p,U_p,length=len)
  count_Se=0
  for (Se in seq_Se){
    count_Se=count_Se+1
    count_Sp=0
    for (Sp in seq_Sp){
      count_Sp=count_Sp+1
      count_p=0
      for (p in seq_p){
        count_p=count_p+1
        record_array[count_Se,count_Sp,count_p]=FUN(c(p,Se,Sp))
      }
      #uniroot(LR_CI,Se=Se,Sp=Sp,c(res$estimate[1],U_p))$root
    }
  }
  record_array[record_array>=0]=0
  record_array[record_array<0]=1
  
  M=melt(record_array)
  M=M[M$value==1,]
  plot3d(x=seq_Se[M$X1],y=seq_Sp[M$X2],z=seq_p[M$X3],col=rgb(red=1,green=0,blue=0),
        xlim=c(seq_Se[1],seq_Se[len]), ylim=c(seq_Sp[1],seq_Sp[len]),
        zlim=c(seq_p[1],seq_p[len]),xlab='Se',ylab='Sp',zlab='p',main=main_name)
  return (length(M$value)/len^3*(U_p-L_p)*(U_Se-L_Se)*(U_Sp-L_Sp))
}


#test the numeric singularity of a matrix
test.singular <- function(m) 
  class(try(solve(m),silent=T))=="matrix"


#this function returns the best group size described JRSS-B paper (Masterpool)
Masterpool.get_intermediate_Doptimality<-function(theta,xL,xU){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  c=Se/(Se+Sp-1)/(1-p)^xL
  delta=(1-Se)/Se
  r=(1-p)^(xU-xL)
  Delta0=r*log(r)/(1-r)
  asdf<-function(a){
    2/a*(1+(1+Delta0/a)/(log(a)-Delta0*(1/a-1)))-1/(delta*c+a)+1/(c-a)
  }
  a=uniroot(asdf,c(r+0.0001,0.9999),tol=0.000001)$root
  return (round(xL+log(a)/log(1-p)))
}

#
Dorfman.equivalent.grpnums.Masterpool<-function(theta,Masterpool.n,Dorfman_c){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  return (floor(Masterpool.n/(1+((1-p)^Dorfman_c*(1-Sp)+(1-(1-p)^Dorfman_c)*Se)*Dorfman_c)))
}



Dorfman.lessfavorable.equivalent.grpnums.Masterpool<-function(theta,Masterpool.n,Dorfman_c){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  c=Dorfman_c
  n=Masterpool.n
  
  pi_x = (1-p)^c*(1-Sp-Se)+Se
  d=3
  delta = sqrt(d^2*c^2*pi_x*(1-pi_x)+4*n*(1+pi_x*c))
  root= (-d*c*sqrt(pi_x*(1-pi_x))+delta)/2/(1+pi_x*c)
  if (root < 0) warning('negative n')
  return (floor(root^2))
  
  
  
}




#this function returns the best group size for Dorfman testing when keep the same expected number of tests as masterpool testing
#criteria could be 'D' or 'Ds'
Dorfman.best_grpsize<-function(theta,Masterpool.n,criteria){
  record_hessian_Dorfman=rep(0,xU-1)
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  count=0
  for (c in 2:xU){
    count=count+1
    #n=round(n_Masterpool/(1+((1-p)^c*(1-Sp)+(1-(1-p)^c)*Se)*c))
    n=Dorfman.equivalent.grpnums.Masterpool(theta,Masterpool.n,c)
    if (criteria == 'D'){
      record_hessian_Dorfman[count]=det(Dorfman.IM.new(theta,n,c))
    }else if (criteria == 'Ds'){
      record_hessian_Dorfman[count]=-solve(Dorfman.IM.new(theta,n,c))[1]
    }
  }
  c=(2:xU)[which(record_hessian_Dorfman==max(record_hessian_Dorfman))]
  return (c)
}


#this function returns the best group size for Dorfman testing when keep the same expected number of tests as masterpool testing
#only change to the above function: n=Dorfman.equivalent.grpnums.Masterpool(theta,Masterpool.n,c) 
#                                    -> n=Dorfman.lessfavourable.equivalent.grpnums.Masterpool(theta,Masterpool.n,c)
#criteria could be 'D' or 'Ds'
Dorfman.lessfavorable.best_grpsize<-function(theta,Masterpool.n,criteria){
  record_hessian_Dorfman=rep(0,xU-1)
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  count=0
  for (c in 2:xU){
    count=count+1
    #n=round(n_Masterpool/(1+((1-p)^c*(1-Sp)+(1-(1-p)^c)*Se)*c))
    n=Dorfman.lessfavorable.equivalent.grpnums.Masterpool(theta,Masterpool.n,c)
    if (criteria == 'D'){
      record_hessian_Dorfman[count]=det(Dorfman.IM.new(theta,n,c))
    }else if (criteria == 'Ds'){
      record_hessian_Dorfman[count]=-solve(Dorfman.IM.new(theta,n,c))[1]
    }
  }
  c=(2:xU)[which(record_hessian_Dorfman==max(record_hessian_Dorfman))]
  return (c)
}






#JRSS-B paper requires 3 different group sizes. THEORETCAL
Masterpool.IM.new<-function(theta,Masterpool.grp_sizes,Masterpool.grp_nums){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  
  xL=Masterpool.grp_sizes[1]
  c_intermediate=Masterpool.grp_sizes[2]
  xU=Masterpool.grp_sizes[3]
  
  hessian_masterpool=0
  
  count=0
  for (c in c(xL,c_intermediate,xU)){
    count=count+1
    hessian1=matrix(sapply(Masterpool_loglikelihood_DD_part1, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)), length(vars))
    hessian2=matrix(sapply(Masterpool_loglikelihood_DD_part2, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)), length(vars))
    hessian_masterpool = hessian_masterpool + hessian1*(Masterpool.grp_nums[count])*(1-(Se-(Se+Sp-1)*(1-p)^c))+hessian2*(Masterpool.grp_nums[count])*(Se-(Se+Sp-1)*(1-p)^c)
  }
  return (-hessian_masterpool)
}


#JRSS-B paper requires 3 different group sizes. THEORETCAL
Masterpool.score.new<-function(theta,Masterpool.grp_sizes,Masterpool.grp_nums){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  
  xL=Masterpool.grp_sizes[1]
  c_intermediate=Masterpool.grp_sizes[2]
  xU=Masterpool.grp_sizes[3]
  
  score_masterpool=0
  
  count=0
  for (c in c(xL,c_intermediate,xU)){
    count=count+1
    score1=matrix(sapply(Masterpool_loglikelihood_D_part1, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)), length(vars))
    score2=matrix(sapply(Masterpool_loglikelihood_D_part2, eval, env=list(Se=Se,Sp=Sp,p=p,c=c)), length(vars))
    score_masterpool = score_masterpool + score1*(Masterpool.grp_nums[count]-Masterpool.pos_grps[count])+score2*Masterpool.pos_grps[count]
  }
  return (score_masterpool)
}





Masterpool.generate_data<-function(theta_true,grp_sizes,grp_nums){
  p_true  = theta_true[1]
  Se_true = theta_true[2]
  Sp_true = theta_true[3]
  pos_grps=rep(0,3) #record the number of postive groups
  for (i in 1:3){
    c=grp_sizes[i]
    n=grp_nums[i]
    N=n*c
    Y_ind_true=rbinom(N,1,p_true) #individual true binary status 
    #group testing stage
    Y_grp_true=ifelse(colSums(matrix(Y_ind_true,ncol=n,nrow=c))>0,1,0) #group binary true status
    Y_grp=test(Y_grp_true,Se_true,Sp_true)
    pos_grps[i]=sum(Y_grp)
  }
  return (pos_grps) 
}




Masterpool.llf.original<-function(theta,grp_nums,grp_sizes,pos_grps){
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  c1=grp_sizes[1]
  c2=grp_sizes[2]
  c3=grp_sizes[3]
  y1=pos_grps[1]
  y2=pos_grps[2]
  y3=pos_grps[3]
  n1=grp_nums[1]
  n2=grp_nums[2]
  n3=grp_nums[3]
  return( y1*log((1-p)^c1*(1-Sp)+(1-(1-p)^c1)*Se)+(n1-y1)*log(1-(1-p)^c1*(1-Sp)-(1-(1-p)^c1)*Se)
          +y2*log((1-p)^c2*(1-Sp)+(1-(1-p)^c2)*Se)+(n2-y2)*log(1-(1-p)^c2*(1-Sp)-(1-(1-p)^c2)*Se)
          +y3*log((1-p)^c3*(1-Sp)+(1-(1-p)^c3)*Se)+(n3-y3)*log(1-(1-p)^c3*(1-Sp)-(1-(1-p)^c3)*Se) )
  
}





#theta=c(p,Se,Sp)
Masterpool.LR.CI<-function(theta){ #likelihood ratio confidence interval
  return (-2*(Masterpool.llf.new(theta)-res_mle$maximum)-qchisq(0.95,3))
}



#wald confidence interval
Masterpool.W.CI<-function(theta){ #likelihood ratio confidence interval
  return (t(theta-res_mle$estimate)%*%Masterpool_IM_mle%*%(theta-res_mle$estimate)-qchisq(0.95,3))
}

#score confidence interval
Masterpool.S.CI<-function(theta){
  IM = Masterpool.IM.new(theta,Masterpool.grp_sizes,Masterpool.grp_nums)
  Score = Masterpool.score.new(theta,Masterpool.grp_sizes,Masterpool.grp_nums)
  return (t(Score)%*%solve(IM)%*%Score-qchisq(0.95,3))
}




Estep1 <- function(theta){
  
  ###############################
  ############# asdf
  ###############################
  
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  
  
  res = (
    #for n-y negative pool, i.e., Y_i=0:
    (n-y)*
    ( 
      c*(E_tildeZ_ij_given_Yi0*log(p) + (1-E_tildeZ_ij_given_Yi0)*log(1-p)) 
      + E_IsigmatildeZ_ij0_given_Yi0*log(Sp) + (1-E_IsigmatildeZ_ij0_given_Yi0)*log(1-Se) 
    )
  #
  +
    #for y positive pool, i.e., Y_i=1, there are Z_category positive individuals:
    sum(
      Z_freq*
        (
          # if Z_ij = 1
          (Z_category) * (E_tildeZ_ij_given_Yi1Zij1*log(p) + (1-E_tildeZ_ij_given_Yi1Zij1)*log(1-p)) +
            # if Z_ij=0
            (c-Z_category) * (E_tildeZ_ij_given_Yi1Zij0*log(p) + (1-E_tildeZ_ij_given_Yi1Zij0)*log(1-p)) +
            
            E_I_Zij0_given_Yi1*log(1-Sp) + (1-E_I_Zij0_given_Yi1)*log(Se) +
            
            #if Z_ij = 1
            (Z_category) * (E_tildeZ_ij_given_Yi1Zij1*log(Se) + (1-E_tildeZ_ij_given_Yi1Zij1)*log(1-Sp)) + 
            #if Z_ij = 0
            (c-Z_category) *  (E_tildeZ_ij_given_Yi1Zij0*log(1-Se) + (1-E_tildeZ_ij_given_Yi1Zij0)*log(Sp))
        )
    )
  )
  return (res)
  
}

Estep2 <- function(theta){
  
  p=theta[1]
  Se=theta[2]
  Sp=theta[3]
  
  
  res=sum(
    rowSums(E_Z_mat_tilde)*log(p)+(c-rowSums(E_Z_mat_tilde))*log(1-p) +
      E_I0_tilde*(1-Y_grp)*log(Sp) + (1-E_I0_tilde)*Y_grp*log(Se) + E_I0_tilde*Y_grp*log(1-Sp) +(1-E_I0_tilde)*(1-Y_grp)*log(1-Se) +
      rowSums((1-E_Z_mat_tilde)*(1-Z_mat)*log(Sp)+E_Z_mat_tilde*Z_mat*log(Se)+(1-E_Z_mat_tilde)*Z_mat*log(1-Sp)+E_Z_mat_tilde*(1-Z_mat)*log(1-Se))*Y_grp
  )
  return(res)
}


interative.input <-function(){
  raw_data=list()
  flag='N'
  c_choice=c()
  while (flag != 'Y'){
    cat('##################################\n','######Please input your data######\n','##################################',sep='')
    c <- readline(prompt='Enter pool size: ')
    c <- as.numeric(unlist(strsplit(c, ",")))
    
    n <- readline(prompt='Enter total number of pools: ')
    n <- as.numeric(unlist(strsplit(n, ",")))
    
    y <- readline(prompt='Enter number of postive pools: ')
    y <- as.numeric(unlist(strsplit(y, ",")))
    if (y!=0){
      Z <- readline(prompt='Enter number of positive individuals in each postive pool: ')
      Z <- as.numeric(unlist(strsplit(Z, ",")))
    }else{
      Z <- c()
    }
    if (y!=length(Z)){
      print('Error. Start again.')
      next
    }
    if (c %in% c_choice){
      start.index=(which(c_choice==c)-1)*5+1
      raw_data[[start.index+1]]=raw_data[[start.index+1]]+n
      raw_data[[start.index+2]]=raw_data[[start.index+2]]+y
      Z=c(raw_data[[start.index+3]],Z)
      Z_category = as.numeric(levels(as.data.frame(table(Z))$Z))
      Z_freq     = as.numeric(as.data.frame(table(Z))$Freq)
      raw_data[[start.index+3]]=Z
      raw_data[[start.index+4]]=Z_category
      raw_data[[start.index+5]]=Z_freq
    }else{
      c_choice=c(c_choice,c)
      Z_category = as.numeric(levels(as.data.frame(table(Z))$Z))
      Z_freq     = as.numeric(as.data.frame(table(Z))$Freq)
      raw_data=append(raw_data,list(c,n,y,Z,Z_category,Z_freq))
    }
    flag <- readline(prompt='Is this the end of your data (Y/N)? ')
  }
  #delete the Z in the raw_data
  raw_data=raw_data[-(seq(4,length(raw_data),6))]
  return (raw_data)
}

############################move functions from main into src 20200427
Dorfman.llf.all<-function(theta){
  num_data=length(my.data)/5
  res=0
  for (i in 1:num_data){
    res=res+Dorfman.llf.each_grpsize(theta,my.data[[(5*(i-1)+1)]],my.data[(5*(i-1)+2):(5*(i-1)+5)])
  }
  return (res)
}




Dorfman.Score.all<-function(theta){
  num_data=length(my.data)/5
  res=0
  for (i in 1:num_data){
    res=res+Dorfman.Score.each_grpsize(theta,my.data[[(5*(i-1)+1)]],my.data[(5*(i-1)+2):(5*(i-1)+5)],Dorfman_loglikelihood_D_part1,Dorfman_loglikelihood_D_part2)
  }
  return (res)
}


Dorfman.IM.all<-function(theta){
  num_data=length(my.data)/5
  res=0
  for (i in 1:num_data){
    res=res+Dorfman.IM.each_grpsize(theta,my.data[[(5*(i-1)+1)]],my.data[(5*(i-1)+2):(5*(i-1)+5)],Dorfman_loglikelihood_DD_part1,Dorfman_loglikelihood_DD_part2)
  }
  return (res)
}


DT.summary<-function(my.data){
  ##calculate mle and related variables
  #A and B construct inequality for 0<p<1, 0.5<Se,Sp<1
  A=matrix(nrow=6,ncol=3,c( 1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1))
  B=c(0,1,-0.5,1,-0.5,1)
  res_mle=maxNM(Dorfman.llf.all,start=c(0.1,0.9,0.9),constraints=list(ineqA=A, ineqB=B),control=list(printLevel=0))
  theta_mle = res_mle$estimate
  theta_mle
  
  Dorfman_IM_mle=Dorfman.IM.all(theta_mle) #information matrix of all data
  
  
  Sigma=sqrt(diag(solve(Dorfman_IM_mle))) #the standard deviation calculated from the information matrix
  
  #Wald CI for each of p, Se, and Sp
  Wald_CI=
    rbind(
      c(res_mle$estimate[1]-qnorm(0.975)*Sigma[1], res_mle$estimate[1]+qnorm(0.975)*Sigma[1]),
      c(res_mle$estimate[2]-qnorm(0.975)*Sigma[2], res_mle$estimate[2]+qnorm(0.975)*Sigma[2]),
      c(res_mle$estimate[3]-qnorm(0.975)*Sigma[3], res_mle$estimate[3]+qnorm(0.975)*Sigma[3])
    )
  row.names(Wald_CI)=c('p','Se','Sp')
  
  
  #Profile CI for each of p, Se, and Sp
  Profile_CI=Dorfman.profileLR.CI.EP()
  row.names(Profile_CI)=c('p','Se','Sp')
  Profile_CI
  
  #summary the output with respect to p, Se, and Sp
  output_pSeSp=cbind(theta_mle,Sigma,Wald_CI,Profile_CI)
  colnames(output_pSeSp)=c('est','std err','95% Wald L','95% Wald U','95% Profile L','95% Profile U')
  return (output_pSeSp)
  
}


DT.plot<-function(my.data,type='Wald'){
  ##calculate mle and related variables
  #A and B construct inequality for 0<p<1, 0.5<Se,Sp<1
  A=matrix(nrow=6,ncol=3,c( 1,-1,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1,-1))
  B=c(0,1,-0.5,1,-0.5,1)
  res_mle=maxNM(Dorfman.llf.all,start=c(0.1,0.9,0.9),constraints=list(ineqA=A, ineqB=B),control=list(printLevel=0))
  theta_mle = res_mle$estimate
  theta_mle
  
  Dorfman_IM_mle=Dorfman.IM.all(theta_mle) #information matrix of all data
  
  
  Sigma=sqrt(diag(solve(Dorfman_IM_mle))) #the standard deviation calculated from the information matrix
  
  #Approach I: plot 
  len=100 # number of points on one dimension of the cube
  
  if (type=='Wald'){
    # #Wald Confidence Region
    L = res_mle$estimate-3*Sigma
    U = res_mle$estimate+3*Sigma
    UL_limit = as.vector(matrix(rbind(L,U),nrow=2,ncol=3))
    my.plot3D(Dorfman.W.CI,UL_limit,len,main_name = 'Joint Wald Region')
  }
  
  if (type=='LR'){
    # #LR
    L = res_mle$estimate-3*Sigma
    U = res_mle$estimate+3*Sigma
    UL_limit = as.vector(matrix(rbind(L,U),nrow=2,ncol=3))
    UL_limit = ifelse(UL_limit>1,1,UL_limit) #truncate the plot if any of p, Se, and Sp is larger than 1.
    my.plot3D(Dorfman.LR.CI,UL_limit,len,main_name = 'Joint LR Region' )
  }
  # if (type='Score'){
  #   # # #Score too slow
  #   # tt=proc.time()
  #   # L = res_mle$estimate-3*Sigma
  #   # U = res_mle$estimate+3*Sigma
  #   # UL_limit = as.vector(matrix(rbind(L,U),nrow=2,ncol=3))
  #   # UL_limit = ifelse(UL_limit>1,1,UL_limit)
  #   # my.plot3D(Dorfman.S.CI,UL_limit,len)
  #   # proc.time()-tt
  # }
}

DT.chr<-function(my.data){
  
  #expressions of characteristic in Dorfman testing
  EC_expression  = expression((1-p)^c*(1-Sp)*(Se+Sp-1) + (1-p)*(1-Se+Se*Sp-Se^2) + Se^2)
  ET_expression  = expression(Se-(1-p)^c*(Se+Sp-1)+1/c)
  PPV_expression = expression(p*Se^2/((1-p)*(1-Sp)*(Se+(1-Se-Sp)*(1-p)^(c-1))+p*Se^2))
  NPV_expression = expression((1-p)*(1-(1-Sp)*(Se+(1-Se-Sp)*(1-p)^(c-1)))/(p*(1-Se^2)+(1-p)*(1-(1-Sp)*(Se+(1-Se-Sp)*(1-p)^(c-1)))))
  
  
  
  ##calculate the Wald CI's of characteristics using delta method
  final_res=list()
  for (i in 1:(length(my.data)/5)){
    grp_size=my.data[[sort(seq(1,length(my.data),5))[i]]]
    character_interval_ET = get_characteristic_at_mle_realdata(ET_expression,grp_size)
    character_interval_EC = get_characteristic_at_mle_realdata(EC_expression,grp_size)
    character_interval_PPV = get_characteristic_at_mle_realdata(PPV_expression,grp_size)
    character_interval_NPV = get_characteristic_at_mle_realdata(NPV_expression,grp_size)
    res=rbind(character_interval_ET,character_interval_EC,character_interval_PPV,character_interval_NPV)
    row.names(res)=c('E(T)','E(C)','PPV','NPV')
    colnames(res)=c('L','U')
    cat('The Wald CI of group testing characteristics of pool size ',grp_size,' is :\n',sep='')
    print(res)
    cat('\n')
    final_res[[i]]=rbind(final_res,res)
    names(final_res)[i]=paste('pool size ',grp_size,sep='')
  }
  return(final_res)
}


###################################################################
#symbolic expressions that can speed up the calculation
###################################################################
#expression of Dorfman likelihood, diff likelihood and diffdiff likelihood
vars=c('p','Se','Sp')
Dorfman_loglikelihood_part1 = expression(log(1-(Se-(Se+Sp-1)*(1-p)^c)))
Dorfman_loglikelihood_part2 = expression(log((1-Sp-Se)*(1-p)^c*(1-Sp)^Z_category*Sp^(c-Z_category)+Se*(p*Se+(1-p)*(1-Sp))^Z_category*((1-p)*Sp+p*(1-Se))^(c-Z_category)))
Dorfman_loglikelihood_D_part1 = get.gradient(Dorfman_loglikelihood_part1,vars)
Dorfman_loglikelihood_D_part2 = get.gradient(Dorfman_loglikelihood_part2,vars)
Dorfman_loglikelihood_DD_part1 = get.hessian(Dorfman_loglikelihood_part1,vars)
Dorfman_loglikelihood_DD_part2 = get.hessian(Dorfman_loglikelihood_part2,vars)





