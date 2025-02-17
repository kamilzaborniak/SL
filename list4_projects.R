## theoretical

#zad 3
rss<- c(1731, 730, 49, 38.9, 32, 29, 28.5, 27.8, 27.6, 26.6)
k<-1:10
p<-10
n<-100
s<-1
aic<- rss+2*k*s
bic <- rss+k*log(n)
ric<- rss+2*k*log(p)
which.min(ric)


### PROJECT 1 ####

## zad 1
load("Lab3.Rdata")
View(xx)
nam.row<- row.names(xx)
nam.col<- colnames(xx)
X<- apply(xx, 2, as.numeric)
row.names(X)<- nam.row
#class(X)
#class(xx)

X_scaled<-  apply(X, 2, function(i) (i-mean(i))/sd(i) + mean(i)) #scale(X) - srednio działa
X_sub10<- X_scaled-10

# MLE
mu.mle<- colMeans(X_sub10[1:5,])

# JS estimation - shrunk towards zero
# warianccja jest równa 1/5 bo 1/n*I   
p<- 300
c_JS <- 1 - (0.2* (p-2) / sum(mu.mle^2))
mu.JS_to0<- c_JS * mu.mle

# JS estimation - shrunk towards common mean
d_JS<- ((p-3)*0.2)/((p-1)* var(mu.mle))
mu.JS_tocm<- (1-d_JS)*mu.mle + d_JS*mean(mu.mle)


X_sub10.rem210.avg<- colMeans(X_sub10[6:210,])
df.means<- data.frame(avg=X_sub10.rem210.avg, MLE=mu.mle, JS_c=mu.JS_to0, JS_d=mu.JS_tocm)

df.means<- as.data.table(df.means) %>%  melt(id.vars="avg", variable.name= "Estimator", value.name = "Value")

ggplot(df.means, aes(x=avg, y=Value, col= Estimator))+
  geom_point(shape=1, size=2)

sum((X_sub10.rem210.avg-mu.mle)^2)
sum((X_sub10.rem210.avg-mu.JS_to0 )^2)
sum((X_sub10.rem210.avg-mu.JS_tocm)^2)

## zad 2

n<- 1000
p<- 950
sigma<- 1000^-0.5
m<-5
proj1_exc2<- function(m){
#set.seed(410)
  X<- replicate(p, rnorm(n,0,sigma))

  # beta p-wymiarowa
  #  k<- 5
  Beta<- numeric(p)
  Beta[1:5]<- 3

  #set.seed(140)
  epsilon<- rnorm(n)

  # Y= X*B + eps
  Y<- X %*% Beta + epsilon

  model<- lm(Y~X[,1:m]-1)
  Beta.est <- coef(model)

  RSS<- sum((Y- X[,1:m]%*% Beta.est)^2)
  #set.seed(14)
  #epsilon_1<- rnorm(n)
  PE.theoretical<- sum( (X[,1:m] %*% (Beta[1:m]-Beta.est))^2)+ n*1 # 2. wzór ppkt. a)
  PE.1<- RSS + 2*1*m # ppkt. b) sigma =1
  PE.2<- RSS + 2*RSS/(n-m) *m # (sigma est.= RSS/n-p)    
  M<- X[,1:m] %*% solve(t(X[,1:m]) %*% X[,1:m]) %*% t(X[,1:m])
  PE.CV <- sum( ( (Y- X[,1:m] %*% Beta.est) / (1-diag(M)) )^2 ) # leave-one-out cv

  result<-data.frame(m=m, RSS=RSS, PE.theoretical=PE.theoretical, PE_1=PE.1, PE_2=PE.2, PE_3=PE.CV)
  return( result)
}
proj1_exc2(5)
sim_k<-replicate(3, proj1_exc2(5), simplify = T) %>% t() %>% as.data.frame()
i<-30
k<- c(2,5,10,100,500,950)
df<-rbind(
  rdply(i, proj1_exc2(2))  %>% as.data.frame(),
  rdply(i, proj1_exc2(5))  %>% as.data.frame(),
  rdply(i, proj1_exc2(10)) %>% as.data.frame(),
  rdply(i, proj1_exc2(100)) %>% as.data.frame(),
  rdply(i, proj1_exc2(500))  %>% as.data.frame(),
  rdply(i, proj1_exc2(950))  %>% as.data.frame()
)
j=500
df[df$m==j,] %>% as.data.frame() -> ddd
colMeans(ddd)
dfdiff<- data.frame(m=df$m, PE_1=df$PE_1-df$PE.theoretical,
                 PE_2=df$PE_2-df$PE.theoretical,
                 PE_3=df$PE_3-df$PE.theoretical) %>% as.data.table() %>% melt(id.vars="m", variable.name= "Estimator", value.name = "Value")
dfdiff<- as.data.frame()
dfdiff$m %<>% as.factor()
  
plt1<-ggplot(dfdiff[1:510, ], aes(x=m, y=Value, col=Estimator))+
  geom_boxplot()
plt12<-ggplot(dfdiff[511:540, ], aes(x=m, y=Value, col=Estimator))+
  geom_boxplot()

