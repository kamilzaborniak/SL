k<- 5
n<- 500
p<- 450
sd <- 1/sqrt(n)

simulation<- function(k, n=500, p=450){
  sd <- 1/sqrt(n)
  beta<- numeric(p)
  epsilon <- 2* rnorm(n)
  beta[1:k]<- 10
  
  X<- replicate(p, rnorm(n, sd=1/sqrt(n)))
  Y<- X%*% beta + epsilon
  
  model_LS<- lm(Y~X-1)
  B_LS <- c(coef(model_LS))
  MSE_LS<- mean((beta-B_LS)^2)
  mu_LS<- mean((X%*%beta- X%*%B_LS)^2)
    
  # Ridge with CV
  
  model_RRcv <- cv.glmnet(X, Y, alpha=0,
                          intercept=FALSE,
                          standardize =FALSE)
  B_RRcv <- coef(model_RRcv, s="lambda.min") %>% c()
  B_RRcv <- B_RRcv[[1]]@x
  MSE_RRcv<- mean((beta-B_RRcv)^2)
  mu_RRcv<- mean((X%*%beta- X%*%B_RRcv)^2)
  # LASSO with CV
  
  model_LASSOcv<- cv.glmnet(X,Y, alpha=1,
                            intercept=FALSE,
                            standardize =FALSE)
  B_LASSOcv<- coef(model_LASSOcv, s="lambda.min")%>% c()
  id<- B_LASSOcv[[1]]@i 
  val<- B_LASSOcv[[1]]
  val<- val@x
  B_LASSOcv<- numeric(p)
  B_LASSOcv[id]<- val
  FDP_LASSO <- sum(B_LASSOcv[(k+1):p]!=0)/max(sum(B_LASSOcv!=0),1)
  power_LASSO<- sum(B_LASSOcv[1:k]!=0)/k
  MSE_LASSOcv<- mean((beta-B_LASSOcv)^2)
  mu_LASSOcv<- mean((X%*%beta- X%*%B_LASSOcv)^2)
  
  # Knockoffs with RR
  
  mu<- rep(0,p)
  Sigma <- diag(1/n,p)
  
  knockoff <- function(X) create.gaussian(X, mu, Sigma)
  foo<- stat.glmnet_coefdiff
  k_stat<- function(X, X_k, y) foo(X,
                                   X_k,
                                   y,
                                   family = "gaussian",
                                   alpha=0)
  selected_RRko<- knockoff.filter(X,
                                  Y,
                                  knockoffs = knockoff,
                                  fdr=0.2,
                                  statistic = k_stat)
  selected_RRko<- selected_RRko[["selected"]]
  FDP_RR_ko<- sum(selected_RRko>k) / max(1, sum(selected_RRko>0))
  power_RR_ko<- sum(selected_RRko<=k)/k
  # Knockoffs with LASSO
  
  k_stat<- function(X, X_k, y) foo(X,
                                   X_k,
                                   y,
                                   alpha=1)
  selected_LASSOko<- knockoff.filter(X,
                                  Y,
                                  knockoffs = knockoff,
                                  fdr=0.2,
                                  statistic = k_stat)
  selected_LASSOko<-selected_LASSOko[["selected"]]
  FDP_LASSO_ko<- sum(selected_LASSOko>k) / max(1, sum(selected_LASSOko>0))
  power_LASSO_ko<- sum(selected_LASSOko<=k)/k
  
  result<- c(FDP_LASSO, FDP_RR_ko, FDP_LASSO_ko,
             power_LASSO, power_RR_ko, power_LASSO_ko,
             MSE_LS, MSE_RRcv, MSE_LASSOcv,
             mu_LS, mu_RRcv, mu_LASSOcv)
  return(result)
  }

simulation(k=5)
replicate(100, simulation(k=5)) -> df_k5
replicate(100, simulation(k=20)) -> df_k20
replicate(100, simulation(k=50))-> df_k50

df_k5 %>% t() %>% as.data.frame() %>% colMeans()
df_k20 %>% t() %>% as.data.frame() %>% colMeans()
df_k50 %>% t() %>% as.data.frame() %>% colMeans()












