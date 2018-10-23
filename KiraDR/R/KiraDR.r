# PCA ---------------------------------------------------------------------
my_pca <- function(data, cor = F)
{
  data <- scale(data,center = T,scale = F) 
  if(cor == F)
  {
    sigma <- cov(data)
    eigenvector.sigma <-eigen(sigma)$vectors
    pca.sigma <- as.matrix(data) %*% eigenvector.sigma
    return(pca.sigma)
  }
  if(cor == T)
  {
    error <- F
    error <- tryCatch(r <- cor(data),
                      warning = function(e)
                      {
                        message("Warning:The standard deviation is zero, and the correlation coefficient matrix can not be calculated, so use the covariance matrix to solve the principal component.")
                        return(T)
                      })
    if(error==F)
    {    
      eigenvector.r <- eigen(r)$vectors
      pca.r <- as.matrix(data) %*% eigenvector.r
      return(pca.r)
    }
    else
    {
      sigma <- cov(data)
      eigenvector.sigma <-eigen(sigma)$vectors
      pca.sigma <- as.matrix(data) %*% eigenvector.sigma
      return(pca.sigma)
    }
  }
}

my_screeplot <- function(pca , n, title = 'Screeplot')
{
  data_var <- apply(pca[,1:n],2,var)
  data_plot <- data.frame(comp = c(1:n), var = data_var)
  ggplot(data_plot,aes(x = comp, y = var, group = 1)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    theme(axis.title.x = element_text(size = 25, face = "bold.italic", vjust = 0.5, hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 25, face = "bold.italic", vjust = 0.5, hjust = 0.5)) + 
    theme(axis.text.x = element_text(size = 15, vjust = 0.5, hjust = 0.5)) + 
    theme(axis.text.y = element_text(size = 15, vjust = 0.5, hjust = 0.5)) + 
    scale_x_continuous(breaks=seq(0, n, n%/%5)) + 
    ggtitle(title) + theme(plot.title=element_text(size = 30))
}



# EFA ---------------------------------------------------------------------
Kmo.test <- function(corr)
{
  ifsolve <- T
  ifsolve <- tryCatch({r <- solve(corr);ifsolve <- T},
                      error = function(e){
                        massage("Matrix is not invertible.")
                        return(F)
                      })
  if(ifsolve == T)
  {
    r <- cov2cor(r)
    diag(corr) <- 0
    diag(r) <- 0
    sumr2 <- sum(r^2)
    sumcorr2 <- sum(corr^2)
    Kmo <- sumcorr2/(sumcorr2 + sumr2)
    cat('Kmo:',Kmo,sep = '')
    return(Kmo)
  }
  else
    return()
}

Bartlett.test <- function (corr, n) 
{
  p <- nrow(corr)
  detr <- det(corr)
  statistic <- -log(detr) * (n - 1 - (2 * p + 5)/6)
  df <- p * (p - 1)/2
  pvalue <- pchisq(statistic, df, lower.tail = FALSE)
  bartlett <- list(chisq = statistic, p.value = pvalue, df = df)
  return(bartlett)
}

my_fa_pa <- function(corr, nfactor)
{
  p <- nrow(corr)
  diag_corr <- diag(corr)
  sum_rank <- sum(diag_corr)
  rowname <- colnames(corr)
  colname <- paste("PA", 1:nfactor, sep="")
  loading <- matrix(0, nrow = p, ncol = nfactor, dimnames = list(rowname, colname))
  k <- 1
  kmax <- 50 # maximum number of iterations
  if(qr(corr)$rank == p) # iterative initial value
    sigma <- diag(1/solve(corr))
  else
    sigma <- 1 - apply(abs(corr - diag(1, nrow = p)),2, max)
  h <- diag_corr - sigma
  repeat
  {
    diag(corr) <- h
    h1 <- h
    eig <- eigen(corr)
    for (i in 1:nfactor)
      loading[,i] <- sqrt(eig$values[i]) * eig$vectors[,i]
    h <- diag(loading %*% t(loading))
    if ((sqrt(sum((h - h1)^2)) < 0.0001)|k == kmax) # reach steady state
      break()
    k <- k+1
  }
  rowname <- c("SS loadings", "Proportion Var", "Cumulative Var")
  var_sum <- matrix(0, nrow = 3, ncol = nfactor, dimnames = list(rowname, colname))
  for (i in 1:nfactor)
  {
    var_sum[1,i]<-sum(loading[,i]^2)
    var_sum[2,i]<-var_sum[1,i]/sum_rank
    var_sum[3,i]<-sum(var_sum[1,1:i])/sum_rank
  }
  method <- c("Principal Factor Method")
  list(method = method, loadings = loading, var_sum = var_sum, sp = diag_corr - h, iterative = k) 
}

my_fa_ml <- function(corr, nfactor)
{
  p <- nrow(corr)
  diag_corr <- diag(corr)
  sum_rank <- sum(diag_corr)
  rowname <- colnames(corr)
  colname <- paste("ML", 1:nfactor, sep="")
  loading <- matrix(0, nrow = p, ncol = nfactor, dimnames = list(rowname, colname))
  k <- 1
  kmax <- 50 # maximum number of iterations
  if(qr(corr)$rank == p) # iterative initial value
    sigma <- diag(1/solve(corr))
  else
    sigma <- 1 - apply(abs(corr - diag(1, nrow = p)),2, max)
  repeat
  {
    d1 <- sigma
    d2 <- 1/sqrt(sigma)
    eig <- eigen(corr * (d2 %o% d2))
    for (i in 1:nfactor)
      loading[,i] <- sqrt(eig$values[i] - 1) * eig$vectors[,i]
    loading <- diag(sqrt(sigma)) %*% loading
    sigma <- diag(corr - loading%*%t(loading))
    if ((sqrt(sum((sigma - d1)^2)) < 0.0001)|k == kmax)
      break()
    k <- k+1
  }
  rowname <- c("SS loadings","Proportion Var","Cumulative Var")
  var_sum <- matrix(0, nrow=3, ncol=nfactor, dimnames=list(rowname, colname))
  for (i in 1:nfactor)
  {
    var_sum[1,i] <- sum(loading[,i]^2)
    var_sum[2,i] <- var_sum[1,i]/sum_rank
    var_sum[3,i] <- sum(var_sum[1,1:i])/sum_rank
  }
  method<-c("Maximum Likelihood Method")
  list(method = method, loadings = loading, var_sum = var_sum, sp = sigma, iterative = k)
}

