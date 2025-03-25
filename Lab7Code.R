################################################################################
# Lab 7 Code
# Avery Johnson
################################################################################

library(tidyverse)

################################################################################
# Task 1: Describe the Population Distribution
################################################################################


population.distribution <- function(alpha, beta) {
  beta.dist <- data.frame(
    alpha = alpha,
    beta = beta,
    mean = alpha / (alpha + beta),
    var = (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1)),
    skew = (2 * (beta - alpha) * sqrt(alpha + beta + 1)) / 
      ((alpha + beta + 2) * sqrt(alpha * beta)),
    kurtosis = (6*((alpha-beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2))) /
      (alpha*beta*(alpha+beta+2)*(alpha+beta+3))
  )
  return(beta.dist)
}

population.table <- bind_rows(
  population.distribution(2,5),
  population.distribution(5,5),
  population.distribution(5,2),
  population.distribution(0.50, 0.50)
)

View(population.table)

# plot these four distributions 
fig.1.data <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, 2, 5))                          # compute the beta PDF

plot1 <- ggplot(data= fig.1.data)+                                           # specify data
  geom_line(aes(x=x, y=beta.pdf, color="Beta(2,5)")) +          # plot beta dist
  geom_hline(yintercept=0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "grey"))+                 # change colors
  theme(legend.position = "bottom")                                    # move legend to bottom


fig.2.data <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   
  mutate(beta.pdf = dbeta(x, 5, 5))                          
         
plot2 <- ggplot(data= fig.2.data)+                                          
  geom_line(aes(x=x, y=beta.pdf, color="Beta(5,5)")) +          
  geom_hline(yintercept=0)+                                            
  theme_bw()+                                                          
  xlab("x")+                                                          
  ylab("Density")+                                                     
  scale_color_manual("", values = c("black", "grey"))+                
  theme(legend.position = "bottom")

fig.3.data <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   
  mutate(beta.pdf = dbeta(x, 5, 2))                          

plot3 <- ggplot(data= fig.3.data)+                                          
  geom_line(aes(x=x, y=beta.pdf, color="Beta(5,2)")) +          
  geom_hline(yintercept=0)+                                            
  theme_bw()+                                                          
  xlab("x")+                                                          
  ylab("Density")+                                                     
  scale_color_manual("", values = c("black", "grey"))+                
  theme(legend.position = "bottom") 

fig.4.data <- tibble(x = seq(-0.25, 1.25, length.out=1000))|>   
  mutate(beta.pdf = dbeta(x, 0.5, 0.5))                          

plot4 <- ggplot(data= fig.4.data)+                                          
  geom_line(aes(x=x, y=beta.pdf, color="Beta(0.5,0.5)")) +          
  geom_hline(yintercept=0)+                                            
  theme_bw()+                                                          
  xlab("x")+                                                          
  ylab("Density")+                                                     
  scale_color_manual("", values = c("black", "grey"))+                
  theme(legend.position = "bottom") 

library(patchwork)
plot1 + plot2 + plot3 + plot4

################################################################################
# Task 2: Compute the Moments
################################################################################

beta.moment <- function(alpha, beta, k, centered){
  mean <- alpha / (alpha + beta)
  if(centered==TRUE) {
    centered.moment <- integrate(function(x) ( (x-mean)^k * dbeta(x, alpha, beta) ), 0, 1)$value
    return(centered.moment)
  } else{
    uncentered.moment <- integrate(function(x) ((x^k) * dbeta(x, alpha, beta)), 0, 1)$value
    return(uncentered.moment)
  }
}

# Test sample values
(ex_mean <- beta.moment(2, 5, 1, centered=FALSE))
(ex_var <- beta.moment(2, 5, 2, centered=TRUE))
(ex_skew <- (beta.moment(2, 5, 3, centered=TRUE)) / ((ex_var)^(3/2)))
(ex_kurt <- ((beta.moment(2, 5, 4, centered=TRUE)) / ((ex_var)^2)) - 3)

(ex_mean2 <- beta.moment(5, 5, 1, centered=FALSE))
(ex_var2 <- beta.moment(5, 5, 2, centered=TRUE))
(ex_skew2 <- (beta.moment(5, 5, 3, centered=TRUE)) / ((ex_var2)^(3/2)))
(ex_kurt2 <- ((beta.moment(5, 5, 4, centered=TRUE)) / ((ex_var2)^2)) - 3)

(ex_mean3 <- beta.moment(5, 2, 1, centered=FALSE))
(ex_var3 <- beta.moment(5, 2, 2, centered=TRUE))
(ex_skew3 <- (beta.moment(5, 2, 3, centered=TRUE)) / ((ex_var3)^(3/2)))
(ex_kurt3 <- ((beta.moment(5, 2, 4, centered=TRUE)) / ((ex_var3)^2)) - 3)

(ex_mean4 <- beta.moment(0.5, 0.5, 1, centered=FALSE))
(ex_var4 <- beta.moment(0.5, 0.5, 2, centered=TRUE))
(ex_skew4 <- (beta.moment(0.5, 0.5, 3, centered=TRUE)) / ((ex_var4)^(3/2)))
(ex_kurt4 <- ((beta.moment(0.5, 0.5, 4, centered=TRUE)) / ((ex_var4)^2)) - 3)


################################################################################
# Task 3: Do Data Summaries Help?
################################################################################

# for beta(2,5) dist
library(e1071)

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 2
beta <- 5
beta.sample1 <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter

beta.sample1 <- tibble(beta.sample1)

sample.stats1 <- beta.sample1 |>
  summarize(
    alpha=alpha,
    beta=beta,
    mean = mean(beta.sample1),
    variance = var(beta.sample1),
    skewness = e1071::skewness(beta.sample1),
    excess_kurt = e1071::kurtosis(beta.sample1)
  )

hist1 <- ggplot() +
  geom_histogram(data=beta.sample1,
                 aes(x=beta.sample1, y=after_stat(density)),
                 breaks=seq(0,1,0.1),
                 fill="grey30",
                 color="lightgrey")+
  stat_density(data=beta.sample1,
               aes(x=beta.sample1, color="Density"), geom="line") +
  geom_line(data=fig.1.data,
            aes(x=x, y=beta.pdf, color="Population Beta Dist")) +
  geom_hline(yintercept = 0) +
  theme_bw()+
  xlab("x")+
  ylab("Density") +
  labs(color="Line") +
  ggtitle("Beta (2,5) Distribution")

# for beta(5,5) dist

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 5
beta <- 5
beta.sample2 <- rbeta(n = sample.size,  # sample size
                      shape1 = alpha,   # alpha parameter
                      shape2 = beta)    # beta parameter

beta.sample2 <- tibble(beta.sample2)

sample.stats2 <- beta.sample2 |>
  summarize(
    alpha=alpha,
    beta=beta,
    mean = mean(beta.sample2),
    variance = var(beta.sample2),
    skewness = e1071::skewness(beta.sample2),
    excess_kurt = e1071::kurtosis(beta.sample2)
  )

hist2 <- ggplot() +
  geom_histogram(data=beta.sample2,
                 aes(x=beta.sample2, y=after_stat(density)),
                 breaks=seq(0,1,0.1),
                 fill="grey30",
                 color="lightgrey")+
  stat_density(data=beta.sample2,
               aes(x=beta.sample2, color="Density"), geom="line") +
  geom_line(data=fig.2.data,
            aes(x=x, y=beta.pdf, color="Population Beta Dist")) +
  geom_hline(yintercept = 0) +
  theme_bw()+
  xlab("x")+
  ylab("Density") +
  labs(color="Line") +
  ggtitle("Beta (5,5) Distribution")

# for beta (5,2) distribution

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 5
beta <- 2
beta.sample3 <- rbeta(n = sample.size,  # sample size
                      shape1 = alpha,   # alpha parameter
                      shape2 = beta)    # beta parameter

beta.sample3 <- tibble(beta.sample3)

sample.stats3 <- beta.sample3 |>
  summarize(
    alpha=alpha,
    beta=beta,
    mean = mean(beta.sample3),
    variance = var(beta.sample3),
    skewness = e1071::skewness(beta.sample3),
    excess_kurt = e1071::kurtosis(beta.sample3)
  )

hist3 <- ggplot() +
  geom_histogram(data=beta.sample3,
                 aes(x=beta.sample3, y=after_stat(density)),
                 breaks=seq(0,1,0.1),
                 fill="grey30",
                 color="lightgrey")+
  stat_density(data=beta.sample3,
               aes(x=beta.sample3, color="Density"), geom="line") +
  geom_line(data=fig.3.data,
            aes(x=x, y=beta.pdf, color="Population Beta Dist")) +
  geom_hline(yintercept = 0) +
  theme_bw()+
  xlab("x")+
  ylab("Density") +
  labs(color="Line") +
  ggtitle("Beta (5,2) Distribution")

# for beta(0.5, 0.5 dist)

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 0.5
beta <- 0.5
beta.sample4 <- rbeta(n = sample.size,  # sample size
                      shape1 = alpha,   # alpha parameter
                      shape2 = beta)    # beta parameter

beta.sample4 <- tibble(beta.sample4)

sample.stats4 <- beta.sample4 |>
  summarize(
    alpha=alpha,
    beta=beta,
    mean = mean(beta.sample4),
    variance = var(beta.sample4),
    skewness = e1071::skewness(beta.sample4),
    excess_kurt = e1071::kurtosis(beta.sample4)
  )

hist4 <- ggplot() +
  geom_histogram(data=beta.sample4,
                 aes(x=beta.sample4, y=after_stat(density)),
                 breaks=seq(0,1,0.1),
                 fill="grey30",
                 color="lightgrey")+
  stat_density(data=beta.sample4,
               aes(x=beta.sample4, color="Density"), geom="line") +
  geom_line(data=fig.4.data,
            aes(x=x, y=beta.pdf, color="Population Beta Dist")) +
  geom_hline(yintercept = 0) +
  theme_bw()+
  xlab("x")+
  ylab("Density") +
  labs(color="Line") +
  ggtitle("Beta (0.5,0.5) Distribution")

hist1 + hist2 + hist3 + hist4 + plot_layout(guides = 'collect')

sample.stats <- bind_rows(
  sample.stats1,
  sample.stats2,
  sample.stats3,
  sample.stats4
)

################################################################################
# Task 4: Is Sample Size Important?
################################################################################

# for the beta(2,5) dist.
library("cumstats")

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details
alpha <- 2
beta <- 5
beta.sample1 <- rbeta(n = sample.size,  # sample size
                      shape1 = alpha,   # alpha parameter
                      shape2 = beta)    # beta parameter

cumstats.data <- data.frame(
  alpha=2,
  beta=5,
  mean = cummean(beta.sample1),
  var = cumvar(beta.sample1),
  skew = cumskew(beta.sample1),
  kurtosis = cumkurt(beta.sample1)
)

cumstats.data <- cumstats.data |>
  mutate(observation = 1:n())

pop.mean <- alpha / (alpha + beta)
pop.var <- (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
pop.skew <- (2 * (beta - alpha) * sqrt(alpha + beta + 1)) / 
  ((alpha + beta + 2) * sqrt(alpha * beta))
pop.kurt <- (6*((alpha-beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2))) /
  (alpha*beta*(alpha+beta+2)*(alpha+beta+3))
  


cumstats.mean <- ggplot(data=cumstats.data) +
  geom_line(aes(x=observation, y=mean, color="Cumulative Mean"), show.legend=F) +
  geom_hline(yintercept=pop.mean)+  
  xlab("observation") +
  ylab("mean") +
  ggtitle("Cumulative Statistics Mean") +
  theme_bw()

cumstats.var <- ggplot(data=cumstats.data) +
  geom_line(aes(x=observation, y=var, color="Cumulative Variance"), show.legend=F) +
  geom_hline(yintercept = pop.var) +  
  xlab("observation") +
  ylab("variance") +
  ggtitle("Cumulative Statistics Variance") +
  theme_bw()

cumstats.skew <- ggplot(data=cumstats.data) +
  geom_line(aes(x=observation, y=skew, color="Cumulative Skewness"), show.legend=F) +
  geom_hline(yintercept = pop.skew) +
  xlab("observation") +
  ylab("skewness") +
  ggtitle("Cumulative Statistics Skewness") +
  theme_bw()

cumstats.kurt <- ggplot(data=cumstats.data) +
  geom_line(aes(x=observation, y=kurtosis-3, color="Cumulative Kurtosis"), show.legend=F) +
  geom_hline(yintercept = pop.kurt) + 
  xlab("observation") +
  ylab("skewness") +
  ggtitle("Cumulative Statistics Kurtosis") +
  theme_bw()

library(patchwork)
cumulative.stats <-  cumstats.mean + cumstats.var + cumstats.skew + cumstats.kurt

# now do this in a for loop for 2:50
for (i in 2:50){
  set.seed(7272 + i)
  sample.size <- 500 # Specify sample details
  alpha <- 2
  beta <- 5
  beta.sample <- rbeta(n = sample.size,  # sample size
                        shape1 = alpha,   # alpha parameter
                        shape2 = beta)    # beta parameter
  
  new.cumstats.data <- data.frame(
    alpha=2,
    beta=5,
    mean = cummean(beta.sample),
    var = cumvar(beta.sample),
    skew = cumskew(beta.sample),
    kurtosis = cumkurt(beta.sample)-3
  )
  
  new.cumstats.data <- new.cumstats.data |>
    mutate(observation = 1:n())
  
  cumstats.mean <- cumstats.mean +
    geom_line(data=new.cumstats.data, aes(x=observation, y=mean), color=i)
  
  cumstats.var <- cumstats.var + 
    geom_line(data=new.cumstats.data, aes(x=observation, y=var), color=i)
  
  cumstats.skew <- cumstats.skew +
    geom_line(data=new.cumstats.data, aes(x=observation, y=skew), color=i)
  
  cumstats.kurt <- cumstats.kurt +
    geom_line(data=new.cumstats.data, aes(x=observation, y=kurtosis), color=i)
  
  cumulative.stats.new <- cumstats.mean + cumstats.var + cumstats.skew + cumstats.kurt
  }

cumulative.stats.new

################################################################################
# Task 5: How can we model the variation? 
################################################################################

new.stats <- data.frame() 

for (i in 1:1000){
  set.seed(7272 + i)
  sample.size <- 500 # Specify sample details
  alpha <- 2
  beta <- 5
  beta.sample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter
  stats <- data.frame(
    alpha = alpha,
    beta = beta,
    mean = mean(beta.sample),
    variance = var(beta.sample),
    skewness = e1071::skewness(beta.sample),
    excess_kurt = e1071::kurtosis(beta.sample)
  )
  
  new.stats <- bind_rows(new.stats, stats)
}

view(new.stats)

# Plot Stats

# Mean Plot
mean.plot <- ggplot() +
  geom_histogram(data = new.stats, aes(x=mean, y=after_stat(density)), 
                 bins=40,
                 fill="grey30",
                 color="lightgray")+
  geom_density(data=new.stats, aes(x=mean, color="Density")) +
  theme_bw() +
  ggtitle("Histogram of Means") +
  xlab("Mean") +
  ylab("Density") +
  labs(color="Line")

# Variance Plot
var.plot <- ggplot() +
  geom_histogram(data = new.stats, aes(x=variance, y=after_stat(density)), 
                 bins=40,
                 fill="grey30",
                 color="lightgray")+
  geom_density(data=new.stats, aes(x=variance, color="Density")) +
  theme_bw() +
  ggtitle("Histogram of Variances") +
  xlab("Variance") +
  ylab("Density") +
  labs(color="Line")

# Skew Plot
skew.plot <- ggplot() +
  geom_histogram(data = new.stats, aes(x=skewness, y=after_stat(density)), 
                 bins=40,
                 fill="grey30",
                 color="lightgray")+
  geom_density(data=new.stats, aes(x=skewness, color="Density")) +
  theme_bw() +
  ggtitle("Histogram of Skewnesses") +
  xlab("Skewness") +
  ylab("Density") +
  labs(color="Line")

# Excess Kurt Plot
kurt.plot <- ggplot() +
  geom_histogram(data = new.stats, aes(x=excess_kurt, y=after_stat(density)), 
                 bins=40,
                 fill="grey30",
                 color="lightgray")+
  geom_density(data=new.stats, aes(x=excess_kurt, color="Density")) +
  theme_bw() +
  ggtitle("Histogram of Excess Kurtosis") +
  xlab("Excess Kurtosis") +
  ylab("Density") +
  labs(color="Line")

stats.plots <- mean.plot + var.plot + skew.plot + kurt.plot + plot_layout(guides = "collect")
stats.plots

################################################################################
# Task 6: Collect and Clean Data
################################################################################
library(tidyverse)
death.rates <- read_csv("countrydeathrates.csv", skip=3)

death.data <- death.rates |>
  select("Country Name", "Country Code", `2022`) |> #keep only relevant columns
  filter(!is.na(`2022`)) |>
  mutate(`2022` = `2022` / 1000) |> #convert to rate
  rename("death rate" = `2022` )

View(death.data)  

################################################################################
# Task 7: What are alpha and beta?
################################################################################
library(nleqslv)

# MOM
MOM.beta <- function(data, par){
  alpha <- par[1]
  beta <- par[2]
  
  EX <- alpha / (alpha + beta)
  EX2 <- ((alpha + 1) * alpha) / ((alpha + beta + 1)*(alpha + beta))
  m1 <- mean(data)
  m2 <- mean(data^2)
  
  return(c((EX-m1), (EX2-m2))) # goal: find alpha and beta so these are 0
}

mom.estimates <- nleqslv(x = c(1, 1), # guess
        fn = MOM.beta,
        data=death.data$`death rate`)

mom.alpha <- mom.estimates$x[1]
mom.beta <- mom.estimates$x[2]

# MLE
llbeta <- function(data, par, neg=FALSE){
  alpha <- par[1]
  beta <- par[2]
  loglik <- sum(log(dbeta(x=data, alpha, beta)))
  
  return(ifelse(neg, -loglik, loglik))
}

mle.estimates <- optim(par = c(1,1),
      fn = llbeta,
      data=death.data$`death rate`,
      neg = T)

mle.alpha <- mle.estimates$par[1]
mle.beta <- mle.estimates$par[2]

density.data <- tibble(x=seq(0,0.023, length.out=1000)) |>
  mutate(mom.density = dbeta(x, mom.alpha, mom.beta),
         mle.density = dbeta(x, mle.alpha, mle.beta))

# histograms 
ggplot() + 
  geom_histogram(data=death.data,
           aes(x=`death rate`,
               y=after_stat(density)),
           fill="grey")+
  geom_line(data = density.data, aes(x=x, y=mom.density, color="MOM Estimate")) +
  geom_line(data = density.data, aes(x=x, y=mle.density, color="MLE Estimate")) +
  theme_bw() +
  xlab("Death Rate") +
  ylab("Density") +
  ggtitle("Death Rates with Estimated Beta Distributions") 

################################################################################
# Task 8: Which estimators should we use
################################################################################

results <- data.frame(mom.alpha = numeric(1000), 
                      mom.beta = numeric(1000), 
                      mle.alpha = numeric(1000), 
                      mle.beta = numeric(1000))

for (i in 1:1000){
  set.seed(7272 + i)
  sample.size <- 266 # Specify sample details
  alpha <- 8
  beta <- 950
  beta.sample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter

  # MOM ESTIMATES
  mom.estimates <- nleqslv(x = c(1, 1), # guess
                           fn = MOM.beta,
                           data=beta.sample)
  
  mom.alpha <- mom.estimates$x[1]
  mom.beta <- mom.estimates$x[2]
  
  # MLE ESTIMATES
  mle.estimates <- optim(par = c(1,1),
                         fn = llbeta,
                         data=beta.sample,
                         neg = T)
  
  mle.alpha <- mle.estimates$par[1]
  mle.beta <- mle.estimates$par[2]
  
  results[i, ] <- c(mom.alpha, mom.beta, mle.alpha, mle.beta)
}
view(results)

# true parameters
alpha_true <- 8
beta_true <- 950

#compute bias, precision, MLE
bias_mom_alpha <- mean(results$mom.alpha) - alpha_true
bias_mom_beta <- mean(results$mom.beta) - beta_true
bias_mle_alpha <- mean(results$mle.alpha) - alpha_true
bias_mle_beta <- mean(results$mle.beta) - beta_true

precision_mom_alpha <- 1 / var(results$mom.alpha)
precision_mom_beta <- 1 / var(results$mom.beta)
precision_mle_alpha <- 1 / var(results$mle.alpha)
precision_mle_beta <- 1 / var(results$mle.beta)

mse_mom_alpha <- var(results$mom.alpha) + bias_mom_alpha^2
mse_mom_beta <- var(results$mom.beta) + bias_mom_beta^2
mse_mle_alpha <- var(results$mle.alpha) + bias_mle_alpha^2
mse_mle_beta <- var(results$mle.beta) + bias_mle_beta^2

summary.table <- data.frame(
  Parameter = c("Alpha", "Beta"),
  Bias_MOM = c(bias_mom_alpha, bias_mom_beta),
  Bias_MLE = c(bias_mle_alpha, bias_mle_beta),
  Precision_MOM = c(precision_mom_alpha, precision_mom_beta),
  Precision_MLE = c(precision_mle_alpha, precision_mle_beta),
  MSE_MOM = c(mse_mom_alpha, mse_mom_beta),
  MSE_MLE = c(mse_mle_alpha, mse_mle_beta)
)

view(summary.table)

# create density plots
p1 <- ggplot() +
  geom_density(data=results, aes(x=mom.alpha, color="Density")) +
  ggtitle("MOM Alpha Density") +
  xlab("MOM Density") +
  ylab("Density")

p2<- ggplot() + 
  geom_density(data=results, aes(x=mom.beta, color="Density")) +
  ggtitle("MOM Beta Density") +
  xlab("MOM Density") +
  ylab("Density")

p3 <- ggplot() +
  geom_density(data=results, aes(x=mle.alpha, color="Density")) +
  ggtitle("MLE Alpha Density") +
  xlab("MLE Density") +
  ylab("Density")

p4<- ggplot() + 
  geom_density(data=results, aes(x=mle.beta, color="Density")) +
  ggtitle("MLE Beta Density") +
  xlab("MLE Density") +
  ylab("Density")

plot.grid <- p1 + p2 + p3 + p4

p.alpha <- ggplot() +
  geom_density(data=results, aes(x=mom.alpha, color="MOM")) +
  geom_density(data=results, aes(x=mle.alpha, color="MLE")) +
  theme_bw() +
  ggtitle("Alpha Density") +
  xlab("Alpha") +
  ylab("Density") +
  scale_color_manual(name = "Estimator", values = c("MOM" = "blue", "MLE" = "red"))
  
p.beta <- ggplot() +
  geom_density(data=results, aes(x=mom.beta, color="MOM")) +
  geom_density(data=results, aes(x=mle.beta, color="MLE")) +
  theme_bw() +
  ggtitle("Beta Density") +
  xlab("Beta") +
  ylab("Density") +
  scale_color_manual(name = "Estimator", values = c("MOM" = "blue", "MLE" = "red"))

p.alpha + p.beta
