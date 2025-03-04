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
    skewness = skewness(beta.sample1),
    excess_kurt = kurtosis(beta.sample1)
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
    skewness = skewness(beta.sample2),
    excess_kurt = kurtosis(beta.sample2)
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
    skewness = skewness(beta.sample3),
    excess_kurt = kurtosis(beta.sample3)
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
    skewness = skewness(beta.sample4),
    excess_kurt = kurtosis(beta.sample4)
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


