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
