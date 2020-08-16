## Consider a mixture of two normal distributions with equal weights (w1 = w2 = 1/2)
## Component 1 has mean 0 and standard deviation 1
## Component 2 has mean 1 and standard deviation 1
## The observation is x = 50
## What is Pr(c = 1 | x)?
dnorm(50, 0, 1)
dnorm(50, 1, 1)
dnorm(50, 0, 1)/(dnorm(50, 0, 1) + dnorm(50, 1, 1))


## What if x=3?  Two ways to do the calculation
## One way:  Direct calculation
z1 = dnorm(3, 0, 1)
z2 = dnorm(3, 1, 1)
z1/(z1+z2)
## A second way:  Compute in the logarithm scale, add b 
## to all values, and then exponentiate before standardizing
lz1 = dnorm(3, 0, 1, log=T)
lz2 = dnorm(3, 1, 1, log=T)
b = 3
exp(lz1+b)/(exp(lz1+b) + exp(lz2+b))

## Going back to the case x - 50:
## Wrong
lz1 = log(dnorm(50, 0, 1))
lz2 = log(dnorm(50, 1, 1))
b = max(lz1, lz2)
exp(lz1-b)/(exp(lz1-b) + exp(lz2-b))
## Wrong
lz1 = log(exp(-0.5*50^2)/sqrt(2*pi))
lz2 = log(exp(-0.5*49^2)/sqrt(2*pi))
b = max(lz1, lz2)
exp(lz1-b)/(exp(lz1-b) + exp(lz2-b))
## Right
lz1 = dnorm(50, 0, 1, log=TRUE)
lz2 = dnorm(50, 1, 1, log=TRUE)
b = max(lz1, lz2)
exp(lz1-b)/(exp(lz1-b) + exp(lz2-b))
## Also right (just more cumbersome)
lz1 = -0.5*log(2*pi) - 0.5*50^2
lz2 = -0.5*log(2*pi) - 0.5*49^2
b = max(lz1, lz2)
exp(lz1-b)/(exp(lz1-b) + exp(lz2-b))