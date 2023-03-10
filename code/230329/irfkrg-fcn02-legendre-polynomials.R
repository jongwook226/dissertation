# irfkrg-fcn02-legendre-polynomials.R

## Legendre Polynomial Functions

P.0 <- function(t) 1
P.1 <- function(t) t
P.2 <- function(t) (1/2)*(3*t^2 - 1)
P.3 <- function(t) (1/2)*(5*t^3 - 3*t)
P.4 <- function(t) (1/8)*(35*t^4 - 30*t^2 + 3)
P.5 <- function(t) (1/8)*(63*t^5 - 70*t^3 + 15*t)
P.6 <- function(t) (1/16)*(231*t^6 - 315*t^4 + 105*t^2 - 5)
P.7 <- function(t) (1/16)*(429*t^7 - 693*t^5 + 315*t^3 - 35*t)
P.8 <- function(t) (1/128)*(6435*t^8 - 12012*t^6 + 6930*t^4 - 1260*t^2 + 35)