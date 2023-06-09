% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Laguerre density functions.R
\name{laguerre_density}
\alias{laguerre_density}
\title{Computes f_{theta,theta_tilde}(z), i.e.,  the Laguerre density as specified
in the paper.
Input:
z           - Vector of points at which the density shall be evaluated
tau,theta,  - see likelihood
theta_tilde
Output: A vector of the same length as z which contains the corresponding
values of the density with quantile tau and parameters theta and
theta_tilde}
\usage{
laguerre_density(z, tau, theta, theta_tilde)
}
\description{
Computes f_{theta,theta_tilde}(z), i.e.,  the Laguerre density as specified
in the paper.
Input:
z           - Vector of points at which the density shall be evaluated
tau,theta,  - see likelihood
theta_tilde
Output: A vector of the same length as z which contains the corresponding
values of the density with quantile tau and parameters theta and
theta_tilde
}
