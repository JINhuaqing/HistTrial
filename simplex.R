euclidean_proj_simplex <- function(v, s=1){
# """ Compute the Euclidean projection on a positive simplex

# Solves the optimisation problem (using the algorithm from [1]):

#     min_w 0.5 * || w - v ||_2^2 , s.t. \sum_i w_i = s, w_i >= 0 

# Parameters
# ----------
# v: (n,) numpy array,
#    n-dimensional vector to project

# s: int, optional, default: 1,
#    radius of the simplex

# Returns
# -------
# w: (n,) numpy array,
#    Euclidean projection of v on the simplex

# Notes
# -----
# The complexity of this algorithm is in O(n log(n)) as it involves sorting v.
# Better alternatives exist for high-dimensional sparse vectors (cf. [1])
# However, this implementation still easily scales to millions of dimensions.

# References
# ----------
# [1] Efficient Projections onto the .1-Ball for Learning in High Dimensions
#     John Duchi, Shai Shalev-Shwartz, Yoram Singer, and Tushar Chandra.
#     International Conference on Machine Learning (ICML 2008)
#     http://www.cs.berkeley.edu/~jduchi/projects/DuchiSiShCh08.pdf
# """
    
    n <- length(v)
    if ((sum(v) == s) & (sum(v>=0)==n)){
        return(v)
    }
    
    
    u <- sort(v, decreasing=TRUE)
    cssv <- cumsum(u)
    rhos <-  (u * (1:n)) > (cssv - s)
    rho <- which(rhos)
         
    if (length(rho)==0){
        return(0)
    }
    
    rho <- rho[length(rho)]
    theta <- (cssv[rho]-s)/rho
    w <- v-theta
    w[w<0] <- 0
    return(w)
 
}


euclidean_proj_l1ball <- function(v, s=1){
#    """ Compute the Euclidean projection on a L1-ball
#
#    Solves the optimisation problem (using the algorithm from [1]):
#
#        min_w 0.5 * || w - v ||_2^2 , s.t. || w ||_1 <= s
#
#    Parameters
#    ----------
#    v: (n,) numpy array,
#       n-dimensional vector to project
#
#    s: int, optional, default: 1,
#       radius of the L1-ball
#
#    Returns
#    -------
#    w: (n,) numpy array,
#       Euclidean projection of v on the L1-ball of radius s
#
#    Notes
#    -----
#    Solves the problem by a reduction to the positive simplex case
#
#    See also
#    --------
#    euclidean_proj_simplex
#    """
    n <- length(v)
    u <- abs(v)
    if (sum(u) < s){
        return(v)
    }
    
    w <- euclidean_proj_simplex(u, s=s)
    
    w <- w*sign(v)
    
    return(w)
}

#w <- euclidean_proj_l1ball(c(-1, 0.2, 0.3, 0.4 ), s=10);w
#sum(w)
