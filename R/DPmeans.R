#' Penalized kmeans algorithm
#' Reference:
#' Brian Kulis & Michael I. Jordan, Revisiting k-means: New Algorithms via Bayesian Nonparametrics, 2012
DPmeans <- function(x, lambda, tol = 1e-3){
  k <- 1
  centers <- mean(x)
  z <- rep(1, length(x))
  f_val <- sum((x - centers)^2) + lambda * length(centers)

  changes <- 1
  while(changes > tol){
    for(i in 1:length(x)){   # O(n)
      # cat("x=", x, "\n")
      # cat("centers=", centers, "\n")
      min_d <- min((x[i] - centers)^2)      # O(k)
      # cat("min_d=", min_d, "\n")

      if(min_d > lambda){
        k <- k+1
        z[i] <- k

        # cat("k=", k, "\n")
        # cat("x[i]=", x[i], "\n")
        centers[k] <- x[i]
      } else {
        z[i] <- which.min((x[i] - centers)^2) # O(k)
      }
    }

    new_centers <- rep(Inf, k)
    group_mean <- aggregate(x, by = list(z), FUN = mean) # O(n)
    new_centers[group_mean$Group.1] <- group_mean$x
    # cat("z=", z, "\n")
    # cat("new_centers=", new_centers, "\n")

    new_f_val <- sum(sapply(x, function(xi){
      min((xi - new_centers)^2)
    })) + lambda * length(new_centers)

    changes <- f_val - new_f_val
    centers <- new_centers
    f_val <- new_f_val

    # cat("centers = ", centers, '\n')
    # cat("f value = ", f_val, '\n')
  }

  sort(unique(centers[centers < Inf]))
}
