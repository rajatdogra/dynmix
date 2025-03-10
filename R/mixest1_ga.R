### GA-enhanced mixture estimation with optimized component selection
### Tryig to reduces computational load by evolving subsets of models for dynmix

mixest1_ga <- function(y, x, 
                      k_models = 10,          # Number of models per GA individual
                      population_size = 20,    # GA population size
                      max_generations = 50,    # Max GA generations
                      mutation_rate = 0.1,     # Mutation probability per bit
                      crossover_rate = 0.7,    # Crossover probability
                      tournament_size = 3,     # Tournament selection size
                      convergence_threshold = 1e-5, # Convergence criterion
                      ftype = NULL,            # Original dynmix parameters
                      lambda = NULL,
                      kappa = NULL,
                      V = NULL,
                      W = NULL,
                      atype = NULL) {
  
  # Ensure x is a matrix and get dimensions
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  colnames(x) <- colnames(x, do.NULL = FALSE, prefix = "X")
  
  # Initialize population of model subsets
  initialize_population <- function() {
    lapply(1:population_size, function(i) {
      matrix(sample(0:1, k_models * p, replace = TRUE), 
             nrow = k_models, ncol = p)
    })
  }
  
  # Crossover operator: uniform row crossover
  crossover <- function(parent1, parent2) {
    crossover_rows <- sample(1:k_models, floor(k_models/2))
    child <- parent1
    child[crossover_rows, ] <- parent2[crossover_rows, ]
    return(child)
  }
  
  # Mutation operator: bit-flip with repair
  mutate <- function(individual) {
    mut_mask <- matrix(runif(k_models * p) < mutation_rate, 
                       nrow = k_models, ncol = p)
    individual[mut_mask] <- 1 - individual[mut_mask]
    # Ensure at least one active variable per model
    apply(individual, 1, function(row) {
      if (sum(row) == 0) row[sample(p, 1)] <- 1
      return(row)
    }) |> t()
  }
  
  # Fitness evaluation: MSE of dynmix with current models
  evaluate_fitness <- function(individual) {
    mods <- cbind(1, individual)  # Add intercept
    colnames(mods) <- c("const", colnames(x))
    tryCatch({
      model <- mixest1(y = y, x = x, mods = mods, ftype = ftype,
                      lambda = lambda, kappa = kappa, V = V, W = W, atype = atype)
      mse <- mean((y - model$y.hat)^2, na.rm = TRUE)
      return(list(mse = mse, model = model))
    }, error = function(e) return(list(mse = Inf, model = NULL)))
  }
  
  # Initialize population and trackers
  population <- initialize_population()
  best_mse <- Inf
  best_model <- NULL
  convergence_count <- 0
  fitness_history <- numeric(max_generations)
  
  # Main GA loop
  for (gen in 1:max_generations) {
    # Evaluate population
    evaluations <- lapply(population, evaluate_fitness)
    mses <- sapply(evaluations, function(e) e$mse)
    valid <- which(mses < Inf)
    
    # Update best solution
    current_best <- which.min(mses)
    if (length(current_best) > 0 && mses[current_best] < best_mse) {
      best_mse <- mses[current_best]
      best_model <- evaluations[[current_best]]$model
      best_model$ga_mods <- cbind(1, population[[current_best]])
      convergence_count <- 0
    } else {
      convergence_count <- convergence_count + 1
    }
    fitness_history[gen] <- best_mse
    
    # Check convergence
    if (convergence_count >= 5 && gen > 10) {
      if (sd(fitness_history[(gen-5):gen]) < convergence_threshold) break
    }
    
    # Selection: tournament selection on valid individuals
    parents <- lapply(1:population_size, function(i) {
      contestants <- sample(valid, tournament_size, replace = TRUE)
      population[[contestants[which.min(mses[contestants])]]]
    })
    
    # Reproduction: crossover and mutation
    new_population <- lapply(1:population_size, function(i) {
      if (runif(1) < crossover_rate && i < population_size) {
        child <- crossover(parents[[i]], parents[[sample(valid, 1)]])
      } else {
        child <- parents[[i]]
      }
      mutate(child)
    })
    population <- new_population
  }
  
  # Attach GA metadata to best model
  if (!is.null(best_model)) {
    best_model$ga_info <- list(
      best_mse = best_mse,
      generations = gen,
      fitness_history = fitness_history[1:gen],
      parameters = list(
        k_models = k_models,
        population_size = population_size,
        mutation_rate = mutation_rate
      )
    )
    class(best_model) <- c("mixest_ga", class(best_model))
  }
  
  return(best_model %||% stop("GA failed to find valid model"))
}

# S3 method for summarizing GA-enhanced model
summary.mixest_ga <- function(object) {
  cat("GA-enhanced Dynamic Mixture Model\n")
  cat("Best MSE:", object$ga_info$best_mse, "\n")
  cat("Selected models:", nrow(object$ga_mods), "\n")
  print(summary.mixest(object))
}