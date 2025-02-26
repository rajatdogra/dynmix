### GA-enhanced mixture estimation
### Extension of Nagy & Suzdaleva (2013) with genetic algorithms
### for component selection and sparsity-inducing priors

source("/Users/rajat/Documents/UW-Application/UW-Studies/research/code-repo/dynmix/IDUB/R/mixest1.R")  # Import original mixest1 function

mixest1_ga <- function(y, x, 
                      population_size = 30,
                      max_generations = 50,
                      mutation_rate = 0.1,
                      crossover_rate = 0.7,
                      tournament_size = 3,
                      convergence_threshold = 1e-6,
                      ftype = NULL,
                      lambda = NULL,
                      kappa = NULL,
                      V = NULL,
                      W = NULL,
                      atype = NULL,
                      max_vars = NULL) {
  
  # Initialize tracking variables
  best_fitness <- Inf
  best_solution <- NULL
  best_model <- NULL
  n_vars <- ncol(x)
  generation_history <- list()
  fitness_history <- numeric(max_generations)
  
  # Initialize population
  population <- initialize_population(n_vars, population_size, max_vars)
  
  # Ensure we have a valid initial solution
  initial_valid <- FALSE
  attempts <- 0
  while (!initial_valid && attempts < 10) {
    tryCatch({
      # Try to get one valid model
      mods <- genome_to_components(population[[1]], n_vars)
      model <- mixest1(y = y, x = x, mods = mods, 
                      ftype = ftype, lambda = lambda, kappa = kappa,
                      V = V, W = W, atype = atype)
      best_model <- model
      best_solution <- population[[1]]
      best_fitness <- mean((y - model$y.hat)^2, na.rm = TRUE)
      initial_valid <- TRUE
    }, error = function(e) {
      population <- initialize_population(n_vars, population_size, max_vars)
      attempts <- attempts + 1
    })
  }
  
  if (!initial_valid) {
    stop("Unable to find valid initial solution. Please check input data and parameters.")
  }
  
  # Main GA loop
  for (gen in 1:max_generations) {
    # Evaluate fitness for each individual in population
    fitness_scores <- numeric(population_size)
    models <- list()
    valid_solutions <- logical(population_size)
    
    for (i in 1:population_size) {
      mods <- genome_to_components(population[[i]], n_vars)
      
      tryCatch({
        model <- mixest1(y = y, x = x, mods = mods, 
                        ftype = ftype, lambda = lambda, kappa = kappa,
                        V = V, W = W, atype = atype)
        
        models[[i]] <- model
        fitness_scores[i] <- -mean((y - model$y.hat)^2, na.rm = TRUE)
        valid_solutions[i] <- TRUE
      }, error = function(e) {
        fitness_scores[i] <- -Inf
        models[[i]] <- NULL
        valid_solutions[i] <- FALSE
      })
    }
    
    # Update best solution if we found any valid solutions
    if (any(valid_solutions)) {
      best_gen_idx <- which.max(fitness_scores)
      if (fitness_scores[best_gen_idx] > -best_fitness) {
        best_fitness <- -fitness_scores[best_gen_idx]
        best_solution <- population[[best_gen_idx]]
        best_model <- models[[best_gen_idx]]
      }
    }
    
    # Store generation history
    generation_history[[gen]] <- list(
      population = population,
      fitness_scores = fitness_scores,
      best_fitness = best_fitness,
      valid_solutions = valid_solutions
    )
    
    fitness_history[gen] <- best_fitness
    
    # Check convergence
    if (gen > 5) {
      recent_improvement <- abs(fitness_history[gen] - fitness_history[gen-5])
      if (recent_improvement < convergence_threshold) {
        break
      }
    }
    
    # Selection and reproduction
    new_population <- list()
    
    for (i in 1:population_size) {
      # Tournament selection from valid solutions if possible
      if (any(valid_solutions)) {
        parent1_idx <- tournament_select(fitness_scores, tournament_size, valid_solutions)
        parent2_idx <- tournament_select(fitness_scores, tournament_size, valid_solutions)
      } else {
        # If no valid solutions in current generation, use random selection
        parent1_idx <- sample(population_size, 1)
        parent2_idx <- sample(population_size, 1)
      }
      
      if (runif(1) < crossover_rate) {
        child <- crossover(population[[parent1_idx]], 
                         population[[parent2_idx]])
      } else {
        child <- population[[parent1_idx]]
      }
      
      child <- mutate(child, mutation_rate)
      new_population[[i]] <- child
    }
    
    population <- new_population
  }
  
  # Return best model found
  if (is.null(best_model)) {
    stop("No valid solution found during optimization")
  }
  
  # Add GA results to output
  best_model$ga_results <- list(
    fitness_history = fitness_history[1:gen],
    generation_history = generation_history[1:gen],
    best_fitness = best_fitness,
    generations_run = gen,
    final_components = genome_to_components(best_solution, n_vars),
    ga_parameters = list(
      population_size = population_size,
      max_generations = max_generations,
      mutation_rate = mutation_rate,
      crossover_rate = crossover_rate,
      tournament_size = tournament_size
    )
  )
  
  return(best_model)
}

# Updated helper functions

tournament_select <- function(fitness_scores, tournament_size, valid_solutions = NULL) {
  if (!is.null(valid_solutions)) {
    valid_indices <- which(valid_solutions)
    if (length(valid_indices) > 0) {
      contestants <- sample(valid_indices, min(tournament_size, length(valid_indices)))
    } else {
      contestants <- sample(length(fitness_scores), tournament_size)
    }
  } else {
    contestants <- sample(length(fitness_scores), tournament_size)
  }
  return(contestants[which.max(fitness_scores[contestants])])
}

# Rest of the helper functions remain the same
initialize_population <- function(n_vars, population_size, max_vars) {
  population <- vector("list", population_size)
  
  for (i in 1:population_size) {
    if (!is.null(max_vars)) {
      n_active <- sample(1:min(max_vars, n_vars), 1)
      genome <- rep(0, n_vars)
      genome[sample(1:n_vars, n_active)] <- 1
    } else {
      genome <- sample(c(0,1), n_vars, replace = TRUE)
    }
    if (sum(genome) == 0) {
      genome[sample(1:n_vars, 1)] <- 1
    }
    population[[i]] <- genome
  }
  
  return(population)
}

crossover <- function(parent1, parent2) {
  crossover_point <- sample(1:(length(parent1)-1), 1)
  child <- c(parent1[1:crossover_point], parent2[(crossover_point+1):length(parent2)])
  return(child)
}

mutate <- function(genome, mutation_rate) {
  for (i in 1:length(genome)) {
    if (runif(1) < mutation_rate) {
      genome[i] <- 1 - genome[i]
    }
  }
  if (sum(genome) == 0) {
    genome[sample(1:length(genome), 1)] <- 1
  }
  return(genome)
}

genome_to_components <- function(genome, n_vars) {
  components <- matrix(c(1, genome), nrow = 1)
  colnames(components) <- c("const", paste0("X", 1:n_vars))
  return(components)
}