
# load Rcpp files ---------------------------------------------------------

library("Rcpp")
sourceCpp("hierLDAem.cpp")


# main --------------------------------------------------------------------

hlda_em <- function(setting, corpus, tree_structure){
  
  K <- setting$num_topics
  num_cats <- length(tree_structure)+1
  # two parameter lists:
  lda_model <- list(gamma = setting$initial_gamma, 
                      eta = setting$initial_eta, 
                      alpha = array(setting$initial_alpha, num_cats) )
  var_model <- list(kappa = matrix(1/K, length(tree_structure)+1, K), 
                tau = array(K, length(tree_structure)+1), 
                nu = matrix(0, corpus$num_sessions, K), 
                rho = rep(list(matrix(), corpus$num_sessions)), 
                lambda = matrix(0, K, corpus$num_apps) )
  for (i in 1:corpus$num_sessions) {
    var_model$rho[[i]] = matrix(0, corpus$sessions[[i]]$length, num_topics)
  }
  
  #initializations:
  suffstats <- matrix(0, K, corpus$num_apps)
  digamma_lambda <- matrix(0, K, corpus$num_apps)
  digamma_lambda_sum <- array(0, K)
  lgamma_lambda_sum <- array(0, K)
  oneoverk <- 1/K
  etaoverv <- lda_model$eta / corpus$num_apps
  #initialize lambda
  corpus_initialize_ss(suffstats, K, corpus$num_apps, corpus$num_sessions, setting$num_sessions_for_init, corpus$sessions)
  opt_lambda(suffstats, K, corpus$num_apps, lgamma_lambda_sum, etaoverv, digamma_lambda_sum, digamma_lambda, var_model$lambda) 
  
  digamma_nu <- matrix(0, corpus$num_sessions, K)
  digamma_nu_sum <- array(0, corpus$num_sessions)
  digamma_sum_over_children <- matrix(0, num_cats, K)
  digamma_sum_over_children_for_kappa <- matrix(0, num_cats, K)
  
  dirichlet_prior <- matrix( (setting$initial_alpha) * oneoverk , num_cats, K)
  dirichlet_prior_root <- array( lda_model$gamma/K , K)
  
  #----EM function----
  em_step <- array(0,num_cats)
  returned_likelihood_from_user <- array(0,num_cats-1)
  
  #initialize values in root loop
  whole_likelihood_old <- 0
  old_node_likelihood <- array(0,num_cats)
  while(1){
    em_step[1] <- em_step[1] + 1
    #loop through each user
    for(node_index in 2:num_cats){
      num_children <- length(tree_structure[[node_index-1]])
      while(1){
        em_step[node_index] <- em_step[node_index] + 1
        digamma_sum_over_children_for_kappa[node_index,] <- 0
        digamma_sum_over_children[node_index,] <- 0
        children_sum <- 0
        # E-STEP: update nu, rho for each session of current user
        for(session_index in tree_structure[[node_index-1]]){
          nu <- var_model$nu[session_index,]
          rho <- var_model$rho[[session_index]]
          digamma_nu_session <- digamma_nu[session_index,]
          children_sum <- children_sum + session_e_step(corpus$sessions[[session_index]], dirichlet_prior[node_index,],
                                                        nu, digamma_lambda, digamma_lambda_sum, setting, session_index-1,
                                                        rho, digamma_nu_session, digamma_nu_sum, oneoverk)
          digamma_nu[session_index,] <- digamma_nu_session
          var_model$nu[session_index,] <- nu
          var_model$rho[[session_index]] <- rho
        }
        # collect outputs to be used to update other parameters
        loop1(var_model$tau, var_model$kappa, K, node_index-1, tree_structure[[node_index-1]],
              user_size = 0, length(tree_structure[[node_index-1]]),
              digamma_sum_over_children_for_kappa, digamma_sum_over_children, digamma_nu,
              digamma_nu_sum)

        tau <- var_model$tau[node_index]
        kappa <- var_model$kappa[node_index,]
        alpha <- lda_model$alpha[node_index]
        # E-STEP: update tau, kappa, alpha for current user
        node_likelihood <- cat_e_step(tau, kappa, alpha, K, num_children, setting,
                                      children_sum, dirichlet_prior[1,], lda_model$alpha[1],
                                      digamma_sum_over_children[node_index,],
                                      digamma_sum_over_children_for_kappa[node_index,], node_index-1,
                                      oneoverk)
        var_model$tau[node_index] <- tau
        var_model$kappa[node_index,] <- kappa
        lda_model$alpha[node_index] <- alpha
        
        # M-STEP: check user-likelihood convergence
        node_converged <- (old_node_likelihood[node_index] - node_likelihood) / old_node_likelihood[node_index]
        
        if ( (em_step[node_index] > 1) && (node_converged < setting$cat_converged) ) {
          returned_likelihood_from_user[node_index-1] <- node_likelihood
          
          break  #start next user loop(or exit to the root loop)
        } else if (em_step[node_index] > setting$cat_max_iter) {
          print("max iteration reached")
          break
        } else { # retry children
          old_node_likelihood[node_index] <- node_likelihood
          
          # update dirichlet_prior
          dirichlet_prior[node_index,] <- lda_model$alpha[node_index] * var_model$kappa[node_index,]
        }
      }
    }
    num_children <- num_cats-1
    
    children_sum <- sum(returned_likelihood_from_user)
    digamma_sum_over_children_for_kappa[1,] <- 0
    digamma_sum_over_children[1,] <- 0
    loop1(var_model$tau, var_model$kappa, K, 0, sessionids = 0,
          user_size = num_cats-1, session_size = 0, digamma_sum_over_children_for_kappa,
          digamma_sum_over_children, digamma_nu, digamma_nu_sum)
    
    tau <- var_model$tau[1]
    kappa <- var_model$kappa[1,]
    alpha <- lda_model$alpha[1]
    # update global alpha,tau,kappa
    root_likelihood <- cat_e_step(tau,kappa,alpha,K, num_children, setting, children_sum,
                                  dirichlet_prior_root, lda_model$gamma, digamma_sum_over_children[1,],
                                  digamma_sum_over_children_for_kappa[1,], node_index = 0, oneoverk)
    var_model$tau[1] <- tau
    var_model$kappa[1,] <- kappa
    lda_model$alpha[1] <- alpha
    #compute a statistic for updating lambda
    collect_lambda_ss(suffstats, var_model, corpus, K)
    whole_likelihood <- root_likelihood + opt_lambda(suffstats, K, corpus$num_apps, lgamma_lambda_sum, etaoverv, digamma_lambda_sum, digamma_lambda,var_model$lambda)
    whole_converged <- (whole_likelihood_old - whole_likelihood) / whole_likelihood_old
    
    # check global convergence
    if ( (em_step[1] > 1) && (whole_converged < setting$em_converged) ) {
      break #converged
    } else if (em_step[1] > setting$em_max_iter) {
      print("max iteration reached")
      break
    } else { # retry children
      whole_likelihood_old <- whole_likelihood
      
      # update dirichlet_prior
      dirichlet_prior[1,] <- lda_model$alpha[1] * var_model$kappa[1,]
    }
  }#end of root while loop
  par_list <- list(lambda = var_model$lambda, kappa = var_model$kappa, tau = var_model$tau,
                   alpha = lda_model$alpha, nu = var_model$nu)
  return(par_list)
}
