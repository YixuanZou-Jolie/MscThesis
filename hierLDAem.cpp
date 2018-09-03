#include <Rcpp.h>
#include <assert.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
void corpus_initialize_ss(NumericMatrix ss, const int num_topics, const int num_apps,
                          const int num_sessions, const int num_sessions_for_init, const List c_s)
{ 
  int i, k, d, n;
  List session;
  
  for (k = 0; k < num_topics; k++){
    for (i = 0; i < num_sessions_for_init; i++){
      d = floor(R::runif(0,1) * num_sessions);//randomly choose a session
      session = as<List>(c_s[d]);
      
      for (n = 0; n < int(session["length"]); n++){ 
        ss( k,int(as<List>(session["apps"])[n])-1 ) += int(as<List>(session["counts"])[n]);
      }
    }
    for (n = 0; n < num_apps; ++n) {
      ss(k,n) += 1; //avoid zero occurance
    }
  }
}
//--------------------------------------------------------------------------
// [[Rcpp::export]]
void loop1(const NumericVector tau, const NumericMatrix kappa,
           const int K, const int node_index,
           const NumericVector sessionids,
           const int user_size, const int session_size,
           NumericMatrix digamma_sum_over_children_for_kappa,
           NumericMatrix digamma_sum_over_children,
           const NumericMatrix digamma_nu,
           const NumericVector digamma_nu_sum)
{ 
  if (0 < node_index){ //user node
    for (unsigned int j = 0; j < session_size; ++j) {
      const int& d = sessionids[j]-1;
      for (int i = 0; i < K; ++i) {
        digamma_sum_over_children_for_kappa(node_index,i) += digamma_nu(d,i);
        digamma_sum_over_children(node_index,i) += digamma_nu(d,i) - digamma_nu_sum[d];
      }
    }
  } else {  //root node
    for (unsigned int j = 0; j < user_size; ++j) {
      const int& c = j;
      const double digamma_tau = R::digamma(tau[c+1]);
      
      for (int i = 0; i < K; ++i) {
        const double digamma_taukappai = R::digamma(tau[c+1] * kappa(c+1,i));
        digamma_sum_over_children_for_kappa(0,i) += digamma_taukappai;
        digamma_sum_over_children(0,i) += digamma_taukappai - digamma_tau;
      }
    }
  }
}

//-----------------------------------------------------------------------------------
// [[Rcpp::export]]
void collect_lambda_ss(NumericMatrix ss, const List var_model, const List c, const int num_topics)
{
  for (int d = 0; d < int(c["num_sessions"]); ++d) {
    const List session = as<List>(c["sessions"])[d];
    for (int l = 0; l < int(session["length"]); ++l) {
      for (int i = 0; i < num_topics; ++i) {
        ss(i,int(as<List>(session["apps"])[l])-1) +=
          as<double>(as<List>(session["counts"])[l]) *
          as<NumericMatrix>(as<List>(var_model["rho"])[d])(l,i);
      }
    }
  }
}

//-----------------------------------------------------------------------
// [[Rcpp::export]]
double opt_lambda(const NumericMatrix ss, const int num_topics, const int num_apps, NumericVector lgamma_lambda_sum,
                  const double etaoverv, NumericVector digamma_lambda_sum,
                  NumericMatrix digamma_lambda, NumericMatrix lambda)
{
  double lambda_likelihood = 0;
  
  for (int i = 0; i < num_topics; ++i) {
    double lambda_sum = 0.0;
    lambda_likelihood -= lgamma_lambda_sum[i];
    
    for (int v = 0; v < num_apps; ++v) {
      // compute likelihood before updating lambda
      lambda_likelihood += lgamma(lambda(i,v));
      lambda_likelihood += (etaoverv - lambda(i,v)) * (digamma_lambda(i,v) - digamma_lambda_sum[i]);
      
      lambda(i,v) = etaoverv + ss(i,v);
      digamma_lambda(i,v) = R::digamma(lambda(i,v));
      lambda_sum += lambda(i,v);
    }
    digamma_lambda_sum[i] = R::digamma(lambda_sum);
    lgamma_lambda_sum[i] = lgamma(lambda_sum);
  }
  return lambda_likelihood;
}
//----------------------------------------------------------------------
// [[Rcpp::export]]
double log_sum(double log_a, double log_b)
{
  double v;
  
  if (log_a < log_b)
  {
    v = log_b+log(1 + exp(log_a-log_b));
  }
  else
  {
    v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}
//-------------------------------------------------------------------------
#define SESSION_DECREASE_ALLOWANCE 1e-6
//[[Rcpp::export]]
double session_e_step(const List session, const NumericVector dirichlet_prior, NumericVector nu,
                      const NumericMatrix digamma_lambda, const NumericVector digamma_lambda_sum, const List setting,
                      const int session_id, NumericMatrix rho,
                      NumericVector digamma_nu, NumericVector digamma_nu_sum, const double oneoverk)
{ 
  const int& numtopics = setting["num_topics"];
  const int& sessionlength = session["length"];
  const double sessiontotaloverk = (double) (session["total"]) / (double) (numtopics);
  NumericVector old_rho(numtopics);
  
  // Initialize rho, nu
  for (int i = 0; i < numtopics; ++i) {
    for (int l = 0; l < sessionlength; ++l) {
      rho(l,i) = oneoverk;
    }
    nu[i] = dirichlet_prior[i] + sessiontotaloverk;
    digamma_nu[i] = R::digamma(nu[i]);
  }
  int session_loop = 0;
  double session_likelihood = 0;
  double session_likelihood_old = 0;
  double session_converged = 1;
  double nu_sum = 0;
  double indep_part_likelihood = 0;
  double dep_part_likelihood = 0;
  
  while ((session_loop < 2)
           || ((session_converged > int(setting["session_converged"])) && (session_loop < int(setting["session_max_iter"])))) {
    session_loop += 1;
    for (int l = 0; l < sessionlength; ++l) {
      double rhosum = 0;
      for (int i = 0; i < numtopics; ++i) {
        old_rho[i] = rho(l,i);
        rho(l,i) = digamma_nu[i] + digamma_lambda(i,int(as<List>(session["apps"])[l])-1) - digamma_lambda_sum[i];
        assert(rho(l,i) != 0);
        assert(!std::isnan(rho(l,i)));
        if (i > 0) {
          rhosum = log_sum(rhosum, rho(l,i));
        } else {
          rhosum = rho(l,i);
        }
        assert(!std::isnan(rhosum));
      }
      for (int i = 0; i < numtopics; ++i) {
        rho(l,i) = exp(rho(l,i) - rhosum);
        nu[i] = nu[i] + int(as<List>(session["counts"])[l]) * (rho(l,i) - old_rho[i]);
        digamma_nu[i] = R::digamma(nu[i]);
        assert(!std::isnan(digamma_nu[i]));
      }
    }
    nu_sum = 0;
    for (int i = 0; i < numtopics; ++i) {
      nu_sum += nu[i];
    }
    digamma_nu_sum[session_id] = R::digamma(nu_sum);
    
    indep_part_likelihood = -lgamma(nu_sum);
    dep_part_likelihood = 0;
    for (int i = 0; i < numtopics; ++i) {
      double delta = (digamma_nu[i] - digamma_nu_sum[session_id]);
      indep_part_likelihood += lgamma(nu[i]) - delta * nu[i];
      dep_part_likelihood += delta * dirichlet_prior[i];
      for (int l = 0; l < sessionlength; ++l) {
        if (rho(l,i) > 0) {
          indep_part_likelihood += rho(l,i) * int(as<List>(session["counts"])[l])
          * (delta + digamma_lambda(i,int(as<List>(session["apps"])[l])-1) - digamma_lambda_sum[i] - log(rho(l,i)));
        }
      }
    }
    assert(!std::isnan(indep_part_likelihood));
    assert(!std::isnan(dep_part_likelihood));
    
    session_likelihood = indep_part_likelihood + dep_part_likelihood;
    session_converged = (session_likelihood_old - session_likelihood) / session_likelihood_old;
    session_likelihood_old = session_likelihood;
  }
  if (session_loop >= int(setting["session_max_iter"])) {
    printf("session loop max reached %d\n", session_id);
  }
  return indep_part_likelihood;
}
//-------------------------------------------------------------------------
// [[Rcpp::export]]
double tetragamma(double x)
{
  double p;
  int i;
  
  x=x+6;
  p=1/(x*x);
  p=(((((0.3 - 0.833333333333333 * p) * p - 0.166666666666666) * p + 0.166666666666666) * p - 0.5) * p - 1/x - 1) * p;
  for (i=0; i<6 ;i++)
  {
    x=x-1;
    p = p - 2 / (x*x*x);
  }
  return(p);
}
//------------------------------------------------------------------------------------
#define KAPPA_NEWTON_THRESH 1e-6
#define KAPPA_MAX_ITER 5000
#define LIKELIHOOD_DECREASE_ALLOWANCE 1e-5
// [[Rcpp::export]]
double opt_kappa(NumericVector kappa, const int ntopics, const int nchildren,
                 const NumericVector dirichlet_prior,
                 const double alpha, const NumericVector tau,
                 const NumericVector digamma_sum_over_children,
                 const int node_index, const double oneoverk)
{
  NumericVector g(ntopics);
  NumericVector h(ntopics);
  NumericVector delta_kappa(ntopics);
  NumericVector new_kappa(ntopics);
  for (int i = 0; i < ntopics; ++i) {
    kappa[i] = oneoverk;
  }
  double	invhsum = 0;
  double	goverhsum = 0;
  double	coefficient = 0;
  double	old_likelihood = 0;
  double	sqr_newton_decrement = 0;
  double	step_size;
  double	indep_new_likelihood = 0;
  double	dep_new_likelihood = 0;
  double	new_likelihood;
  double	expected_increase;
  int		iter = 0;
  
  for (int i = 0; i < ntopics; ++i) {
    double const taukappai = tau[0] * kappa[i];
    double const alphakappai = alpha * kappa[i];
    double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
    double const digammataukappai = R::digamma(taukappai);
    double const logkappai = log(kappa[i]);
    dep_new_likelihood += digammataukappai * common;
    indep_new_likelihood -= nchildren * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
    indep_new_likelihood += alphakappai * digamma_sum_over_children[i];
    dep_new_likelihood += lgamma(taukappai);
  }
  new_likelihood = indep_new_likelihood + dep_new_likelihood;
  do {
    iter++;
    invhsum = 0;
    goverhsum = 0;
    coefficient = 0;
    for (int i = 0; i < ntopics; ++i) {
      double const taukappai = tau[0] * kappa[i];
      double const alphakappai = alpha * kappa[i];
      double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
      double const digammataukappai = R::digamma(taukappai);
      double const trigammataukappai = R::trigamma(taukappai);
      double const logkappai = log(kappa[i]);
      
      g[i] = tau[0] * trigammataukappai * common
        - nchildren * alpha * (R::digamma(alphakappai) - logkappai + digammataukappai - 1)
        - nchildren / kappa[i]
        + alpha * digamma_sum_over_children[i];
        
        h[i] = tau[0] * tau[0] * tetragamma(taukappai) * common
          - tau[0] * trigammataukappai * (tau[0] + 2 * alpha * nchildren)
          - alpha * alpha * R::trigamma(alphakappai) * nchildren
          + alpha * nchildren / kappa[i]
          + nchildren / (kappa[i] * kappa[i]);
          
          invhsum += 1 / h[i];
          goverhsum += g[i] / h[i];
    }
    
    old_likelihood = new_likelihood;
    
    coefficient = goverhsum / invhsum;
    sqr_newton_decrement = 0;
    expected_increase = 0;
    step_size = 1;

    for (int i = 0; i < ntopics; ++i) {
      delta_kappa[i] = (coefficient - g[i]) / h[i];
      sqr_newton_decrement -= h[i] * delta_kappa[i] * delta_kappa[i]; // this one is maximization
      expected_increase += g[i] * delta_kappa[i];
      if (delta_kappa[i] < 0) {
        double limit = (kappa[i] - 1e-10) / -(delta_kappa[i]);
        if (step_size > limit) {
          step_size = limit;
        }
      }
    }
    if (sqr_newton_decrement < KAPPA_NEWTON_THRESH * 2 || step_size < 1e-8 ) {
      break;
    }
    // backtracking line search
    while(1) {
      indep_new_likelihood = 0.0;
      dep_new_likelihood = 0.0;
      for (int i = 0; i < ntopics; ++i) {
        new_kappa[i] = kappa[i] + step_size * delta_kappa[i];
        
        double const taukappai = tau[0] * new_kappa[i];
        double const alphakappai = alpha * new_kappa[i];
        double const common = dirichlet_prior[i] + nchildren * (1 - alphakappai) - taukappai;
        double const logkappai = log(new_kappa[i]);
        
        dep_new_likelihood += R::digamma(taukappai) * common;
        indep_new_likelihood -= nchildren * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
        indep_new_likelihood += alphakappai * digamma_sum_over_children[i];
        dep_new_likelihood += lgamma(taukappai);
      }
      new_likelihood = indep_new_likelihood + dep_new_likelihood;
      if (new_likelihood > old_likelihood + 0.4 * step_size * expected_increase) {
        break;
      }
      step_size *= 0.9;
      if (step_size < 1e-8) break;
    }
    if (step_size < 1e-8) break;
    
    for (int i = 0; i < ntopics; ++i) {
      kappa[i] = new_kappa[i];
      assert(!std::isnan(kappa[i]));
      assert(kappa[i] > 0);
    }
  } while (iter < KAPPA_MAX_ITER);
  if (iter >= KAPPA_MAX_ITER) {
    exit(-1);
  }
  return new_likelihood;
}
//--------------------------------------------------------------------------
#define TAU_NEWTON_THRESH  1e-6
#define TAU_MAX_ITER  5000
// [[Rcpp::export]]
double opt_tau(NumericVector tau_input, const NumericVector kappa,
               const int ntopics, const int nchildren,
               const NumericVector dirichlet_prior,
               const double alpha,
               const int node_index)
{
  double tau;
  int		iter = 0;
  double	d1 = 0;
  double	d2 = 0;
  double	init_tau = 100;
  double	log_tau = 0;
  double	likelihood = 0;
  double	old_likelihood = 0;
  
  tau = init_tau;
  log_tau = log(tau);
  do {
    iter++;
    d1 = 0;
    d2 = 0;
    likelihood = 0;
    
    double const common2 = nchildren * alpha * (ntopics - 1) / tau;
    double const trigammatau = R::trigamma(tau);
    double const digammatau = R::digamma(tau);
    for (int i = 0; i < ntopics; ++i) {
      double const taukappai = tau * kappa[i];
      double const trigammataukappai = R::trigamma(taukappai);
      double const common = dirichlet_prior[i] - taukappai + nchildren * (1 - alpha * kappa[i]);
      
      d1 += (trigammataukappai * kappa[i] - trigammatau) * common;
      d2 += kappa[i] * kappa[i] * ( tetragamma(taukappai) * common - trigammataukappai);
      d2 -= tetragamma(tau) * common;
      likelihood += (R::digamma(taukappai) - digammatau) * common
        + lgamma(taukappai);
    }
    d1 += common2 / tau;
    d2 += trigammatau - 2 * common2 / tau / tau;
    
    likelihood -= lgamma(tau);
    likelihood -= common2;
    
    assert(!std::isnan(d1));
    assert(!std::isnan(d2));
    assert(!std::isnan(likelihood));
    assert( (old_likelihood == 0) || (likelihood >= old_likelihood) );
    
    old_likelihood = likelihood;
    if (fabs(d1) < TAU_NEWTON_THRESH) {
      break;
    }
    log_tau = log_tau - d1 / (d2 * tau + d1);
    tau = exp(log_tau);
    if (std::isnan(tau) || tau < 1e-10) {
      init_tau = init_tau * 10;
      printf("warning: tau is nan; new init = %5.5f\n", init_tau);
      tau = init_tau;
      log_tau = log(tau);
      old_likelihood = 0;
    }
  } while (iter < TAU_MAX_ITER);
  
  if (iter >= TAU_MAX_ITER) {
    printf("tau iter max reached\n");
    exit(-1);
  }
  tau_input[0] = tau;
  return likelihood;
}
//--------------------------------------------------------------------------
#define ALPHA_NEWTON_THRESH  1e-6
#define ALPHA_MAX_ITER  5000
#define ALPHA_DECREASE_ALLOWANCE  1e-8
// [[Rcpp::export]]
double opt_alpha(NumericVector alpha_old, const int ntopics, const int nchildren,
                 const NumericVector tau_input,
                 const NumericVector kappa, 
                 const NumericVector digamma_sum_over_children,
                 const int node_index)
{ 
  double alpha;
  double tau = tau_input[0];
  
  int		iter = 0;
  double	d1 = 0;
  double	d2 = 0;
  double	likelihood = 0;
  double	old_likelihood = 0;
  double	init_alpha = 100;
  double	log_alpha = 0;
  double precompute = 0;
  
  precompute = nchildren * (R::digamma(tau) - (ntopics - 1) / tau);
  for (int i = 0; i < ntopics; ++i) {
    precompute += nchildren * kappa[i] * (log(kappa[i]) - R::digamma(tau * kappa[i]));
    precompute += kappa[i] * digamma_sum_over_children[i];
  }
  
  alpha = init_alpha;
  log_alpha = log(alpha);
  do {
    iter++;
    likelihood = 0;
    d1 = 0;
    d2 = 0;
    
    for (int i = 0; i < ntopics; ++i) {
      const double alphakappai = alpha * kappa[i];
      
      likelihood -= lgamma(alphakappai);
      d1 -= kappa[i] * R::digamma(alphakappai);
      d2 -= kappa[i] * kappa[i] * R::trigamma(alphakappai);
    }
    likelihood = nchildren * (likelihood +  lgamma(alpha))+ alpha * precompute;
    d1 = nchildren * (d1 + R::digamma(alpha)) + precompute;
    d2 = nchildren * (d2 + R::trigamma(alpha));
    
    assert(!std::isnan(likelihood));
    assert(!std::isnan(d1));
    assert(!std::isnan(d2));
    
    old_likelihood = likelihood;
    if (fabs(d1) < ALPHA_NEWTON_THRESH) {
      break;
    }
    
    log_alpha = log_alpha - d1 / (d2 * alpha + d1);
    alpha = exp(log_alpha);
    if (std::isnan(alpha) || alpha < 1e-10) {
      init_alpha = init_alpha * 10;
      printf("warning: alpha is nan; new init = %5.5f\n", init_alpha);
      alpha = init_alpha;
      log_alpha = log(alpha);
      old_likelihood = 0;
    }
  } while (iter < ALPHA_MAX_ITER);
  
  if (iter >= ALPHA_MAX_ITER) {
    printf("alpha iter max reached\n");
    exit(-1);
  }
  alpha_old[0] = alpha;
  return likelihood;
}
//-------------------------------cat_e_step.cpp-------------------------------------------
// [[Rcpp::export]]
double cat_e_step(NumericVector tau, NumericVector kappa, NumericVector alpha,
                  const int K, const int num_children, const List setting,
                  const double children_sum,
                  const NumericVector dirichlet_prior,
                  const NumericVector alpha_pi,
                  const NumericVector digamma_sum_over_children,
                  const NumericVector digamma_sum_over_children_for_kappa,
                  const int node_index, const double oneoverk)
{ double dep_likelihood = 0.0;
  double indep_likelihood = 0.0;
  
  double kappa_tau_likelihood = 0.0;
  double kappa_tau_likelihood_old = 0.0;
  double kappa_tau_converged;
  int kappa_tau_loop = 0;
  
  while ((kappa_tau_loop < 2) ||
         ((kappa_tau_converged > int(setting["kappa_tau_converged"])) && (kappa_tau_loop < int(setting["kappa_tau_max_iter"])))) {
    kappa_tau_loop += 1;
    kappa_tau_likelihood = 0.0;
    opt_tau(tau,
            kappa,
            K, num_children,
            dirichlet_prior,
            alpha[0],
                 node_index);
    kappa_tau_likelihood += opt_kappa(kappa, K, num_children,
                                      dirichlet_prior,
                                      alpha[0], tau,
                                      digamma_sum_over_children_for_kappa, node_index, oneoverk);
    for (int i = 0; i < K; ++i) {
      kappa_tau_likelihood += R::digamma(tau[0]) * (-dirichlet_prior[i] + tau[0] * kappa[i] + (alpha[0] * kappa[i] - 1) * num_children);
    }
    kappa_tau_likelihood -= lgamma(tau[0]);
    kappa_tau_likelihood -= alpha[0] * num_children * (K - 1) / tau[0];

    kappa_tau_converged = (kappa_tau_likelihood_old - kappa_tau_likelihood) / kappa_tau_likelihood_old;
    kappa_tau_likelihood_old = kappa_tau_likelihood;
  }
  
  if (kappa_tau_loop >= int(setting["kappa_tau_max_iter"])) {
    printf("kappa_tau_loop max reached\n");
    exit(-1);
  }
  
  if (setting["estimate_alpha"]) {
    indep_likelihood += opt_alpha(alpha, K, num_children, tau, kappa, digamma_sum_over_children, node_index);
  } else {
    double precompute = 0.0;
    double alpha_likelihood = 0.0;
    
    precompute = num_children * (R::digamma(tau[0]) - (K - 1) / tau[0]);
    for (int i = 0; i < K; ++i) {
      const double alphakappai = alpha[0] * kappa[i];
      
      alpha_likelihood -= lgamma(alphakappai);
      precompute += num_children * kappa[i] * (log(kappa[i]) - R::digamma(tau[0] * kappa[i]));
      precompute += kappa[i] * digamma_sum_over_children[i];
    }
    alpha_likelihood = num_children * (alpha_likelihood + lgamma(alpha[0])) + alpha[0] * precompute;
    indep_likelihood += alpha_likelihood;
  }
  
  indep_likelihood += children_sum;
  const double digamma_tau = R::digamma(tau[0]);
  indep_likelihood -= num_children * K * digamma_tau;
  for (int i = 0; i < K; ++i) {
    const double& kappai = kappa[i];
    const double taukappai = tau[0] * kappai;
    const double digammataukappai = R::digamma(taukappai);
    const double common = (digammataukappai - digamma_tau);
    
    indep_likelihood -= num_children * (log(kappai) - digammataukappai);
    indep_likelihood += lgamma(taukappai);
    indep_likelihood -= taukappai * common;
    dep_likelihood += dirichlet_prior[i] * common;
  }
  indep_likelihood -= lgamma(tau[0]);
  
  assert(!std::isnan(indep_likelihood));
  assert(!std::isnan(dep_likelihood));
  
  double cat_likelihood = indep_likelihood+dep_likelihood;
  return cat_likelihood ;
}