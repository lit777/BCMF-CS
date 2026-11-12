#include <Rcpp.h>

#include "decision_tree.h"
#include "MCMC_utils_test.h"

using namespace std;
using namespace R;
using namespace Rcpp;
using namespace sugar;

// [[Rcpp::export]]
List MCMC(
    const NumericMatrix& Xpred,
    const NumericMatrix& Xpred1,
    const NumericMatrix& Xpred2,
    const NumericMatrix& Xpred3,
    const NumericMatrix& Xpred_mult,
    const NumericMatrix& XpredM_mult,
    const NumericVector& Y_trt,
    NumericVector& M_out,
    NumericVector& Y_out,
    const double p_grow,   // Prob. of GROW
    const double p_prune,  // Prob. of PRUNE
    const double p_change, // Prob. of CHANGE
    const int m_ps,          // Num. of Trees: PS
    const int m_med1,          // Num. of Trees: Prognostic Func in Med
    const int m_med2,          // Num. of Trees: Modifier Func in Med
    const int m_out1,          // Num. of Trees: Prognostic Func in Out
    const int m_out2,          // Num. of Trees: Modifier Func for treatment effect
    const int m_out3,          // Num. of Trees: Modifier Func for mediator effect
    const int nu,
    double lambda,
    double lambda_m, 
    double lambda_y,
    double dir_alpha, 
    double dir_alpha_m, 
    double dir_alpha_y1,
    double dir_alpha_y2,
    double alpha, 
    double alpha_modifier, 
    double beta, 
    const int n_iter,
    const bool verbose = false
) {

    // Data preparation
    const int P  = Xpred.ncol(); // number of covariates
    const int P1  = Xpred1.ncol(); // number of covariates
    const int P2  = Xpred2.ncol(); // number of covariates
    const int P3  = Xpred3.ncol(); // number of covariates
    const int n  = Xpred.nrow(); // number of observations

    NumericVector M1(n);
    NumericVector M0(n);
    NumericVector Y11(n), Y00(n), Y10(n);
    
    const double shift = mean(Y_out);
    const double ysd = sd(Y_out);
    const double mshift = mean(M_out);
    const double msd = sd(M_out);
    
    Y_out = (Y_out - shift) / ysd; // scaling
    M_out = (M_out - mshift) / msd; // scaling

    NumericVector Xcut[P]; // e.g. unique value of potential confounders
    for (int j = 0; j < P; j++)
    {
        NumericVector temp;
        temp = unique(Xpred(_, j));
        temp.sort();
        Xcut[j] = clone(temp);
    }
    
    NumericVector Xcut1[P1]; // e.g. unique value of potential confounders
    for (int j = 0; j < P1; j++)
    {
        NumericVector temp;
        temp = unique(Xpred1(_, j));
        temp.sort();
        Xcut1[j] = clone(temp);
    }
    
    NumericVector Xcut2[P2]; // e.g. unique value of potential confounders
    for (int j = 0; j < P2; j++)
    {
        NumericVector temp;
        temp = unique(Xpred2(_, j));
        temp.sort();
        Xcut2[j] = clone(temp);
    }
    
    NumericVector Xcut3[P3]; // e.g. unique value of potential confounders
    for (int j = 0; j < P3; j++)
    {
      NumericVector temp;
      temp = unique(Xpred3(_, j));
      temp.sort();
      Xcut3[j] = clone(temp);
    }
    

    // Initial Setup
    // Priors, initial values and hyper-parameters
    NumericVector Z = Rcpp::rnorm(n, R::qnorm(mean(Y_trt), 0, 1, true, false), 1); // latent variable
    NumericVector prob = {p_grow, p_prune, p_change};

    double sigma2 = 1;
    NumericVector sigma2_m       (n_iter + 1); // create placeholder for sigma2_m
    sigma2_m(0)       = 1;
    NumericVector sigma2_y       (n_iter + 1); // create placeholder for sigma2_y
    sigma2_y(0)       = 1;
    
    
    // sigma_mu based on min/max of Z, M (A=1) and Y (A=0)
    double sigma_mu   = std::max(pow(min(Z)  / (-2 * sqrt(m_ps)), 2), pow(max(Z)  / (2 * sqrt(m_ps)), 2));
    double sigma_mu_m_mu = pow((3 / 2), 2) / m_med1;
    double sigma_mu_m_tau = pow((3 / 2), 2) / m_med2;
    double sigma_mu_y_mu = pow((3 / 2), 2) / m_out1;
    double sigma_mu_y_tau1 = pow((3 / 2), 2) / m_out2;
    double sigma_mu_y_tau2 = pow((3 / 2), 2) / m_out3;

    
    // Initial values of R
    NumericVector R  = clone(Z);
    NumericVector R_M = Rcpp::rnorm(n, 0, 1);
    NumericVector R_Y = Rcpp::rnorm(n, 0, 1);
    
    // Initial values for the selection probabilities 
    NumericVector post_dir_alpha  = rep(1.0, P+3);
    NumericVector post_dir_alpha3  = rep(1.0, P);
    NumericVector post_dir_alpha5 = rep(1.0, P2);
    NumericVector post_dir_alpha6 = rep(1.0, P2);
    
    int thin       = 10;    //NumericVector post_dir_alpha  = rep(1.0, P+1);
    int burn_in    = n_iter / 2;
    int n_post     = (n_iter - burn_in) / thin; // number of post sample
    int thin_count = 1;

    
    IntegerMatrix ind (n_post, P);
    ind(0, _) = rep(0, P);
    IntegerMatrix ind1 (n_post, P+1);
    ind1(0, _) = rep(0, P+1);
    IntegerMatrix ind2 (n_post, P+3);
    ind2(0, _) = rep(0, P+3);
    int ind_idx = 0;


    NumericVector Effect (n_post);
    NumericMatrix predicted_M1 (n, n_post);
    NumericMatrix predicted_M0 (n, n_post);
    NumericMatrix predicted_Y11 (n, n_post);
    NumericMatrix predicted_Y00 (n, n_post);
    NumericMatrix predicted_Y10 (n, n_post);
    NumericMatrix predicted_Y (n, n_post);
    NumericMatrix predicted_zeta (Xpred_mult.nrow(), n_post);
    NumericMatrix predicted_d (Xpred_mult.nrow(), n_post);
    NumericMatrix predicted_tau (Xpred_mult.nrow(), n_post);
    

    int post_sample_idx = 0;

    IntegerMatrix Obs1_list (n, m_ps);
    IntegerMatrix Obs2_list (n, m_med1);
    IntegerMatrix Obs3_list (n, m_med2);
    IntegerMatrix Obs4_list (n, m_out1);
    IntegerMatrix Obs5_list (n, m_out2);
    IntegerMatrix Obs6_list (n, m_out3);
    

    // Place-holder for the posterior samples
    NumericMatrix Tree1  (n,    m_ps);
    NumericMatrix Tree2  (n,    m_med1);
    NumericMatrix Tree3  (n,    m_med2);
    NumericMatrix Tree3_pred  (n,    m_med2);
    NumericMatrix Tree4  (n,    m_out1);
    NumericMatrix Tree5  (n,    m_out2);
    NumericMatrix Tree6  (n,    m_out3);
    NumericMatrix Tree5_pred  (n,    m_out2);
    NumericMatrix Tree6_pred  (n,    m_out3);
    
    DecisionTree dt1_list[m_ps];
    DecisionTree dt2_list[m_med1];
    DecisionTree dt3_list[m_med2];
    DecisionTree dt4_list[m_out1];
    DecisionTree dt5_list[m_out2];
    DecisionTree dt6_list[m_out3];

    NumericVector dir_alpha_hist = rep(dir_alpha, n_iter);
    NumericVector dir_alpha_m_hist = rep(dir_alpha_m, n_iter);
    NumericVector dir_alpha_y1_hist = rep(dir_alpha_y1, n_iter);
    NumericVector dir_alpha_y2_hist = rep(dir_alpha_y2, n_iter);
    

    for (int t = 0; t < m_ps; t++)
    {
        dt1_list[t] = DecisionTree(n, t);
    }
    for (int t = 0; t < m_med1; t++)
    {
        dt2_list[t] = DecisionTree(n, t);
    }
    for (int t = 0; t < m_med2; t++)
    {
      dt3_list[t] = DecisionTree(n, t);
    }
    for (int t = 0; t < m_out1; t++)
    {
      dt4_list[t] = DecisionTree(n, t);
    }
    for (int t = 0; t < m_out2; t++)
    {
      dt5_list[t] = DecisionTree(n, t);
    }
    for (int t = 0; t < m_out3; t++)
    {
      dt6_list[t] = DecisionTree(n, t);
    }

    // Obtaining namespace of MCMCpack package
    Environment MCMCpack = Environment::namespace_env("MCMCpack");

    // Picking up rinvgamma() and rdirichlet() function from MCMCpack package
    Function rinvgamma  = MCMCpack["rinvgamma"];
    Function rdirichlet = MCMCpack["rdirichlet"];

    NumericVector prop_prob = rdirichlet(1, rep(1, P+3));
    NumericVector prop_prob3 = rdirichlet(1, rep(1, P));
    NumericVector prop_prob5 = rdirichlet(1, rep(1, P2));
    NumericVector prop_prob6 = rdirichlet(1, rep(1, P2));

  
    
    ////////////////////////////////////////
    //////////   Run main MCMC    //////////
    ////////////////////////////////////////

    
    


    
    for (int iter = 1; iter <= n_iter; iter++)
    {

        Rcout << "Rcpp iter : " << iter << " of " << n_iter << std::endl;
        
        update_Z(Z, Y_trt, Tree1);
        
        // ------ Exposure Model

        for (int t=0; t<m_ps; t++) {
          // decision trees
          update_R(R, Z, Tree1, t);
          
            NumericVector prop_prob_exp = prop_prob[Rcpp::Range(3, P+2)];
            prop_prob_exp = prop_prob_exp / (sum(prop_prob) - prop_prob(0) - prop_prob(1) - prop_prob(2));

          // create new tree instance
          if (dt1_list[t].length()==1) {     // tree has no node yet
            dt1_list[t].GROW_first(
                Xpred, Xcut, sigma2, sigma_mu, R, Obs1_list,
                p_prune, p_grow, alpha, beta, prop_prob_exp
            );
          } else {
            int step = sample(3, 1, false, prob)(0);
            switch (step) {
            case 1:   // GROW step
              dt1_list[t].GROW(
                  Xpred, Xcut, sigma2, sigma_mu, R, Obs1_list,
                  p_prune, p_grow, alpha, beta, prop_prob_exp
              );
              break;
              
            case 2:   // PRUNE step
              dt1_list[t].PRUNE(
                  Xpred, Xcut, sigma2, sigma_mu, R, Obs1_list,
                  p_prune, p_grow, alpha, beta, prop_prob_exp
              );
              break;
              
            case 3:   // CHANGE step
              dt1_list[t].CHANGE(
                  Xpred, Xcut, sigma2, sigma_mu, R, Obs1_list,
                  p_prune, p_grow, alpha, beta, prop_prob_exp
              );
              break;
              
            default: {};
            } // end of switch
          } // end of tree instance
          dt1_list[t].Mean_Parameter(Tree1, sigma2, sigma_mu, R, Obs1_list);
        } // end of Exposure Model

 
        // ------ Mediator Model (mu function)
        for (int t = 0; t < m_med1; t++)
        {
            update_R_mu(R_M, M_out, Tree2, Tree3, Y_trt, t);

          NumericVector prop_prob_m = prop_prob[Rcpp::Range(2, P+2)];
          prop_prob_m = prop_prob_m / (sum(prop_prob) - prop_prob(0) - prop_prob(1));
          
          
            if (dt2_list[t].length() == 1)
            {
                // tree has no node yet
                dt2_list[t].GROW_first(
                    Xpred1, Xcut1, sigma2_m(iter - 1), sigma_mu_m_mu, R_M, Obs2_list,
                    p_prune, p_grow, alpha, beta, prop_prob_m);
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);

                switch (step)
                {
                    case 1: // GROW step
                        dt2_list[t].GROW(
                            Xpred1, Xcut1, sigma2_m(iter - 1), sigma_mu_m_mu, R_M, Obs2_list,
                           p_prune, p_grow, alpha, beta, prop_prob_m
                        );
                        break;

                    case 2: // PRUNE step
                        dt2_list[t].PRUNE(
                            Xpred1, Xcut1, sigma2_m(iter - 1), sigma_mu_m_mu, R_M, Obs2_list,
                            p_prune, p_grow, alpha, beta, prop_prob_m
                        );
                        break;

                    case 3: // CHANGE step
                        dt2_list[t].CHANGE(
                            Xpred1, Xcut1, sigma2_m(iter - 1), sigma_mu_m_mu, R_M, Obs2_list,
                            p_prune, p_grow, alpha, beta, prop_prob_m
                        );
                        break;

                    default: {};
                } // end of switch
            }     // end of tree instance
            
            dt2_list[t].Mean_Parameter(Tree2, sigma2_m(iter - 1), sigma_mu_m_mu, R_M, Obs2_list);
        }
        

        NumericVector sigma2_m_vec(n); 
        for (int i = 0; i < n; i++)
        {
          sigma2_m_vec(i) = sigma2_m(iter - 1) / pow(Y_trt(i)-0.5, 2);
        }
        
        // ------ Mediator Model (tau function)
        for (int t = 0; t < m_med2; t++)
        {
            update_R_tau(R_M, M_out, Tree2, Tree3, Y_trt, t);
            
            
            if (dt3_list[t].length() == 1)
            {
                // tree has no node yet
                dt3_list[t].GROW_first_weight(
                        Xpred, Xcut, sigma2_m_vec, sigma_mu_m_tau, R_M / (Y_trt-0.5), Obs3_list,
                        p_prune, p_grow, 0.25, 3, prop_prob3);
            }
            else
            {
                int step = sample(3, 1, false, prob)(0);
                //   Rcout << "Rcpp step start: " << step <<  std::endl;
                
                switch (step)
                {
                case 1: // GROW step
                    dt3_list[t].GROW_weight(
                            Xpred, Xcut, sigma2_m_vec, sigma_mu_m_tau, R_M / (Y_trt-0.5), Obs3_list,
                            p_prune, p_grow, 0.25, 3, prop_prob3
                    );
                    break;
                    
                case 2: // PRUNE step
                    dt3_list[t].PRUNE_weight(
                            Xpred, Xcut, sigma2_m_vec, sigma_mu_m_tau, R_M / (Y_trt-0.5), Obs3_list,
                            p_prune, p_grow, 0.25, 3, prop_prob3
                    );
                    break;
                    
                case 3: // CHANGE step
                    dt3_list[t].CHANGE_weight(
                            Xpred, Xcut, sigma2_m_vec, sigma_mu_m_tau, R_M / (Y_trt-0.5), Obs3_list,
                            p_prune, p_grow, 0.25, 3, prop_prob3
                    );
                    break;
                    
                default: {};
                } // end of switch
            }     // end of tree instance
            
            dt3_list[t].Mean_Parameter_weight(Tree3, sigma2_m_vec, sigma_mu_m_tau, R_M / (Y_trt-0.5), Obs3_list);
            dt3_list[t].Predict(Tree3_pred, Xcut, Xpred_mult, Xpred_mult.nrow());
        }
        
        
        NumericVector m_new = clone(rowSums(Tree2)) + clone(rowSums(Tree3)) * (Y_trt - 0.5);
        
        
        if (iter > 1)
        {
          for (int i = 0; i < n; i++)
          {
            M0(i) = sum(clone(Tree2)(i,_))+sum(clone(Tree3)(i,_))*(-0.5);
            M1(i) = sum(clone(Tree2)(i,_))+sum(clone(Tree3)(i,_))*(0.5);
          }
          
        }
    
        
        //  Sample variance parameter
        {
          NumericVector sigma2_m_temp = rinvgamma(1, nu / 2 + n / 2, nu * lambda_m / 2 + sum(pow(M_out - m_new, 2)) / 2);
           sigma2_m(iter) = sigma2_m_temp(0);
        }
        

        // ------ Outcome Model (mu function)
        for (int t = 0; t < m_out1; t++)
        {
          update_R_mu1(R_Y, Y_out, Tree4, Tree5, Tree6, Y_trt, M_out, t);
          
          if (dt4_list[t].length() == 1)
          {
            // tree has no node yet
            dt4_list[t].GROW_first(
                Xpred3, Xcut3, sigma2_y(iter - 1), sigma_mu_y_mu, R_Y, Obs4_list,
                p_prune, p_grow, alpha, beta, prop_prob);
          }
          else
          {
            int step = sample(3, 1, false, prob)(0);
            
            switch (step)
            {
            case 1: // GROW step
              dt4_list[t].GROW(
                  Xpred3, Xcut3, sigma2_y(iter - 1), sigma_mu_y_mu, R_Y, Obs4_list,
                  p_prune, p_grow, alpha, beta, prop_prob
              );
              break;
              
            case 2: // PRUNE step
              dt4_list[t].PRUNE(
                  Xpred3, Xcut3, sigma2_y(iter - 1), sigma_mu_y_mu, R_Y, Obs4_list,
                  p_prune, p_grow, alpha, beta, prop_prob
              );
              break;
              
            case 3: // CHANGE step
              dt4_list[t].CHANGE(
                  Xpred3, Xcut3, sigma2_y(iter - 1), sigma_mu_y_mu, R_Y, Obs4_list,
                  p_prune, p_grow, alpha, beta, prop_prob
              );
              break;
              
            default: {};
            } // end of switch
          }     // end of tree instance
          
          dt4_list[t].Mean_Parameter(Tree4, sigma2_y(iter - 1), sigma_mu_y_mu, R_Y, Obs4_list);
        }


         NumericVector sigma2_y_vec(n); 
          for (int i = 0; i < n; i++)
          {
            sigma2_y_vec(i) = sigma2_y(iter - 1) / pow(Y_trt(i)-0.5, 2);
          }
        
        
        // ------ Outcome Model (tau1 function)
        for (int t = 0; t < m_out2; t++)
        {
          update_R_tau1(R_Y, Y_out, Tree4, Tree5, Tree6, Y_trt, M_out, t);
          
          if (dt5_list[t].length() == 1)
          {
            // tree has no node yet
            dt5_list[t].GROW_first_weight(
                Xpred2, Xcut2, sigma2_y_vec, sigma_mu_y_tau1, R_Y / (Y_trt-0.5), Obs5_list,
                p_prune, p_grow, 0.25, 3, prop_prob5);
          }
          else
          {
            int step = sample(3, 1, false, prob)(0);
            
            switch (step)
            {
            case 1: // GROW step
              dt5_list[t].GROW_weight(
                  Xpred2, Xcut2, sigma2_y_vec, sigma_mu_y_tau1, R_Y / (Y_trt-0.5), Obs5_list,
                  p_prune, p_grow, 0.25, 3, prop_prob5
              );
              break;
              
            case 2: // PRUNE step
              dt5_list[t].PRUNE_weight(
                  Xpred2, Xcut2, sigma2_y_vec, sigma_mu_y_tau1, R_Y / (Y_trt-0.5), Obs5_list,
                  p_prune, p_grow, 0.25, 3, prop_prob5
              );
              break;
              
            case 3: // CHANGE step
              dt5_list[t].CHANGE_weight(
                  Xpred2, Xcut2, sigma2_y_vec, sigma_mu_y_tau1, R_Y / (Y_trt-0.5), Obs5_list,
                  p_prune, p_grow, 0.25, 3, prop_prob5
              );
              break;
              
            default: {};
            } // end of switch
          }     // end of tree instance
          
          dt5_list[t].Mean_Parameter_weight(Tree5, sigma2_y_vec, sigma_mu_y_tau1, R_Y / (Y_trt-0.5), Obs5_list);
          dt5_list[t].Predict(Tree5_pred, Xcut2, XpredM_mult, XpredM_mult.nrow());
        }


      
          NumericVector sigma2_y_vec1(n); 
          for (int i = 0; i < n; i++)
          {
            sigma2_y_vec1(i) = sigma2_y(iter - 1) / pow(M_out(i), 2);
          }
          
          
          // ------ Outcome Model (tau2 function)
          for (int t = 0; t < m_out3; t++)
          {
            update_R_tau2(R_Y, Y_out, Tree4, Tree5, Tree6, Y_trt, M_out, t);
            
            if (dt6_list[t].length() == 1)
            {
              // tree has no node yet
              dt6_list[t].GROW_first_weight(
                  Xpred2, Xcut2, sigma2_y_vec1, sigma_mu_y_tau2, R_Y / M_out, Obs6_list,
                  p_prune, p_grow, 0.25, 3, prop_prob6);
            }
            else
            {
              int step = sample(3, 1, false, prob)(0);
              
              switch (step)
              {
              case 1: // GROW step
                dt6_list[t].GROW_weight(
                    Xpred2, Xcut2, sigma2_y_vec1, sigma_mu_y_tau2, R_Y / M_out, Obs6_list,
                    p_prune, p_grow, 0.25, 3, prop_prob6
                );
                break;
                
              case 2: // PRUNE step
                dt6_list[t].PRUNE_weight(
                    Xpred2, Xcut2, sigma2_y_vec1, sigma_mu_y_tau2, R_Y / M_out, Obs6_list,
                    p_prune, p_grow, 0.25, 3, prop_prob6
                );
                break;
                
              case 3: // CHANGE step
                dt6_list[t].CHANGE_weight(
                    Xpred2, Xcut2, sigma2_y_vec1, sigma_mu_y_tau2, R_Y / M_out, Obs6_list,
                    p_prune, p_grow, 0.25, 3, prop_prob6
                );
                break;
                
              default: {};
              } // end of switch
            }     // end of tree instance
            
            dt6_list[t].Mean_Parameter_weight(Tree6, sigma2_y_vec1, sigma_mu_y_tau2, R_Y / M_out, Obs6_list);
            dt6_list[t].Predict(Tree6_pred, Xcut2, XpredM_mult, XpredM_mult.nrow());
          }
                    
        NumericVector y_new = clone(rowSums(Tree4)) + clone(rowSums(Tree5)) * (Y_trt - 0.5) + clone(rowSums(Tree6)) * (M_out);


        //  Sample variance parameter
        {
          NumericVector sigma2_y_temp = rinvgamma(1, nu / 2 + n / 2, nu * lambda_y / 2 + sum(pow(Y_out - y_new, 2)) / 2);
          sigma2_y(iter) = sigma2_y_temp(0);
        }
        

          for (int i = 0; i < n; i++)
          {
            Y11(i) = (clone(rowSums(Tree4))(i) + clone(rowSums(Tree5))(i) * (0.5) + clone(rowSums(Tree6))(i) * (M1(i)))*ysd + shift;
            Y00(i) = (clone(rowSums(Tree4))(i) + clone(rowSums(Tree5))(i) * (-0.5) + clone(rowSums(Tree6))(i) * (M0(i)))*ysd + shift;
            Y10(i) = (clone(rowSums(Tree4))(i) + clone(rowSums(Tree5))(i) * (0.5) + clone(rowSums(Tree6))(i) * (M0(i)))*ysd + shift;
          }
          
          
        // Num. of inclusion of each potential confounder
        NumericVector add1(P), add2(P1), add3(P), add4(P3), add5(P2), add6(P2);

        
        for (int t = 0; t < m_ps; t++)
        {
          add1 += dt1_list[t].num_included(P);
        }
        
        for (int t = 0; t < m_med1; t++)
        {
          add2 += dt2_list[t].num_included(P1);
        }
        
        for (int t = 0; t < m_med2; t++)
        {
          add3 += dt3_list[t].num_included(P);
        }
        for (int t = 0; t < m_out1; t++)
        {
          add4 += dt4_list[t].num_included(P3);
        }
        for (int t = 0; t < m_out2; t++)
        {
          add5 += dt5_list[t].num_included(P2);
        }
        for (int t = 0; t < m_out3; t++)
        {
          add6 += dt6_list[t].num_included(P2);
        }
        
       
        if (iter < n_iter/10) {
          post_dir_alpha3 = rep(1.0, P) + add3;
        } else {
          double p_dir_alpha_m = max(rnorm(dir_alpha_m, 0.1), pow(0.1, 10));
          
          NumericVector SumS(P);
          log_with_LB(SumS, prop_prob3);
          
          double dir_lik_p, dir_lik, ratio;
          
          dir_lik_p =
            sum(     SumS* (rep(p_dir_alpha_m/P, P)-1))
            + lgamma(sum(   rep(p_dir_alpha_m/P, P)))
            - sum(   lgamma(rep(p_dir_alpha_m/P, P)));
          
          dir_lik =
            sum(     SumS*( rep(dir_alpha_m/P, P)-1))
            + lgamma(sum(   rep(dir_alpha_m/P, P)))
            - sum(   lgamma(rep(dir_alpha_m/P, P)));
          
          ratio =
            dir_lik_p
            + log(pow(p_dir_alpha_m/(p_dir_alpha_m+P), 0.5-1)
              * pow(P/(p_dir_alpha_m+P), 1-1)
              * abs(1 / (p_dir_alpha_m+P)
              - p_dir_alpha_m/pow(p_dir_alpha_m+P,2)))
            + dnorm(dir_alpha_m, p_dir_alpha_m, 0.1, true)
            - dir_lik
            - log(pow(dir_alpha_m/(dir_alpha_m+P), 0.5-1)
              * pow(P/(dir_alpha_m+P), 1-1)
              * abs(1/(dir_alpha_m+P) - dir_alpha_m/pow(dir_alpha_m+P, 2)))
            - dnorm(p_dir_alpha_m, dir_alpha_m, 0.1, true);
          
          if (ratio > log(R::runif(0,1))) {
            dir_alpha_m = p_dir_alpha_m;
          }
          
          post_dir_alpha3 = rep(dir_alpha_m/P, P) + add3;
        } // end of M.H. algorithm

        prop_prob3 = rdirichlet(1, post_dir_alpha3);



        if (iter < n_iter/10) {
          post_dir_alpha5 = rep(1.0, P2) + add5;
        } else {
          double p_dir_alpha_y1 = max(rnorm(dir_alpha_y1, 0.1), pow(0.1, 10));
          
          NumericVector SumS(P2);
          log_with_LB(SumS, prop_prob5);
          
          double dir_lik_p, dir_lik, ratio;
          
          dir_lik_p =
            sum(     SumS* (rep(p_dir_alpha_y1/P2, P2)-1))
            + lgamma(sum(   rep(p_dir_alpha_y1/P2, P2)))
            - sum(   lgamma(rep(p_dir_alpha_y1/P2, P2)));
          
          dir_lik =
            sum(     SumS*( rep(dir_alpha_y1/P2, P2)-1))
            + lgamma(sum(   rep(dir_alpha_y1/P2, P2)))
            - sum(   lgamma(rep(dir_alpha_y1/P2, P2)));
          
          ratio =
            dir_lik_p
            + log(pow(p_dir_alpha_y1/(p_dir_alpha_y1+P2), 0.5-1)
              * pow(P2/(p_dir_alpha_y1+P2), 1-1)
              * abs(1 / (p_dir_alpha_y1+P2)
              - p_dir_alpha_y1/pow(p_dir_alpha_y1+P2,2)))
            + dnorm(dir_alpha_y1, p_dir_alpha_y1, 0.1, true)
            - dir_lik
            - log(pow(dir_alpha_y1/(dir_alpha_y1+P2), 0.5-1)
              * pow(P2/(dir_alpha_y1+P2), 1-1)
              * abs(1/(dir_alpha_y1+P2) - dir_alpha_y1/pow(dir_alpha_y1+P2, 2)))
            - dnorm(p_dir_alpha_y1, dir_alpha_y1, 0.1, true);
          
          if (ratio > log(R::runif(0,1))) {
            dir_alpha_y1 = p_dir_alpha_y1;
          }
          
          post_dir_alpha5 = rep(dir_alpha_y1/P2, P2) + add5;
        } // end of M.H. algorithm

        prop_prob5 = rdirichlet(1, post_dir_alpha5);



if (iter < n_iter/10) {
  post_dir_alpha6 = rep(1.0, P2) + add6;
} else {
  double p_dir_alpha_y2 = max(rnorm(dir_alpha_y2, 0.1), pow(0.1, 10));
  
  NumericVector SumS(P2);
  log_with_LB(SumS, prop_prob6);
  
  double dir_lik_p, dir_lik, ratio;
  
  dir_lik_p =
    sum(     SumS* (rep(p_dir_alpha_y2/P2, P2)-1))
    + lgamma(sum(   rep(p_dir_alpha_y2/P2, P2)))
    - sum(   lgamma(rep(p_dir_alpha_y2/P2, P2)));
    
    dir_lik =
    sum(     SumS*( rep(dir_alpha_y2/P2, P2)-1))
      + lgamma(sum(   rep(dir_alpha_y2/P2, P2)))
      - sum(   lgamma(rep(dir_alpha_y2/P2, P2)));
      
      ratio =
      dir_lik_p
      + log(pow(p_dir_alpha_y2/(p_dir_alpha_y2+P2), 0.5-1)
              * pow(P2/(p_dir_alpha_y2+P2), 1-1)
              * abs(1 / (p_dir_alpha_y2+P2)
              - p_dir_alpha_y2/pow(p_dir_alpha_y2+P2,2)))
              + dnorm(dir_alpha_y2, p_dir_alpha_y2, 0.1, true)
              - dir_lik
              - log(pow(dir_alpha_y2/(dir_alpha_y2+P2), 0.5-1)
              * pow(P2/(dir_alpha_y2+P2), 1-1)
              * abs(1/(dir_alpha_y2+P2) - dir_alpha_y2/pow(dir_alpha_y2+P2, 2)))
              - dnorm(p_dir_alpha_y2, dir_alpha_y2, 0.1, true);
              
              if (ratio > log(R::runif(0,1))) {
                dir_alpha_y2 = p_dir_alpha_y2;
              }
              
              post_dir_alpha6 = rep(dir_alpha_y2/P2, P2) + add6;
} // end of M.H. algorithm

prop_prob6 = rdirichlet(1, post_dir_alpha6);


     
        
        if (iter < n_iter/10) {
          post_dir_alpha = rep(1.0, P3);
        } else {
          double p_dir_alpha = std::max(rnorm(dir_alpha, 0.1), pow(0.1, 10));
          
          NumericVector SumS(P3);
          log_with_LB(SumS, prop_prob);
          
          double dir_lik_p, dir_lik, ratio;
          
          dir_lik_p =
            sum(     SumS* (rep(p_dir_alpha/(P3), (P3))-1))
            + lgamma(sum(   rep(p_dir_alpha/(P3), (P3))))
            - sum(   lgamma(rep(p_dir_alpha/(P3), (P3))));
            
            dir_lik =
            sum(     SumS*( rep(dir_alpha/(P3), (P3))-1))
              + lgamma(sum(   rep(dir_alpha/(P3), (P3))))
              - sum(   lgamma(rep(dir_alpha/(P3), (P3))));
              
              ratio =
              dir_lik_p
              + log(pow(p_dir_alpha/(p_dir_alpha+(P3)), 0.5-1)
                      * abs(1 / (p_dir_alpha+(P3))
                      - p_dir_alpha/pow(p_dir_alpha+P3,2)))
                      + dnorm(dir_alpha, p_dir_alpha, 0.1, true)
                      - dir_lik
                      - log(pow(dir_alpha/(dir_alpha+(P3)), 0.5-1)
                      * abs(1/(dir_alpha+(P3)) - dir_alpha/pow(dir_alpha+(P3), 2)))
                      - dnorm(p_dir_alpha, dir_alpha, 0.1, true);
                      
                      if (ratio > log(R::runif(0,1))) {
                        dir_alpha = p_dir_alpha;
                      }
                      
                      post_dir_alpha = rep(dir_alpha/(P3), (P3));
        } // end of M.H. algorithm

        
        // M.H. algorithm for the inclusion probabilities
        {
          double dir_lik_p, dir_lik, ratio;
          
          NumericVector add11(add1.length() + 3);
          add11(0) = add4(0);
          add11(1) = add4(1);
          add11(2) = max(add4(2), add2(0));
          add11[Rcpp::Range(3, add1.length()+2)] = add1;

          NumericVector add22(add2.length() + 2);
          add22(0) = add4(0);
          add22(1) = add4(1);
          add22[Rcpp::Range(2, add1.length()+1)] = add2;
          
          NumericVector p_prop_prob = rdirichlet(1, add4 + post_dir_alpha + add11 + add22);
          NumericVector log_p_prop_prob(P+3), log_prop_prob(P+3);
          log_with_LB(log_p_prop_prob, p_prop_prob);
          log_with_LB(log_prop_prob,   prop_prob);
          
            dir_lik_p = sum(add1) * log(1/(1-p_prop_prob(0)-p_prop_prob(1)-p_prop_prob(2)))
              + sum(add2) * log(1/(1-p_prop_prob(0)-p_prop_prob(1)))
              + (add4(0) + post_dir_alpha(0) - 1.0) * log_p_prop_prob(0)
              + (add4(1) + post_dir_alpha(1) - 1.0) * log_p_prop_prob(1)
              + (add4(2) + add2(0) + post_dir_alpha(2) - 1.0) * log_p_prop_prob(2)
              + sum((add2[Rcpp::Range(1,P)] + add4[Rcpp::Range(3,P + 2)] + add1 + post_dir_alpha[Rcpp::Range(3,P + 2)] - 1.0) * log_p_prop_prob[Rcpp::Range(3,P + 2)]);
              
              
              dir_lik = sum(add1) * log(1/(1-prop_prob(0)-prop_prob(1)-prop_prob(2)))
                + sum(add2) * log(1/(1-prop_prob(0)-prop_prob(1)))
                + (add4(0) + post_dir_alpha(0) - 1.0) * log_prop_prob(0)
                + (add4(1) + post_dir_alpha(1) - 1.0) * log_prop_prob(1)
                + (add4(2) + add2(0) + post_dir_alpha(2) - 1.0) * log_prop_prob(2)
                + sum((add2[Rcpp::Range(1,P)] + add4[Rcpp::Range(3,P + 2)] + add1 + post_dir_alpha[Rcpp::Range(3,P + 2)] - 1.0) * log_prop_prob[Rcpp::Range(3,P + 2)]);
                

          ratio = dir_lik_p
            + sum((add4 + post_dir_alpha + add11 + add22 - 1.0) * log_prop_prob)
            - dir_lik
            - sum((add4 + post_dir_alpha + add11 + add22 - 1.0) * log_p_prop_prob);
                
          if (ratio > log(R::runif(0,1))) {
            prop_prob = clone(p_prop_prob);
          }
        }
        
        
        
        // Sampling E[Y(1)-Y(0)]
        if (iter > burn_in)
        {
            if (thin_count < thin)
            {
                thin_count++;
            }
            else
            {
                thin_count = 1;
                predicted_M1 (_, post_sample_idx) = M1*msd+mshift;
                predicted_M0 (_, post_sample_idx) = M0*msd+mshift;
                predicted_Y11 (_, post_sample_idx) = Y11;
                predicted_Y00 (_, post_sample_idx) = Y00;
                predicted_Y10 (_, post_sample_idx) = Y10;
                predicted_zeta(_, post_sample_idx) = rowSums(clone(Tree5_pred))*ysd;
                predicted_d(_, post_sample_idx) = rowSums(clone(Tree6_pred))*ysd/msd;
                predicted_tau(_, post_sample_idx) = rowSums(clone(Tree3_pred))*msd;
                post_sample_idx++;
                
                IntegerVector ind_temp = ifelse(add1 > 0.0, 1, 0);
                IntegerVector ind_temp1 = ifelse(add2 > 0.0, 1, 0);
                IntegerVector ind_temp2 = ifelse(add4 > 0.0, 1, 0);

                ind(ind_idx, _) = ind_temp; // indicator of whether confounders are included
                ind1(ind_idx, _) = ind_temp1; // indicator of whether confounders are included

                ind2(ind_idx, _) = ind_temp2; // indicator of whether confounders are included

                ind_idx++;

            }
        }

        Rcpp::checkUserInterrupt(); // check for break in R
    } // end of MCMC iterations

    List L = List::create(
        Named("predicted_Y00") = predicted_Y00,
        Named("predicted_Y11") = predicted_Y11,
        Named("predicted_Y10") = predicted_Y10,
        Named("predicted_zeta") = predicted_zeta,
        Named("predicted_d") = predicted_d,
        Named("predicted_tau") = predicted_tau,
        Named("ind") = ind,
        Named("ind1") = ind1,
        Named("ind2") = ind2,
        Named("predicted_M0") = predicted_M0,
        Named("predicted_M1") = predicted_M1
    );

    return L;
}
