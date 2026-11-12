// GROW_first_weight.h
#pragma once

#include "decision_tree.h"
#include "functions.h"
#include "subset.h"

using namespace Rcpp;

// Fun. of grow (first) alteration
void DecisionTree::GROW_first_weight(
    const NumericMatrix& xpred,
    const NumericVector xcut[],
    const NumericVector& sigma2, double sigma_mu,
    const NumericVector& R,
    IntegerMatrix& Obs_list,

    double p_prune, double p_grow,
    double alpha, double beta,
    const NumericVector& prop_prob
) {
    // GROW FIRST
    const int n = xpred.nrow();
    const int P = xpred.ncol();
    const IntegerVector row_idx = Rcpp::Range(0, n - 1); // Create vector s.t. [0, 1, ..., n - 1]
    int prop_pred, prop_rule;
    double value;
    prop_pred = sample(P, 1, false, prop_prob)(0) - 1;   // pick a predictor and match index with cpp
    prop_rule = sample(xcut[prop_pred].length() - 1, 1)(0); // sample from U(2, length(xcut(prop_pred)))
    value = xcut[prop_pred](prop_rule);                     // value for separation
    
    IntegerVector R_L = row_idx[xpred(_, prop_pred) <  value];
    IntegerVector R_R = row_idx[xpred(_, prop_pred) >= value];

    NumericVector xpred_prop_pred = xpred(_, prop_pred);

    // Transition ratio (log scale)
    double TRANS = log(p_prune) - log(std::max(prop_prob(prop_pred), 0.0)) + log(xcut[prop_pred].length() - 1) - log(p_grow);

    // Likelihood ratio (log scale)
    int nlL = R_L.length();
    int nlR = R_R.length();
    double sum_R_L = sum_by_idx_weight(R, sigma2, R_L);
    double sum_R_R = sum_by_idx_weight(R, sigma2, R_R);
    double sum_R = sum_by_idx_weight(R, sigma2, row_idx);
    
    double var_R_L = var_by_idx_weight(sigma2, R_L, sigma_mu);
    double var_R_R = var_by_idx_weight(sigma2, R_R, sigma_mu);
    double var_R = var_by_idx_weight(sigma2, row_idx, sigma_mu);
    
    double LH = 0.5 * (1 / var_R_L) * pow(sum_R_L, 2) + 0.5 * (1 / var_R_R) * pow(sum_R_R, 2) - 0.5 * (1 / var_R) * pow(sum_R, 2) + 0.5 * log(var_R) - 0.5 * log(var_R_L) - 0.5 * log(var_R_R); 

    // Structure ratio (log scale)
    int d = 0;
    double STR = log(alpha) + 2 * log((1 - alpha / pow(2 + d, beta))) - log(pow(1 + d, beta) - alpha) + log(std::max(prop_prob(prop_pred), 0.0)) - log(xcut[prop_pred].length() - 1);

    double r = TRANS + LH + STR;

    if (r > log(R::runif(0, 1)))
    {
        // New tree structure
        this->Split    = prop_pred;
        this->Terminal = rep(0, 1);
        this->Value    = prop_rule;
        this->position = append(this->position, 2, 3);
        this->parent   = append(this->parent, 1, 1);
        this->Terminal = append(this->Terminal, 1, 1);
        this->Split    = append(this->Split, NA_INTEGER, NA_INTEGER);
        this->Value    = append(this->Value, NA_INTEGER, NA_INTEGER);
        this->MU       = append(this->MU, NA_REAL, NA_REAL);
        this->begin    = append(this->begin, 0, nlL);
        this->end      = append(this->end, nlL - 1, n - 1);

        // update obs
        R_L.sort();
        R_R.sort();
        for (int i = 0; i < nlL; i++)
            Obs_list(i, this->id) = R_L(i);
        for (int i = 0; i < nlR; i++)
            Obs_list(i + nlL, this->id) = R_R(i);
    }
} // end of GROW_first()
