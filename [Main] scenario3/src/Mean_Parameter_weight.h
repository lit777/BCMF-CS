// Mean_Parameter_weight.h
#pragma once

#include "decision_tree.h"
#include "functions.h"
#include "subset.h"

using namespace Rcpp;

// Sampling mean parameters
void DecisionTree::Mean_Parameter_weight(
    NumericMatrix& Tree,
    const NumericVector& sigma2,
    double sigma_mu,
    const NumericVector& R,
    const IntegerMatrix& Obs_list
) {

    NumericVector T(this->n);
    IntegerVector tnode_idx = this->terminal_nodes();
    for (int i = 0; i < tnode_idx.length(); i++)
    {
        int begin = this->begin(tnode_idx(i)); // match index
        int end   = this->end(tnode_idx(i));

        IntegerVector Obs_ind(end - begin + 1);
        if (begin == end)
            Obs_ind = Obs_list(begin, this->id);
        else
            subset_by_range(Obs_ind, Obs_list, begin, end, this->id); // Obs_ind = Obs_list[begin:end, ]

        double var_R = var_by_idx_weight(sigma2, Obs_ind, sigma_mu);
        double sum_R = sum_by_idx_weight(R, sigma2, Obs_ind);
        double Var  = 1 / var_R;
        double Mean = Var * sum_R;

        this->MU(tnode_idx(i)) = R::rnorm(Mean, sqrt(Var));
        T[Obs_ind] = this->MU(tnode_idx(i));
    }
    Tree(_, this->id) = clone(T);
}
