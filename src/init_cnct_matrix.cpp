// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
#define ARMA_64BIT_WORD 1
#include "progress.hpp"
using namespace Rcpp;
using namespace arma;

// Developer notes...
// - These functions uses index by row number extensively, as such be careful with the R vs C++ index start differences of 1 vs 0.
// - It is important to include explicit namespace (i.e. Rcpp:: or you may get some weird behaviour and not know what the bug is...)
// - Converting from Rcpp and arma can be a pain, try to stick to one data type when possible...
// - In pkg, the [[Rcpp::depends(RcppArmadillo)]] and [[Rcpp::depends(RcppProgress)]] not used, not necessary once you have them as part of the Imports...
// - The defined ARMA must come after the include... when in a package

// Internal function: *which_cpp_and*
// Cpp version of a which() operation to compare T & T logic of two logical vectors
NumericVector which_cpp_and(Rcpp::LogicalVector x, Rcpp::LogicalVector y, Rcpp::LogicalVector z, bool r_style = false) {

  // Determine length of which vector
  int lengthV = 0;
  int n = x.size();

  for(int i = 0; i < n; i++) {
    if(x[i] == true && y[i] == true && z[i] == true) {
      lengthV += 1;
    }
  }

  // Use length to assign out...
  NumericVector out(lengthV);

  // Populate the out vector based on && logic
  for(int i = 0, j = 0; i < n; i++) {
    if(x[i] == true && y[i] == true && z[i] == true) {
      out[j++] = i;
    }
  }

  // Determine style/format of output indexing
  if(r_style == true) {
    return out + 1;
  } else {
    return out;
  }
}

// Internal function: *matchCpp*
// Cpp version of using matching for extracting prob vector
NumericVector matchCpp(Rcpp::IntegerVector x, Rcpp::IntegerVector y, Rcpp::NumericVector prob) {
  int n = x.size();
  Rcpp::NumericVector match_index(n);

  Rcpp::IntegerVector matches = match(x, y)-1; // Find index matches from left within right
  for(int i = 0; i < n; i++) {
    if(IntegerVector::is_na(matches[i])) {
      match_index[i] = NA_REAL;
    } else {
      match_index[i] = prob[matches[i]];
    }
  }
  return match_index;
}

//' Create adjacency matrix for contacts
//'
//' Fill sparse adjacency matrix based upon a specific contact probability structure and contact limits.
//'
//' Used primarily as internal component to generate the population S4 object's contact structure. C++ was used to improve efficiency; works well for
//' up to 10-100k agents before noticeably slow. Further performance gains could be provided if agents are split into blocks (e.g. 100k rows) and computer in parallel; however,
//' this will have  a large assumption that each block will not have contacts between each other unless additional operations are done to add contacts after each block completes.
//'
//' @param adj_matrix Adjacency matrix, typically of type SparseMatrix; most likely will be 'blank' to be initialized/filled by this function.
//' @param contact_matrix Matrix of probabilities for each type of contact (e.g. probabilities of 3 ages groups and 2 genders interacting, for a total of 36 cross combinations).
//' @param category_id Unique numeric value for the type of contact group (e.g. Female of age 20).
//' @param num_contacts Integer for number of contacts for particular agent.
//' @param sample_order Integer vector to determine order to fill the adjacency matrix (e.g. by default will go in order of agent list, recommended to provide a random order to ensure contacts fill without preference to a categorical sorting).
//' @param display_progress Boolean value to determine if a progress bar is displayed (useful for when creating large networks, such as those >10,000 agents).
//'
//' @export
// [[Rcpp::export]]
arma::sp_mat init_adj_matrix_cpp(arma::sp_mat adj_matrix,
                                 Rcpp::NumericMatrix contact_matrix,
                                 Rcpp::IntegerVector category_id,
                                 Rcpp::NumericVector num_contacts,
                                 Rcpp::IntegerVector sample_order,
                                 bool display_progress = true) {

  // Initialize data objects from inputs
  Rcpp::IntegerVector category_id_cpp = category_id - 1;
  Rcpp::NumericVector index(num_contacts.length()); // Should equate to number of inputs
  Rcpp::IntegerVector sample_order_c = clone(sample_order);
  sample_order_c = sample_order_c - 1; // Was having issues doing sampling directly here from seq...

  // Progress bar
  Progress p(sample_order_c.length(), display_progress);

  // ##########  MAIN LOOP ##########
  // For every step of index (i.e. entries) update the adj matrix
  for(int i = 0; i < sample_order_c.length(); i++) {

    // Check if progress aborted... returns blank matrix if so... (remove if really slow?)
    if (Progress::check_abort())  return arma::sp_mat(sample_order_c.length(),sample_order_c.length());

    // Calculate number of contacts among all in adj matrix
    NumericVector tmp_compare = as<NumericVector>(Rcpp::wrap(arma::rowvec(arma::sum(adj_matrix,0))));

    // Determine remaining entities with vacant contact slots to fill
    int contacts2fill = num_contacts[sample_order_c[i]] - tmp_compare[sample_order_c[i]];

    // Conditional changes based on filling orders (skips all calculations if no other contacts needed based on prior steps)
    if (contacts2fill > 0) {

      // Set data objects for each index step (to join in probabilities by ID for person)
      Rcpp::NumericVector tmp_prob = contact_matrix(_ , category_id_cpp[sample_order_c[i]]);
      Rcpp::NumericVector contact_prob_v = matchCpp(category_id_cpp, Rcpp::seq(1, contact_matrix.nrow())-1, tmp_prob);
      Rcpp::IntegerVector prob_index = Rcpp::seq_along(category_id);
      prob_index = prob_index - 1; // -1 to refer to C++ version of rows

      // Determine which contacts are valid to select from (improve by is ALL probs are 0, then skip?)
      index = which_cpp_and( (num_contacts - tmp_compare) > 0, prob_index != sample_order_c[i], contact_prob_v > 0, false); // 3 conditions (ensure contacts left, not itself, and non zero prob)

      if(index.length() > 0){
        if (index.length() == 1) {
          prob_index = as<IntegerVector>(prob_index[index]);
        } else {
          if(contacts2fill > index.length()) {contacts2fill = index.length();}
          Rcpp::NumericVector sample_temp = Rcpp::sample(Rcpp::as<NumericVector>(prob_index[index]), contacts2fill, false, Rcpp::as<NumericVector>(contact_prob_v[index])); // Need to coerce to a new numeric inline to avoid issues with the subset index
          prob_index = Rcpp::as<IntegerVector>(prob_index[sample_temp]);
        }

        // Update the adj matrix
        for(int j = 0; j < prob_index.length(); j++) {
          adj_matrix(prob_index[j], sample_order_c[i]) = 1;
          adj_matrix(sample_order_c[i], prob_index[j]) = 1;
        }
      }
    }
    p.increment(); // update progress
  }
  return adj_matrix;
}


