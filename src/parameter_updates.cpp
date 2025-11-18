#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat update_Ag(List inits, List priors) {
  mat Ag = inits["Ag"];
  cube G = inits["G"];
  mat Ag_scale = priors["Ag_scale"];
  double G_df = priors["G_df"];
  double Ag_df = priors["Ag_df"];
  uword num_time  = G.n_slices;
  uword num_group = G.n_rows;
  mat Ag_covar(num_group, num_group, fill::zeros);
  Ag_covar += inv(Ag_scale);
  for (uword time = 0; time < num_time; time++) {
    Ag_covar += inv(G.slice(time));
  }
  Ag = rwish(num_time * G_df + Ag_df, inv(Ag_covar));
  return Ag;
}

//[[Rcpp::export]]
arma::cube update_beta_ucar(List inits, List spatial_data) {
  cube beta = inits["beta"];
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  mat tau2 = inits["tau2"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = theta - Z;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region(isl).n_elem;
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sd_beta = sqrt(tau2(grp, time) / num_island_region);
        double mean_beta = mean(get_subregs(tmZ, island_region(isl), grp, time));
        beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
      }
    }
  }
  return beta;
}

//[[Rcpp::export]]
arma::cube update_beta_mstcar(List inits, List spatial_data) {
  cube beta = inits["beta"];
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  vec tau2 = inits["tau2"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = theta - Z;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region(isl).n_elem;
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sd_beta = sqrt(tau2(grp) / num_island_region);
        double mean_beta = mean(get_subregs(tmZ, island_region(isl), grp, time));
        beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
      }
    }
  }
  return beta;
}

//[[Rcpp::export]]
arma::cube update_beta_ucar_restricted(List inits, List spatial_data, List params) {
  cube beta = inits["beta"];
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  mat tau2 = inits["tau2"];
  mat sig2 = inits["sig2"];
  field<uvec> island_region = spatial_data["island_region"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = theta - Z;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region(isl).n_elem;
    mat var_t = tau2 + (tau2 + sig2) / m0;
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sd_beta = sqrt(tau2(grp, time) / num_island_region);
        double mean_beta = mean(get_subregs(tmZ, island_region(isl), grp, time));
        if (method == "binomial") {
          double pi_beta = pow(A(grp, time) - 1, 2) + 4 * (A(grp, time) - 1 / var_t(grp, time));
          double beta_thres = ((1 - A(grp, time)) + sqrt(pi_beta)) / 2;
          beta_thres = log(beta_thres / (1 - beta_thres));
          beta_thres = (beta_thres < 0) ? 0 : beta_thres;
          double beta_max = R::pnorm(beta_thres, mean_beta, sd_beta, true, false);
          if (beta_max > 0) {
            double u = R::runif(0, beta_max);
            beta(isl, grp, time) = R::qnorm(u, mean_beta, sd_beta, true, false);
          }
        } else if (method == "poisson") {
          beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
        }
      }
    }
    
  }
  return beta;
}

//[[Rcpp::export]]
arma::cube update_G_mcar(List inits, List spatial_data, List priors) {
  cube G = inits["G"];
  cube Z = inits["Z"];
  double G_df = priors["G_df"];
  mat G_scale = priors["G_scale"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_time = Z.n_slices;
  double df_G = num_region - num_island + G_df;
  for (uword time = 0; time < num_time; time++) {
    mat scale_G = G_scale;
    for (uword reg = 0; reg < num_region; reg++) {
      vec Zit = get_grp(Z, reg, time);
      vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      scale_G += adjacency(reg).n_elem * Zit * Zit.t() - sum_adj * Zit.t();
    }
    G.slice(time) = riwish(df_G, scale_G);
  }
  return G;
}

//[[Rcpp::export]]
arma::cube update_G_mstcar(List inits, List spatial_data, List priors) {
  cube G = inits["G"];
  cube Z = inits["Z"];
  mat Ag = inits["Ag"];
  vec rho = inits["rho"];
  double G_df = priors["G_df"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  cube Ags(num_group, num_group, num_time, fill::zeros);
  Ags.each_slice() += Ag;
  vec r  = rho;
  vec sr = sqrt(1 - pow(rho, 2));
  for (uword reg = 0; reg < num_region; reg++) {
    double num_adj = adjacency[reg].n_elem;
    mat Zmikt = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
    rowvec Zt = get_grp(Z, reg, 0).t();
    Ags.slice(0) += num_adj * Zmikt.col(0) * Zt;
    for (uword time = 1; time < num_time; time++) {
      vec Zt  = 1 / sr % get_grp(Z, reg, time);
      vec Ztl = r / sr % get_grp(Z, reg, time - 1);
      vec Zm  = 1 / sr % Zmikt.col(time);
      vec Zml = r / sr % Zmikt.col(time - 1);
      Ags.slice(time) += num_adj * ((Zm - Zml) * (Zt - Ztl).t());
    }
  }
  for (uword time = 0; time < num_time; time++) {
    G.slice(time) = riwish((num_region - num_island) + G_df, Ags.slice(time));
  }
  return G;
}

//[[Rcpp::export]]
arma::rowvec update_rho(List inits, List spatial_data, List priors, arma::vec& r_accept) {
  vec rho = inits["rho"];
  cube G = inits["G"];
  cube Z = inits["Z"];
  double rho_a = priors["rho_a"];
  double rho_b = priors["rho_b"];
  vec rho_sd = priors["rho_sd"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  vec logit_rho    = log(rho / (1 - rho));
  vec rand         = Rcpp::rnorm(num_group, 0, 1);
  vec expit_rho    = rand % rho_sd + logit_rho;
  vec rho_star_0   = 1 / (1 + exp(-expit_rho));
  vec r(num_group, fill::zeros);
  vec ra(num_group, fill::zeros);
  vec rb(num_group, fill::zeros);
  vec rc(num_group, fill::zeros);
  cube Zm(num_region, num_group, num_time);
  for (uword reg = 0; reg < num_region; reg++) {
    Zm.row(reg) = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
  }
  for (uword grp = 0; grp < num_group; grp++) {
    vec rho_star = rho;
    rho_star[grp] = rho_star_0[grp];
    ra[grp] = (1 - pow(rho[grp], 2)) / (1 - pow(rho_star[grp], 2));
    field<mat> Sein_rho = Sig_eta_i(G, rho);
    field<mat> Sein_rho_star = Sig_eta_i(G, rho_star);
    field<mat> Sein_diff(num_time, num_time);    
    for (uword time1 = 0; time1 < num_time; time1++) {
      uword time1_l = (time1 == 0) ? 0 : time1 - 1;
      uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
      for (uword time2 = time1_l; time2 < time2_u; time2++) {
        Sein_diff(time2, time1) = Sein_rho_star(time2, time1) - Sein_rho(time2, time1);
      }
    }
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword time1 = 0; time1 < num_time; time1++) {
        uword time2_l = (time1 == 0) ? 0 : time1 - 1;
        uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
        for (uword time2 = time2_l; time2 < time2_u; time2++) {
          mat rb_ik = get_grp(Z, reg, time2).t() * Sein_diff(time2, time1) * get_grp(Zm, reg, time1);
          rb[grp] += adjacency[reg].n_elem * rb_ik[0] / 2;
        }
      }
    }
    rc[grp] = pow(rho_star[grp] / rho[grp], rho_a) * pow((1 - rho_star[grp]) / (1 - rho[grp]), rho_b);
    r[grp] = exp((num_region - num_island) * (num_time - 1) / 2 * log(ra[grp]) - rb[grp]) * rc[grp];
    if (r[grp] >= R::runif(0, 1)) {
      r_accept[grp]++;
      rho[grp] = rho_star[grp];
    }
  }
  return rho.t();
}

//[[Rcpp::export]]
arma::mat update_sig2_ucar(List inits, List spatial_data, List priors) {
  mat sig2 = inits["sig2"];
  cube Z = inits["Z"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < num_region; reg++) {
        sum_adj += Z(reg, grp, time) * sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (num_region - num_island) / 2 + sig_a;
      double sig_scale = 1 / ((sum(pow(get_row(Z, grp, time), 2) % num_adj) - sum_adj) / 2 + sig_b);
      sig2(grp, time) = 1 / R::rgamma(sig_shape, sig_scale);
    }
  }
  return sig2;
}

//[[Rcpp::export]]
arma::mat update_sig2_ucar_restricted(List inits, List spatial_data, List priors, List params) {
  mat sig2 = inits["sig2"];
  cube Z = inits["Z"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec num_island_region = spatial_data["num_island_region"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < num_region; reg++) {
        sum_adj += Z(reg, grp, time) * sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (num_region - num_island) / 2 + sig_a;
      double sig_scale = 1 / ((sum(pow(get_row(Z, grp, time), 2) % num_adj) - sum_adj) / 2 + sig_b);
      double sig_thres = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % num_island_region / num_region);
        pi = exp(pi) / (1 + exp(pi));
        sig_thres = (1 / ((A(grp, time) + pi) * (1 - pi)) - tau2(grp, time) * (1 + 1 / m0)) * m0;
      } else if (method == "poisson") {
        sig_thres = (log(1 / A(grp, time) + 1) - tau2(grp, time) * (1 + 1 / m0)) * m0;
      }
      sig_thres = (sig_thres < 0) ? 0 : sig_thres;
      double u = R::runif(0, R::pgamma(1 / sig_thres, sig_shape, sig_scale, true, false));
      sig2(grp, time) = 1 / R::qgamma(u, sig_shape, sig_scale, true, false);
    }
  }
  return sig2;
}

//[[Rcpp::export]]
arma::mat update_tau2_ucar(List inits, List spatial_data, List priors) {
  mat tau2 = inits["tau2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  cube Z = inits["Z"];
  uvec island_id = spatial_data["island_id"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2);
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double tau_scale = 1 / (sum(get_row(square_resid, grp, time)) / 2 + tau_b);
      tau2(grp, time) = 1 / R::rgamma(tau_shape, tau_scale);
    }
  }
  return tau2;
}

//[[Rcpp::export]]
arma::mat update_tau2_ucar_restricted(List inits, List spatial_data, List priors, List params) {
  mat tau2 = inits["tau2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  cube Z = inits["Z"];
  mat sig2 = inits["sig2"];
  uvec num_island_region = spatial_data["num_island_region"];
  uvec island_id = spatial_data["island_id"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2);
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double tau_scale = 1 / (sum(get_row(square_resid, grp, time)) / 2 + tau_b);
      double tau_thres = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % num_island_region / num_region);
        pi = exp(pi) / (1 + exp(pi));
        tau_thres = (1 / ((A(grp, time) + pi) * (1 - pi)) - sig2(grp, time) / m0) / (1 + 1 / m0);
      } else if (method == "poisson") {
        tau_thres = (log(1 / A(grp, time) + 1) - sig2(grp, time) / m0) / (1 + 1 / m0);
      }
      tau_thres = (tau_thres < 0) ? 0 : tau_thres;
      double u = R::runif(0, R::pgamma(1 / tau_thres, tau_shape, tau_scale, true, false));
      tau2(grp, time) = 1 / R::qgamma(u, tau_shape, tau_scale, true, false);
    }
  }
  return tau2;
}

//[[Rcpp::export]]
arma::mat update_tau2_mstcar(List inits, List spatial_data, List priors) {
  mat tau2 = inits["tau2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  cube Z = inits["Z"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  double tau_shape = num_region * num_time / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2) / 2;
  for (uword grp = 0; grp < num_group; grp++) {
    double tau_scale = 1 / (accu(square_resid.col(grp)) + tau_b);
    tau2[grp] = 1 / R::rgamma(tau_shape, tau_scale);
  }
  return tau2;
}

//[[Rcpp::export]]
arma::cube update_theta(List inits, List spatial_data, List priors, List params, List data, arma::cube& t_accept) {
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  cube Y = data["Y"];
  cube n = data["n"];
  cube theta_sd = priors["theta_sd"];
  uvec island_id = spatial_data["island_id"];
  String method = params["method"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  for (uword time = 0; time < num_time; time++) {
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword grp = 0; grp < num_group; grp++) {
        double theta_star = R::rnorm(theta(reg, grp, time), theta_sd(reg, grp, time));
        double rk1 = Y(reg, grp, time) * (theta_star - theta(reg, grp, time));
        double rk2 = 0;
        if (method == "binomial") {
          rk2 = n(reg, grp, time) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg, grp, time))));
        } 
        if (method == "poisson") {
          rk2 = n(reg, grp, time) * (exp(theta_star) - exp(theta(reg, grp, time)));
        }
        double rk3a = pow(theta_star            - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk3b = pow(theta(reg, grp, time) - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk   = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
        if (rk >= R::runif(0, 1)) {
          t_accept(reg, grp, time)++;
          theta(reg, grp, time) = theta_star;
        }
      }
    }
  }
  return theta;
}

//[[Rcpp::export]]
arma::cube update_Z_ucar(List inits, List spatial_data) {
  cube Z = inits["Z"];
  mat sig2 = inits["sig2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sum_adj = sum(get_subregs(Z, adjacency(reg), grp, time));
        double rate_diff = (theta(reg, grp, time) - beta(island_id(reg), grp, time));
        double var_Z  = 1 / (1 / tau2(grp, time) + num_adj(reg) / sig2(grp, time));
        double mean_Z = var_Z * (rate_diff / tau2(grp, time) + sum_adj / sig2(grp, time));
        Z(reg, grp, time) = R::rnorm(mean_Z, sqrt(var_Z));
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region(isl)), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}

//[[Rcpp::export]]
arma::cube update_Z_mcar(List inits, List spatial_data) {
  cube Z = inits["Z"];
  cube G = inits["G"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  field<mat> Z_cov(max(num_adj) + 1, num_time);
  field<mat> Z_coveig(max(num_adj) + 1, num_time);
  cube rate_diff = theta - get_regs(beta, island_id);
  vec unique_num_adj = unique(num_adj);
  for (uword time = 0; time < num_time; time++) {
    vec taut = tau2.col(time);
    mat Gt = G.slice(time);
    for (uword count : unique_num_adj) {
      Z_cov(count, time) = inv(diagmat(1 / taut) + count * Gt);
      Z_coveig(count, time) = geteig(Z_cov(count, time));
    }
    for (uword reg = 0; reg < num_region; reg++) {
      vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      vec Z_mean = Z_cov(num_adj(reg), time) * (get_grp(rate_diff, reg, time) / taut + Gt * sum_adj);
      vec Z_new  = cpp_rmvnorm(Z_mean, Z_coveig(num_adj(reg), time));
      for (uword grp = 0; grp < num_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region[isl]), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}

//[[Rcpp::export]]
arma::cube update_Z_mstcar(List inits, List spatial_data) {
  cube Z = inits["Z"];
  cube G = inits["G"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  vec rho = inits["rho"];
  vec tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region  = Z.n_rows;
  uword num_group   = Z.n_cols;
  uword num_time    = Z.n_slices;
  uword num_island  = island_region.n_elem;
  field<mat> Sein   = Sig_eta_i(G, rho);
  field<mat> SeSein = Sig_eta(Sein);
  field<mat> Z_cov(num_time, max(num_adj) + 1);
  field<mat> Z_coveig(num_time, max(num_adj) + 1);
  vec unique_num_adj = unique(num_adj);
  for (uword time = 0; time < num_time; time++) {
    for (uword count : unique_num_adj) {
      Z_cov   (time, count) = inv(diagmat(1 / tau2) + count * Sein(time, time));
      Z_coveig(time, count) = geteig(Z_cov(time, count));
    }
  }
  cube rate_diff = theta - get_regs(beta, island_id);
  for (uword reg = 0; reg < num_region; reg++) {
    mat nZm = mean(get_regs(Z, adjacency[reg]), 0);
    for (uword time = 0; time < num_time; time++) {
      vec muZp = nZm.col(time);
      if (time > 0) {
        muZp += SeSein(time, time - 1) * (nZm.col(time - 1) - get_grp(Z, reg, time - 1));
      }
      if (time < num_time - 1) {
        muZp += SeSein(time, time + 1) * (nZm.col(time + 1) - get_grp(Z, reg, time + 1));
      }
      mat Z_mean = Z_cov(time, num_adj[reg]) * (get_grp(rate_diff, reg, time) / tau2 + (num_adj[reg] * Sein(time, time) * muZp));
      vec Z_new  = cpp_rmvnorm(Z_mean, Z_coveig(time, num_adj[reg]));
      for (uword grp = 0; grp < num_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region[isl]), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}