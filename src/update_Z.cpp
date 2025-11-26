#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_Z_ucar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube Z = current_sample["Z"];
  mat sig2 = current_sample["sig2"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  mat tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
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
        double rate_diff = (lambda(reg, grp, time) - beta(island_id(reg), grp, time));
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
  current_sample["Z"] = Z;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_Z_mcar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube Z = current_sample["Z"];
  cube G = current_sample["G"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  mat tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
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
  cube rate_diff = lambda - get_regs(beta, island_id);
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
  current_sample["Z"] = Z;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_Z_mstcar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube Z = current_sample["Z"];
  cube G = current_sample["G"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  vec rho = current_sample["rho"];
  vec tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
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
  cube rate_diff = lambda - get_regs(beta, island_id);
  for (uword reg = 0; reg < num_region; reg++) {
    mat nZm = mean(get_regs(Z, adjacency[reg]), 0);
    if (num_time == 1) nZm = nZm.t();
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
  current_sample["Z"] = Z;
  RSTr_obj["current_sample"] = current_sample;
}