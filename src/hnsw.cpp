//  RcppHNSW -- Rcpp bindings to HNSW library for Approximate Nearest Neighbours
//
//  Copyright (C) 2018  James Melville
//
//  This file is part of RcppHNSW
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <algorithm>
#include <iostream>
#include "hnswlib.h"
#include <Rcpp.h>


template<typename dist_t, typename Distance>
class Hnsw {
public:

  // dim - length of the vectors being added
  // maxElements - size of the data being added
  // efConstruction - controls the quality of the graph. Higher values lead to improved recall at the
  //  expense of longer build time. Suggested range: 100-2000 (default: 200).
  // M - Controls maximum number of neighbors in the zero and above-zero layers. Higher values lead to
  //  better recall and shorter retrieval times, at the expense of longer indexing time.
  //  Suggested range: 5-100 (default: 16).
  Hnsw(const int dim, const size_t maxElements, const size_t M = 16, const size_t efConstruction = 200) :
    dim(dim), cur_l(0)
  {
    l2space = new Distance(dim);
    appr_alg = new hnswlib::HierarchicalNSW<dist_t>(l2space, maxElements, M, efConstruction);
  }

  Hnsw(const int dim, const std::string path_to_index) :
  dim(dim), cur_l(0)
  {
    l2space = new Distance(dim);
    appr_alg = new hnswlib::HierarchicalNSW<dist_t>(l2space, path_to_index);
    cur_l = appr_alg->cur_element_count;
  }

  void addItem(const std::vector<dist_t> dv)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());
    appr_alg->addPoint(&fv[0], (size_t) cur_l);
    ++cur_l;
  }

  std::vector<hnswlib::labeltype> getNNs(const std::vector<dist_t> dv, size_t k)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype >> result =
      appr_alg->searchKnn(&fv[0], k);

    if (result.size() != k) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    std::vector<hnswlib::labeltype> items;
    items.reserve(k);
    for (size_t i = 0; i < k; i++) {
      auto &result_tuple = result.top();
      items.push_back(result_tuple.second);
      result.pop();
    }
    std::reverse(items.begin(), items.end());

    return items;
  }

  Rcpp::List getNNsList(const std::vector<dist_t> dv, size_t k,
                        bool include_distances)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
      appr_alg->searchKnn(&fv[0], k);

    if (result.size() != k) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    std::vector<hnswlib::labeltype> items;
    items.reserve(k);

    if (include_distances) {
      std::vector<dist_t> distances;
      distances.reserve(k);

      for (size_t i = 0; i < k; i++) {
        auto &result_tuple = result.top();
        distances.push_back(result_tuple.first);
        items.push_back(result_tuple.second);
        result.pop();
      }

      std::reverse(distances.begin(), distances.end());
      std::reverse(items.begin(), items.end());

      return Rcpp::List::create(
        Rcpp::Named("item") = items,
        Rcpp::Named("distance") = distances);
    }
    else {
      for (size_t i = 0; i < k; i++) {
        auto &result_tuple = result.top();
        items.push_back(result_tuple.second);
        result.pop();
      }

      std::reverse(items.begin(), items.end());

      return Rcpp::List::create(
        Rcpp::Named("item") = items);
    }
  }

  void callSave(const std::string path_to_index) {
    appr_alg->saveIndex(path_to_index);
  }

  ~Hnsw()
  {
    delete l2space;
    delete appr_alg;
  }

  int dim;
  hnswlib::labeltype cur_l;
  hnswlib::HierarchicalNSW<dist_t> *appr_alg;
  hnswlib::SpaceInterface<float> *l2space;
};

typedef Hnsw<float, hnswlib::L2Space> HnswL2;

RCPP_EXPOSED_CLASS_NODECL(HnswL2)
RCPP_MODULE(HnswL2) {
  Rcpp::class_<HnswL2>("HnswL2")
  .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, ef, M")
  .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
  .method("addItem",    &HnswL2::addItem,    "add item")
  .method("save",       &HnswL2::callSave,   "save index to file")
  .method("getNNs",     &HnswL2::getNNs,     "retrieve Nearest Neigbours given vector")
  .method("getNNsList", &HnswL2::getNNsList, "retrieve Nearest Neigbours given vector")
  ;
}
