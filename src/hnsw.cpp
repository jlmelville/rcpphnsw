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

template <typename dist_t, bool DoNormalize = false>
struct Normalizer {
  static void normalize(std::vector<dist_t>& v) {}
};

template <typename dist_t>
struct Normalizer<dist_t, true> {
  static void normalize(std::vector<dist_t>& v) {
    const size_t dim = v.size();
    float norm = 0.0f;
    for (size_t i = 0; i < dim; i++) {
      norm += v[i] * v[i];
    }
    norm = 1.0 / std::sqrt(norm);

    for (size_t i = 0; i < dim; i++) {
      v[i] *= norm;
    }
  }
};

template<typename dist_t, typename Distance, bool DoNormalize = false>
class Hnsw {
public:

  // dim - length of the vectors being added
  // maxElements - size of the data being added
  // efConstruction - controls the quality of the graph. Higher values lead to improved recall at the
  //  expense of longer build time. Suggested range: 100-2000 (default: 200).
  // M - Controls maximum number of neighbors in the zero and above-zero layers. Higher values lead to
  //  better recall and shorter retrieval times, at the expense of longer indexing time.
  //  Suggested range: 5-100 (default: 16).
  Hnsw(const int dim, const size_t maxElements, const size_t M = 16,
       const size_t efConstruction = 200) :
    dim(dim), cur_l(0)
  {
    space = new Distance(dim);
    appr_alg = new hnswlib::HierarchicalNSW<dist_t>(space, maxElements, M,
                                                    efConstruction);
  }

  Hnsw(const int dim, const std::string path_to_index) :
  dim(dim), cur_l(0)
  {
    space = new Distance(dim);
    appr_alg = new hnswlib::HierarchicalNSW<dist_t>(space, path_to_index);
    cur_l = appr_alg->cur_element_count;
  }

  void addItem(const std::vector<dist_t> dv)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    Normalizer<dist_t, DoNormalize>::normalize(fv);

    appr_alg->addPoint(&fv[0], (size_t) cur_l);
    ++cur_l;
  }

  std::vector<hnswlib::labeltype> getNNs(const std::vector<dist_t> dv, size_t k)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    Normalizer<dist_t, DoNormalize>::normalize(fv);

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

    Normalizer<dist_t, DoNormalize>::normalize(fv);

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
    delete space;
    delete appr_alg;
  }

  int dim;
  bool normalize;
  hnswlib::labeltype cur_l;
  hnswlib::HierarchicalNSW<dist_t> *appr_alg;
  hnswlib::SpaceInterface<float> *space;
};

typedef Hnsw<float, hnswlib::L2Space, false> HnswL2;
typedef Hnsw<float, hnswlib::InnerProductSpace, true> HnswCosine;
typedef Hnsw<float, hnswlib::InnerProductSpace, false> HnswIp;

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

RCPP_EXPOSED_CLASS_NODECL(HnswCosine)
RCPP_MODULE(HnswCosine) {
  Rcpp::class_<HnswCosine>("HnswCosine")
  .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, ef, M")
  .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
  .method("addItem",    &HnswCosine::addItem,    "add item")
  .method("save",       &HnswCosine::callSave,   "save index to file")
  .method("getNNs",     &HnswCosine::getNNs,     "retrieve Nearest Neigbours given vector")
  .method("getNNsList", &HnswCosine::getNNsList, "retrieve Nearest Neigbours given vector")
  ;
}

RCPP_EXPOSED_CLASS_NODECL(HnswIp)
  RCPP_MODULE(HnswIp) {
    Rcpp::class_<HnswIp>("HnswIp")
    .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, ef, M")
    .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
    .method("addItem",    &HnswIp::addItem,    "add item")
    .method("save",       &HnswIp::callSave,   "save index to file")
    .method("getNNs",     &HnswIp::getNNs,     "retrieve Nearest Neigbours given vector")
    .method("getNNsList", &HnswIp::getNNsList, "retrieve Nearest Neigbours given vector")
    ;
  }
