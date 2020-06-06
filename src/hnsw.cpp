//  RcppHNSW -- Rcpp bindings to hnswlib library for Approximate Nearest Neighbours
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
#include <limits>
#include <memory>
#include <Rcpp.h>

#include "hnswlib.h"

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
    norm = 1.0f / (std::sqrt(norm) + 1e-30f);

    for (size_t i = 0; i < dim; i++) {
      v[i] *= norm;
    }
  }
};

template<typename dist_t, typename Distance, bool DoNormalize = false>
class Hnsw {
public:

  // dim - length of the vectors being added
  // max_elements - size of the data being added
  // M - Controls maximum number of neighbors in the zero and above-zero
  //  layers. Higher values lead t better recall and shorter retrieval times,
  //  at the expense of longer indexing time. Suggested range: 5-100
  //  (default: 16).
  // ef_construction - controls the quality of the graph. Higher values lead to
  //  improved recall at the expense of longer build time. Suggested range:
  //  100-2000 (default: 200).
  Hnsw(const int dim, const size_t max_elements, const size_t M = 16,
       const size_t ef_construction = 200) :
    dim(dim), cur_l(0),
    space(std::unique_ptr<Distance>(new Distance(dim))),
    appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
      new hnswlib::HierarchicalNSW<dist_t>(space.get(), max_elements, M,
                                           ef_construction)))
  { }

  Hnsw(const int dim, const std::string path_to_index) :
  dim(dim), cur_l(0),
  space(std::unique_ptr<Distance>(new Distance(dim))),
  appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
    new hnswlib::HierarchicalNSW<dist_t>(space.get(), path_to_index)))
  {
    cur_l = appr_alg->cur_element_count;
  }

  Hnsw(const int dim, const std::string path_to_index,
       const size_t max_elements) :
  dim(dim), cur_l(0),
  space(std::unique_ptr<Distance>(new Distance(dim))),
  appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
    new hnswlib::HierarchicalNSW<dist_t>(space.get(), path_to_index, false,
                                         max_elements)))
  {
    cur_l = appr_alg->cur_element_count;
  }

  void setEf(const size_t ef) {
    appr_alg->ef_ = ef;
  }

  void addItem(Rcpp::NumericVector dv)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    addItemNoCopy(fv);
  }

  void addItemNoCopy(std::vector<dist_t>& dv)
  {
    Normalizer<dist_t, DoNormalize>::normalize(dv);

    appr_alg->addPoint(dv.data(), static_cast<size_t>(cur_l));
    ++cur_l;
  }

  void addItems(Rcpp::NumericMatrix items) {
    for (int i = 0; i < items.nrow(); i++) {
      Rcpp::NumericMatrix::Row row = items.row(i);
      std::vector<dist_t> dv(row.size());
      std::copy(row.begin(), row.end(), dv.begin());
      addItemNoCopy(dv);
    }
  }

  std::vector<hnswlib::labeltype> getNNs(const std::vector<dist_t>& dv, size_t k)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    return getNNsNoCopy(fv, k);
  }

  std::vector<hnswlib::labeltype> getNNsNoCopy(std::vector<dist_t>& fv, size_t k)
  {
    Normalizer<dist_t, DoNormalize>::normalize(fv);

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
      appr_alg->searchKnn(fv.data(), k);

    if (result.size() != k) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    std::vector<hnswlib::labeltype> items;
    items.reserve(k);
    for (size_t i = 0; i < k; i++) {
      auto &result_tuple = result.top();
      items.push_back(result_tuple.second + 1);
      result.pop();
    }
    std::reverse(items.begin(), items.end());

    return items;
  }

  Rcpp::List getNNsList(const std::vector<dist_t>& dv, size_t k,
                        bool include_distances)
  {
    std::vector<dist_t> fv(dv.size());
    std::copy(dv.begin(), dv.end(), fv.begin());

    return getNNsListNoCopy(fv, k, include_distances);
  }

  Rcpp::List getNNsListNoCopy(std::vector<dist_t>& fv, size_t k,
                        bool include_distances)
  {
    Normalizer<dist_t, DoNormalize>::normalize(fv);

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
      appr_alg->searchKnn(fv.data(), k);

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
        items.push_back(result_tuple.second + 1);
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
        items.push_back(result_tuple.second + 1);
        result.pop();
      }

      std::reverse(items.begin(), items.end());

      return Rcpp::List::create(
        Rcpp::Named("item") = items);
    }
  }

  Rcpp::List getAllNNsList(Rcpp::NumericMatrix fm, size_t k,
                           bool include_distances = true)
  {
    Rcpp::IntegerMatrix allItems(fm.nrow(), k);
    Rcpp::NumericMatrix allDistances(fm.nrow(), k);

    for (int i = 0; i < fm.nrow(); i++) {
      Rcpp::NumericMatrix::Row row = fm.row(i);
      std::vector<dist_t> dv(row.size());
      std::copy(row.begin(), row.end(), dv.begin());
      Rcpp::List result = getNNsListNoCopy(dv, k, include_distances);
      std::vector<hnswlib::labeltype> items = result["item"];

      if (include_distances) {
        std::vector<dist_t> dist = result["distance"];

        for (size_t k = 0; k < items.size(); k++) {
          allItems(i, k) = items[k];
          allDistances(i, k) = dist[k];
        }
      }
      else {
        for (size_t k = 0; k < items.size(); k++) {
          allItems(i, k) = items[k];
        }
      }
    }

    if (include_distances) {
      return Rcpp::List::create(
        Rcpp::Named("item") = allItems,
        Rcpp::Named("distance") = allDistances);
    }
    else {
      return Rcpp::List::create(
        Rcpp::Named("item") = allItems);
    }
  }

  Rcpp::IntegerMatrix getAllNNs(Rcpp::NumericMatrix fm, size_t k)
  {
    Rcpp::IntegerMatrix output(fm.nrow(), k);

    for (int i = 0; i < fm.nrow(); i++) {
      Rcpp::NumericMatrix::Row row = fm.row(i);
      std::vector<dist_t> dv(row.size());
      std::copy(row.begin(), row.end(), dv.begin());
      std::vector<hnswlib::labeltype> result = getNNsNoCopy(dv, k);

      for (size_t k = 0; k < result.size(); k++) {
        output(i, k) = result[k];
      }
    }

    return output;
  }

  void callSave(const std::string path_to_index) {
    appr_alg->saveIndex(path_to_index);
  }

  std::size_t size() const {
    return appr_alg->cur_element_count;
  }

  void markDeleted(std::size_t label) {
    if (label < 1 || label > size()) {
      Rcpp::stop("Bad label");
    }
    // internally labels are zero-indexed
    appr_alg->markDelete(label - 1);
  }

  void resizeIndex(std::size_t new_size) {
    appr_alg->resizeIndex(new_size);
  }

  int dim;
  bool normalize;
  hnswlib::labeltype cur_l;
  std::unique_ptr<Distance> space;
  std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>> appr_alg;
};

typedef Hnsw<float, hnswlib::L2Space, false> HnswL2;
typedef Hnsw<float, hnswlib::InnerProductSpace, true> HnswCosine;
typedef Hnsw<float, hnswlib::InnerProductSpace, false> HnswIp;

RCPP_EXPOSED_CLASS_NODECL(HnswL2)
RCPP_MODULE(HnswL2) {
  Rcpp::class_<HnswL2>("HnswL2")
  .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, M, ef")
  .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
  .constructor<int32_t, std::string, size_t>("constructor with dimension, loading from filename, number of items")
  .method("setEf",      &HnswL2::setEf,      "set ef value")
  .method("addItem",    &HnswL2::addItem,    "add item")
  .method("addItems",   &HnswL2::addItems,   "add items")
  .method("save",       &HnswL2::callSave,   "save index to file")
  .method("getNNs",     &HnswL2::getNNs,     "retrieve Nearest Neigbours given vector")
  .method("getNNsList", &HnswL2::getNNsList, "retrieve Nearest Neigbours given vector")
  .method("getAllNNs",  &HnswL2::getAllNNs,  "retrieve Nearest Neigbours given matrix")
  .method("getAllNNsList",  &HnswL2::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
  .method("size",       &HnswL2::size,        "number of items added to the index")
  .method("markDeleted",  &HnswL2::markDeleted, "remove the item with the specified label from the index")
  .method("resizeIndex",  &HnswL2::resizeIndex, "resize the index to use this number of items")
  ;
}

RCPP_EXPOSED_CLASS_NODECL(HnswCosine)
RCPP_MODULE(HnswCosine) {
  Rcpp::class_<HnswCosine>("HnswCosine")
  .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, M, ef")
  .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
  .constructor<int32_t, std::string, size_t>("constructor with dimension, loading from filename, number of items")
  .method("setEf",      &HnswCosine::setEf,      "set ef value")
  .method("addItem",    &HnswCosine::addItem,    "add item")
  .method("addItems",   &HnswCosine::addItems,   "add items")
  .method("save",       &HnswCosine::callSave,   "save index to file")
  .method("getNNs",     &HnswCosine::getNNs,     "retrieve Nearest Neigbours given vector")
  .method("getNNsList", &HnswCosine::getNNsList, "retrieve Nearest Neigbours given vector")
  .method("getAllNNs",  &HnswCosine::getAllNNs,  "retrieve Nearest Neigbours given matrix")
  .method("getAllNNsList",  &HnswCosine::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
  .method("size",       &HnswCosine::size,       "number of items added to the index")
  .method("markDeleted",  &HnswCosine::markDeleted, "remove the item with the specified label from the index")
  .method("resizeIndex",  &HnswCosine::resizeIndex, "resize the index to use this number of items")
  ;
}

RCPP_EXPOSED_CLASS_NODECL(HnswIp)
  RCPP_MODULE(HnswIp) {
    Rcpp::class_<HnswIp>("HnswIp")
    .constructor<int32_t, size_t, size_t, size_t>("constructor with dimension, number of items, M, ef")
    .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
    .constructor<int32_t, std::string, size_t>("constructor with dimension, loading from filename, number of items")
    .method("setEf",      &HnswIp::setEf,      "set ef value")
    .method("addItem",    &HnswIp::addItem,    "add item")
    .method("addItems",   &HnswIp::addItems,   "add items")
    .method("save",       &HnswIp::callSave,   "save index to file")
    .method("getNNs",     &HnswIp::getNNs,     "retrieve Nearest Neigbours given vector")
    .method("getNNsList", &HnswIp::getNNsList, "retrieve Nearest Neigbours given vector")
    .method("getAllNNs",  &HnswIp::getAllNNs,  "retrieve Nearest Neigbours given matrix")
    .method("getAllNNsList",  &HnswIp::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
    .method("size",       &HnswIp::size,       "number of items added to the index")
    .method("markDeleted",  &HnswIp::markDeleted, "remove the item with the specified label from the index")
    .method("resizeIndex",  &HnswIp::resizeIndex, "resize the index to use this number of items")
    ;
  }
