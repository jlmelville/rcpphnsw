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
#include <thread>

#include <Rcpp.h>

#include "hnswlib.h"

#include "RcppPerpendicular/RcppPerpendicular.h"


template <typename dist_t, bool DoNormalize = false>
struct Normalizer {
  static void normalize(std::vector<dist_t>& vec) {}
};

template <typename dist_t>
struct Normalizer<dist_t, true> {
  static const constexpr float FLOAT_MIN = 1e-30F;

  static void normalize(std::vector<dist_t>& vec) {
    const std::size_t dim = vec.size();
    float norm = 0.0F;
    for (std::size_t i = 0; i < dim; i++) {
      norm += vec[i] * vec[i];
    }
    norm = 1.0F / (std::sqrt(norm) + FLOAT_MIN);

    for (std::size_t i = 0; i < dim; i++) {
      vec[i] *= norm;
    }
  }
};

template<typename dist_t, typename Distance, bool DoNormalize = false>
class Hnsw {
  static const constexpr std::size_t M_DEFAULT = 16;
  static const constexpr std::size_t EF_CONSTRUCTION_DEFAULT = 200;
public:
  // dim - length of the vectors being added
  // max_elements - size of the data being added
  // M - Controls maximum number of neighbors in the zero and above-zero
  //  layers. Higher values lead to better recall and shorter retrieval times,
  //  at the expense of longer indexing time. Suggested range: 5-100
  //  (default: 16).
  // ef_construction - controls the quality of the graph. Higher values lead to
  //  improved recall at the expense of longer build time. Suggested range:
  //  100-2000 (default: 200).
  Hnsw(int dim, std::size_t max_elements, std::size_t M = M_DEFAULT,
       std::size_t ef_construction = EF_CONSTRUCTION_DEFAULT) :
  dim(dim), normalize(false), cur_l(0), numThreads(0), grainSize(1),
  space(std::unique_ptr<Distance>(new Distance(dim))),
  appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
      new hnswlib::HierarchicalNSW<dist_t>(space.get(), max_elements, M,
                                           ef_construction)))
  { }

  Hnsw(int dim, const std::string &path_to_index) :
  dim(dim), normalize(false), cur_l(0), numThreads(0), grainSize(1),
  space(std::unique_ptr<Distance>(new Distance(dim))),
  appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
      new hnswlib::HierarchicalNSW<dist_t>(space.get(), path_to_index)))
  {
    cur_l = appr_alg->cur_element_count;
  }

  Hnsw(int dim, const std::string &path_to_index, std::size_t max_elements) :
  dim(dim), normalize(false), cur_l(0), numThreads(0), grainSize(1),
  space(std::unique_ptr<Distance>(new Distance(dim))),
  appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
      new hnswlib::HierarchicalNSW<dist_t>(space.get(), path_to_index, false,
                                           max_elements)))
  {
    cur_l = appr_alg->cur_element_count;
  }

  void setEf(std::size_t ef) {
    appr_alg->ef_ = ef;
  }

  void addItem(Rcpp::NumericVector item)
  {
    std::vector<dist_t> item_copy(item.size());
    std::copy(item.begin(), item.end(), item_copy.begin());

    addItemImpl(item_copy, cur_l);
  }

  void addItemImpl(std::vector<dist_t>& item, std::size_t label)
  {
    Normalizer<dist_t, DoNormalize>::normalize(item);

    appr_alg->addPoint(item.data(), label);
    ++cur_l;
  }

  class AddWorker {
    Hnsw<dist_t, Distance, DoNormalize> &hnsw;
    const std::vector<double> &data;
    std::size_t nrow;
    std::size_t ncol;
    std::size_t index_start;

    public:
        AddWorker(Hnsw<dist_t, Distance, DoNormalize> &hnsw,
              const std::vector<double> &data,
              std::size_t nrow,
              std::size_t ncol,
              std::size_t index_start) :
      hnsw(hnsw), data(data), nrow(nrow), ncol(ncol), index_start(index_start)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      std::vector<dist_t> item_copy(ncol);

      for (std::size_t i = begin; i < end; i++) {
        for (std::size_t j = 0; j < ncol; j++) {
          item_copy[j] =  data[nrow * j + i];
        }
        hnsw.addItemImpl(item_copy, index_start + i);
      }
    }
  };

  void addItems(const Rcpp::NumericMatrix &items) {
    const std::size_t nrow = items.nrow();
    const std::size_t ncol = items.ncol();
    auto data = Rcpp::as<std::vector<double>>(items);

    AddWorker worker(*this, data, nrow, ncol, cur_l);
    RcppPerpendicular::parallel_for(nrow, worker, numThreads);
    cur_l = size();
  }


  auto getNNs(const std::vector<dist_t>& item, std::size_t nnbrs) -> std::vector<hnswlib::labeltype> 
  {
    std::vector<dist_t> item_copy(item);

    bool found_all = true;
    std::vector<hnswlib::labeltype> nbr_labels = getNNsImpl(item_copy, nnbrs, found_all);
    if (!found_all) {
      Rcpp::stop("Unable to find nnbrs results. Probably ef or M is too small");
    }

    return nbr_labels;
  }

  auto getNNsList(const std::vector<dist_t>& item, std::size_t nnbrs,
                        bool include_distances) -> Rcpp::List
  {
    std::vector<dist_t> item_copy(item);

    bool found_all = true;
    std::vector<dist_t> distances(0);
    std::vector<hnswlib::labeltype> nbr_labels =
      getNNsImpl(item_copy, nnbrs, include_distances, distances, found_all);
    if (!found_all) {
      Rcpp::stop("Unable to find nnbrs results. Probably ef or M is too small");
    }

    auto nbr_list = Rcpp::List::create(Rcpp::Named("item") = nbr_labels);
    if (include_distances) {
      nbr_list["distance"] = distances;
    }
    return nbr_list;
  }

  auto getNNsImpl(
      std::vector<dist_t>& item, std::size_t nnbrs, bool include_distances,
      std::vector<dist_t>& distances, bool& found_all) -> std::vector<hnswlib::labeltype>
  {
    found_all = true;
    Normalizer<dist_t, DoNormalize>::normalize(item);

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
      appr_alg->searchKnn(item.data(), nnbrs);

    const std::size_t nresults = result.size();
    if (nresults != nnbrs) {
      found_all = false;
    }

    std::vector<hnswlib::labeltype> items;
    items.reserve(nnbrs);

    if (include_distances) {
      distances.reserve(nnbrs);
      distances.clear();

      for (std::size_t i = 0; i < nresults; i++) {
        auto &result_tuple = result.top();
        distances.push_back(result_tuple.first);
        items.push_back(result_tuple.second + 1);
        result.pop();
      }
      if (!found_all) {
        for (std::size_t i = 0; i != nnbrs - nresults; i++) {
          distances.push_back((std::numeric_limits<dist_t>::max)());
          items.push_back(-1);
        }
      }

      std::reverse(distances.begin(), distances.end());
      std::reverse(items.begin(), items.end());
    }
    else {
      for (std::size_t i = 0; i < nresults; i++) {
        auto &result_tuple = result.top();
        items.push_back(result_tuple.second + 1);
        result.pop();
      }
      if (!found_all) {
        for (std::size_t i = 0; i != nnbrs - nresults; i++) {
          items.push_back(-1);
        }
      }

      std::reverse(items.begin(), items.end());
    }

    return items;
  }

  auto getNNsImpl(std::vector<dist_t>& item,
                                             std::size_t nnbrs, bool& found_all) -> std::vector<hnswlib::labeltype>
  {
    bool include_distances = false;
    std::vector<dist_t> distances(0);
    return getNNsImpl(item, nnbrs, include_distances, distances, found_all);
  }

  class SearchListWorker {
    Hnsw<dist_t, Distance, DoNormalize> &hnsw;
    const std::vector<double> &data;

    const std::size_t nitems;
    const std::size_t ndim;
    const std::size_t nnbrs;
    bool include_distances;

public:
    std::vector<hnswlib::labeltype> idx_vec;
    std::vector<dist_t> dist_vec;
    bool ok{true};

    SearchListWorker(Hnsw<dist_t, Distance, DoNormalize> &hnsw,
                     const std::vector<double> &data, std::size_t nitems,
                     std::size_t ndim, std::size_t nnbrs, bool include_distances) :
      hnsw(hnsw), data(data), nitems(nitems), ndim(ndim), nnbrs(nnbrs),
      include_distances(include_distances), idx_vec(nitems * nnbrs),
      dist_vec(nitems * nnbrs)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      std::vector<dist_t> item_copy(ndim);
      std::vector<dist_t> distances(0);

      for (std::size_t i = begin; i < end; i++) {
        bool ok_row = true;
        for (std::size_t j = 0; j < ndim; j++) {
          item_copy[j] = data[j * nitems + i];
        }

        std::vector<hnswlib::labeltype> items =
          hnsw.getNNsImpl(item_copy, nnbrs, include_distances, distances, ok_row);
        if (!ok_row) {
          ok = false;
          break;
        }

        if (include_distances) {
          for (std::size_t k = 0; k < nnbrs; k++) {
            idx_vec[k * nitems + i] = items[k];
            dist_vec[k * nitems + i] = distances[k];
          }
        }
        else {
          for (std::size_t k = 0; k < nnbrs; k++) {
            idx_vec[k * nitems + i] = items[k];
          }
        }
      }
    }
  };

  auto getAllNNsList(const Rcpp::NumericMatrix &items, std::size_t nnbrs,
                           bool include_distances = true) -> Rcpp::List
  {
    const std::size_t nrow = items.nrow();
    const std::size_t ncol = items.ncol();
    auto data = Rcpp::as<std::vector<double>>(items);

    SearchListWorker worker(*this, data, nrow, ncol, nnbrs, include_distances);
    RcppPerpendicular::parallel_for(nrow, worker, numThreads);
    if (!worker.ok) {
      Rcpp::stop("Unable to find nnbrs results. Probably ef or M is too small");
    }

    auto result = Rcpp::List::create(
      Rcpp::Named("item") = Rcpp::IntegerMatrix(nrow, nnbrs, worker.idx_vec.begin())
    );
    if (include_distances) {
      result["distance"] = Rcpp::NumericMatrix(nrow, nnbrs, worker.dist_vec.begin());
    }
    return result;
  }


  class SearchWorker {
    Hnsw<dist_t, Distance, DoNormalize> &hnsw;
    const std::vector<double> &data;

    const std::size_t nitems;
    const std::size_t ndim;
    const std::size_t nnbrs;

    bool include_distances{false};
    std::vector<dist_t> distances;

public:
    std::vector<hnswlib::labeltype> idx_vec;
    bool ok{true};

    SearchWorker(Hnsw<dist_t, Distance, DoNormalize> &hnsw,
                 const std::vector<double> &data, std::size_t nitems,
                 std::size_t ndim, std::size_t nnbrs) :
      hnsw(hnsw), data(data), nitems(nitems), ndim(ndim), nnbrs(nnbrs), 
      distances(0), idx_vec(nitems * nnbrs) {}

    void operator()(std::size_t begin, std::size_t end) {
      std::vector<dist_t> item_copy(ndim);
      for (std::size_t i = begin; i < end; i++) {
        for (std::size_t j = 0; j < ndim; j++) {
          item_copy[j] = data[j * nitems + i];
        }

        bool ok_row = true;
        std::vector<hnswlib::labeltype> items =
          hnsw.getNNsImpl(item_copy, nnbrs, include_distances, distances, ok_row);
        if (!ok_row) {
          ok = false;
          break;
        }

        for (std::size_t k = 0; k < nnbrs; k++) {
          idx_vec[k * nitems + i] = items[k];
        }
      }
    }
  };

  auto getAllNNs(const Rcpp::NumericMatrix &items, std::size_t nnbrs) -> Rcpp::IntegerMatrix
  {
    const std::size_t nrow = items.nrow();
    const std::size_t ncol = items.ncol();

    auto data = Rcpp::as<std::vector<double>>(items);

    SearchWorker worker(*this, data, nrow, ncol, nnbrs);

    RcppPerpendicular::parallel_for(nrow, worker, numThreads);
    if (!worker.ok) {
      Rcpp::stop("Unable to find nnbrs results. Probably ef or M is too small");
    }

    Rcpp::IntegerMatrix idx(nrow, nnbrs, worker.idx_vec.begin());
    return idx;
  }

  void callSave(const std::string &path_to_index) {
    appr_alg->saveIndex(path_to_index);
  }

  auto size() const -> std::size_t {
    return appr_alg->cur_element_count;
  }

  void setNumThreads(std::size_t numThreads) {
    this->numThreads = numThreads;
  }

  void setGrainSize(std::size_t grainSize) {
    this->grainSize = grainSize;
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

private:
  int dim;
  bool normalize;
  hnswlib::labeltype cur_l;
  std::size_t numThreads;
  std::size_t grainSize;
  std::unique_ptr<Distance> space;
  std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>> appr_alg;
};

using HnswL2 = Hnsw<float, hnswlib::L2Space, false>;
using HnswCosine = Hnsw<float, hnswlib::InnerProductSpace, true>;
using HnswIp = Hnsw<float, hnswlib::InnerProductSpace, false>;

RCPP_EXPOSED_CLASS_NODECL(HnswL2)
  RCPP_MODULE(HnswL2) {
    Rcpp::class_<HnswL2>("HnswL2")
    .constructor<int32_t, std::size_t, std::size_t, std::size_t>("constructor with dimension, number of items, M, ef")
    .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
    .constructor<int32_t, std::string, std::size_t>("constructor with dimension, loading from filename, number of items")
    .method("setEf",      &HnswL2::setEf,      "set ef value")
    .method("addItem",    &HnswL2::addItem,    "add item")
    .method("addItems",   &HnswL2::addItems,   "add items")
    .method("save",       &HnswL2::callSave,   "save index to file")
    .method("getNNs",     &HnswL2::getNNs,     "retrieve Nearest Neigbours given vector")
    .method("getNNsList", &HnswL2::getNNsList, "retrieve Nearest Neigbours given vector")
    .method("getAllNNs",  &HnswL2::getAllNNs,  "retrieve Nearest Neigbours given matrix")
    .method("getAllNNsList",  &HnswL2::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
    .method("size",       &HnswL2::size,        "number of items added to the index")
    .method("setNumThreads",  &HnswL2::setNumThreads, "set the number of threads to use")
    .method("setGrainSize",  &HnswL2::setGrainSize, "set minimum grain size for using multiple threads")
    .method("markDeleted",  &HnswL2::markDeleted, "remove the item with the specified label from the index")
    .method("resizeIndex",  &HnswL2::resizeIndex, "resize the index to use this number of items")
    ;
  }

RCPP_EXPOSED_CLASS_NODECL(HnswCosine)
  RCPP_MODULE(HnswCosine) {
    Rcpp::class_<HnswCosine>("HnswCosine")
    .constructor<int32_t, std::size_t, std::size_t, std::size_t>("constructor with dimension, number of items, M, ef")
    .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
    .constructor<int32_t, std::string, std::size_t>("constructor with dimension, loading from filename, number of items")
    .method("setEf",      &HnswCosine::setEf,      "set ef value")
    .method("addItem",    &HnswCosine::addItem,    "add item")
    .method("addItems",   &HnswCosine::addItems,   "add items")
    .method("save",       &HnswCosine::callSave,   "save index to file")
    .method("getNNs",     &HnswCosine::getNNs,     "retrieve Nearest Neigbours given vector")
    .method("getNNsList", &HnswCosine::getNNsList, "retrieve Nearest Neigbours given vector")
    .method("getAllNNs",  &HnswCosine::getAllNNs,  "retrieve Nearest Neigbours given matrix")
    .method("getAllNNsList",  &HnswCosine::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
    .method("size",       &HnswCosine::size,       "number of items added to the index")
    .method("setNumThreads",  &HnswCosine::setNumThreads, "set the number of threads to use")
    .method("setGrainSize",  &HnswCosine::setGrainSize, "set minimum grain size for using multiple threads")
    .method("markDeleted",  &HnswCosine::markDeleted, "remove the item with the specified label from the index")
    .method("resizeIndex",  &HnswCosine::resizeIndex, "resize the index to use this number of items")
    ;
  }

RCPP_EXPOSED_CLASS_NODECL(HnswIp)
  RCPP_MODULE(HnswIp) {
    Rcpp::class_<HnswIp>("HnswIp")
    .constructor<int32_t, std::size_t, std::size_t, std::size_t>("constructor with dimension, number of items, M, ef")
    .constructor<int32_t, std::string>("constructor with dimension, loading from filename")
    .constructor<int32_t, std::string, std::size_t>("constructor with dimension, loading from filename, number of items")
    .method("setEf",      &HnswIp::setEf,      "set ef value")
    .method("addItem",    &HnswIp::addItem,    "add item")
    .method("save",       &HnswIp::callSave,   "save index to file")
    .method("getNNs",     &HnswIp::getNNs,     "retrieve Nearest Neigbours given vector")
    .method("getNNsList", &HnswIp::getNNsList, "retrieve Nearest Neigbours given vector")
    .method("getAllNNs",  &HnswIp::getAllNNs,  "retrieve Nearest Neigbours given matrix")
    .method("getAllNNsList",  &HnswIp::getAllNNsList,  "retrieve Nearest Neigbours given matrix")
    .method("size",       &HnswIp::size,       "number of items added to the index")
    .method("setNumThreads",  &HnswIp::setNumThreads, "set the number of threads to use")
    .method("setGrainSize",  &HnswIp::setGrainSize, "set minimum grain size for using multiple threads")
    .method("markDeleted",  &HnswIp::markDeleted, "remove the item with the specified label from the index")
    .method("resizeIndex",  &HnswIp::resizeIndex, "resize the index to use this number of items")
    ;
  }
