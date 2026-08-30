//  RcppHNSW -- Rcpp bindings to hnswlib library for Approximate Nearest
//  Neighbors
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
#include <atomic>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <Rcpp.h>

#include "rcpphnsw.h"

#include "pforr/pforr.h"

namespace {

constexpr int R_INTEGER_MAX = (std::numeric_limits<int>::max)();

auto is_numeric_value(SEXP value) -> bool {
  const int type = TYPEOF(value);
  return (type == INTSXP || type == REALSXP) && !Rf_inherits(value, "factor");
}

auto numeric_value_at(SEXP value, R_xlen_t index) -> double {
  if (TYPEOF(value) == INTSXP) {
    const int result = INTEGER_ELT(value, index);
    return result == NA_INTEGER ? NA_REAL : static_cast<double>(result);
  }
  return REAL_ELT(value, index);
}

auto checked_whole_number_at(SEXP value, R_xlen_t index, const char *name,
                             int lower, int upper = R_INTEGER_MAX)
    -> std::size_t {
  const double result = numeric_value_at(value, index);
  if (!std::isfinite(result) || result != std::floor(result) ||
      result < static_cast<double>(lower) ||
      result > static_cast<double>(upper)) {
    Rcpp::stop("%s cannot be outside the whole-number range %d to %d", name,
               lower, upper);
  }
  return static_cast<std::size_t>(result);
}

auto check_whole_number(SEXP value, const char *name, int lower,
                        int upper = R_INTEGER_MAX) -> std::size_t {
  if (!is_numeric_value(value) || Rf_xlength(value) != 1) {
    Rcpp::stop("%s cannot be outside the whole-number range %d to %d", name,
               lower, upper);
  }
  return checked_whole_number_at(value, 0, name, lower, upper);
}

auto check_item_identifier_at(SEXP value, R_xlen_t index, int upper)
    -> std::size_t {
  const double result = numeric_value_at(value, index);
  if (!std::isfinite(result) || result != std::floor(result) || result < 1.0 ||
      result > static_cast<double>(upper)) {
    Rcpp::stop("Invalid index requested: item identifiers must be whole "
               "numbers from 1 to %d",
               upper);
  }
  return static_cast<std::size_t>(result);
}

auto check_logical(SEXP value, const char *name) -> bool {
  if (TYPEOF(value) != LGLSXP || Rf_xlength(value) != 1 ||
      LOGICAL_ELT(value, 0) == NA_LOGICAL) {
    Rcpp::stop("%s must be TRUE or FALSE", name);
  }
  return LOGICAL_ELT(value, 0) == TRUE;
}

auto check_path(SEXP value, const char *name) -> std::string {
  if (TYPEOF(value) != STRSXP || Rf_xlength(value) != 1 ||
      STRING_ELT(value, 0) == NA_STRING ||
      Rf_length(STRING_ELT(value, 0)) == 0) {
    Rcpp::stop("%s must be one non-empty, non-missing string", name);
  }
  return Rcpp::as<std::string>(value);
}

auto checked_product(std::size_t left, std::size_t right,
                     const char *description) -> std::size_t {
  if (left != 0 && right > (std::numeric_limits<std::size_t>::max)() / left) {
    Rcpp::stop("%s is too large", description);
  }
  const std::size_t result = left * right;
  if (result > static_cast<std::size_t>(R_XLEN_T_MAX)) {
    Rcpp::stop("%s is too large for an R vector", description);
  }
  return result;
}

auto checked_float_at(SEXP values, R_xlen_t index, const char *name) -> float {
  const double value = numeric_value_at(values, index);
  const double float_max =
      static_cast<double>((std::numeric_limits<float>::max)());
  if (!std::isfinite(value) || value < -float_max || value > float_max) {
    Rcpp::stop(
        "%s must contain only finite values representable as single-precision "
        "floats",
        name);
  }
  const float converted = static_cast<float>(value);
  if (!std::isfinite(converted)) {
    Rcpp::stop(
        "%s must contain only finite values representable as single-precision "
        "floats",
        name);
  }
  return converted;
}

} // namespace

template <typename dist_t, bool DoNormalize = false> struct Normalizer {
  static void normalize(dist_t *, std::size_t) {}
};

template <typename dist_t> struct Normalizer<dist_t, true> {
  static void normalize(dist_t *values, std::size_t size) {
    double squared_norm = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
      const double value = static_cast<double>(values[i]);
      squared_norm += value * value;
    }
    if (!(squared_norm > 0.0) || !std::isfinite(squared_norm)) {
      Rcpp::stop("Cosine vectors must have a positive finite norm after "
                 "conversion to single precision");
    }

    const double inverse_norm = 1.0 / std::sqrt(squared_norm);
    for (std::size_t i = 0; i < size; ++i) {
      values[i] = static_cast<dist_t>(values[i] * inverse_norm);
    }
  }
};

struct NoDistanceProcess {
  template <typename dist_t>
  static void process_distances(std::vector<dist_t> &) {}
};

struct SquareRootDistanceProcess {
  template <typename dist_t>
  static void process_distances(std::vector<dist_t> &vec) {
    for (std::size_t i = 0; i < vec.size(); i++) {
      vec[i] = std::sqrt(vec[i]);
    }
  }
};

template <typename dist_t, typename Distance, bool DoNormalize,
          typename DistanceProcess>
class Hnsw {
  struct NewIndexConfig {
    int dim;
    std::size_t max_elements;
    std::size_t M;
    std::size_t ef_construction;
    std::size_t random_seed;
  };

  struct LoadIndexConfig {
    int dim;
    std::string path;
    std::size_t max_elements;
  };

  struct CheckedItems {
    int nitems;
    std::vector<dist_t> data;
  };

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
  Hnsw(SEXP dim, SEXP max_elements, SEXP M, SEXP ef_construction)
      : Hnsw(validateNewIndex(dim, max_elements, M, ef_construction, R_NilValue,
                              false)) {}

  Hnsw(SEXP dim, SEXP max_elements, SEXP M, SEXP ef_construction,
       SEXP random_seed)
      : Hnsw(validateNewIndex(dim, max_elements, M, ef_construction,
                              random_seed, true)) {}

  Hnsw(SEXP dim, SEXP path_to_index)
      : Hnsw(validateLoadIndex(dim, path_to_index, R_NilValue, false)) {}

  Hnsw(SEXP dim, SEXP path_to_index, SEXP max_elements)
      : Hnsw(validateLoadIndex(dim, path_to_index, max_elements, true)) {}

  void setEf(SEXP ef) {
    ensureUsable();
    appr_alg->ef_ = check_whole_number(ef, "ef", 1);
  }

  void addItem(SEXP item) {
    ensureUsable();
    std::vector<dist_t> item_copy = checkedItem(item, "Item to add");
    ensureCapacityFor(1);
    const std::size_t label = currentSize();
    try {
      addItemImpl(item_copy, label);
    } catch (...) {
      usable = false;
      throw;
    }
  }

  void addItemsCol(SEXP items) { addItemsAdapter(items, false); }

  void addItems(SEXP items) { addItemsAdapter(items, true); }

  auto getNNs(SEXP item, SEXP nnbrs) -> std::vector<hnswlib::labeltype> {
    ensureUsable();
    std::vector<dist_t> item_copy = checkedItem(item, "Query item");
    const std::size_t checked_nnbrs = checkK(nnbrs);

    bool found_all = true;
    std::vector<hnswlib::labeltype> nbr_labels =
        getNNsImpl(item_copy, checked_nnbrs, found_all);
    if (!found_all) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    return nbr_labels;
  }

  auto getNNsList(SEXP item, SEXP nnbrs, SEXP include_distances) -> Rcpp::List {
    ensureUsable();
    std::vector<dist_t> item_copy = checkedItem(item, "Query item");
    const std::size_t checked_nnbrs = checkK(nnbrs);
    const bool checked_include_distances =
        check_logical(include_distances, "include_distances");

    bool found_all = true;
    std::vector<dist_t> distances;
    std::vector<hnswlib::labeltype> nbr_labels =
        getNNsImpl(item_copy, checked_nnbrs, checked_include_distances,
                   distances, found_all);
    if (!found_all) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    auto nbr_list = Rcpp::List::create(Rcpp::Named("item") = nbr_labels);
    if (checked_include_distances) {
      DistanceProcess::process_distances(distances);
      nbr_list["distance"] = distances;
    }
    return nbr_list;
  }

  auto getNNsImpl(std::vector<dist_t> &item, std::size_t nnbrs,
                  bool include_distances, std::vector<dist_t> &distances,
                  bool &found_all) -> std::vector<hnswlib::labeltype> {
    found_all = true;

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> result =
        appr_alg->searchKnn(item.data(), nnbrs);

    const std::size_t nresults = result.size();
    if (nresults != nnbrs) {
      found_all = false;
    }

    std::vector<hnswlib::labeltype> result_items;
    result_items.reserve(nnbrs);

    if (include_distances) {
      distances.reserve(nnbrs);
      distances.clear();

      for (std::size_t i = 0; i < nresults; i++) {
        auto &result_tuple = result.top();
        distances.push_back(result_tuple.first);
        result_items.push_back(result_tuple.second + 1);
        result.pop();
      }
      std::reverse(distances.begin(), distances.end());
      std::reverse(result_items.begin(), result_items.end());
    } else {
      for (std::size_t i = 0; i < nresults; i++) {
        auto &result_tuple = result.top();
        result_items.push_back(result_tuple.second + 1);
        result.pop();
      }
      std::reverse(result_items.begin(), result_items.end());
    }

    return result_items;
  }

  auto getNNsImpl(std::vector<dist_t> &item, std::size_t nnbrs, bool &found_all)
      -> std::vector<hnswlib::labeltype> {
    std::vector<dist_t> distances;
    return getNNsImpl(item, nnbrs, false, distances, found_all);
  }

  auto getAllNNsList(SEXP items, SEXP nnbrs, SEXP include_distances)
      -> Rcpp::List {
    return getAllNNsListAdapter(items, nnbrs, include_distances, true);
  }

  auto getAllNNs(SEXP items, SEXP nnbrs) -> Rcpp::IntegerMatrix {
    return getAllNNsAdapter(items, nnbrs, true);
  }

  auto getAllNNsListCol(SEXP items, SEXP nnbrs, SEXP include_distances)
      -> Rcpp::List {
    return getAllNNsListAdapter(items, nnbrs, include_distances, false);
  }

  auto getAllNNsCol(SEXP items, SEXP nnbrs) -> Rcpp::IntegerMatrix {
    return getAllNNsAdapter(items, nnbrs, false);
  }

  auto getItems(SEXP ids) -> Rcpp::NumericMatrix {
    ensureUsable();
    if (!is_numeric_value(ids) ||
        Rf_getAttrib(ids, R_DimSymbol) != R_NilValue) {
      Rcpp::stop("ids must be a numeric vector of item identifiers");
    }
    const R_xlen_t nitems_xlen = Rf_xlength(ids);
    if (nitems_xlen > R_INTEGER_MAX) {
      Rcpp::stop("ids cannot contain more than INT_MAX values");
    }
    const auto nitems = static_cast<std::size_t>(nitems_xlen);
    checked_product(static_cast<std::size_t>(dim), nitems,
                    "Requested item result");

    const std::size_t total_count = currentSize();
    if (total_count == 0 && nitems != 0) {
      Rcpp::stop("Invalid index requested: index has no items");
    }

    std::vector<hnswlib::labeltype> checked_ids(nitems);
    for (R_xlen_t i = 0; i < nitems_xlen; ++i) {
      const std::size_t id =
          check_item_identifier_at(ids, i, static_cast<int>(total_count));
      checked_ids[static_cast<std::size_t>(i)] = id - 1;
    }

    std::vector<dist_t> data = getItemsImpl(checked_ids);
    Rcpp::NumericMatrix result(dim, static_cast<int>(nitems));
    std::copy(data.begin(), data.end(), result.begin());
    return Rcpp::transpose(result);
  }

  void callSave(SEXP path_to_index) {
    ensureUsable();
    appr_alg->saveIndex(check_path(path_to_index, "path_to_index"));
  }

  auto size() const -> std::size_t {
    ensureUsable();
    return currentSize();
  }

  void setNumThreads(SEXP value) {
    ensureUsable();
    numThreads = check_whole_number(value, "n_threads", 0);
  }

  void setGrainSize(SEXP value) {
    ensureUsable();
    const std::size_t checked = check_whole_number(value, "grain_size", 0);
    grainSize = (std::max)(checked, static_cast<std::size_t>(1));
  }

  void markDeleted(SEXP value) {
    ensureUsable();
    const std::size_t total_count = currentSize();
    if (total_count == 0) {
      Rcpp::stop("Bad label: index has no items");
    }
    const std::size_t label =
        check_whole_number(value, "label", 1, static_cast<int>(total_count));
    appr_alg->markDelete(label - 1);
  }

  void resizeIndex(SEXP value) {
    ensureUsable();
    const std::size_t total_count = currentSize();
    const std::size_t new_size =
        check_whole_number(value, "new_size", static_cast<int>(total_count));
    appr_alg->resizeIndex(new_size);
  }

private:
  explicit Hnsw(const NewIndexConfig &config)
      : dim(config.dim), usable(true), numThreads(0), grainSize(1),
        space(std::unique_ptr<Distance>(new Distance(config.dim))),
        appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
            new hnswlib::HierarchicalNSW<dist_t>(
                space.get(), config.max_elements, config.M,
                config.ef_construction, config.random_seed))) {}

  explicit Hnsw(const LoadIndexConfig &config)
      : dim(config.dim), usable(true), numThreads(0), grainSize(1),
        space(std::unique_ptr<Distance>(new Distance(config.dim))),
        appr_alg(std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>>(
            new hnswlib::HierarchicalNSW<dist_t>(space.get(), config.path,
                                                 false, config.max_elements))) {
    if (currentSize() > static_cast<std::size_t>(R_INTEGER_MAX) ||
        appr_alg->getMaxElements() > static_cast<std::size_t>(R_INTEGER_MAX)) {
      Rcpp::stop("Loaded index capacity cannot be larger than INT_MAX");
    }
  }

  void ensureUsable() const {
    if (!usable) {
      Rcpp::stop("Index is unusable after a failed insertion; discard it and "
                 "rebuild or reload the index");
    }
  }

  auto currentSize() const -> std::size_t {
    return appr_alg->cur_element_count;
  }

  void addItemImpl(std::vector<dist_t> &item, std::size_t label) {
    appr_alg->addPoint(item.data(), label);
  }

  static auto validateNewIndex(SEXP dim, SEXP max_elements, SEXP M,
                               SEXP ef_construction, SEXP random_seed,
                               bool has_random_seed) -> NewIndexConfig {
    NewIndexConfig config{
        static_cast<int>(check_whole_number(dim, "dimension", 1)),
        check_whole_number(max_elements, "capacity", 1),
        check_whole_number(M, "M", 2, 10000),
        check_whole_number(ef_construction, "ef_construction", 1), 100};
    if (has_random_seed) {
      config.random_seed = check_whole_number(random_seed, "random_seed", 0);
    }
    return config;
  }

  static auto validateLoadIndex(SEXP dim, SEXP path_to_index, SEXP max_elements,
                                bool has_max_elements) -> LoadIndexConfig {
    LoadIndexConfig config{
        static_cast<int>(check_whole_number(dim, "dimension", 1)),
        check_path(path_to_index, "path_to_index"), 0};
    if (has_max_elements) {
      config.max_elements = check_whole_number(max_elements, "capacity", 1);
    }
    return config;
  }

  auto checkedItem(SEXP item, const char *name) const -> std::vector<dist_t> {
    if (!is_numeric_value(item) ||
        Rf_getAttrib(item, R_DimSymbol) != R_NilValue ||
        Rf_xlength(item) != dim) {
      Rcpp::stop("%s has incorrect dimensions: expected a numeric vector of "
                 "length %d",
                 name, dim);
    }

    std::vector<dist_t> result(static_cast<std::size_t>(dim));
    for (int i = 0; i < dim; ++i) {
      result[static_cast<std::size_t>(i)] =
          static_cast<dist_t>(checked_float_at(item, i, name));
    }
    Normalizer<dist_t, DoNormalize>::normalize(result.data(), result.size());
    return result;
  }

  auto checkedItems(SEXP items, bool by_row, const char *name) const
      -> CheckedItems {
    if (!is_numeric_value(items) || !Rf_isMatrix(items)) {
      Rcpp::stop("%s must be a numeric matrix", name);
    }
    SEXP dimensions = Rf_getAttrib(items, R_DimSymbol);
    const int nrow = INTEGER_ELT(dimensions, 0);
    const int ncol = INTEGER_ELT(dimensions, 1);
    const int observed_dim = by_row ? ncol : nrow;
    const int nitems = by_row ? nrow : ncol;
    if (observed_dim != dim) {
      Rcpp::stop("%s have incorrect dimensions: expected dimension %d", name,
                 dim);
    }

    const std::size_t total =
        checked_product(static_cast<std::size_t>(nitems),
                        static_cast<std::size_t>(dim), "Input matrix");
    std::vector<dist_t> result(total);
    for (int item = 0; item < nitems; ++item) {
      const std::size_t output_offset =
          static_cast<std::size_t>(item) * static_cast<std::size_t>(dim);
      for (int coordinate = 0; coordinate < dim; ++coordinate) {
        const R_xlen_t input_offset =
            by_row ? item + static_cast<R_xlen_t>(nitems) * coordinate
                   : coordinate + static_cast<R_xlen_t>(dim) * item;
        result[output_offset + static_cast<std::size_t>(coordinate)] =
            static_cast<dist_t>(checked_float_at(items, input_offset, name));
      }
      Normalizer<dist_t, DoNormalize>::normalize(result.data() + output_offset,
                                                 dim);
    }
    return {nitems, std::move(result)};
  }

  void ensureCapacityFor(std::size_t nitems) const {
    const std::size_t index_start = currentSize();
    const std::size_t capacity = appr_alg->max_elements_;
    if (index_start > capacity || nitems > capacity - index_start) {
      Rcpp::stop("Index is too small to contain all items");
    }
  }

  void addItemsAdapter(SEXP items, bool by_row) {
    ensureUsable();
    CheckedItems checked = checkedItems(items, by_row, "Items to add");
    const std::size_t nitems = static_cast<std::size_t>(checked.nitems);
    ensureCapacityFor(nitems);
    const std::size_t index_start = currentSize();
    auto data_begin = checked.data.cbegin();
    std::atomic<bool> insertion_started{false};
    auto worker = [&](std::size_t begin, std::size_t end) {
      for (auto i = begin; i < end; i++) {
        auto first = data_begin + static_cast<std::size_t>(dim) * i;
        std::vector<dist_t> item_copy(first, first + dim);
        insertion_started.store(true, std::memory_order_relaxed);
        addItemImpl(item_copy, index_start + i);
      }
    };
    std::size_t parallel_begin = 0;
    try {
      if (index_start == 0 && nitems != 0) {
        worker(0, 1);
        parallel_begin = 1;
      }
      pforr::parallel_for(parallel_begin, nitems, worker, numThreads,
                          grainSize);
    } catch (...) {
      if (insertion_started.load(std::memory_order_relaxed)) {
        usable = false;
      }
      throw;
    }
  }

  auto checkK(SEXP value) const -> std::size_t {
    const std::size_t nnbrs = check_whole_number(value, "k", 1);
    const std::size_t total_count = currentSize();
    const std::size_t deleted_count = appr_alg->getDeletedCount();
    if (deleted_count > total_count) {
      Rcpp::stop("Index has an invalid deleted-item count");
    }
    const std::size_t active_count = total_count - deleted_count;
    if (nnbrs > active_count) {
      Rcpp::stop("k cannot be larger than active item count %llu; unable to "
                 "find requested results",
                 static_cast<unsigned long long>(active_count));
    }
    return nnbrs;
  }

  auto getAllNNsImpl(const std::vector<dist_t> &data, std::size_t nitems,
                     std::size_t nnbrs, bool by_row, bool include_distances,
                     std::vector<hnswlib::labeltype> &idx_vec,
                     std::vector<dist_t> &dist_vec) -> bool {
    std::vector<unsigned char> chunk_status(nitems, 1);
    auto data_begin = data.cbegin();
    auto worker = [&](std::size_t begin, std::size_t end,
                      std::size_t chunk_id) {
      std::vector<dist_t> distances;
      for (auto i = begin; i < end; i++) {
        auto first = data_begin + static_cast<std::size_t>(dim) * i;
        std::vector<dist_t> item_copy(first, first + dim);
        bool ok_row = true;
        std::vector<hnswlib::labeltype> nbr_labels =
            getNNsImpl(item_copy, nnbrs, include_distances, distances, ok_row);
        if (!ok_row) {
          chunk_status[chunk_id] = 0;
          break;
        }
        for (std::size_t k = 0; k < nnbrs; ++k) {
          const std::size_t output_offset =
              by_row ? k * nitems + i : nnbrs * i + k;
          idx_vec[output_offset] = nbr_labels[k];
          if (include_distances) {
            dist_vec[output_offset] = distances[k];
          }
        }
      }
    };
    pforr::parallel_for_indexed(0, nitems, worker, numThreads, grainSize);
    return std::all_of(chunk_status.cbegin(), chunk_status.cend(),
                       [](unsigned char status) { return status != 0; });
  }

  auto getAllNNsListAdapter(SEXP items, SEXP nnbrs, SEXP include_distances,
                            bool by_row) -> Rcpp::List {
    ensureUsable();
    CheckedItems checked = checkedItems(items, by_row, "Query items");
    const std::size_t checked_nnbrs = checkK(nnbrs);
    const bool checked_include_distances =
        check_logical(include_distances, "include_distances");
    const std::size_t nitems = static_cast<std::size_t>(checked.nitems);
    const std::size_t result_size =
        checked_product(nitems, checked_nnbrs, "Nearest-neighbor result");

    std::vector<hnswlib::labeltype> idx_vec(result_size);
    std::vector<dist_t> dist_vec(checked_include_distances ? result_size : 0);
    bool found_all =
        getAllNNsImpl(checked.data, nitems, checked_nnbrs, by_row,
                      checked_include_distances, idx_vec, dist_vec);
    if (!found_all) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }

    Rcpp::IntegerMatrix idx_result = makeIndexMatrix(
        idx_vec, checked.nitems, static_cast<int>(checked_nnbrs), by_row);
    auto result = Rcpp::List::create(Rcpp::Named("item") = idx_result);
    if (checked_include_distances) {
      DistanceProcess::process_distances(dist_vec);
      result["distance"] = makeDistanceMatrix(
          dist_vec, checked.nitems, static_cast<int>(checked_nnbrs), by_row);
    }
    return result;
  }

  auto getAllNNsAdapter(SEXP items, SEXP nnbrs, bool by_row)
      -> Rcpp::IntegerMatrix {
    ensureUsable();
    CheckedItems checked = checkedItems(items, by_row, "Query items");
    const std::size_t checked_nnbrs = checkK(nnbrs);
    const std::size_t nitems = static_cast<std::size_t>(checked.nitems);
    const std::size_t result_size =
        checked_product(nitems, checked_nnbrs, "Nearest-neighbor result");

    std::vector<hnswlib::labeltype> idx_vec(result_size);
    std::vector<dist_t> dist_vec;
    bool found_all = getAllNNsImpl(checked.data, nitems, checked_nnbrs, by_row,
                                   false, idx_vec, dist_vec);
    if (!found_all) {
      Rcpp::stop("Unable to find k results. Probably ef or M is too small");
    }
    return makeIndexMatrix(idx_vec, checked.nitems,
                           static_cast<int>(checked_nnbrs), by_row);
  }

  static auto makeIndexMatrix(const std::vector<hnswlib::labeltype> &values,
                              int nitems, int nnbrs, bool by_row)
      -> Rcpp::IntegerMatrix {
    Rcpp::IntegerMatrix result(by_row ? nitems : nnbrs,
                               by_row ? nnbrs : nitems);
    std::transform(
        values.begin(), values.end(), result.begin(),
        [](hnswlib::labeltype value) { return static_cast<int>(value); });
    return result;
  }

  static auto makeDistanceMatrix(const std::vector<dist_t> &values, int nitems,
                                 int nnbrs, bool by_row)
      -> Rcpp::NumericMatrix {
    Rcpp::NumericMatrix result(by_row ? nitems : nnbrs,
                               by_row ? nnbrs : nitems);
    std::copy(values.begin(), values.end(), result.begin());
    return result;
  }

  auto getItemsImpl(const std::vector<hnswlib::labeltype> &ids)
      -> std::vector<dist_t> {
    // this method assumes all the ids are valid
    const std::size_t nitems = ids.size();
    std::vector<dist_t> data(dim * nitems);

    auto worker = [&](std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i != end; i++) {
        auto obs = appr_alg->template getDataByLabel<dist_t>(ids[i]);
        std::copy(obs.begin(), obs.end(), data.begin() + i * dim);
      }
    };

    pforr::parallel_for(0, nitems, worker, numThreads, grainSize);

    return data;
  }
  int dim;
  bool usable;
  std::size_t numThreads;
  std::size_t grainSize;
  std::unique_ptr<Distance> space;
  std::unique_ptr<hnswlib::HierarchicalNSW<dist_t>> appr_alg;
};

using HnswL2 = Hnsw<float, hnswlib::L2Space, false, NoDistanceProcess>;
using HnswCosine =
    Hnsw<float, hnswlib::InnerProductSpace, true, NoDistanceProcess>;
using HnswIp =
    Hnsw<float, hnswlib::InnerProductSpace, false, NoDistanceProcess>;
using HnswEuclidean =
    Hnsw<float, hnswlib::L2Space, false, SquareRootDistanceProcess>;

RCPP_EXPOSED_CLASS_NODECL(HnswL2)
RCPP_MODULE(HnswL2) {
  Rcpp::class_<HnswL2>("HnswL2")
      .constructor<SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef")
      .constructor<SEXP, SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef, random seed")
      .constructor<SEXP, SEXP>("load raw index: dimension, filename")
      .constructor<SEXP, SEXP, SEXP>(
          "load raw index: dimension, filename, capacity")
      .method("setEf", &HnswL2::setEf, "set search ef")
      .method("addItem", &HnswL2::addItem, "add item")
      .method("addItems", &HnswL2::addItems,
              "add items where each item is stored row-wise")
      .method("addItemsCol", &HnswL2::addItemsCol,
              "add items where each item is stored column-wise")
      .method("getItems", &HnswL2::getItems,
              "retrieve vectors for one-based item identifiers")
      .method("save", &HnswL2::callSave, "save raw index to file")
      .method("getNNs", &HnswL2::getNNs,
              "retrieve nearest neighbors for one vector")
      .method("getNNsList", &HnswL2::getNNsList,
              "retrieve nearest neighbors and optional distances for one "
              "vector")
      .method("getAllNNs", &HnswL2::getAllNNs,
              "retrieve nearest neighbors for row-wise matrix items")
      .method("getAllNNsList", &HnswL2::getAllNNsList,
              "retrieve nearest neighbors and optional distances for row-wise "
              "matrix items")
      .method("getAllNNsCol", &HnswL2::getAllNNsCol,
              "retrieve nearest neighbors for column-wise matrix items; "
              "results are column-wise")
      .method("getAllNNsListCol", &HnswL2::getAllNNsListCol,
              "retrieve nearest neighbors and optional distances for "
              "column-wise matrix items; results are column-wise")
      .method("size", &HnswL2::size,
              "total items added, including deleted items")
      .method("setNumThreads", &HnswL2::setNumThreads,
              "set maximum threads; zero or one is serial")
      .method("setGrainSize", &HnswL2::setGrainSize,
              "set minimum items per thread; zero is treated as one")
      .method("markDeleted", &HnswL2::markDeleted,
              "mark a label deleted without reclaiming capacity")
      .method("resizeIndex", &HnswL2::resizeIndex,
              "change maximum index capacity");
}

RCPP_EXPOSED_CLASS_NODECL(HnswCosine)
RCPP_MODULE(HnswCosine) {
  Rcpp::class_<HnswCosine>("HnswCosine")
      .constructor<SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef")
      .constructor<SEXP, SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef, random seed")
      .constructor<SEXP, SEXP>("load raw index: dimension, filename")
      .constructor<SEXP, SEXP, SEXP>(
          "load raw index: dimension, filename, capacity")
      .method("setEf", &HnswCosine::setEf, "set search ef")
      .method("addItem", &HnswCosine::addItem, "add item")
      .method("addItems", &HnswCosine::addItems,
              "add items where each item is stored row-wise")
      .method("addItemsCol", &HnswCosine::addItemsCol,
              "add items where each item is stored column-wise")
      .method("getItems", &HnswCosine::getItems,
              "retrieve normalized vectors for one-based item identifiers")
      .method("save", &HnswCosine::callSave, "save raw index to file")
      .method("getNNs", &HnswCosine::getNNs,
              "retrieve nearest neighbors for one vector")
      .method("getNNsList", &HnswCosine::getNNsList,
              "retrieve nearest neighbors and optional distances for one "
              "vector")
      .method("getAllNNs", &HnswCosine::getAllNNs,
              "retrieve nearest neighbors for row-wise matrix items")
      .method("getAllNNsList", &HnswCosine::getAllNNsList,
              "retrieve nearest neighbors and optional distances for row-wise "
              "matrix items")
      .method("getAllNNsCol", &HnswCosine::getAllNNsCol,
              "retrieve nearest neighbors for column-wise matrix items; "
              "results are column-wise")
      .method("getAllNNsListCol", &HnswCosine::getAllNNsListCol,
              "retrieve nearest neighbors and optional distances for "
              "column-wise matrix items; results are column-wise")
      .method("size", &HnswCosine::size,
              "total items added, including deleted items")
      .method("setNumThreads", &HnswCosine::setNumThreads,
              "set maximum threads; zero or one is serial")
      .method("setGrainSize", &HnswCosine::setGrainSize,
              "set minimum items per thread; zero is treated as one")
      .method("markDeleted", &HnswCosine::markDeleted,
              "mark a label deleted without reclaiming capacity")
      .method("resizeIndex", &HnswCosine::resizeIndex,
              "change maximum index capacity");
}

RCPP_EXPOSED_CLASS_NODECL(HnswIp)
RCPP_MODULE(HnswIp) {
  Rcpp::class_<HnswIp>("HnswIp")
      .constructor<SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef")
      .constructor<SEXP, SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef, random seed")
      .constructor<SEXP, SEXP>("load raw index: dimension, filename")
      .constructor<SEXP, SEXP, SEXP>(
          "load raw index: dimension, filename, capacity")
      .method("setEf", &HnswIp::setEf, "set search ef")
      .method("addItem", &HnswIp::addItem, "add item")
      .method("addItems", &HnswIp::addItems,
              "add items where each item is stored row-wise")
      .method("addItemsCol", &HnswIp::addItemsCol,
              "add items where each item is stored column-wise")
      .method("getItems", &HnswIp::getItems,
              "retrieve vectors for one-based item identifiers")
      .method("save", &HnswIp::callSave, "save raw index to file")
      .method("getNNs", &HnswIp::getNNs,
              "retrieve nearest neighbors for one vector")
      .method("getNNsList", &HnswIp::getNNsList,
              "retrieve nearest neighbors and optional distances for one "
              "vector")
      .method("getAllNNs", &HnswIp::getAllNNs,
              "retrieve nearest neighbors for row-wise matrix items")
      .method("getAllNNsList", &HnswIp::getAllNNsList,
              "retrieve nearest neighbors and optional distances for row-wise "
              "matrix items")
      .method("getAllNNsCol", &HnswIp::getAllNNsCol,
              "retrieve nearest neighbors for column-wise matrix items; "
              "results are column-wise")
      .method("getAllNNsListCol", &HnswIp::getAllNNsListCol,
              "retrieve nearest neighbors and optional distances for "
              "column-wise matrix items; results are column-wise")
      .method("size", &HnswIp::size,
              "total items added, including deleted items")
      .method("setNumThreads", &HnswIp::setNumThreads,
              "set maximum threads; zero or one is serial")
      .method("setGrainSize", &HnswIp::setGrainSize,
              "set minimum items per thread; zero is treated as one")
      .method("markDeleted", &HnswIp::markDeleted,
              "mark a label deleted without reclaiming capacity")
      .method("resizeIndex", &HnswIp::resizeIndex,
              "change maximum index capacity");
}

RCPP_EXPOSED_CLASS_NODECL(HnswEuclidean)
RCPP_MODULE(HnswEuclidean) {
  Rcpp::class_<HnswEuclidean>("HnswEuclidean")
      .constructor<SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef")
      .constructor<SEXP, SEXP, SEXP, SEXP, SEXP>(
          "create index: dimension, capacity, M, construction ef, random seed")
      .constructor<SEXP, SEXP>("load raw index: dimension, filename")
      .constructor<SEXP, SEXP, SEXP>(
          "load raw index: dimension, filename, capacity")
      .method("setEf", &HnswEuclidean::setEf, "set search ef")
      .method("addItem", &HnswEuclidean::addItem, "add item")
      .method("addItems", &HnswEuclidean::addItems,
              "add items where each item is stored row-wise")
      .method("addItemsCol", &HnswEuclidean::addItemsCol,
              "add items where each item is stored column-wise")
      .method("getItems", &HnswEuclidean::getItems,
              "retrieve vectors for one-based item identifiers")
      .method("save", &HnswEuclidean::callSave, "save raw index to file")
      .method("getNNs", &HnswEuclidean::getNNs,
              "retrieve nearest neighbors for one vector")
      .method("getNNsList", &HnswEuclidean::getNNsList,
              "retrieve nearest neighbors and optional distances for one "
              "vector")
      .method("getAllNNs", &HnswEuclidean::getAllNNs,
              "retrieve nearest neighbors for row-wise matrix items")
      .method("getAllNNsList", &HnswEuclidean::getAllNNsList,
              "retrieve nearest neighbors and optional distances for row-wise "
              "matrix items")
      .method("getAllNNsCol", &HnswEuclidean::getAllNNsCol,
              "retrieve nearest neighbors for column-wise matrix items; "
              "results are column-wise")
      .method("getAllNNsListCol", &HnswEuclidean::getAllNNsListCol,
              "retrieve nearest neighbors and optional distances for "
              "column-wise matrix items; results are column-wise")
      .method("size", &HnswEuclidean::size,
              "total items added, including deleted items")
      .method("setNumThreads", &HnswEuclidean::setNumThreads,
              "set maximum threads; zero or one is serial")
      .method("setGrainSize", &HnswEuclidean::setGrainSize,
              "set minimum items per thread; zero is treated as one")
      .method("markDeleted", &HnswEuclidean::markDeleted,
              "mark a label deleted without reclaiming capacity")
      .method("resizeIndex", &HnswEuclidean::resizeIndex,
              "change maximum index capacity");
}
