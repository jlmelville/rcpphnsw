// Taken from RcppParallel.h and then modified slightly to rename header guards
// and namespaces to avoid any potential clashes. RcppParallel is licensed under
// GPLv2 or later:

// pfor.h a version of parallel for based on RcppParallel
// Copyright (C) 2020 James Melville
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#ifndef PFORR
#define PFORR

#include <thread>
#include <utility>
#include <vector>

namespace pforr {

using IndexRange = std::pair<std::size_t, std::size_t>;

template <typename Worker>
auto worker_thread(Worker &worker, const IndexRange &range) -> void {
  try {
    worker(range.first, range.second);
  } catch (...) {
  }
}

// Function to calculate the ranges for a given input
inline auto split_input_range(const IndexRange &range, std::size_t n_threads,
                              std::size_t grain_size)
  -> std::vector<IndexRange> {

    // compute grain_size (including enforcing requested minimum)
    std::size_t length = range.second - range.first;
    if (n_threads == 1)
      grain_size = length;
    else if ((length % n_threads) == 0) // perfect division
      grain_size = (std::max)(length / n_threads, grain_size);
    else // imperfect division, divide by threads - 1
      grain_size = (std::max)(length / (n_threads - 1), grain_size);

    // allocate ranges
    std::vector<IndexRange> ranges;
    std::size_t begin = range.first;
    while (begin < range.second) {
      std::size_t end = (std::min)(begin + grain_size, range.second);
      ranges.emplace_back(begin, end);
      begin = end;
    }

    return ranges;
  }

// Execute the Worker over the IndexRange in parallel
template <typename Worker>
inline void parallel_for(std::size_t begin, std::size_t end, Worker &worker,
                         std::size_t n_threads, std::size_t grain_size = 1) {
  if (n_threads == 0) {
    worker(begin, end);
    return;
  }
  // split the work
  IndexRange input_range(begin, end);
  std::vector<IndexRange> ranges =
    split_input_range(input_range, n_threads, grain_size);

  std::vector<std::thread> threads;
  threads.reserve(ranges.size());
  for (auto &range : ranges) {
    threads.push_back(
      std::thread(&worker_thread<Worker>, std::ref(worker), range));
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return;
}

} // namespace pforr

#endif // PFORR