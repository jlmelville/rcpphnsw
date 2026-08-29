#include "hnswlib.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr std::size_t dimension = 8;
constexpr std::size_t default_rounds = 4;
constexpr std::size_t default_threads = 8;
constexpr std::size_t default_points_per_thread = 128;
constexpr std::size_t default_m = 8;
constexpr std::size_t default_seed = 100;
constexpr std::size_t max_rounds = 100;
constexpr std::size_t max_threads = 64;
constexpr std::size_t max_points_per_thread = 4096;
constexpr std::size_t max_m = 256;

auto parse_count(const char *text, const char *name, std::size_t minimum,
                 std::size_t maximum) -> std::size_t {
  const std::string value(text);
  if (value.empty()) {
    throw std::invalid_argument(std::string(name) + " must not be empty");
  }
  for (const char character : value) {
    if (character < '0' || character > '9') {
      throw std::invalid_argument(std::string(name) +
                                  " must be a whole number");
    }
  }

  std::size_t consumed = 0;
  unsigned long long parsed = 0;
  try {
    parsed = std::stoull(value, &consumed);
  } catch (const std::exception &) {
    throw std::invalid_argument(std::string(name) + " is out of range");
  }
  if (consumed != value.size() || parsed < minimum || parsed > maximum) {
    throw std::invalid_argument(std::string(name) + " must be in " +
                                std::to_string(minimum) + ".." +
                                std::to_string(maximum));
  }
  return static_cast<std::size_t>(parsed);
}

auto make_point(std::size_t label) -> std::array<float, dimension> {
  std::array<float, dimension> point{};
  for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
    point[coordinate] = static_cast<float>((label + 1) * (coordinate + 3)) /
                        static_cast<float>(dimension + coordinate + 1);
  }
  return point;
}

void run_round(std::size_t thread_count, std::size_t points_per_thread,
               std::size_t m, std::size_t seed) {
  if (thread_count >
      (std::numeric_limits<std::size_t>::max() - 1) / points_per_thread) {
    throw std::invalid_argument("requested index capacity is out of range");
  }
  const std::size_t capacity = 1 + thread_count * points_per_thread;
  hnswlib::L2Space space(dimension);
  hnswlib::HierarchicalNSW<float> index(&space, capacity, m, 32, seed);

  const auto first_point = make_point(0);
  index.addPoint(first_point.data(), 0);
  const int initial_max_level = index.maxlevel_;

  std::atomic<std::size_t> ready{0};
  std::atomic<bool> start{false};
  std::exception_ptr worker_error = nullptr;
  std::mutex worker_error_mutex;
  std::vector<std::thread> workers;
  workers.reserve(thread_count);

  try {
    for (std::size_t thread_id = 0; thread_id < thread_count; ++thread_id) {
      workers.emplace_back([&, thread_id]() {
        ready.fetch_add(1, std::memory_order_release);
        while (!start.load(std::memory_order_acquire)) {
          std::this_thread::yield();
        }
        try {
          for (std::size_t point_id = 0; point_id < points_per_thread;
               ++point_id) {
            const std::size_t label =
                1 + thread_id * points_per_thread + point_id;
            const auto point = make_point(label);
            index.addPoint(point.data(), label);
          }
        } catch (...) {
          std::lock_guard<std::mutex> lock(worker_error_mutex);
          if (worker_error == nullptr) {
            worker_error = std::current_exception();
          }
        }
      });
    }
  } catch (...) {
    start.store(true, std::memory_order_release);
    for (auto &worker : workers) {
      worker.join();
    }
    throw;
  }

  while (ready.load(std::memory_order_acquire) != thread_count) {
    std::this_thread::yield();
  }
  start.store(true, std::memory_order_release);
  for (auto &worker : workers) {
    worker.join();
  }
  if (worker_error != nullptr) {
    std::rethrow_exception(worker_error);
  }
  if (index.getCurrentElementCount() != capacity) {
    throw std::runtime_error("concurrent insertion count did not match");
  }
  if (index.maxlevel_ <= initial_max_level) {
    throw std::runtime_error("workload did not exercise maximum-level growth");
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc > 6) {
    std::cerr << "usage: " << argv[0]
              << " [rounds [threads [points-per-thread [M [seed]]]]]"
              << std::endl;
    return 2;
  }

  try {
    const std::size_t rounds =
        argc > 1 ? parse_count(argv[1], "rounds", 1, max_rounds)
                 : default_rounds;
    const std::size_t threads =
        argc > 2 ? parse_count(argv[2], "threads", 2, max_threads)
                 : default_threads;
    const std::size_t points_per_thread =
        argc > 3 ? parse_count(argv[3], "points-per-thread", 1,
                               max_points_per_thread)
                 : default_points_per_thread;
    const std::size_t m =
        argc > 4 ? parse_count(argv[4], "M", 2, max_m) : default_m;
    const std::size_t seed =
        argc > 5 ? parse_count(argv[5], "seed", 0,
                               std::numeric_limits<unsigned int>::max())
                 : default_seed;

    for (std::size_t round = 0; round < rounds; ++round) {
      run_round(threads, points_per_thread, m, seed);
    }
    std::cout << "completed " << rounds << " round(s) with " << threads
              << " concurrent inserters and " << points_per_thread
              << " points per inserter; M = " << m << ", seed = " << seed
              << "; maximum-level growth observed in every round" << std::endl;
  } catch (const std::exception &error) {
    std::cerr << "hnswlib RNG TSan diagnostic: " << error.what() << std::endl;
    return 2;
  }
  return 0;
}
