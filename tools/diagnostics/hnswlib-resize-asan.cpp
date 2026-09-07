#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>

namespace allocation_probe {

bool enabled = false;
std::size_t allocation_count = 0;
std::size_t failure_index = (std::numeric_limits<std::size_t>::max)();
std::ptrdiff_t active_allocations = 0;

auto should_fail() noexcept -> bool {
  if (!enabled) {
    return false;
  }
  ++allocation_count;
  if (allocation_count == failure_index) {
    enabled = false;
    return true;
  }
  return false;
}

void begin(std::size_t fail_at) noexcept {
  allocation_count = 0;
  failure_index = fail_at;
  enabled = true;
}

void end() noexcept { enabled = false; }

void record_allocation() noexcept { ++active_allocations; }

void release(void *memory) noexcept {
  if (memory != nullptr) {
    --active_allocations;
  }
  std::free(memory);
}

} // namespace allocation_probe

void *operator new(std::size_t size) {
  if (allocation_probe::should_fail()) {
    throw std::bad_alloc();
  }
  if (void *memory = std::malloc(size)) {
    allocation_probe::record_allocation();
    return memory;
  }
  throw std::bad_alloc();
}

void *operator new[](std::size_t size) {
  if (allocation_probe::should_fail()) {
    throw std::bad_alloc();
  }
  if (void *memory = std::malloc(size)) {
    allocation_probe::record_allocation();
    return memory;
  }
  throw std::bad_alloc();
}

void operator delete(void *memory) noexcept {
  allocation_probe::release(memory);
}
void operator delete[](void *memory) noexcept {
  allocation_probe::release(memory);
}
void operator delete(void *memory, std::size_t) noexcept {
  allocation_probe::release(memory);
}
void operator delete[](void *memory, std::size_t) noexcept {
  allocation_probe::release(memory);
}

extern "C" void *__real_realloc(void *memory, std::size_t size);

extern "C" void *__wrap_realloc(void *memory, std::size_t size) {
  if (allocation_probe::should_fail()) {
    return nullptr;
  }
  return __real_realloc(memory, size);
}

#include "hnswlib.h"

namespace {

constexpr std::size_t original_capacity = 8;
constexpr std::size_t stored_items = 2;

auto parse_failure_index(const char *text) -> std::size_t {
  const std::string value(text);
  if (value.empty()) {
    throw std::invalid_argument("failure index must not be empty");
  }
  std::size_t result = 0;
  for (const char character : value) {
    if (character < '0' || character > '9') {
      throw std::invalid_argument("failure index must be a positive integer");
    }
    const std::size_t digit = static_cast<std::size_t>(character - '0');
    if (result > ((std::numeric_limits<std::size_t>::max)() - digit) / 10) {
      throw std::invalid_argument("failure index is too large");
    }
    result = result * 10 + digit;
  }
  if (result == 0) {
    throw std::invalid_argument("failure index must be positive");
  }
  return result;
}

void verify_stored_state(hnswlib::HierarchicalNSW<float> &index,
                         const float *query) {
  if (index.getCurrentElementCount() != stored_items ||
      index.element_levels_.size() < stored_items ||
      index.link_list_locks_.size() < stored_items ||
      index.data_level0_memory_ == nullptr || index.linkLists_ == nullptr) {
    throw std::runtime_error("resize failure left stored items unsafe");
  }
  if (index.searchKnn(query, 1).empty()) {
    throw std::runtime_error("resize failure left stored items unsearchable");
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2 || argc > 3) {
      throw std::invalid_argument(
          "usage: hnswlib-resize-asan zero|growth|shrink [failure-index]");
    }

    const std::string mode(argv[1]);
    const bool zero_mode = mode == "zero";
    if (zero_mode && argc != 2) {
      throw std::invalid_argument("zero mode does not take a failure index");
    }
    if (!zero_mode && mode != "growth" && mode != "shrink") {
      throw std::invalid_argument("mode must be zero, growth, or shrink");
    }
    const bool inject_failure = argc == 3;
    const std::size_t failure_index =
        inject_failure ? parse_failure_index(argv[2])
                       : (std::numeric_limits<std::size_t>::max)();
    const std::ptrdiff_t allocation_baseline =
        allocation_probe::active_allocations;
    std::size_t observed_allocation_count = 0;

    {
      hnswlib::L2Space space(2);
      hnswlib::HierarchicalNSW<float> index(&space, original_capacity);
      const float first[2] = {1.0F, 0.0F};
      const float second[2] = {0.0F, 1.0F};
      index.addPoint(first, 0);
      index.addPoint(second, 1);

      if (zero_mode) {
        bool rejected = false;
        try {
          index.resizeIndex(0);
        } catch (const std::runtime_error &) {
          rejected = true;
        }
        if (!rejected || index.getMaxElements() != original_capacity) {
          throw std::runtime_error(
              "zero-capacity resize was not rejected safely");
        }
        verify_stored_state(index, first);
      } else {
        const std::size_t target = mode == "growth" ? 16 : stored_items;

        allocation_probe::begin(failure_index);
        bool failed = false;
        try {
          index.resizeIndex(target);
        } catch (const std::bad_alloc &) {
          failed = true;
        } catch (const std::runtime_error &) {
          failed = true;
        }
        allocation_probe::end();
        observed_allocation_count = allocation_probe::allocation_count;

        if (!inject_failure) {
          if (failed || index.getMaxElements() != target) {
            throw std::runtime_error("uninjected resize did not succeed");
          }
        } else {
          if (!failed) {
            throw std::runtime_error(
                "requested allocation failure was not reached");
          }
          verify_stored_state(index, first);
        }
      }
    }

    if (allocation_probe::active_allocations != allocation_baseline) {
      throw std::runtime_error("resize lifecycle leaked C++ allocations");
    }
    if (zero_mode) {
      std::cout << "zero-capacity resize rejected before allocation\n";
    } else if (!inject_failure) {
      std::cout << mode << " allocations=" << observed_allocation_count << '\n';
    } else {
      std::cout << mode << " failure=" << failure_index
                << " safely destructible\n";
    }
    return 0;
  } catch (const std::exception &error) {
    allocation_probe::end();
    std::cerr << "hnswlib-resize-asan: " << error.what() << '\n';
    return 2;
  }
}
