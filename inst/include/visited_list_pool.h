/*
 * Modified by the RcppHNSW project on 2026-09-06.
 * Changes make initial visited-list pool construction exception-safe. See the
 * RcppHNSW COPYRIGHTS file and the source package's tools/vendor/patches
 * directory for exact changes.
 */
#pragma once

#include <memory>
#include <mutex>
#include <string.h>
#include <deque>

namespace hnswlib {
typedef unsigned short int vl_type;

class VisitedList {
 public:
    vl_type curV;
    vl_type *mass;
    unsigned int numelements;

    VisitedList(int numelements1) {
        curV = -1;
        numelements = numelements1;
        mass = new vl_type[numelements];
    }

    void reset() {
        curV++;
        if (curV == 0) {
            memset(mass, 0, sizeof(vl_type) * numelements);
            curV++;
        }
    }

    ~VisitedList() { delete[] mass; }
};
///////////////////////////////////////////////////////////
//
// Class for multi-threaded pool-management of VisitedLists
//
/////////////////////////////////////////////////////////

class VisitedListPool {
    std::deque<VisitedList *> pool;
    std::mutex poolguard;
    int numelements;

    void clearPool() noexcept {
        while (pool.size()) {
            VisitedList *rez = pool.front();
            pool.pop_front();
            delete rez;
        }
    }

 public:
    VisitedListPool(int initmaxpools, int numelements1) {
        numelements = numelements1;
        try {
            for (int i = 0; i < initmaxpools; i++) {
                std::unique_ptr<VisitedList> visited_list(new VisitedList(numelements));
                pool.push_front(visited_list.get());
                visited_list.release();
            }
        } catch (...) {
            clearPool();
            throw;
        }
    }

    VisitedList *getFreeVisitedList() {
        VisitedList *rez;
        {
            std::unique_lock <std::mutex> lock(poolguard);
            if (pool.size() > 0) {
                rez = pool.front();
                pool.pop_front();
            } else {
                rez = new VisitedList(numelements);
            }
        }
        rez->reset();
        return rez;
    }

    void releaseVisitedList(VisitedList *vl) {
        std::unique_lock <std::mutex> lock(poolguard);
        pool.push_front(vl);
    }

    ~VisitedListPool() { clearPool(); }
};
}  // namespace hnswlib
