#pragma once

// Standard
#include <algorithm>
#include <execution>
#include <vector>

/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 */
template <typename S, typename T>  // "S" is a scalar type; "T" is the type of data associated with
                                   // each interval
class IITree {
   public:
    void add(const S& start, const S& end, const T& data) { a.emplace_back(start, end, data); }
    void remove(size_t i) { a.erase(a.begin() + i); }
    void remove(std::vector<size_t>& indices) {
        std::sort(indices.begin(), indices.end(), std::greater<size_t>());
        for (size_t i : indices) a.erase(a.begin() + i);
    }

    void index() {
        std::sort(std::execution::par, a.begin(), a.end(), IntervalLess());
        max_level = index_core(a);
    }

    void indexNoSort() { max_level = index_core(a); }

    bool overlap(const S& start, const S& end, std::vector<size_t>& out) const {
        int t = 0;
        StackCell stack[64];
        out.clear();
        if (max_level < 0) return false;
        stack[t++] = StackCell(max_level, (1LL << max_level) - 1, 0);
        while (t) {
            StackCell z = stack[--t];
            if (z.k <= 3) {
                size_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL << (z.k + 1)) - 1;
                if (i1 >= a.size()) i1 = a.size();
                for (i = i0; i < i1 && a[i].start < end; ++i)
                    if (start < a[i].end) out.push_back(i);
            } else if (z.w == 0) {
                size_t y = z.x - (1LL << (z.k - 1));
                stack[t++] = StackCell(z.k, z.x, 1);
                if (y >= a.size() || a[y].max > start) stack[t++] = StackCell(z.k - 1, y, 0);
            } else if (z.x < a.size() && a[z.x].start < end) {
                if (start < a[z.x].end) out.push_back(z.x);
                stack[t++] = StackCell(z.k - 1, z.x + (1LL << (z.k - 1)), 0);
            }
        }
        return !out.empty();
    }

    size_t size() const { return a.size(); }
    const S& start(size_t i) const { return a[i].start; }
    void setStart(size_t i, const S& s) { a[i].start = s; }
    const S& end(size_t i) const { return a[i].end; }
    void setEnd(size_t i, const S& e) { a[i].end = e; }
    const T& data(size_t i) const { return a[i].data; }
    T& data(size_t i) { return a[i].data; }

   private:
    struct StackCell {
        size_t x{};
        int k{}, w{};
        StackCell() = default;
        StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_) {}
    };
    struct Interval {
        S start, end, max;
        T data;
        Interval(const S& start, const S& end, const T& data)
            : start(start), end(end), max(end), data(data) {}
    };
    struct IntervalLess {
        bool operator()(const Interval& lhs, const Interval& rhs) const {
            return lhs.start < rhs.start;
        }
    };

    std::vector<Interval> a;
    int max_level{};
    int index_core(std::vector<Interval>& intervalls) {
        size_t i{}, last_i{};
        S last{};
        size_t k{};
        if (intervalls.empty()) return -1;
        for (i = 0; i < intervalls.size(); i += 2)
            last_i = i, last = intervalls[i].max = intervalls[i].end;
        for (k = 1; 1UL << k <= intervalls.size(); ++k) {
            size_t x = 1UL << (k - 1), i0 = (x << 1) - 1, step = x << 2;
            for (i = i0; i < intervalls.size(); i += step) {
                S el = intervalls[i - x].max;
                S er = i + x < intervalls.size() ? intervalls[i + x].max : last;
                S e = intervalls[i].end;
                e = e > el ? e : el;
                e = e > er ? e : er;
                intervalls[i].max = e;
            }
            last_i = last_i >> k & 1 ? last_i - x : last_i + x;
            if (last_i < intervalls.size() && intervalls[last_i].max > last)
                last = intervalls[last_i].max;
        }
        return k - 1;
    }

   public:
    const std::vector<Interval>& intervals() const { return a; }
};
