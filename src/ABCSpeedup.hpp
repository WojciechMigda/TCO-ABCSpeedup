/*******************************************************************************
 * Copyright (c) 2015 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the GNU LGPL v3
 *******************************************************************************
 *
 * Filename: ABCSpeedup.hpp
 *
 * Description:
 *      description
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2015-01-15   wm              Initial version
 *
 ******************************************************************************/

/*
 * IDEAS:
 * 1. naive: return separate results from each batch (testing the baseline)
 * 2. merge clusters from batches, skip those with distance > thr1, merge
 *    those with distance < thr2, for the gray zone: either leave them separate
 *    or repeat the process:)
 * 3. another batch-based approach: align dendrograms from batches, then apply
 *    flat clusters on top of the combined result, (CCLS-11-04 paper ???)
 * 4. probabilistic: create random batch, repeat clustering and collect results,
 *    somehow use the statistics (which?) to retrieve cluster assignments. (too slow)
 */
#ifndef ABCSPEEDUP_HPP_
#define ABCSPEEDUP_HPP_

//#pragma GCC optimize ( "-Ofast" )
#pragma GCC optimize ( "-ffast-math" )
#pragma GCC optimize ( "-fopenmp" )
#pragma GCC target ( "sse2" )
//#pragma GCC target ( "inline-stringops-dynamically" ) // not really
//#pragma GCC target ( "inline-all-stringops" ) // not really
//#pragma GCC target ( "fpmath=sse" )

#include "timestamp.hpp"
#include "json.hpp"
#include "seq.hpp"

#include <vector>
#include <iostream>
#include <string>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include <cassert>
#include <cfenv>
#include <iterator>

typedef float real_t;
//typedef int16_t real_t;

void done(void)
{
    std::cerr << "Done!" << std::endl << std::endl;
}

std::string vs_concat(std::vector<std::string> && vs)
{
    std::string result;

    std::for_each(vs.begin(), vs.end(),
        [&result](std::string & s)
        {
            std::string s_(std::move(s));
            const char * trimpoint = std::find_if(s_.c_str(), s_.c_str() + s_.length(), [](const char ch){return ch != '0';});
            result.append(trimpoint);
        }
    );

    return result;
}

enum {CACHE_LINE_SIZE = 32};

std::size_t fastLD(const std::string& s1, const std::string& s2)
{
    struct
    {
        uint8_t str1[CACHE_LINE_SIZE];
        uint8_t str2[CACHE_LINE_SIZE];
        uint8_t rowA[CACHE_LINE_SIZE * 2];
        uint8_t rowB[CACHE_LINE_SIZE * 2];
    } local __attribute__((aligned(32)));

    uint8_t * v0 = local.rowA;
    uint8_t * v1 = local.rowB;

    std::copy(s1.c_str(), s1.c_str() + s1.length(), local.str1);
    std::copy(s2.c_str(), s2.c_str() + s2.length(), local.str2);
    std::iota(v0, v0 + sizeof (local.rowA), 0);

    for (std::size_t iidx = 0; iidx < s1.length(); ++iidx)
    {
        v1[0] = iidx + 1;

        for (std::size_t jidx = 0; jidx < s2.length(); ++jidx)
        {
            bool cost = local.str1[iidx] != local.str2[jidx];

            const auto temp = std::min(v1[jidx] + 1, v0[jidx + 1] + 1); // dodaj 1 poza min

            v1[jidx + 1] = std::min(temp, v0[jidx] + cost);
        }

        std::swap(v0, v1);
    }

    return v0[s2.length()];
}

std::size_t LevenshteinDistance(const std::string& s1, const std::string& s2)
{
    if (s1.length() <= CACHE_LINE_SIZE && s2.length() <= CACHE_LINE_SIZE)
    {
        return fastLD(s1, s2);
    }

    std::size_t dp[s1.length() + 1][s2.length() + 1];
    for (std::size_t i = 0; i <= s1.length(); i++)
    {
        dp[i][0] = i;
    }
    for (std::size_t i = 0; i <= s2.length(); i++)
    {
        dp[0][i] = i;
    }
    for (std::size_t i = 1; i <= s1.length(); i++)
    {
        for (std::size_t j = 1; j <= s2.length(); j++)
        {
            dp[i][j] = std::min(dp[i - 1][j] + 1, dp[i][j - 1] + 1);
            dp[i][j] = std::min(dp[i][j],
                dp[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1));
        }
    }
    return dp[s1.length()][s2.length()];
}

double GlobalAlignment(
    const std::string& s1,
    const std::string& s2,
    double match,
    double mismatch,
    double gap,
    double extend)
{
    assert(gap == extend);
    double dp[s1.length() + 1][s2.length() + 1];
    for (size_t i = 0; i <= s1.length(); i++)
    {
        dp[i][0] = gap * i;
    }
    for (size_t i = 0; i <= s2.length(); i++)
    {
        dp[0][i] = gap * i;
    }
    for (size_t i = 1; i <= s1.length(); i++)
    {
        for (size_t j = 1; j <= s2.length(); j++)
        {
            dp[i][j] = dp[i - 1][j - 1]
                + (s1[i - 1] == s2[j - 1] ? match : mismatch);
            dp[i][j] = std::max(dp[i][j], dp[i - 1][j] + gap);
            dp[i][j] = std::max(dp[i][j], dp[i][j - 1] + gap);
        }
    }
    return dp[s1.length()][s2.length()];
}

double GetLD(const Seq& s1, const Seq& s2)
{
    if (s1.junc.length() == s2.junc.length())
    {
        return s1.junc.length() - GlobalAlignment(s1.junc, s2.junc, 1, 0, -50, -50);
    }
    else
    {
        return LevenshteinDistance(s1.junc, s2.junc);
    }
}

double vCompare(const Seq& s1, const Seq& s2)
{
    if (0 != s1.v_gene.compare(s2.v_gene))
    {
        return 8;
    }
    if (0 != s1.v_all.compare(s2.v_all))
    {
        return 1;
    }
    return 0;
}

double jCompare(const Seq& s1, const Seq& s2)
{
    if (0 != s1.j_gene.compare(s2.j_gene))
    {
        return 8;
    }
    if (0 != s1.j_all.compare(s2.j_all))
    {
        return 1;
    }
    return 0;
}

double SharedMuts(const Seq& s1, const Seq& s2)
{
    if (0 == s1.id.compare(s2.id))
    {
        return 0.;
    }
    double bonus = 0.;
    size_t p1 = 0, p2 = 0;

    const std::size_t s1_muts_size = s1.muts.size();
    const std::size_t s2_muts_size = s2.muts.size();
    const decltype(&s1.muts[0]) s1_muts_p = s1.muts.data();
    const decltype(&s2.muts[0]) s2_muts_p = s2.muts.data();

    while (p1 < s1_muts_size && p2 < s2_muts_size)
    {
        if (s1_muts_p[p1] < s2_muts_p[p2])
        {
            p1++;
        }
        else if (s2_muts_p[p2] < s1_muts_p[p1])
        {
            p2++;
        }
        else
        {
            p1++;
            bonus += 0.35;
        }
    }
    return bonus;
}

//constexpr real_t CUTOFF = 32;
constexpr real_t CUTOFF = 0.32;

real_t GetScore(const Seq& s1, const Seq& s2)
{
    if (0 == s1.id.compare(s2.id))
    {
        return 0.;
    }

    const double vPenalty = vCompare(s1, s2);
    const double jPenalty = jCompare(s1, s2);
    const double lenPenalty = std::fabs(0. + s1.junc.length() - s2.junc.length()) * 2;
    const double editLength = (double)std::min(s1.junc.length(), s2.junc.length());

//    const double LD = GetLD(s1, s2);
    double LD;
//    const double LD = 0.0;
    if (vPenalty > 7 || jPenalty > 7 || lenPenalty > 7 )
    {
        LD = 6;
//        LD = GetLD(s1, s2);
    }
    else
    {
        LD = GetLD(s1, s2);
    }

    double mutBonus = SharedMuts(s1, s2);
//    double mutBonus = 0.0;

    if (mutBonus > (LD + vPenalty + jPenalty))
    {
        mutBonus = (LD + vPenalty + jPenalty - 0.001);
    }

//    const double score_no_LD = vPenalty + jPenalty + lenPenalty - mutBonus;
//
//    const double score = (score_no_LD < CUTOFF * editLength) ?
//        GetLD(s1, s2) + score_no_LD
//        :
//        score_no_LD;

    const double score = (LD + vPenalty + jPenalty + lenPenalty - mutBonus);
//    std::cerr << score << " " << LD << " " << vPenalty << " " << jPenalty << " " << lenPenalty << " " << mutBonus << std::endl;
    return (score / editLength);
}

std::unique_ptr<real_t[]> &&
BuildCondensedMatrix(
    std::vector<Seq *>::const_iterator cbegin,
    std::vector<Seq *>::const_iterator cend,
    std::unique_ptr<real_t[]> && scores)
{
    std::cerr << "Computing condensed matrix" << std::endl;

    const std::size_t NELEM = std::distance(cbegin, cend);

    std::vector<std::size_t> sumTable;
    sumTable.reserve(NELEM);
    sumTable.push_back(NELEM - 1);
    for (size_t i = 1; i < NELEM; i++)
    {
        sumTable.push_back(sumTable.back() + NELEM - i - 1);
    }

    const size_t nPairs = NELEM * (NELEM - 1) / 2;

#pragma omp parallel for
    for (std::size_t i = 0; i < nPairs; i++)
    {
        std::size_t first = std::lower_bound(sumTable.begin(), sumTable.end(), i + 1) -
            sumTable.begin();
        std::size_t second = NELEM - 1 - (sumTable[first] - (i + 1));
        //fprintf(stderr, "%lu %lu\n", first, second);
        scores[i] = GetScore(**(cbegin + first), **(cbegin + second));
    }

    return std::move(scores);
}

typedef int_fast32_t t_index;
typedef double t_float;

struct node
{
    t_index node1;
    t_index node2;
    t_float dist;

    inline friend bool operator<(const node & a, const node & b)
    {
        return (a.dist < b.dist);
    }
};

class cluster_result
{
private:
    std::unique_ptr<node[]> Z;
    t_index pos;

public:
    explicit cluster_result(const t_index size)
    :
        Z(new node[size]),
        pos(0)
    {
    }

    void append(const t_index node1, const t_index node2, const t_float dist)
    {
        Z[pos].node1 = node1;
        Z[pos].node2 = node2;
        Z[pos].dist = dist;
        ++pos;
    }

    node * operator[](const t_index idx) const
    {
        return &Z[idx];
    }

    /* Define several methods to postprocess the distances. All these functions
     are monotone, so they do not change the sorted order of distances. */

    void sqrt() const
    {
        std::for_each(&Z[0], &Z[pos],
            [](node & value)
            {
                value.dist = ::sqrt(value.dist);
            }
        );
    }

    void sqrt(const t_float) const
    { // ignore the argument
        sqrt();
    }

    void sqrtdouble(const t_float) const
    { // ignore the argument
        std::for_each(&Z[0], &Z[pos],
            [](node & value)
            {
                value.dist = ::sqrt(2 * value.dist);
            }
        );
    }

#ifdef R_pow
#define my_pow R_pow
#else
#define my_pow pow
#endif

    void power(const t_float p) const
    {
        t_float const q = 1 / p;
        std::for_each(&Z[0], &Z[pos],
            [&q](node & value)
            {
                value.dist = my_pow(value.dist, q);
            }
        );
    }

    void plusone(const t_float) const
    { // ignore the argument
        std::for_each(&Z[0], &Z[pos],
            [](node & value)
            {
                value.dist += 1;
            }
        );
    }

    void divide(const t_float denom) const
    {
        std::for_each(&Z[0], &Z[pos],
            [&denom](node & value)
            {
                value.dist /= denom;
            }
        );
    }
};

class nan_error{};
#ifdef FE_INVALID
class fenv_error{};
#endif

class doubly_linked_list
{
    /*
     Class for a doubly linked list. Initially, the list is the integer range
     [0, size]. We provide a forward iterator and a method to delete an index
     from the list.

     Typical use: for (i=L.start; L<size; i=L.succ[I])
     or
     for (i=somevalue; L<size; i=L.succ[I])
     */
public:
    t_index start;
    std::unique_ptr<t_index[]> succ;

private:
    std::unique_ptr<t_index[]> pred;
    // Not necessarily private, we just do not need it in this instance.

public:
    doubly_linked_list(const t_index size)
    // Initialize to the given size.
    :
        start(0),
        succ(new t_index[size + 1]),
        pred(new t_index[size + 1])
    {
        for (t_index i = 0; i < size; ++i)
        {
            pred[i + 1] = i;
            succ[i] = i + 1;
        }
        // pred[0] is never accessed!
        //succ[size] is never accessed!
    }

    ~doubly_linked_list()
    {
    }

    void remove(const t_index idx)
    {
        // Remove an index from the list.
        if (idx == start)
        {
            start = succ[idx];
        }
        else
        {
            succ[pred[idx]] = succ[idx];
            pred[succ[idx]] = pred[idx];
        }
        succ[idx] = 0; // Mark as inactive
    }

    bool is_inactive(t_index idx) const
    {
        return (succ[idx] == 0);
    }
};

/*
  Lookup function for a union-find data structure.

  The function finds the root of idx by going iteratively through all
  parent elements until a root is found. An element i is a root if
  nodes[i] is zero. To make subsequent searches faster, the entry for
  idx and all its parents is updated with the root element.
 */
class union_find
{
private:
    std::unique_ptr<t_index[]> parent;
    t_index nextparent;

public:
    union_find(const t_index size)
    :
        parent(new t_index[size > 0 ? 2 * size - 1 : 0]),
        nextparent(size)
    {
        std::fill(&parent[0], &parent[size > 0 ? 2 * size - 1 : 0], 0);
    }

    t_index Find(t_index idx) const
    {
        if (parent[idx] != 0)
        { // a → b
            t_index p = idx;
            idx = parent[idx];
            if (parent[idx] != 0)
            { // a → b → c
                do
                {
                    idx = parent[idx];
                } while (parent[idx] != 0);
                do
                {
                    t_index tmp = parent[p];
                    parent[p] = idx;
                    p = tmp;
                } while (parent[p] != idx);
            }
        }
        return idx;
    }

    void Union(const t_index node1, const t_index node2)
    {
        parent[node1] = parent[node2] = nextparent++;
    }
};

inline void f_average(
    real_t * const b,
    const real_t a,
    const t_float s,
    const t_float t)
{
    *b = s * a + t * (*b);
#ifndef FE_INVALID
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif
    if (fc_isnan(*b))
    {
        throw(nan_error());
    }
#if HAVE_DIAGNOSTIC
#pragma GCC diagnostic pop
#endif
#endif
}

template <typename t_members>
void NN_chain_core(
    const t_index N,
    std::unique_ptr<real_t[]> && D,
    std::unique_ptr<t_index[]> && members,
    cluster_result & Z2)
{
    /*
    N: integer
    D: condensed distance matrix N*(N-1)/2
    Z2: output data structure

    This is the NN-chain algorithm, described on page 86 in the following book:

    Fionn Murtagh, Multidimensional Clustering Algorithms,
    Vienna, Würzburg: Physica-Verlag, 1985.
     */
    const auto D_ = [&N, &D](const std::size_t r_, const std::size_t c_) -> real_t &
    {
        const std::size_t index = (static_cast<std::ptrdiff_t>(2*N-3-(r_))*(r_)>>1)+(c_)-1;

        return D[index];
    };

    t_index i;

    std::unique_ptr<t_index[]> NN_chain(new t_index[N]);
    t_index NN_chain_tip = 0;

    t_index idx1, idx2;

    t_float size1, size2;
    doubly_linked_list active_nodes(N);

    t_float min;

#ifdef FE_INVALID
    if (feclearexcept(FE_INVALID))
        throw fenv_error();
#endif

    for (t_index j = 0; j < N - 1; ++j)
    {
        if (NN_chain_tip <= 3)
        {
            NN_chain[0] = idx1 = active_nodes.start;
            NN_chain_tip = 1;

            idx2 = active_nodes.succ[idx1];
            min = D_(idx1, idx2);
            for (i = active_nodes.succ[idx2]; i < N; i = active_nodes.succ[i])
            {
                if (D_(idx1, i) < min)
                {
                    min = D_(idx1, i);
                    idx2 = i;
                }
            }
        }  // a: idx1   b: idx2
        else
        {
            NN_chain_tip -= 3;
            idx1 = NN_chain[NN_chain_tip - 1];
            idx2 = NN_chain[NN_chain_tip];
            min = idx1 < idx2 ? D_(idx1, idx2) : D_(idx2, idx1);
        }  // a: idx1   b: idx2

        do
        {
            NN_chain[NN_chain_tip] = idx2;

            for (i = active_nodes.start; i < idx2; i = active_nodes.succ[i])
            {
                if (D_(i, idx2) < min)
                {
                    min = D_(i, idx2);
                    idx1 = i;
                }
            }
            for (i = active_nodes.succ[idx2]; i < N; i = active_nodes.succ[i])
            {
                if (D_(idx2, i) < min)
                {
                    min = D_(idx2, i);
                    idx1 = i;
                }
            }

            idx2 = idx1;
            idx1 = NN_chain[NN_chain_tip++];

        } while (idx2 != NN_chain[NN_chain_tip - 2]);

//        std::cerr << "Z2.append " << idx1 << " " << idx2 << " " << min << std::endl;
        Z2.append(idx1, idx2, min);

        if (idx1 > idx2)
        {
            t_index tmp = idx1;
            idx1 = idx2;
            idx2 = tmp;
        }

        size1 = static_cast<t_float>(members[idx1]);
        size2 = static_cast<t_float>(members[idx2]);
        members[idx2] += members[idx1];

        // Remove the smaller index from the valid indices (active_nodes).
        active_nodes.remove(idx1);

        t_float s = size1 / (size1 + size2);
        t_float t = size2 / (size1 + size2);
        for (i = active_nodes.start; i < idx1; i = active_nodes.succ[i])
            f_average(&D_(i, idx2), D_(i, idx1), s, t);
        // Update the distance matrix in the range (idx1, idx2).
        for (; i < idx2; i = active_nodes.succ[i])
            f_average(&D_(i, idx2), D_(idx1, i), s, t);
        // Update the distance matrix in the range (idx2, N).
        for (i = active_nodes.succ[idx2]; i < N; i = active_nodes.succ[i])
            f_average(&D_(idx2, i), D_(idx1, i), s, t);

    }
#ifdef FE_INVALID
    if (fetestexcept(FE_INVALID))
        throw fenv_error();
#endif
}

class linkage_output
{
private:
    t_float * Z;

public:
    linkage_output(t_float * const Z_) :
        Z(Z_)
    {
    }

    void append(const t_index node1, const t_index node2, const t_float dist,
        const t_float size)
    {
        if (node1 < node2)
        {
            *(Z++) = static_cast<t_float>(node1);
            *(Z++) = static_cast<t_float>(node2);
        }
        else
        {
            *(Z++) = static_cast<t_float>(node2);
            *(Z++) = static_cast<t_float>(node1);
        }
        *(Z++) = dist;
        *(Z++) = size;
    }
};

template<const bool sorted>
std::unique_ptr<double[]> &&
generate_SciPy_dendrogram(
    std::unique_ptr<double[]> && Z,
    cluster_result && Z2,
    const t_index N)
{
    auto Z_ = [&Z](const t_index & _r, const t_index & _c) -> t_index
    {
        return (Z[_r * 4 + _c]);
    };
    auto size_ = [&N, &Z_](const t_index & r_) -> t_index
    {
        return (((r_ < N) ? 1 : Z_(r_ - N, 3)));
    };

    // The array "nodes" is a union-find data structure for the cluster
    // identities (only needed for unsorted cluster_result input).
    union_find nodes(sorted ? 0 : N);
    if (!sorted)
    {
        std::stable_sort(Z2[0], Z2[N - 1]);
    }

    linkage_output output(&Z[0]);
    t_index node1, node2;

    for (node const * NN = Z2[0]; NN != Z2[N - 1]; ++NN)
    {
        // Get two data points whose clusters are merged in step i.
        if (sorted)
        {
            node1 = NN->node1;
            node2 = NN->node2;
        }
        else
        {
            // Find the cluster identifiers for these points.
            node1 = nodes.Find(NN->node1);
            node2 = nodes.Find(NN->node2);
            // Merge the nodes in the union-find data structure by making them
            // children of a new node.
            nodes.Union(node1, node2);
        }
        output.append(node1, node2, NN->dist, size_(node1) + size_(node2));
    }

    return std::move(Z);
}

std::unique_ptr<double[]> &&
linkage(
    const t_index N,
    std::unique_ptr<real_t[]> && matrix,
    std::unique_ptr<double[]> && Z)
{
    std::unique_ptr<t_index[]> members(new t_index[N]);
    std::fill(&members[0], &members[N], 1);
    cluster_result Z2(N - 1);
    NN_chain_core<t_index>(N, std::move(matrix), std::move(members), Z2);
    Z = std::move(generate_SciPy_dendrogram<false>(std::move(Z), std::move(Z2), N));

//    std::cerr << "linkage "; std::copy(&Z[0], &Z[4 * (N - 1)], std::ostream_iterator<double>(std::cerr, " ")); std::cerr << std::endl;

    return std::move(Z);
}

enum
{
    CPY_LIN_LEFT  = 0,
    CPY_LIN_RIGHT = 1,
    CPY_LIN_DIST = 2,
    CPY_LIN_CNT = 3,
    CPY_LIS = 4
};

std::pair<std::unique_ptr<double[]>, std::unique_ptr<double[]>>
get_max_dist_for_each_cluster(
    std::unique_ptr<double[]> && Z,
    std::unique_ptr<double[]> && max_dists,
    const std::size_t N)
{
#define CPY_MAX(_x, _y) ((_x > _y) ? (_x) : (_y))
#define CPY_BITS_PER_CHAR (sizeof(unsigned char) * 8)
#define CPY_FLAG_ARRAY_SIZE_BYTES(num_bits) (CPY_CEIL_DIV((num_bits), \
                                                          CPY_BITS_PER_CHAR))
#define CPY_GET_BIT(_xx, i) (((_xx)[(i) / CPY_BITS_PER_CHAR] >> \
                             ((CPY_BITS_PER_CHAR-1) - \
                              ((i) % CPY_BITS_PER_CHAR))) & 0x1)
#define CPY_SET_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] |= \
                              ((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))
#define CPY_CLEAR_BIT(_xx, i) ((_xx)[(i) / CPY_BITS_PER_CHAR] &= \
                              ~((0x1) << ((CPY_BITS_PER_CHAR-1) \
                                         -((i) % CPY_BITS_PER_CHAR))))

#ifndef CPY_CEIL_DIV
#define CPY_CEIL_DIV(x, y) ((((double)x)/(double)y) == \
                            ((double)((x)/(y))) ? ((x)/(y)) : ((x)/(y) + 1))
#endif
#ifdef CPY_DEBUG
#define CPY_DEBUG_MSG(...) fprintf(stderr, __VA_ARGS__)
#else
#define CPY_DEBUG_MSG(...)
#endif

    std::unique_ptr<int[]> curNode(new int [N]);
    std::size_t ndid, lid, rid; // TODO byly int-y
    int k;

    const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(N);
    std::unique_ptr<unsigned char[]> lvisited(new unsigned char[bff]);
    std::unique_ptr<unsigned char[]> rvisited(new unsigned char[bff]);
    const double * Zrow;
    double max_dist;

    k = 0;
    curNode[k] = (N * 2) - 2;
    std::fill(&lvisited[0], &lvisited[bff], 0);
    std::fill(&rvisited[0], &rvisited[bff], 0);

    while (k >= 0)
    {
        ndid = curNode[k];
        Zrow = &Z[((ndid - N) * CPY_LIS)];
        lid = (int) Zrow[CPY_LIN_LEFT];
        rid = (int) Zrow[CPY_LIN_RIGHT];
        if (lid >= N && !CPY_GET_BIT(lvisited, ndid - N))
        {
            CPY_SET_BIT(lvisited, ndid - N);
            curNode[k + 1] = lid;
            k++;
            continue;
        }
        if (rid >= N && !CPY_GET_BIT(rvisited, ndid - N))
        {
            CPY_SET_BIT(rvisited, ndid - N);
            curNode[k + 1] = rid;
            k++;
            continue;
        }
        max_dist = Zrow[CPY_LIN_DIST];
        if (lid >= N)
        {
            max_dist = CPY_MAX(max_dist, max_dists[lid - N]);
        }
        if (rid >= N)
        {
            max_dist = CPY_MAX(max_dist, max_dists[rid - N]);
        }
        max_dists[ndid - N] = max_dist;
        CPY_DEBUG_MSG("i=%d maxdist[i]=%5.5f verif=%5.5f\n", ndid - N, max_dist,
            max_dists[ndid - N]);
        k--;
    }

    return std::make_pair(std::move(Z), std::move(max_dists));
}

std::unique_ptr<std::size_t[]>
form_flat_clusters_from_monotonic_criterion(
    std::unique_ptr<double[]> && Z,
    std::unique_ptr<double[]> && mono_crit,
    std::unique_ptr<std::size_t[]> && T,
    const real_t CUTOFF,
    const std::size_t N)
{
    std::unique_ptr<int[]> curNode(new int [N]);

    const int bff = CPY_FLAG_ARRAY_SIZE_BYTES(N);
    std::unique_ptr<unsigned char[]> lvisited(new unsigned char[bff]);
    std::unique_ptr<unsigned char[]> rvisited(new unsigned char[bff]);

    std::size_t ndid, lid, rid, nc;
    int k;
    int ms;
    double max_crit;
    const double * Zrow;

    /** number of clusters formed so far. */
    nc = 0;
    /** are we in part of a tree below the cutoff? .*/
    ms = -1;
    k = 0;
    curNode[k] = (N * 2) - 2;
    std::fill(&lvisited[0], &lvisited[bff], 0);
    std::fill(&rvisited[0], &rvisited[bff], 0);

    while (k >= 0)
    {
        ndid = curNode[k];
        Zrow = &Z[((ndid - N) * CPY_LIS)];
        lid = (int) Zrow[CPY_LIN_LEFT];
        rid = (int) Zrow[CPY_LIN_RIGHT];
        max_crit = mono_crit[ndid - N];
        CPY_DEBUG_MSG("cutoff: %5.5f maxc: %5.5f nc: %d\n", cutoff, max_crit, nc);
        if (ms == -1 && max_crit <= CUTOFF)
        {
            CPY_DEBUG_MSG("leader: i=%d\n", ndid);
            ms = k;
            nc++;
        }
        if (lid >= N && !CPY_GET_BIT(lvisited, ndid - N))
        {
            CPY_SET_BIT(lvisited, ndid - N);
            curNode[k + 1] = lid;
            k++;
            continue;
        }
        if (rid >= N && !CPY_GET_BIT(rvisited, ndid - N))
        {
            CPY_SET_BIT(rvisited, ndid - N);
            curNode[k + 1] = rid;
            k++;
            continue;
        }
        if (ndid >= N)
        {
            if (lid < N)
            {
                if (ms == -1)
                {
                    nc++;
                    T[lid] = nc;
                }
                else
                {
                    T[lid] = nc;
                }
            }
            if (rid < N)
            {
                if (ms == -1)
                {
                    nc++;
                    T[rid] = nc;
                }
                else
                {
                    T[rid] = nc;
                }
            }
            if (ms == k)
            {
                ms = -1;
            }
        }
        k--;
    }

    return std::move(T);
}

std::unique_ptr<std::size_t[]>
form_flat_clusters_from_dist(
    std::unique_ptr<double[]> && Z,
    std::unique_ptr<std::size_t[]> && T,
    const real_t CUTOFF,
    const std::size_t N)
{
    std::unique_ptr<double[]> max_dists(new double[N]);

    std::pair<std::unique_ptr<double[]>, std::unique_ptr<double[]>>
        result(get_max_dist_for_each_cluster(std::move(Z), std::move(max_dists), N));

    Z = std::move(result.first);
    max_dists = std::move(result.second);
//    std::cerr << "Z/m " << Z[0] << " " << max_dists[0] << std::endl;

    //CPY_DEBUG_MSG("cupid: n=%d cutoff=%5.5f MD[0]=%5.5f MD[n-1]=%5.5f\n", n, cutoff, max_dists[0], max_dists[n-2]);
    T = form_flat_clusters_from_monotonic_criterion(std::move(Z), std::move(max_dists), std::move(T), CUTOFF, N);

//    std::cerr << "flat_clusters "; std::copy(&T[0], &T[N], std::ostream_iterator<int>(std::cerr, " ")); std::cerr << std::endl;

    return std::move(T);
}

std::vector<std::size_t>
fcluster(std::vector<Seq *>::const_iterator cbegin, std::vector<Seq *>::const_iterator cend)
{
    std::cerr << std::endl << "start " << __FUNCTION__ << std::endl;
    std::cerr << "Processing " << std::distance(cbegin, cend) << " sequences" << std::endl;

    const std::size_t NELEM = std::distance(cbegin, cend);
    std::vector<std::size_t> result(std::distance(cbegin, cend));

    std::unique_ptr<real_t[]> dist_matrix(new real_t[NELEM * (NELEM - 1) / 2]);

    Timestamp then(std::chrono::steady_clock::now());
    dist_matrix = std::move(BuildCondensedMatrix(cbegin, cend, std::move(dist_matrix)));
    std::cerr << "It took " << then.now(std::chrono::steady_clock::now()) << std::endl;

    std::unique_ptr<double[]> linkage_matrix(new double[4 * (NELEM - 1)]);
    linkage_matrix = std::move(linkage(NELEM, std::move(dist_matrix), std::move(linkage_matrix)));

    std::unique_ptr<std::size_t[]> flat_cluster(new std::size_t[NELEM]);

    flat_cluster = form_flat_clusters_from_dist(
        std::move(linkage_matrix),
        std::move(flat_cluster),
        CUTOFF,
        NELEM);

    std::cerr << std::endl << "end " << __FUNCTION__ << std::endl;

    std::copy(&flat_cluster[0], &flat_cluster[NELEM], result.data());

    return result;
}

std::size_t
process_pseqs_in_batches(
    std::vector<Seq *> & pseqs,
    const std::size_t BATCH_SIZE,
    const std::size_t BACKSTEP_SIZE
    )
{
    const std::size_t N = pseqs.size();
    const std::size_t nbatch = (N + BATCH_SIZE - 1) / BATCH_SIZE;
    const std::size_t batch_size = N % nbatch ? 1 + N / nbatch : N / nbatch;

    std::cerr << "Doing " << nbatch << " batch(es), " << batch_size << " sequence(s) each." << std::endl;

    std::vector<Seq *>::const_iterator batch_pos = pseqs.cbegin();

    std::size_t max_cluster_no{0};
    std::size_t backstep{0};

    for (std::size_t bi{0}; batch_pos < pseqs.cend(); ++bi)
    {
        std::cerr << "Now doing batch no. " << (bi + 1) << std::endl;
        std::cerr << "Range: " << std::distance(pseqs.cbegin(), batch_pos) << " " << std::distance(pseqs.cbegin(), batch_pos + batch_size + backstep) << std::endl;

        std::vector<std::size_t> bres =
            std::move(fcluster(batch_pos, std::min(batch_pos + batch_size + backstep, pseqs.cend())));

        // backstep
        if (batch_pos + batch_size + backstep < pseqs.cend())
        {
            auto found = std::find(bres.cend() - BACKSTEP_SIZE, bres.cend(), bres.back());
            backstep = std::distance(found, bres.cend());
            std::cerr << "Backstep is " << backstep << std::endl << std::endl;
            backstep = backstep == BACKSTEP_SIZE ? 0 : backstep;
        }
        else
        {
            // for the last batch we process full results
            backstep = 0;
        }

        auto bres_backed_end = bres.cend() - backstep;

        const std::size_t max_b_cluster_no = *(std::max_element(bres.cbegin(), bres_backed_end));
        std::cerr << "Max received cluster no. " << max_b_cluster_no << std::endl;

        for (std::size_t idx{0}; idx < (std::size_t)std::distance(bres.cbegin(), bres_backed_end); ++idx)
        {
            const std::size_t offset = std::distance(pseqs.cbegin(), batch_pos);
            pseqs[offset + idx]->m_clustno = bres[idx] + max_cluster_no;
        }

        std::cerr << "Did " << std::distance(pseqs.cbegin(), batch_pos) << " assignments so far." << std::endl;
        std::cerr << "In this batch, front " << bres.front() << " to backstep " << *(bres_backed_end - 1) << ", back " << bres.back() << std::endl << std::endl;

        max_cluster_no += max_b_cluster_no;

        batch_pos += bres_backed_end - bres.cbegin();
    }

    return max_cluster_no;
}

void
rearrange_clusters(std::vector<Seq *> & pseqs)
{
    typedef std::pair<std::size_t, std::size_t> cluster_type;

    std::vector<cluster_type> clusters;

    std::size_t curr = 0;
    std::size_t next = 0;

    while (curr != pseqs.size())
    {
        while (next != pseqs.size() && pseqs[next]->m_clustno == pseqs[curr]->m_clustno)
        {
            ++next;
        }
        clusters.emplace_back(curr, next);
//        std::cerr << clusters.back().first << " " << clusters.back().second << std::endl;
        curr = next;
    }

    std::vector<cluster_type>::iterator cit1 = clusters.begin();
    std::vector<cluster_type>::iterator cit2 = std::next(cit1);
    std::vector<cluster_type>::iterator cit3 = std::prev(clusters.end());

    while (cit2 != std::prev(clusters.end()))
    {
        std::vector<cluster_type>::iterator cit3 = std::prev(clusters.end());
        do
        {
            if (GetScore(*pseqs[cit1->first], *pseqs[cit3->first]) < CUTOFF * 4)
            {
                cluster_type moved = *cit3;
                clusters.erase(cit3);
                clusters.insert(cit2, moved);
                ++cit2;
                continue;
            }
            --cit3;
        } while (cit3 != cit2);

//        cit1 = cit2;
        ++cit1;
        ++cit2;
    }

//    bool swapped = false;
//    do
//    {
//        swapped = false;
//        while (cit3 != clusters.end())
//        {
//            if (GetScore(*pseqs[cit1->first], *pseqs[cit2->first]) > GetScore(*pseqs[cit1->first], *pseqs[cit3->first]))
//            {
//                std::swap(*cit2, *cit3);
//            }
//            ++cit1;
//            ++cit2;
//            ++cit3;
//        }
//    } while (swapped);

    for (auto cluster : clusters)
    {
//        std::cerr << cluster.first << " " << cluster.second << std::endl;
    }

    std::vector<Seq *> collage;

    for (auto cluster : clusters)
    {
        std::copy(pseqs.cbegin() + cluster.first, pseqs.cbegin() + cluster.second, std::back_inserter(collage));
    }
    pseqs = collage;
}

std::vector<std::string>
fcluster_1(std::vector<Seq> && seqs)
{
    std::cerr << std::endl << "start " << __FUNCTION__ << std::endl;

    const std::size_t N = seqs.size();

    std::vector<Seq *> pseqs(N);
    std::transform(seqs.begin(), seqs.end(), pseqs.begin(),
        [](Seq & val)
        {
        return &val;
        }
    );

    auto seqp_comparator = [](const void * p, const void * q) -> int
    {
        const Seq * lhs = *(const Seq **)p;
        const Seq * rhs = *(const Seq **)q;

        const std::string & lhs_v = lhs->v_gene;
        const std::string & rhs_v = rhs->v_gene;

        const int v_rank = lhs_v.compare(rhs_v);
        if (v_rank == 0)
        {
            const std::string & lhs_j = lhs->j_gene;
            const std::string & rhs_j = rhs->j_gene;

            const int j_rank = lhs_j.compare(rhs_j);
            if (j_rank == 0)
            {
                return lhs->junc.compare(rhs->junc);
            }
            else
            {
                return j_rank;
            }
        }
        else
        {
            return v_rank;
        }
    };

    std::cerr << "Sorting sequences..." << std::endl;
    std::qsort(pseqs.data(), pseqs.size(), sizeof (*pseqs.data()), seqp_comparator);
    done();

    std::vector<std::string> result;

    std::size_t clustno = process_pseqs_in_batches(pseqs, 3000, 600);

    if (clustno > 2)
    {
        rearrange_clusters(pseqs);
        process_pseqs_in_batches(pseqs, 10000, 2000);
    }

//    std::transform(pseqs.cbegin(), pseqs.cend(), std::ostream_iterator<std::size_t>(std::cerr, "\n"),
//        [](const Seq * pseq)
//        {
//            return pseq->m_clustno;
//        }
//    );

    std::transform(seqs.cbegin(), seqs.cend(), std::back_inserter(result),
        [](const Seq & seq)
        {
            return std::to_string(seq.m_clustno);
        }
    );

    std::cerr << std::endl << "end " << __FUNCTION__ << std::endl;

    return result;
}

std::vector<std::string>
fcluster_2(std::vector<Seq> && seqs)
{
    std::cerr << std::endl << "start " << __FUNCTION__ << std::endl;

    const std::size_t N = seqs.size();

    std::vector<Seq *> pseqs(N);
    std::transform(seqs.begin(), seqs.end(), pseqs.begin(),
        [](Seq & val)
        {
        return &val;
        }
    );

    auto seqp_comparator = [](const void * p, const void * q) -> int
    {
        const Seq * lhs = *(const Seq **)p;
        const Seq * rhs = *(const Seq **)q;

        const std::string & lhs_v = lhs->v_gene;
        const std::string & rhs_v = rhs->v_gene;

        const int v_rank = lhs_v.compare(rhs_v);
        if (v_rank == 0)
        {
            const std::string & lhs_j = lhs->j_gene;
            const std::string & rhs_j = rhs->j_gene;

            const int j_rank = lhs_j.compare(rhs_j);
            if (j_rank == 0)
            {
                return lhs->junc.compare(rhs->junc);
            }
            else
            {
                return j_rank;
            }
        }
        else
        {
            return v_rank;
        }
    };

    std::cerr << "Sorting sequences..." << std::endl;
    std::qsort(pseqs.data(), pseqs.size(), sizeof (*pseqs.data()), seqp_comparator);

    done();

    std::vector<std::string> result;
    std::size_t curr_clust_no{1};

    for (auto it = pseqs.cbegin(); it < (pseqs.end() - 1); ++it)
    {
        (*it)->m_clustno = curr_clust_no;
        if (GetScore(**it, **(it + 1)) > CUTOFF)
        {
            ++curr_clust_no;
        }
//        std::cerr << (*it)->m_clustno << " " << GetScore(**it, **(it + 1)) << std::endl;
    }

    std::transform(seqs.cbegin(), seqs.cend(), std::back_inserter(result),
        [](const Seq & seq)
        {
            return std::to_string(seq.m_clustno);
        }
    );

    std::cerr << std::endl << "end " << __FUNCTION__ << std::endl;

    return result;
}

std::vector<std::string>
fcluster(std::vector<std::string> && json)
{
    std::cerr << "Consuming input collection..." << std::endl;
    Timestamp then1(std::chrono::steady_clock::now());
    std::string json_str{vs_concat(std::move(json))};
    std::cerr << "It took " << then1.now(std::chrono::steady_clock::now()) << std::endl;
    done();

    std::cerr << "Parsing JSON..." << std::endl;
    jsonxx::Array json_array;
    Timestamp then2(std::chrono::steady_clock::now());
    json_array.parse(json_str);
    std::cerr << "It took " << then2.now(std::chrono::steady_clock::now()) << std::endl;
    done();

    const std::size_t n = json_array.size();
    if (json_array.empty())
    {
        std::cerr << "Error - No Sequence in the input file." << std::endl;
        exit(0);
    }

    std::cerr << "Building sequences..." << std::endl;
    std::vector<Seq> seqs;
    for (size_t i = 0; i < json_array.size() && i < size_t(n); i++)
    {
        auto o = json_array.get<jsonxx::Object>(i);
        seqs.emplace_back(i, o, "junc_aa");
    }
    json_array.reset();
    done();
    std::cerr << "Got " << seqs.size() << " sequences." << std::endl;

    // Seqs are ready, do the hard work now

    std::vector<std::string> result = fcluster_1(std::move(seqs));
//    std::vector<std::string> result = fcluster_2(std::move(seqs));

    return result;
}

struct ABCSpeedup
{
    std::vector<std::string> cluster(std::vector<std::string> & json)
    {
        return fcluster(std::move(json));
    }
};

#endif /* ABCSPEEDUP_HPP_ */
