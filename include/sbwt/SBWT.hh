#pragma once

//#include <bit>
#include <bitset>
#include <cstdint>
#include <iostream>

#include <vector>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "NodeBOSSInMemoryConstructor.hh"
#include "throwing_streams.hh"
#include "suffix_group_optimization.hh"
#include "kmc_construct.hh"
#include "globals.hh"
#include "Kmer.hh"
#include "batch_query.hh"
#include <map>
#include <optional>

/*

This file contains a class that implements the SBWT index described in the paper:

Alanko, J. N., Puglisi, S. J., & Vuohtoniemi, J. (2022). Succinct k-mer Set 
Representations Using Subset Rank Queries on the Spectral Burrows-Wheeler 
Transform (SBWT). bioRxiv.

And methods to support batched k-mer lookup queries with the SBWT described in the paper:

Batched k-mer lookup on the Spectral Burrows-Wheeler Transform
*/

using namespace std;

namespace sbwt{

const std::string SBWT_VERSION = "v0.1"; // Update this after breaking changes. This is serialized with the index and checked when loading.

// Assumes that a root node always exists
template <typename subset_rank_t>
class SBWT{

private:

    subset_rank_t subset_rank; // The subset rank query implementation
    sdsl::bit_vector suffix_group_starts; // Marks the first column of every suffix group (see paper)
    vector<int64_t> C; // The array of cumulative character counts

    vector<pair<int64_t,int64_t> > kmer_prefix_precalc; // SBWT intervals for all p-mers with p = precalc_k.
    int64_t precalc_k = 0;

    int64_t n_nodes; // Number of nodes (= columns) in the data structure
    int64_t n_kmers; // Number of k-mers indexed in the data structure
    int64_t k; // The k-mer k

    static constexpr char alphabet[4] = {'A', 'C', 'G', 'T'};

    int64_t get_char_idx(char c) const{
        switch(c){
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    }

public:

    struct BuildConfig{
        vector<string> input_files; /**< List of paths to input filenames. */
        int k = 30; /**< Length of the k-mers. */
        bool build_streaming_support = true; /**< Whether we should build the streaming query support. */
        int n_threads = 1; /**< Number of parallel threads in construction. */
        int min_abundance = 1; /**< k-mers occurring fewer than this many times are discarded. */
        int max_abundance = 1e9; /**< k-mers occurring more than this many times are discarded */
        int ram_gigas = 2; /**< RAM budget in gigabytes. Not strictly enforced. */
        int precalc_k = 0; /**< We will precalculate and store the SBWT intervals of all DNA-strings of this length */
        string temp_dir = "."; /**< Path to the directory for the temporary files. */
    };

    /**
     * @brief Construct an empty SBWT.
     * 
     */
    SBWT() : n_nodes(0), n_kmers(0), k(0){}

    /**
     * @brief Construct SBWT from plain matrix SBWT bit vectors.
     * 
     * @param A_bits Row of character A in the plain matrix SBWT.
     * @param C_bits Row of character C in the plain matrix SBWT.
     * @param G_bits Row of character G in the plain matrix SBWT.
     * @param T_bits Row of character T in the plain matrix SBWT.
     * @param streaming_support The streaming support bit vector. Can be empty.
     * @param k Length of the k-mers.
     * @param number_of_kmers Number of k-mers in the data structure.
     */
    SBWT(const sdsl::bit_vector& A_bits, 
         const sdsl::bit_vector& C_bits, 
         const sdsl::bit_vector& G_bits, 
         const sdsl::bit_vector& T_bits, 
         const sdsl::bit_vector& streaming_support, // Streaming support may be empty
         int64_t k, 
         int64_t number_of_kmers,
         int64_t precalc_k);

    /**
     * @brief Construct SBWT using the KMC-based construction algorithm.
     * 
     * @param config construction paramters.
     */
    SBWT(const BuildConfig& config); 

    // Accessors

    /**
     * @brief Get a const reference to the internal subset rank structure.
     */
    const subset_rank_t& get_subset_rank_structure() const {return subset_rank;}

    /**
     * @brief Get a const reference to the internal streaming support bit vector.
     */
    const sdsl::bit_vector& get_streaming_support() const {return suffix_group_starts;}

    /**
     * @brief Compute and return a bit vector that marks which nodes do not correspond to a full k-mer.
     */
    sdsl::bit_vector compute_dummy_node_marks() const;

    /**
     * @brief Get a const reference to the cumulative character count array.
     */
    const vector<int64_t>& get_C_array() const {return C;}

    /**
     * @brief Get a const reference to the k-mer prefix precalc
     */
    const vector<pair<int64_t, int64_t>>& get_precalc() const {return kmer_prefix_precalc;}

    /**
     * @brief Get the precalc k-mer prefix length
     */
    int64_t get_precalc_k() const {return precalc_k;}


    /**
     * @brief Precalculate all SBWT intervals of all strings of length prefix_length. These will be used in search.
     */
    void do_kmer_prefix_precalc(int64_t prefix_length);

    /**
     * @brief Get the number of subsets in the SBWT (= number of columns in the plain matrix representation).
     */
    int64_t number_of_subsets() const {return n_nodes;}

    /**
     * @brief Get the number of k-mers indexed in the data structure.
     */
    int64_t number_of_kmers() const {return n_kmers;}

    /**
     * @brief Get the length of the k-mers.
     */
    int64_t get_k() const {return k;}

    /**
     * @brief Follow an edge in the de Bruijn graph. If called on a dummy node, follows an edge in the dummy node tree.
     * 
     * @param node The node to move from.
     * @param c The character to follow.
     * @return The node id of the node at the end of the edge from `node` labeled with `c`, or -1 if does not exist.
     */
    int64_t forward(int64_t node, char c) const;

    /**
     * @brief Search for a k-mer as an std::string.
     * 
     * @param kmer The k-mer to search for. If the length of this is longer than k, then only the first k-mer is searched.
     * @return The rank of the k-mer in the data structure, or -1 if the k-mer is not in the index.
     * @see streaming_search()
     */
    int64_t search(const string& kmer) const;

    /**
     * @brief Search for a k-mer as C-string.
     * 
     * @param kmer The k-mer to search for.
     * @return The rank of the k-mer in the data structure, or -1 if the k-mer is not in the index.
     * @see streaming_search()
     */
    int64_t search(const char* kmer) const;

    /**
     * @brief Prepare the queues for batch querying.
     * 
     * @param queries The input query vector
     * @param read_queues The queues to be filled with queries
     * @param results The vector of result indices
     */
    void prepare_queues(std::vector<batch_query>& queries, std::array<std::vector<batch_query>, 4>& read_queues, std::vector<int64_t>& results) const;

    /**
     * @brief The main logic of vertical batch querying.
     * 
     * @param read_queue The queue of queries being read/processed
     * @param results The vector of result indices
     * @param write_queues The four write queues (one for each character)
     * @param qk_idx The index of current k-mer character being processed
     */
    void process_queue(std::vector<batch_query>& read_queue, std::vector<int64_t>& results, std::array<std::vector<batch_query>, 4>& write_queues, int64_t qk_idx) const;


    void process_queue_scanning(std::vector<batch_query>& read_queue, std::vector<int64_t>& results, std::array<std::vector<batch_query>, 4>& write_queues, int64_t qk_idx, const sdsl::bit_vector** DNA_bitvectors, vector<pair<int64_t, int64_t>>& s_pointers, vector<pair<int64_t, int64_t>>& e_pointers) const;

    /**
     * @brief Empty the queue into results once processing has completed.
     * 
     * @param read_queue The queue of queries to transform into results
     * @param results The vector of result indices to fill up based on the query intervals in the read_queue
     */
    void empty_queue(std::vector<batch_query>& read_queue, std::vector<int64_t>& results) const;
 
    /**
     * @brief Search for a batch of k-mers as batch_query types.
     * 
     * @param batch The vector of k-mers to search for.
     * @return A vector of the rank of each k-mer in the data structure, or -1 if the k-mer is not in the index.
     */
    std::vector<int64_t> batch_search(std::vector<batch_query>&) const;

    /**
     * @brief Scan a bit vector from "start" to "end' instead of using sdsl rank operations.
     * 
     * @param start The starting index
     * @param end The ending index
     * @param Bit_v The bit vector
     */
    int64_t scanning(int64_t start, const int64_t end, const sdsl::bit_vector& Bit_v) const;

    /**
     * @brief The main logic of vertical batch querying scanning bitvectors instead of using rank.
     * 
     * @param read_queue The queue of queries being read/processed
     * @param results The vector of result indices
     * @param write_queues The four write queues (one for each character)
     * @param qk_idx The index of current k-mer character being processed
     * @param DNA_bitvectors The SBWT matrix bit vectors
     * @param s_pointers The vector with pairs of starting indices of lexicographic rank interval, and rank values, for each char
     * @param e_pointers The vector with pairs of ending indices of lexicographic rank interval, and rank values, for each char
     */
    void process_queues_scanning(std::vector<batch_query>& read_queue, std::vector<int64_t>& results, std::array<std::vector<batch_query>, 4>& write_queues, int64_t qk_idx, const sdsl::bit_vector** DNA_bitvectors, vector<pair<int64_t, int64_t>>& s_pointers, vector<pair<int64_t, int64_t>>& e_pointers) const;

    /**
     * @brief Search for a batch of k-mers as batch_query types without usign rank.
     * 
     * @param queries The vector of k-mers to search for.
     * @param DNA_bitvectors The SBWT binary matrix bitvcetors
     * @return A vector of the rank of each k-mer in the data structure, or -1 if the k-mer is not in the index.
     */
    std::vector<int64_t> batch_search_scanning(std::vector<batch_query>& queries, const sdsl::bit_vector** DNA_bitvectors) const;

    /**
     * @brief Search for a batch of k-mers using the LCS array 
     * 
     * @param batch The vector of reads to search for.
     * @param lcs The LCS array 
     * @return A vector of the rank of each k-mer in the data structure, or -1 if the k-mer is not in the index.
     */
    std::vector<int64_t> batch_streaming_search_lcs(std::vector<StreamingQuery>&, const sdsl::int_vector<>& lcs) const;

    // Helper functions used in streaming search
    tuple<int64_t, int64_t, int64_t, int64_t, int64_t> expand_interval(int64_t start, int64_t end, char c, int64_t match_len, int64_t rank_start, int64_t rank_past_end, const sdsl::int_vector<>& lcs) const;
    tuple<int64_t, int64_t, int64_t, int64_t, int64_t> expand_interval_alternative(int64_t start, int64_t end, char c, int64_t match_len, int64_t rank_start, int64_t rank_past_end, const sdsl::int_vector<>& lcs) const;

    /**
     * @brief Search for all k-mers the query using the lcs array for streaming
     * 
     * @param query The query string
     * @param query_length Length of the query string
     * @param lcs The LCS array 
     * @return A vector of the rank of each k-mer in the data structure, or -1 if the k-mer is not in the index.
     */
    std::vector<int64_t> streaming_search_lcs(const char* query, int64_t query_length, const sdsl::int_vector<>& lcs) const;


    /**
     * @brief Searches for up to the first len characters of the input. For k-mer lookups, it's
     *        better to use `search` because it uses a precalculated lookup table to speed up the search.
     * @param input The input string.
     * @param len The length of the input string.
     * @return {{l,r}, d}, where [l,r] is the colexicographic interval of the longest
     *         prefix of the input that is found in the index, and len is the length of that prefix. 
     * @see search()
     */
    std::pair<std::pair<int64_t, int64_t>, int64_t> partial_search(const char* input, int64_t len) const;

    /**
     * @brief Searches for up to the first len characters of the input. For k-mer lookups, it's
     *        better to use `search` because it uses a precalculated lookup table to speed up the search.
     * @param input The input string.
     * @return {{l,r}, d}, where [l,r] is the colexicographic interval of the longest
     *         prefix of the input that is found in the index, and len is the length of that prefix. 
     * @see search()
     */
    std::pair<std::pair<int64_t, int64_t>, int64_t> partial_search(const std::string& input) const;

    /**
     * @brief Run SBWT search iterations for characters in S starting from interval I.
     * 
     * @param S The string to search for. Can be of any length.
     * @param I The SBWT interval to start from.
     * @return The updated interval, or {-1,-1} if it was not found.
     */
    std::pair<int64_t,int64_t> update_sbwt_interval(const string& S, std::pair<int64_t,int64_t> I) const;

    /**
     * @brief Run SBWT search iterations for characters in S starting from interval I.
     * 
     * @param S The string to search for.
     * @param S_length The length of S.
     * @param I The SBWT interval to start from.
     * @return The updated interval, or {-1,-1} if it was not found.
     */
    std::pair<int64_t,int64_t> update_sbwt_interval(const char* S, int64_t S_length, std::pair<int64_t,int64_t> I) const;

    /**
     * @brief Query all k-mers of the input std::string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. The result will be the same as if search() was called for each k-mer of the input from left to right in order.
     * @see search()
     */
    vector<int64_t> streaming_search(const string& input) const;

    /**
     * @brief Query all k-mers of the input C-string. Requires that the streaming support had been built.
     * 
     * @throws std::runtime_error If the streaming support has not been built.
     * @param input The input string 
     * @param len Length of the input string
     * @return vector<int64_t> The ranks of the k-mers of the input in the data structure, with -1 for those that are not found in the index. The result will be the same as if search() was called for each k-mer of the input from left to right in order.
     * @see search()
     */
    vector<int64_t> streaming_search(const char* input, int64_t len) const;

    /**
     * @brief Whether streaming support is built for the data structure.
     * 
     * @return true If streaming support has been built.
     * @return false If streaming support has not been built.
     */
    bool has_streaming_query_support() const {return suffix_group_starts.size() > 0;}
    
    /**
     * @brief Write the data structure into the given output stream.
     * 
     * @param out The output stream.
     * @return int64_t Number of bytes written.
     * @see load()
     */
    int64_t serialize(ostream& out) const; // Returns the number of bytes written

    /**
     * @brief Write the data structure into the given file.
     * 
     * @param filename The output file.
     * @return int64_t Number of bytes written.
     * @see load()
     */
    int64_t serialize(const string& filename) const; // Returns the number of bytes written

    /**
     * @brief Load the serialized data structure from an input stream.
     * 
     * @param in The input stream.
     * @see serialize()
     */
    void load(istream& in);

    /**
     * @brief Load the serialized data structure from an input file.
     * 
     * @param filename The input file.
     * @see serialize()
     */
    void load(const string& filename);

    /**
     * @brief Reconstruct all k-mers in the data structure.
     * 
     * @return string The reconstructed k-mers concatenated into a single string in colexicographic order, including dummy k-mers.
     *         The i-th k-mer starts at position i*k in the string.
     */
    string reconstruct_all_kmers() const;

    /**
     * @brief Retrieve the k-mer with the given colexicographic rank (including dummy k-mers). Has time complexity O(k log n).
     * 
     * @param colex_rank The colexicographic rank, between 0 and number_of_subsets() - 1.
     * @param buf The output array where the k-mer will be stored. Must have at least k bytes of space.
     */
    void get_kmer(int64_t colex_rank, char* buf) const;

    /**
     * @brief Retrieve the k-mer with the given colexicographic rank (including dummy k-mers), using a subset select support on the SBWT. Has time complexity O(kt), where t is the time for a select query.
     * 
     * @param colex_rank The colexicographic rank, between 0 and number_of_subsets() - 1.
     * @param buf The output array where the k-mer will be stored. Must have at least k bytes of space.
     * @param ss class with a const member function taking a 1-based rank int64_t r and a char c, returning the smallest position p such that subsetrank(p+1,c) >= r, where subsetrank is on the subset sequence of this SBWT.
     * @see SubsetMatrixSelectSupport
     */

    template<typename subset_select_support_t>
    void get_kmer_fast(int64_t colex_rank, char* buf, const subset_select_support_t& ss) const;
};


template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const sdsl::bit_vector& A_bits, const sdsl::bit_vector& C_bits, const sdsl::bit_vector& G_bits, const sdsl::bit_vector& T_bits, const sdsl::bit_vector& streaming_support, int64_t k, int64_t n_kmers, int64_t precalc_k){
    subset_rank = subset_rank_t(A_bits, C_bits, G_bits, T_bits);

    this->n_nodes = A_bits.size();
    this->k = k;
    this->suffix_group_starts = streaming_support;
    this->n_kmers = n_kmers;

    // Get the C-array
    C.clear(); C.resize(4);
    C[0] = 1; // There is one incoming ghost-dollar to the root node
    C[1] = C[0] + subset_rank.rank(n_nodes, 'A');
    C[2] = C[1] + subset_rank.rank(n_nodes, 'C');
    C[3] = C[2] + subset_rank.rank(n_nodes, 'G');

    do_kmer_prefix_precalc(precalc_k);

}

template <typename subset_rank_t>
SBWT<subset_rank_t>::SBWT(const BuildConfig& config){
    string old_temp_dir = get_temp_file_manager().get_dir();
    get_temp_file_manager().set_dir(config.temp_dir);

    NodeBOSSKMCConstructor<SBWT<subset_rank_t>> builder;
    builder.build(config.input_files, *this, config.k, config.n_threads, config.ram_gigas, config.build_streaming_support, config.min_abundance, config.max_abundance, config.precalc_k);

    get_temp_file_manager().set_dir(old_temp_dir); // Return the old temporary directory

    // Precalc will be done in the other constructor which is called by the builder
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::forward(int64_t node, char c) const{
    if(!has_streaming_query_support())
        throw std::runtime_error("Error: Streaming support required for SBWT::forward");

    // Go to start of the suffix group.
    while(!suffix_group_starts[node]) node--; // Guaranteed to terminate because the first node is always marked

    int64_t r1 = subset_rank.rank(node, c);
    int64_t r2 = subset_rank.rank(node+1, c);
    if(r1 == r2) return -1; // No edge found. TODO: could save one rank query if we had direct access to the SBWT sets

    return C[get_char_idx(c)] + r1;
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const string& kmer) const{
    assert(kmer.size() == k);
    return search(kmer.c_str());
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::search(const char* kmer) const{
    pair<int64_t, int64_t> I;
    
    if(precalc_k > 0){ // Precalc is available
        // Find the index of the k-mer prefix in the precalc table (see do_kmer_prefix_precalc).
        uint64_t precalc_idx = 0;
        for(int64_t i = 0; i < precalc_k; i++){
            int64_t char_idx = DNA_to_char_idx(kmer[precalc_k-1-i]);
            if(char_idx == -1) return -1; // non-ACGT character
            precalc_idx = (precalc_idx << 2) | char_idx; // Add the character
        }

        // Continue search from precalculated interval
        I = update_sbwt_interval(kmer + precalc_k, k - precalc_k, kmer_prefix_precalc[precalc_idx]);
    }

    else // No precalc
        I = update_sbwt_interval(kmer, k, {0, n_nodes-1});

    if(I.first != I.second){
        cerr << "Bug: k-mer search did not give a singleton interval: " << I.first << " " << I.second << endl;
        exit(1);
    }
    return I.first;
}


// Prepare the queues based on the first character of each query
template <typename subset_rank_t>
void SBWT<subset_rank_t>::prepare_queues(std::vector<batch_query>& queries, std::array<std::vector<batch_query>, 4>& read_queues, std::vector<int64_t>& results) const {

    // Get the interval for each character in the alphabet
    std::array<std::pair<int64_t, int64_t>, 4> init_intervals;
    for (int64_t char_idx = 0; char_idx < 3; ++char_idx) {
        int64_t l = C[char_idx]; 
        int64_t r = C[char_idx+1] - 1;
        init_intervals[char_idx] = {l, r};
    }
    init_intervals[3]={C[3], n_nodes - 1};

    // Init the queues according to the first character of each query
    size_t query_idx = 0;
    for (auto &query : queries) {
        char c = query.char_at(0);
        int64_t char_idx = get_char_idx(c);
        query.update_interval(init_intervals[char_idx]);
        if (c != -1) read_queues[char_idx].push_back(query);
        else results[query_idx] = -1; // non-acgt character
        query_idx += 1;
    }
}

// Vertical batch search
template <typename subset_rank_t>
void SBWT<subset_rank_t>::process_queue(std::vector<batch_query>& read_queue, std::vector<int64_t>& results, std::array<std::vector<batch_query>, 4>& write_queues, int64_t qk_idx) const {

    // Temporary interval
    pair<int64_t, int64_t> I;
 
    // Rip thru each queue
    for (auto &query : read_queue) {
        char c = query.char_at(qk_idx);
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) {
          results[query.get_idx()] = -1; // Invalid character
          continue;
        }
 
        I = query.get_interval();
        I.first = C[char_idx] + subset_rank.rank(I.first, c);
        I.second = C[char_idx] + subset_rank.rank(I.second+1, c) - 1;
        // Not found
        if(I.first > I.second) { 
            results[query.get_idx()] = -1;
        } 
        // Still alive... Update the interval and add to new queue
        else {
            query.update_interval(I);
            write_queues[char_idx].push_back(query);
        } 
    }
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::empty_queue(std::vector<batch_query>& read_queue, std::vector<int64_t>& results) const {

    // Temporary interval
    pair<int64_t, int64_t> I;
 
    // Rip thru each queue
    for (const auto &query : read_queue) {
        I = query.get_interval();
        // Bug...
        if(I.first != I.second){
            cerr << "While emptying the queue... Bug: k-mer search did not give a singleton interval: " << I.first << " " << I.second << endl;
            exit(1);
        }
        results[query.get_idx()] = I.first;
    }
}

template <typename subset_rank_t>
std::vector<int64_t> SBWT<subset_rank_t>::batch_search(std::vector<batch_query>& queries) const {

    // Storage for results
    std::vector<int64_t> results(queries.size());

    // Track the current character and effective length of k-mers 
    size_t qk_idx = 0;
    
    // 0 = A, 1 = C, 2 = G, 3 = T
    std::array<std::vector<batch_query>, 4> read_queues;
    std::array<std::vector<batch_query>, 4> write_queues;

    // Pay for the worst-case storage cost up-front
    for (size_t i = 0; i < 4; ++i) {
        read_queues[i].reserve(queries.size());
        write_queues[i].reserve(queries.size());
    }

    // Do not apply precalc to vertical batch
    prepare_queues(queries, read_queues, results);
    qk_idx = 1;
   
    // For each successive character up to k
    for (; qk_idx < k; ++qk_idx) {

      // Process the queues in A>C>G>T order
      for (auto &queue : read_queues) {
          process_queue(queue, results, write_queues, qk_idx);
          queue.clear();
      }
      // Read queues and now write queues, vv.
      std::swap(read_queues, write_queues);
    }

    // Processing is complete. Time to collect answers
    for (auto &queue : read_queues) {
        empty_queue(queue, results);
    }
    
    return results;
}



template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::scanning(int64_t start, const int64_t end, const sdsl::bit_vector& Bit_v) const {

    int64_t count_of_1s = 0;

    int64_t len = end-start;
    if (len < 64){
        count_of_1s += __builtin_popcountll(Bit_v.get_int(start, len));
        return count_of_1s;
    }

    // from start to 64-bit boundary
    if (start % 64 != 0) {
        int64_t max_start = ((start / 64) +1) * 64;
            count_of_1s += __builtin_popcountll(Bit_v.get_int(start,max_start - start)); // do not consider the last char now
            start = max_start;
    }

    // 64-bit blocks
    while (start + 64 <= end) {
        int64_t block_idx = start / 64;
        count_of_1s += __builtin_popcountll(Bit_v.data()[block_idx]);
        start += 64;
    }

    // till the end
    if(start < end){
        count_of_1s += __builtin_popcountll(Bit_v.get_int(start,end-start));
    }
    
    return count_of_1s;
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::process_queue_scanning(std::vector<batch_query>& read_queue, std::vector<int64_t>& results, std::array<std::vector<batch_query>, 4>& write_queues, int64_t qk_idx, const sdsl::bit_vector** DNA_bitvectors, vector<pair<int64_t, int64_t>>& s_pointers, vector<pair<int64_t, int64_t>>& e_pointers) const {
    
    // Temporary interval
    pair<int64_t, int64_t> I;

    // Rip thru each queue
    for (auto &query : read_queue) {
        char c = query.char_at(qk_idx);
        int64_t char_idx = get_char_idx(c);
        if(char_idx == -1) {
          results[query.get_idx()] = -1; // Invalid character
          continue;
        }
        I = query.get_interval();
        const sdsl::bit_vector& Bit_v = *(DNA_bitvectors[char_idx]);
        
        if (I.first > s_pointers[char_idx].first){
            s_pointers[char_idx].second += scanning(s_pointers[char_idx].first,I.first,Bit_v); // scan the missing part, update the number of c found
            s_pointers[char_idx].first = I.first; //update the scanned length

        }
        if (e_pointers[char_idx] < s_pointers[char_idx]){ e_pointers[char_idx] = s_pointers[char_idx];}
        if (I.second+1 > e_pointers[char_idx].first){
            e_pointers[char_idx].second += scanning(e_pointers[char_idx].first,I.second+1,Bit_v); // scan the missing part
            e_pointers[char_idx].first = I.second+1;

        } 
        // We only want to know how many chars are there inside the interval
        I.first = C[char_idx] + s_pointers[char_idx].second;
        I.second = C[char_idx] + e_pointers[char_idx].second - 1;
        
        // Not found
        if(I.first > I.second) { 
            results[query.get_idx()] = -1;
        } 
        // Still alive... Update the interval and add to new queue
        else {
            query.update_interval(I);
            write_queues[char_idx].push_back(query);
        } 
    }
}

// Horizontal batch search scanning keeping track of pointers [plain-matirx only]
template <typename subset_rank_t>
std::vector<int64_t> SBWT<subset_rank_t>::batch_search_scanning(std::vector<batch_query>& queries, const sdsl::bit_vector** DNA_bitvectors) const { // bool pointers,

    std::vector<int64_t> results(queries.size());
    // Track the current character and effective length of k-mers 
    size_t qk_idx = 0;
    
    // 0 = A, 1 = C, 2 = G, 3 = T
    std::array<std::vector<batch_query>, 4> read_queues;
    std::array<std::vector<batch_query>, 4> write_queues;

    // Pay now
    for (size_t i = 0; i < 4; ++i) {
        read_queues[i].reserve(queries.size());
        write_queues[i].reserve(queries.size());
    }

    // Do not apply precalc to vertical batch
    prepare_queues(queries, read_queues, results);
    qk_idx = 1;
 
    // For each successive character up to k
    for (; qk_idx < k; ++qk_idx) {

        // {x,y} y(.second) how many c before x excluded
        vector<pair<int64_t, int64_t>> s_pointers = {{0,0},{0,0},{0,0},{0,0}};
        vector<pair<int64_t, int64_t>> e_pointers = {{0,0},{0,0},{0,0},{0,0}};

        // Rip thru each queue
        for (auto &queue : read_queues) {
            process_queue_scanning(queue, results, write_queues, qk_idx, DNA_bitvectors, s_pointers, e_pointers);
            queue.clear();
            s_pointers = e_pointers;
        }
        std::swap(read_queues, write_queues);

    }

    // Collect answers
    for (auto &queue : read_queues) {
        empty_queue(queue, results);
    }

    return results;
}

template <typename subset_rank_t>
tuple<int64_t, int64_t, int64_t, int64_t, int64_t> SBWT<subset_rank_t>::expand_interval(int64_t start, int64_t end, char c, int64_t match_len, int64_t rank_start, int64_t rank_past_end, const sdsl::int_vector<>& lcs) const {

    int64_t c_idx = from_ACGT_to_0123_lookup_table[c];
    check_true(c_idx >= 0, "Interval expansion called on invalid character");

    int64_t new_match_len = match_len;
    while(rank_start == rank_past_end && new_match_len > 0){
        new_match_len--;

        while(start > 0 && lcs[start] >= new_match_len) {
            start--;
            rank_start -= subset_rank.contains_by_char_idx(start, c_idx);
        }
        while(end+1 < this->number_of_subsets() && lcs[end+1] >= new_match_len) {
            end++;
            rank_past_end += subset_rank.contains_by_char_idx(end, c_idx);
        }
    }
    
    return {start, end, new_match_len, rank_start, rank_past_end};
}

template <typename subset_rank_t>
tuple<int64_t, int64_t, int64_t, int64_t, int64_t> SBWT<subset_rank_t>::expand_interval_alternative(int64_t start, int64_t end, char c, int64_t match_len, int64_t rank_start, int64_t rank_past_end, const sdsl::int_vector<>& lcs) const {

    int64_t new_match_len = match_len;
    while(rank_start == rank_past_end && new_match_len > 0){
        new_match_len--;

        while(start > 0 && lcs[start] >= new_match_len) {
            start--;
        }
        while(end+1 < this->number_of_subsets() && lcs[end+1] >= new_match_len) {
            end++;
        }

        rank_start = subset_rank.rank(start, c);
        rank_past_end = subset_rank.rank(end+1, c);
    }
    
    return {start, end, new_match_len, rank_start, rank_past_end};
}

template <typename subset_rank_t>
std::vector<int64_t> SBWT<subset_rank_t>::batch_streaming_search_lcs(std::vector<StreamingQuery>& queries, const sdsl::int_vector<>& lcs) const {

    check_true(lcs.size() == this->number_of_subsets(), "LCS must have the same length as the SBWT");

    // We assume all queries are of equal length. Verify that.
    check_true(queries.size() > 0, "Query vector is empty");
    int64_t n = queries[0].read.size();
    for(auto& q : queries) check_true(q.read.size() == n, "All queries need to have the same length");
    check_true(n >= get_k(), "Queries need to be longer or equal to length k");

    int64_t total_number_of_kmers = queries.size() * (n - get_k() + 1);

    // Setup output writing pointers for queries
    vector<int64_t> all_answers(total_number_of_kmers);
    int64_t cumul_kmer_count = 0;
    for(StreamingQuery& query : queries){
        query.result_write_ptr = all_answers.data() + cumul_kmer_count;
        cumul_kmer_count += n - get_k() + 1;
    }

    // Five queues: N, A, C, G, T. The N-queue contains everything that is not from ACGT.
    std::array<std::vector<StreamingQuery>, 5> read_queues; 
    std::array<std::vector<StreamingQuery>, 5> write_queues;

    std::swap(read_queues[0], queries); // Put the queries in the first read queue
    // From now on we must not accidentally access the original query vector

    auto reset_interval = [this](StreamingQuery& query, int64_t character_count){
        query.colex_start = 0;
        query.colex_end = this->number_of_subsets() - 1;
        query.match_len = 0;
    };

    for(int64_t round = 0; round < n; round++) {
        for(int64_t queue_idx = 0; queue_idx < 5; queue_idx++) {
            for(StreamingQuery& query : read_queues[queue_idx]) {
                char c = query.read.get(query.pos); 
                int64_t char_idx = get_char_idx(c);

                if(char_idx >= 0) {
                    int64_t rank_start = subset_rank.rank(query.colex_start, c);
                    int64_t rank_past_end = subset_rank.rank(query.colex_end+1, c);

                    // Expand interval
                    tie(query.colex_start, query.colex_end, query.match_len, rank_start, rank_past_end) = expand_interval(query.colex_start, query.colex_end, c, query.match_len, rank_start, rank_past_end, lcs);

                    // Update interval
                    query.colex_start = C[char_idx] + rank_start; 
                    query.colex_end = C[char_idx] + rank_past_end - 1; 

                    if(query.colex_end >= query.colex_start) {
                        // Successful extension
                        query.match_len = min((uint8_t) (query.match_len + 1), (uint8_t) k);
                    } else {
                        // Unsuccessful extension -> back to interval of the empty string
                        reset_interval(query, 0);
                    }
                } else { // Invalid character -> expand all the way to the full interval
                    reset_interval(query, 0);
                }
                
                // Report result
                if(query.pos >= k-1){
                    // We have now read a full k-mer
                    if(query.match_len == this->get_k()) {
                        check_true(query.colex_start == query.colex_end, "BUG: interval not singleton"); // Must be singleton interval, otherwise bug
                        query.write_answer(query.colex_start);
                    } else {
                        query.write_answer(-1);
                    }
                }

                query.pos++;

                // Determine output queue.
                int64_t write_queue_index = char_idx + 1; // -1 0 1 2 3 -> 0 1 2 3 4

                // Write with std::move to avoid heap copies
                write_queues[write_queue_index].push_back(std::move(query));

            }

        }

        // Make write queues into read queues for the next round
        for(int64_t queue_idx = 0; queue_idx < 5; queue_idx++) {
            read_queues[queue_idx].clear();
            std::swap(read_queues[queue_idx], write_queues[queue_idx]);
        }
    }

    return all_answers;

}

template <typename subset_rank_t>
std::vector<int64_t> SBWT<subset_rank_t>::streaming_search_lcs(const char* query, int64_t query_length, const sdsl::int_vector<>& lcs) const {
    if(query_length < this->get_k()) return {};

    int64_t n_kmers = query_length - this->get_k() + 1;
    vector<int64_t> answers;
    answers.reserve(n_kmers);

    // Algorithm state
    int64_t colex_start = 0;
    int64_t colex_end = number_of_subsets();
    int64_t match_len = 0;
    int64_t rank_start, rank_past_end;

    auto reset_interval = [&colex_start, &colex_end, &match_len, &rank_start, &rank_past_end, this](){
        colex_start = 0;
        colex_end = this->number_of_subsets() - 1;
        match_len = 0;
        rank_start = 0;
        rank_past_end = 0;
    };

    for(int64_t i = 0; i < query_length; i++){
        char c = query[i];
        int64_t char_idx = from_ACGT_to_0123_lookup_table[c];

        // Expand interval
        if(char_idx >= 0) { // c is from the alphabet ACGT
            rank_start = subset_rank.rank(colex_start, c);
            rank_past_end = subset_rank.rank(colex_end+1, c);
            tie(colex_start, colex_end, match_len, rank_start, rank_past_end) = expand_interval(colex_start, colex_end, c, match_len, rank_start, rank_past_end, lcs);

            // Update interval
            colex_start = C[char_idx] + rank_start; 
            colex_end = C[char_idx] + rank_past_end - 1; 

            if(colex_end >= colex_start) {
                // Successful extension
                match_len = min(match_len + 1, get_k());
            } else {
                // Unsuccessful extension -> back to interval of the empty string
                reset_interval();
            }
        } else { // Invalid character -> expand all the way to the full interval
            reset_interval();
        }

        // Store result
        if(i >= k-1){
            // We have now read a full k-mer
            if(match_len == get_k()) {
                check_true(colex_start == colex_end, "BUG: interval not singleton"); // Must be singleton interval, otherwise bug
                answers.push_back(colex_start);
            } else {
                answers.push_back(-1);
            }
        }
    }

    return answers;
}


template<typename subset_rank_t>
std::pair<int64_t,int64_t> SBWT<subset_rank_t>::update_sbwt_interval(const string& S, pair<int64_t,int64_t> I) const{
    return update_sbwt_interval(S.c_str(), S.size(), I);
}

template<typename subset_rank_t>
std::pair<int64_t,int64_t> SBWT<subset_rank_t>::update_sbwt_interval(const char* S, int64_t S_length, pair<int64_t,int64_t> I) const{
    if(I.first == -1) return I;
    for(int64_t i = 0; i < S_length; i++){
        char c = toupper(S[i]);
        int64_t char_idx = get_char_idx(S[i]);
        if(char_idx == -1) return {-1,-1}; // Invalid character

        I.first = C[char_idx] + subset_rank.rank(I.first, c);
        I.second = C[char_idx] + subset_rank.rank(I.second+1, c) - 1;

        if(I.first > I.second) return {-1,-1}; // Not found
    }

    return {I.first, I.second};
}


// Utility function: Serialization for a std::vector
// Returns number of bytes written
template<typename T>
int64_t serialize_std_vector(const vector<T>& v, ostream& os){
    // Write C-array
    int64_t n_bytes = sizeof(T) * v.size();
    os.write((char*)&n_bytes, sizeof(n_bytes));
    os.write((char*)v.data(), n_bytes);
    return sizeof(n_bytes) + n_bytes;
}

template<typename T>
vector<T> load_std_vector(istream& is){
    int64_t n_bytes = 0;
    is.read((char*)&n_bytes, sizeof(n_bytes));
    assert(n_bytes % sizeof(T) == 0);
    vector<T> v(n_bytes / sizeof(T));
    is.read((char*)(v.data()), n_bytes);
    return v;
}


template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::serialize(ostream& os) const{
    int64_t written = 0;

    written += serialize_string(SBWT_VERSION, os);

    written += subset_rank.serialize(os);
    written += suffix_group_starts.serialize(os);

    written += serialize_std_vector(C, os);

    // Write precalc
    written += serialize_std_vector(kmer_prefix_precalc, os);
    os.write((char*)&precalc_k, sizeof(precalc_k));
    written += sizeof(precalc_k);

    // Write number of nodes
    os.write((char*)&n_nodes, sizeof(n_nodes));
    written += sizeof(n_nodes);

    // Write number of k-mers
    os.write((char*)&n_kmers, sizeof(n_kmers));
    written += sizeof(n_kmers);

    // Write k
    os.write((char*)&k, sizeof(k));
    written += sizeof(k);

    return written;
}

template <typename subset_rank_t>
int64_t SBWT<subset_rank_t>::serialize(const string& filename) const{
    throwing_ofstream out(filename, ios::binary);
    return serialize(out.stream);
}


template <typename subset_rank_t>
void SBWT<subset_rank_t>::load(istream& is){
    string version = load_string(is);
    if(version != SBWT_VERSION){
        throw std::runtime_error("Error: Corrupt index file, or the index was constructed with an incompatible version of SBWT.");
    }

    subset_rank.load(is);
    suffix_group_starts.load(is);
    C = load_std_vector<int64_t>(is);
    kmer_prefix_precalc = load_std_vector<pair<int64_t, int64_t>>(is);
    is.read((char*)&precalc_k, sizeof(precalc_k));
    is.read((char*)&n_nodes, sizeof(n_nodes));
    is.read((char*)&n_kmers, sizeof(n_kmers));
    is.read((char*)&k, sizeof(k));

}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::load(const string& filename){
    throwing_ifstream in(filename, ios::binary);
    load(in.stream);
}


template <typename subset_rank_t>
std::pair<std::pair<int64_t, int64_t>, int64_t> SBWT<subset_rank_t>::partial_search(const char* input, int64_t len) const{
    int64_t l = 0;
    int64_t r = n_nodes-1;
    for(int64_t i = 0; i < len; i++){
        char c = toupper(input[i]);
        int64_t l_new, r_new;
        std::tie(l_new, r_new) = update_sbwt_interval(&c, 1, {l,r});
        if(l_new == -1) return {{l,r},i};
        l = l_new; r = r_new;
    }
    return {{l,r},len}; // All found
}

template <typename subset_rank_t>
std::pair<std::pair<int64_t, int64_t>, int64_t> SBWT<subset_rank_t>::partial_search(const string& input) const{
    return partial_search(input.c_str(), input.size());
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const char* input, int64_t len) const{
    if(suffix_group_starts.size() == 0)
        throw std::runtime_error("Error: streaming search support not built");
    
    vector<int64_t> ans;
    if(len < k) return ans;

    // Search the first k-mer
    const char* first_kmer_start = input;
    ans.push_back(search(first_kmer_start)); 

    for(int64_t i = 1; i < len - k + 1; i++){
        if(ans.back() == -1){
            // Need to search from scratch
            ans.push_back(search(first_kmer_start + i));
        } else{
            // Got to the start of the suffix group and do one search iteration
            int64_t column = ans.back();
            while(suffix_group_starts[column] == 0) column--; // can not go negative because the first column is always marked

            char c = toupper(input[i+k-1]);
            int64_t char_idx = get_char_idx(c);
        
            if(char_idx == -1) ans.push_back(-1); // Not found
            else{
                int64_t node_left = column;
                int64_t node_right = column;
                node_left = C[char_idx] + subset_rank.rank(node_left, c);
                node_right = C[char_idx] + subset_rank.rank(node_right+1, c) - 1;
                if(node_left == node_right) ans.push_back(node_left);
                else ans.push_back(-1);
                // Todo: could save one subset rank query if we have fast access to the SBWT columns
            }
        }
    }
    return ans;
}

template <typename subset_rank_t>
vector<int64_t> SBWT<subset_rank_t>::streaming_search(const string& input) const{
    return streaming_search(input.c_str(), input.size());
}

template <typename subset_rank_t>
sdsl::bit_vector SBWT<subset_rank_t>::compute_dummy_node_marks() const{
    int64_t count = 0;
    vector<pair<int64_t,int64_t>> dfs_stack; // pairs (node, depth)
    dfs_stack.push_back({0, 0}); // Root node
    // dfs to depth k-1
    // the dummy part is a tree so no visited-list is required

    string ACGT = "ACGT";
    sdsl::bit_vector marks(n_nodes, 0);
    int64_t v,d; // node,depth
    while(!dfs_stack.empty()){
        tie(v,d) = dfs_stack.back();
        dfs_stack.pop_back();
        if(d < k){
            count++;
            marks[v] = 1;
        }
        if(d < k-1){ // Push children
            for(char c : ACGT){
                int64_t u = forward(v, c);
                if(u != -1) dfs_stack.push_back({u,d+1});
            }
        }
    }
    return marks;
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::do_kmer_prefix_precalc(int64_t prefix_length){
    if(prefix_length == 0) return;
    if(prefix_length > 20){
        throw std::runtime_error("Error: Can't precalc longer than 20-mers (would take over 4^20 = 2^40 bytes");
    }

    if(prefix_length > k)
        throw std::runtime_error("Error: Precalc length is longer than k (" + to_string(prefix_length) + " > " + to_string(k) + ")");
    
    uint64_t n_kmers_to_precalc = 1 << (2*prefix_length); // Four to the power prefix_length

    // Initialize member variables
    kmer_prefix_precalc.resize(n_kmers_to_precalc);
    this->precalc_k = prefix_length;

    uint64_t data = 0; // K-mer packed 2 bits per nucleotide
    string prefix(prefix_length, '\0');

    while(n_kmers_to_precalc--){
        for(int64_t i = 0; i < prefix_length; i++){
            char c = char_idx_to_DNA((data >> (2*i)) & 0x3); // Decode the i-th character
            prefix[i] = c;
        }

        kmer_prefix_precalc[data] = update_sbwt_interval(prefix, {0, n_nodes-1});
        data++;
    }

}

template <typename subset_rank_t>
string SBWT<subset_rank_t>::reconstruct_all_kmers() const {

    int64_t n_nodes = this->number_of_subsets(); 
    vector<int64_t> C_array(4);

    vector<char> last; // last[i] = incoming character to node i
    last.push_back('$');

    C_array[0] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(subset_rank.contains(i,'A')) last.push_back('A');

    C_array[1] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(subset_rank.contains(i,'C')) last.push_back('C');

    C_array[2] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(subset_rank.contains(i,'G')) last.push_back('G');
    
    C_array[3] = last.size();
    for(int64_t i = 0; i < n_nodes; i++) if(subset_rank.contains(i,'T')) last.push_back('T');

    if(last.size() != n_nodes){
        cerr << "BUG " << last.size() << " " << n_nodes << endl;
        exit(1);
    }

    string kmers_concat(n_nodes * k, '\0');

    for(int64_t round = 0; round < k; round++){
        //cerr << "round " << round << "/" << k-1 << endl;
        for(int64_t i = 0; i < n_nodes; i++){
            int64_t pos = k-1-round;
            kmers_concat[i*k + pos] = last[i];
        }

        // Propagate the labels one step forward in the graph
        vector<char> propagated(n_nodes, '$');
        int64_t A_ptr = C_array[0];
        int64_t C_ptr = C_array[1];
        int64_t G_ptr = C_array[2];
        int64_t T_ptr = C_array[3];
        for(int64_t i = 0; i < n_nodes; i++){
            if(subset_rank.contains(i,'A')) propagated[A_ptr++] = last[i];
            if(subset_rank.contains(i,'C')) propagated[C_ptr++] = last[i];
            if(subset_rank.contains(i,'G')) propagated[G_ptr++] = last[i];
            if(subset_rank.contains(i,'T')) propagated[T_ptr++] = last[i];
        }
        last = propagated;
    }

    return kmers_concat;
}

template <typename subset_rank_t>
void SBWT<subset_rank_t>::get_kmer(int64_t colex_rank, char* buf) const {
    for(int64_t i = 0; i < this->k; i++){
        if(colex_rank == 0){
            buf[k-1-i] = '$';
        } else{ 
            int64_t char_idx = 0;
            while(char_idx+1 < 4 && colex_rank >= C[char_idx+1]) char_idx++;
            char c = SBWT<subset_rank_t>::alphabet[char_idx];
            buf[k-1-i] = c;

            // Step backward

            int64_t char_rel_rank = colex_rank - C[char_idx];
            // Find the index p of the SBWT subset that contains the occurrence of c with rank char_rel_rank
            int64_t p = 0;
            int64_t step = this->number_of_subsets();
            while(step > 0){
                while(p + step <= this->number_of_subsets() && this->subset_rank.rank(p + step, c) <= char_rel_rank)
                    p += step;
                step /= 2;
            }
            colex_rank = p;
        }
    }
}

template <typename subset_rank_t>
template <typename subset_select_support_t>
void SBWT<subset_rank_t>::get_kmer_fast(int64_t colex_rank, char* buf, const subset_select_support_t& ss) const{
    for(int64_t i = 0; i < this->k; i++){
        if(colex_rank == 0){
            buf[k-1-i] = '$';
        } else{ 
            int64_t char_idx = 0;
            while(char_idx+1 < 4 && colex_rank >= C[char_idx+1]) char_idx++;
            char c = SBWT<subset_rank_t>::alphabet[char_idx];
            buf[k-1-i] = c;

            // Step backward
            // Find the index p of the SBWT subset that contains the occurrence of c with rank char_rel_rank
            int64_t char_rel_rank = colex_rank - C[char_idx];
            char_rel_rank++; // 1-based rank
            colex_rank = ss.select(char_rel_rank, c);
        }
    }

}

} // namespace sbwt
