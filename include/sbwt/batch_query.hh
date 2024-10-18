#pragma once

#include <vector>
#include <set>
#include <utility>
#include <iostream>
#include "Kmer.hh"

namespace sbwt{

using namespace std;

typedef Kmer<MAX_KMER_LENGTH> kmer_t;

class batch_query {

  private:
    kmer_t kmer;
    pair<int64_t, int64_t> interval;
    int64_t payload; // could be something other than i64
    int64_t idx;

  public:
    batch_query() {}
    batch_query(const string s, int64_t n_nodes, int64_t pl, int64_t i) : kmer(s), interval(0, n_nodes), payload(pl), idx(i) {}
    batch_query(const string s, int64_t l, int64_t r, int64_t pl, int64_t i) : kmer(s), interval(l, r), payload(pl), idx(i) {}


    bool operator<(const batch_query &other) {
        return interval.first < other.interval.first;
    }

    void update_interval(const int64_t l, const int64_t r) {
        interval.first = l;
        interval.second = r;
    }

    void update_interval(const pair<int64_t, int64_t> p) {
        interval = p;
    }

    void update_payload(int64_t p) {
        payload = p;
    }

    const kmer_t get_kmer() const {
        return kmer;
    }

    char char_at(int64_t idx) const {
        return kmer.get(idx);
    }

    int64_t get_payload() const {
        return payload;
    }

    std::pair<int64_t, int64_t> get_interval() const {
        return interval;
    }

    int64_t get_idx() const {
        return idx;
    }

};

class StreamingQuery{

    private:

    void flush_result_cache(){
        memcpy(result_write_ptr, result_cache, elements_in_result_cache * sizeof(int64_t));
        result_write_ptr += elements_in_result_cache;
        elements_in_result_cache = 0;
    }

    public:

    static const uint8_t RESULT_CACHE_CAPACITY = 8;

    Kmer<256> read;

    // Results are written to a local result cache before flushed to the "global" result vector at result_write_ptr.
    int64_t result_cache[RESULT_CACHE_CAPACITY];
    uint8_t elements_in_result_cache;
    uint8_t match_len; // Length of the current match, up to k. Declared here for word-alignment purposes
    uint32_t pos; // Position of the next character to be processed in read. Declared here for word-alignment purposes.
    int64_t* result_write_ptr; // Need to set this before calling write_answer !!
    int64_t colex_start, colex_end; // Current interval. End is inclusive.

    StreamingQuery(const char* read, int64_t read_length, int64_t n_nodes_in_sbwt) 
        : read(read, read_length), 
          elements_in_result_cache(0),
          match_len(0),
          pos(0),
          result_write_ptr(nullptr), // Needs to be set later
          colex_start(0),
          colex_end(n_nodes_in_sbwt-1) {}

    // result_write_ptr must be set before calling this!!!
    void write_answer(int64_t answer){
        check_true(this->result_write_ptr, "Result write pointer not set"); // This branch should be never taken so this is probably not too expensive
        if(elements_in_result_cache == RESULT_CACHE_CAPACITY) flush_result_cache();
        result_cache[elements_in_result_cache++] = answer;
        if(pos == read.size()-1) flush_result_cache(); // Done
    }

};

}