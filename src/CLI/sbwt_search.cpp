#include <string>
#include <cstring>
#include "cxxopts.hpp"
#include "globals.hh"
#include "SBWT.hh"
#include "SubsetWT.hh"
#include "stdlib_printing.hh"
#include "SeqIO/SeqIO.hh"
#include "SubsetMatrixRank.hh"
#include "SeqIO/buffered_streams.hh"
#include "variants.hh"
#include "commands.hh"
#include "batch_query.hh"
#include <filesystem>
#include <cstdio>

using namespace std;

using namespace sbwt;

// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector(const vector<int64_t>& v, writer_t& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        int64_t i = 0;
        if(x == -1){
            buffer[0] = '1';
            buffer[1] = '-';
            i = 2;
        } else{
            while(x > 0){
                buffer[i++] = '0' + (x % 10);
                x /= 10;
            }
        }
        std::reverse(buffer, buffer + i);
        buffer[i] = ' ';
        out.write(buffer, i+1);
    }
    out.write(&newline, 1);
}

// Assumes values of v are -1 or larger
template <typename writer_t>
inline void print_vector_with_newlines(const vector<int64_t>& v, writer_t& out){
    // Fast manual integer-to-string conversion
    char buffer[32];
    char newline = '\n';
    for(int64_t x : v){
        int64_t i = 0;
        if(x == -1){
            buffer[0] = '1';
            buffer[1] = '-';
            i = 2;
        } else{
            while(x > 0){
                buffer[i++] = '0' + (x % 10);
                x /= 10;
            }
        }
        std::reverse(buffer, buffer + i);
        buffer[i] = ' ';
        out.write(buffer, i+1);
        out.write(&newline, 1);
    }
}


template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t io_micros = 0;
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    while(true){
        int64_t i0 = cur_time_micros();
        int64_t len = reader.get_next_read_to_buffer();
        io_micros += cur_time_micros() - i0;
        if(len == 0) break;

        int64_t t0 = cur_time_micros();
        vector<int64_t> out_buffer = sbwt.streaming_search(reader.read_buf, len);
        total_micros += cur_time_micros() - t0;

        number_of_queries += out_buffer.size();

        // Write out
        print_vector(out_buffer, writer);
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros+io_micros) / number_of_queries) + " (including I/O etc)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)(total_micros+io_micros)) + " (including I/O etc)", LogLevel::MAJOR);

    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_streaming_lcs(reader_t& reader, writer_t& writer, const sbwt_t& sbwt, const sdsl::int_vector<>& lcs){

    int64_t io_micros = 0;
    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    while(true){
        int64_t i0 = cur_time_micros();
        int64_t len = reader.get_next_read_to_buffer();
        io_micros += cur_time_micros() - i0;
        if(len == 0) break;

        int64_t t0 = cur_time_micros();
        vector<int64_t> out_buffer = sbwt.streaming_search_lcs(reader.read_buf, len, lcs);
        total_micros += cur_time_micros() - t0;

        number_of_queries += out_buffer.size();

        // Write out
        print_vector(out_buffer, writer);
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros+io_micros) / number_of_queries) + " (including I/O etc)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)(total_micros+io_micros)) + " (including I/O etc)", LogLevel::MAJOR);

    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_batch(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros_read = 0;
    int64_t total_micros_query = 0;
    int64_t number_of_queries = 0;
    std::vector<batch_query> queries;

    int64_t t0 = cur_time_micros();
    
    // Step 1: Load the batch
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        queries.emplace_back(std::string(reader.read_buf, len), 0, sbwt.number_of_subsets(), number_of_queries);
        number_of_queries++;
    }

    total_micros_read += cur_time_micros() - t0;
    write_log("read: " + to_string(number_of_queries) + " queries in " + to_string(total_micros_read) + " microsecs", LogLevel::MAJOR);
    
    // Step 2: Run the queries
    t0 = cur_time_micros();
    std::vector<int64_t> ans = sbwt.batch_search(queries);
    total_micros_query += cur_time_micros() - t0;
    
    print_vector_with_newlines(ans, writer);

    write_log("us/query: " + to_string((double)total_micros_query / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros_query+total_micros_read) / number_of_queries) + " (including I/O)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)total_micros_query+total_micros_read) + " (including I/O etc)", LogLevel::MAJOR);
    
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_batch_scanning(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;
    const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};

    int64_t total_micros_read = 0;
    int64_t total_micros_query = 0;
    int64_t number_of_queries = 0;
    std::vector<batch_query> queries;

    int64_t t0 = cur_time_micros();
    
    // Step 1: Load the batch
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        queries.emplace_back(std::string(reader.read_buf, len), 0, sbwt.number_of_subsets(), number_of_queries);
        number_of_queries++;
    }

    total_micros_read += cur_time_micros() - t0;
    write_log("read: " + to_string(number_of_queries) + " queries in " + to_string(total_micros_read) + " microsecs", LogLevel::MAJOR);
    
    // Step 2: Run the queries
    t0 = cur_time_micros();
    std::vector<int64_t> ans = sbwt.batch_search_scanning(queries, DNA_bitvectors);

    total_micros_query += cur_time_micros() - t0;
    
    print_vector_with_newlines(ans, writer);

    write_log("us/query: " + to_string((double)total_micros_query / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros_query+total_micros_read) / number_of_queries) + " (including I/O)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)total_micros_query+total_micros_read) + " (including I/O etc)", LogLevel::MAJOR);
    
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_batch_scanning_rank(reader_t& reader, writer_t& writer, const sbwt_t& sbwt, const int64_t r){

    const sdsl::bit_vector& A_bits = sbwt.get_subset_rank_structure().A_bits;
    const sdsl::bit_vector& C_bits = sbwt.get_subset_rank_structure().C_bits;
    const sdsl::bit_vector& G_bits = sbwt.get_subset_rank_structure().G_bits;
    const sdsl::bit_vector& T_bits = sbwt.get_subset_rank_structure().T_bits;
    const sdsl::bit_vector* DNA_bitvectors[4] = {&A_bits, &C_bits, &G_bits, &T_bits};

    int64_t total_micros_read = 0;
    int64_t total_micros_query = 0;
    int64_t number_of_queries = 0;
    std::vector<batch_query> queries;

    int64_t t0 = cur_time_micros();
    
    // Step 1: Load the batch
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        // XXX this needs fixing with a read_id
        queries.emplace_back(std::string(reader.read_buf, len), 0, sbwt.number_of_subsets(), number_of_queries);
        number_of_queries++;
    }

    total_micros_read += cur_time_micros() - t0;
    write_log("read: " + to_string(number_of_queries) + " queries in " + to_string(total_micros_read) + " microsecs", LogLevel::MAJOR);
    
    // Step 2: Run the queries
    t0 = cur_time_micros();

    std::vector<int64_t> ans = sbwt.batch_search_scanning_rank(queries, DNA_bitvectors,r);


    total_micros_query += cur_time_micros() - t0;
    
    print_vector_with_newlines(ans, writer);

    write_log("us/query: " + to_string((double)total_micros_query / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros_query+total_micros_read) / number_of_queries) + " (including I/O)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)total_micros_query+total_micros_read) + " (including I/O etc)", LogLevel::MAJOR);
    
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_batch_streaming_lcs(reader_t& reader, writer_t& writer, const sbwt_t& sbwt, const sdsl::int_vector<>& lcs){

    int64_t total_micros_read = 0;
    int64_t total_micros_query = 0;
    int64_t number_of_queries = 0;
    std::vector<StreamingQuery> queries;

    int64_t t0 = cur_time_micros();

    // Get a C array that also includes the sum of counts of all characters in the end
    vector<int64_t> C_array = sbwt.get_C_array();
    C_array.push_back(C_array.back() + sbwt.get_subset_rank_structure().rank(sbwt.number_of_subsets(), 'T'));
    
    // Step 1: Load the reads
    int64_t query_idx = 0;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        if(len < sbwt.get_k()) continue;

        char c = reader.read_buf[0]; // Ok since len > 0 here
        char c_idx = from_ACGT_to_0123_lookup_table[c];
        int64_t c_count = (c_idx == -1 ? 0 : C_array[c_idx+1] - C_array[c_idx]);

        queries.emplace_back(reader.read_buf, len, sbwt.number_of_subsets());

        number_of_queries += len - sbwt.get_k() + 1; // Is ok since len >= k here
    }

    total_micros_read += cur_time_micros() - t0;
    write_log("read: " + to_string(number_of_queries) + " queries in " + to_string(total_micros_read) + " microsecs", LogLevel::MAJOR);
    
    // Step 2: Run the queries
    t0 = cur_time_micros();

    std::vector<int64_t> ans = sbwt.batch_streaming_search_lcs(queries, lcs);

    total_micros_query += cur_time_micros() - t0;
    
    print_vector_with_newlines(ans, writer);

    write_log("us/query: " + to_string((double)total_micros_query / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros_query+total_micros_read) / number_of_queries) + " (including I/O)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)total_micros_query+total_micros_read) + " (including I/O etc)", LogLevel::MAJOR);
    
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_batch_pointers(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros_read = 0;
    int64_t total_micros_query = 0;
    int64_t number_of_queries = 0;
    std::vector<batch_query> queries;

    int64_t t0 = cur_time_micros();
    
    // Step 1: Load the batch
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        // XXX this needs fixing with a read_id
        queries.emplace_back(std::string(reader.read_buf, len), 0, sbwt.number_of_subsets(), number_of_queries);
        number_of_queries++;
    }

    total_micros_read += cur_time_micros() - t0;
    write_log("read: " + to_string(number_of_queries) + " queries in " + to_string(total_micros_read) + " microsecs", LogLevel::MAJOR);
    
    // Step 2: Run the queries
    t0 = cur_time_micros();

    std::vector<int64_t> ans = sbwt.batch_search_pointers(queries);


    total_micros_query += cur_time_micros() - t0;
    
    print_vector_with_newlines(ans, writer);

    write_log("us/query: " + to_string((double)total_micros_query / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    write_log("us/query: " + to_string((double)(total_micros_query+total_micros_read) / number_of_queries) + " (including I/O)", LogLevel::MAJOR);
    write_log("us total: " + to_string((double)total_micros_query+total_micros_read) + " (including I/O etc)", LogLevel::MAJOR);
    
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_queries_not_streaming(reader_t& reader, writer_t& writer, const sbwt_t& sbwt){

    int64_t total_micros = 0;
    int64_t number_of_queries = 0;
    int64_t k = sbwt.get_k();
    vector<int64_t> out_buffer;
    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;

        for(int64_t i = 0; i < len - k + 1; i++){
            int64_t t0 = cur_time_micros();
            int64_t ans = sbwt.search(reader.read_buf + i);
            total_micros += cur_time_micros() - t0;
            number_of_queries++;
            out_buffer.push_back(ans);
        }

        print_vector(out_buffer, writer);
        out_buffer.clear();
    }
    write_log("us/query: " + to_string((double)total_micros / number_of_queries) + " (excluding I/O etc)", LogLevel::MAJOR);
    return number_of_queries;
}

template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file(const string& infile, const string& outfile, const sbwt_t& sbwt, std::optional<sdsl::int_vector<>>& lcs, bool batch_processing){
    reader_t reader(infile);
    writer_t writer(outfile);

    if(lcs.has_value()){
        if(batch_processing) {
            write_log("Running BATCH queries with LCS STREAMING from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
            return run_queries_batch_streaming_lcs<sbwt_t, reader_t, writer_t>(reader, writer, sbwt, lcs.value());
        } else {
            write_log("Running naive queries with LCS STREAMING from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
            return run_queries_streaming_lcs<sbwt_t, reader_t, writer_t>(reader, writer, sbwt, lcs.value());
        }
    } 

    if(batch_processing) {
        write_log("Running batch queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_batch<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
    else if(sbwt.has_streaming_query_support()){
        write_log("Running non-batched bit vector streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Overload [plain matrix]
template<typename sbwt_t, typename reader_t, typename writer_t>
int64_t run_file(const string& infile, const string& outfile, const sbwt_t& sbwt, std::optional<sdsl::int_vector<>>& lcs, bool batch_processing, bool pointers){
    reader_t reader(infile);
    writer_t writer(outfile);

    if(lcs.has_value()){
        if(batch_processing) {
            write_log("Running BATCH queries with LCS STREAMING from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
            return run_queries_batch_streaming_lcs<sbwt_t, reader_t, writer_t>(reader, writer, sbwt, lcs.value());
        } else {
            write_log("Running naive queries with LCS STREAMING from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
            return run_queries_streaming_lcs<sbwt_t, reader_t, writer_t>(reader, writer, sbwt, lcs.value());
        }
    } 

    if(batch_processing) {
        if(pointers){
            write_log("Running BATCHED queries with POINTERS (scanning) from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
            return run_queries_batch_scanning<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
        }
        else { 
          write_log("Running batched queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
          return run_queries_batch<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
        }
    }
    else if(sbwt.has_streaming_query_support()){
        write_log("Running streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
    else{
        write_log("Running non-streaming queries from input file " + infile + " to output file " + outfile , LogLevel::MAJOR);
        return run_queries_not_streaming<sbwt_t, reader_t, writer_t>(reader, writer, sbwt);
    }
}

// Returns number of queries executed
template<typename sbwt_t>
int64_t run_queries(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, std::optional<sdsl::int_vector<>>& lcs, bool gzip_output, bool batch_processing){

    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>> in_gzip;
    typedef seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef seq_io::Buffered_ofstream<seq_io::zstr::ofstream> out_gzip;
    typedef seq_io::Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = seq_io::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing);
        }
    }
    return n_queries_run;

}

// Returns number of queries executed
// Overload [plain-matrix]
template<typename sbwt_t>
int64_t run_queries(const vector<string>& infiles, const vector<string>& outfiles, const sbwt_t& sbwt, std::optional<sdsl::int_vector<>>& lcs, bool gzip_output, bool batch_processing, bool pointers){
    if(infiles.size() != outfiles.size()){
        string count1 = to_string(infiles.size());
        string count2 = to_string(outfiles.size());
        throw std::runtime_error("Number of input and output files does not match (" + count1 + " vs " + count2 + ")");
    }

    typedef seq_io::Reader<seq_io::Buffered_ifstream<seq_io::zstr::ifstream>> in_gzip;
    typedef seq_io::Reader<seq_io::Buffered_ifstream<std::ifstream>> in_no_gzip;

    typedef seq_io::Buffered_ofstream<seq_io::zstr::ofstream> out_gzip;
    typedef seq_io::Buffered_ofstream<std::ofstream> out_no_gzip;

    int64_t n_queries_run = 0;
    for(int64_t i = 0; i < infiles.size(); i++){
        bool gzip_input = seq_io::figure_out_file_format(infiles[i]).gzipped;
        if(gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing, pointers);
        }
        if(gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing, pointers);
        }
        if(!gzip_input && gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing, pointers);
        }
        if(!gzip_input && !gzip_output){
            n_queries_run += run_file<sbwt_t, in_no_gzip, out_no_gzip>(infiles[i], outfiles[i], sbwt, lcs, batch_processing, pointers);
        }
    }
    return n_queries_run;

}


int search_main(int argc, char** argv){

    int64_t micros_start = cur_time_micros();

    set_log_level(LogLevel::MINOR);

    cxxopts::Options options(argv[0], "Query all k-mers of all input reads.");

    options.add_options()
        ("o,out-file", "Output filename.", cxxopts::value<string>())
        ("i,index-file", "Index input file.", cxxopts::value<string>())
        ("l,lcs-file", "LCS array input file for batched streaming search.", cxxopts::value<string>())
        ("q,query-file", "The query in FASTA or FASTQ format, possibly gzipped. Multi-line FASTQ is not supported. If the file extension is .txt, this is interpreted as a list of query files, one per line. In this case, --out-file is also interpreted as a list of output files in the same manner, one line for each input file.", cxxopts::value<string>())
        ("z,gzip-output", "Writes output in gzipped form. This can shrink the output files by an order of magnitude.", cxxopts::value<bool>()->default_value("false"))
        ("b,batch", "Apply batch processing instead of one-by-one.", cxxopts::value<bool>()->default_value("false"))
        ("p,pointers", "If batch processing is applied, use pointers to scan bitvectors insted of using rank.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    int64_t old_argc = argc; // Must store this because the parser modifies it
    auto opts = options.parse(argc, argv);

    if (old_argc == 1 || opts.count("help")){
        std::cerr << options.help() << std::endl;
        exit(1);
    }

    string indexfile = opts["index-file"].as<string>();
    check_readable(indexfile);

    // Interpret input file
    string queryfile = opts["query-file"].as<string>();
    vector<string> input_files;
    bool multi_file = queryfile.size() >= 4 && queryfile.substr(queryfile.size() - 4) == ".txt";
    if(multi_file){
        input_files = readlines(queryfile);
    } else{
        input_files = {queryfile};
    }
    for(string file : input_files) check_readable(file);

    // Interpret mode
    bool batch_processing = opts["batch"].as<bool>();
    bool pointers = opts["pointers"].as<bool>();

    // Interpret output file
    string outfile = opts["out-file"].as<string>();
    bool gzip_output = opts["gzip-output"].as<bool>();
    vector<string> output_files;
    if(multi_file){
        output_files = readlines(outfile);
    } else{
        output_files = {outfile};
    }
    for(string file : output_files) check_writable(file);

    vector<string> variants = get_available_variants();

    throwing_ifstream in(indexfile, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(std::find(variants.begin(), variants.end(), variant) == variants.end()){
        cerr << "Error loading index from file: unrecognized variant specified in the file" << endl;
        return 1;
    }

    // Load LCS file if given
    std::optional<sdsl::int_vector<>> lcs;
    try { 
        string lcs_filename = opts["lcs-file"].as<string>(); // Throws if not present
        lcs = sdsl::int_vector<>();
        sdsl::load_from_file(*lcs, lcs_filename);
        write_log("Loaded LCS array of length " + to_string(lcs->size()), LogLevel::MAJOR);
    } catch(cxxopts::option_has_no_value_exception& e){
        write_log("No LCS array given", LogLevel::MAJOR);
        // LCS filename not given. That's ok.
    }


    write_log("Loading the index variant " + variant, LogLevel::MAJOR);
    int64_t number_of_queries = 0;
    
    // We only support Plain-Matrix or EF-Split indexes for batch processing
    if (variant == "plain-matrix"){
        plain_matrix_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, lcs, gzip_output, batch_processing, pointers);
    }
    else if (variant == "mef-split"){
        mef_split_sbwt_t sbwt;
        sbwt.load(in.stream);
        number_of_queries += run_queries(input_files, output_files, sbwt, lcs, gzip_output, batch_processing);
    } else {
        cerr << "UNSUPPORTED VARIANT: " << variant << endl;
        exit(1);
    }

    int64_t total_micros = cur_time_micros() - micros_start;
    write_log("us/query end-to-end: " + to_string((double)total_micros / number_of_queries), LogLevel::MAJOR);

    return 0;

}
