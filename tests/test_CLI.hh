#pragma once

#include "setup_tests.hh"
#include "buffered_streams.hh"
#include "globals.hh"
#include "commands.hh"
#include "variants.hh"
#include <gtest/gtest.h>

string read_gzipped_file(string filename){
    zstr::ifstream in(filename);
    string S;
    char c;
    while(in.get(c)) S += c;
    return S;
}

TEST(CLI, end_to_end_build_and_query){
    vector<string> seqs1 = {"ACTAGTGTAGCTACAAA","ATGTGCTGATGCTAGCATTTTTTT"};
    vector<string> seqs2 = {"GTGTACTAGTGTGTAGTCGAT"};
    string seqfile1 = get_temp_file_manager().create_filename("",".fna.gz");
    string seqfile2 = get_temp_file_manager().create_filename("",".fna.gz");

    {
        SeqIO::Writer<zstr::ofstream> writer(seqfile1);
        for(string seq : seqs1) writer.write_sequence(seq.c_str(), seq.size());
    } // End of scope flushes stream (would not flush properly with flush() because of how zstr works)

    {
        SeqIO::Writer<zstr::ofstream> writer(seqfile2);
        for(string seq : seqs2) writer.write_sequence(seq.c_str(), seq.size());
    } // End of scope flushes stream (would not flush properly with flush() because of how zstr works)
    
    string seqfile_list_file = get_temp_file_manager().create_filename("",".txt");
    write_to_file(seqfile1 +"\n" + seqfile2 + "\n", seqfile_list_file);

    // Construct the index
    string indexfile = get_temp_file_manager().create_filename("",".sbwt");
    plain_matrix_sbwt_t sbwt;
    string tempdir = get_temp_file_manager().get_dir();
    vector<string> build_args = {"build","-i",seqfile_list_file,"-o",indexfile,"-k","6","--add-reverse-complements","--temp-dir",tempdir};
    Argv build_argv(build_args);
    build_main(build_argv.size, build_argv.array);

    // Queries

    vector<string> queries = {"GGAGAACTAGTGTAGCTACAAAGAGAG", "AGTGTGTAGCAAAATGTGCTGATGCTAGCAAAAAAAA", "CTCTACACACTTC"}; 

    string q1 = get_temp_file_manager().create_filename("",".fq");
    string q2 = get_temp_file_manager().create_filename("",".fna");
    string q3 = get_temp_file_manager().create_filename("",".fq.gz");
    string q4 = get_temp_file_manager().create_filename("",".fna.gz");

    {
        // Artifical scope to flush the streams in the scope at the end of the scope.
        // It has to be done this way because the flush method of zstr::ofstream does
        // not actually flush anything.
        
        SeqIO::Writer<std::ofstream> w1(q1);
        SeqIO::Writer<std::ofstream> w2(q2);
        SeqIO::Writer<zstr::ofstream> w3(q3);
        SeqIO::Writer<zstr::ofstream> w4(q4);

        for(string S : queries){
            w1.write_sequence(S.c_str(), S.size());
            w2.write_sequence(S.c_str(), S.size());
            w3.write_sequence(S.c_str(), S.size());
            w4.write_sequence(S.c_str(), S.size());
        }
    }

    string o1 = get_temp_file_manager().create_filename("",".txt");
    string o2 = get_temp_file_manager().create_filename("",".txt");
    string o3 = get_temp_file_manager().create_filename("",".txt");
    string o4 = get_temp_file_manager().create_filename("",".txt");

    string input_file_list = get_temp_file_manager().create_filename("",".txt");
    write_to_file(q1 + "\n" + q2 + "\n" + q3 + "\n" + q4 + "\n", input_file_list);

    string output_file_list = get_temp_file_manager().create_filename("",".txt");
    write_to_file(o1 + "\n" + o2 + "\n" + o3 + "\n" + o4 + "\n", output_file_list);

    vector<string> args = {"search", "-o", output_file_list, "-i", indexfile, "-q", input_file_list}; 
    Argv ARGS(args);

    search_main(ARGS.size, ARGS.array);

    string correct_answer = "-1 -1 -1 -1 -1 14 43 72 21 61 80 56 69 18 51 39 64 -1 -1 -1 -1 -1 \n22 62 83 61 80 -1 -1 -1 -1 -1 6 -1 -1 27 82 59 78 53 44 76 47 26 77 52 41 -1 -1 -1 -1 4 4 4 \n-1 -1 40 65 10 30 -1 -1 \n";
    string answer_file = get_temp_file_manager().create_filename("",".txt");
    write_to_file(correct_answer, answer_file);

    ASSERT_TRUE(files_are_equal(answer_file, o1));
    ASSERT_TRUE(files_are_equal(o1, o2));
    ASSERT_TRUE(files_are_equal(o2, o3));
    ASSERT_TRUE(files_are_equal(o3, o4));

    string output_file_list_gz = get_temp_file_manager().create_filename("",".txt");
    write_to_file(o1 + ".gz\n" + o2 + ".gz\n" + o3 + ".gz\n" + o4 + ".gz\n", output_file_list_gz);

    vector<string> args_gz = {"search", "-o", output_file_list_gz, "-i", indexfile, "-q", input_file_list, "--gzip-output"}; 
    Argv ARGS_gz(args_gz);

    search_main(ARGS_gz.size, ARGS_gz.array);

    ASSERT_EQ(read_gzipped_file(o1 + ".gz"), correct_answer);
    ASSERT_EQ(read_gzipped_file(o2 + ".gz"), correct_answer);
    ASSERT_EQ(read_gzipped_file(o3 + ".gz"), correct_answer);
    ASSERT_EQ(read_gzipped_file(o4 + ".gz"), correct_answer);


}

