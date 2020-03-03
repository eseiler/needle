#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "ibf.h"
#include "minimizer.h"
#include "search.h"

#ifndef DATA_DIR
#  define DATA_DIR @DATA_DIR@
#endif

void initialization_args(arguments & args)
{
    args.compressed = false;
    args.k = 4;
    args.window_size = 4;
    args.seed = 0;
}

void initialization_ibf_args(ibf_arguments & args)
{
    args.expression_levels = {1};
    args.bin_size = {1000};
    args.path_out = DATA_DIR;
}

TEST(minimizer, small_example)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,2,2,1,2,2,2,1,1,1,1};

    auto hash_table = compute_occurrences(seq, 4, 4, 0, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 4, 4, 0, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected,  minimizer_occurences);
}

TEST(minimizer, small_example_k_w_unequal)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,1,2,1,1};

    auto hash_table = compute_occurrences(seq, 4, 8, 0, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 4, 8, 0, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected, minimizer_occurences);
}

TEST(minimizer, small_example_gaps)
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> seq = {"ACGTCGACGTTTAG"_dna4};
    std::vector<uint32_t> expected{2,1,2,1,1};

    auto hash_table = compute_occurrences(seq, 5, 8, 0b10101, 0);

    std::vector<uint32_t> minimizer_occurences;
    for (auto & minHash : compute_minimizer(seq[0], 5, 8, 0b10101, 0))
    {
        minimizer_occurences.push_back(hash_table[minHash]);
    }

    EXPECT_EQ(expected, minimizer_occurences);
}

TEST(ibf, median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, mean)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};
    ibf_args.normalization_method = "mean";

    std::vector<uint32_t> expected{3};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(ibf, random)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};
    ibf_args.normalization_method = "random";
    ibf_args.random = 40;

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_TRUE((3 == medians[0]));
}

TEST(ibf, genom_median)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};
    ibf_args.genome_file = std::string(DATA_DIR) + "/mini_genom.fasta";

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

// Test a genome which has so different minimizers that some minimizers in seq can not be found
TEST(ibf, genom_median_no_match)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example2.fasta"};
    ibf_args.genome_file = std::string(DATA_DIR) + "/mini_genom.fasta";

    std::vector<uint32_t> expected{1};

    std::vector<uint32_t> medians = ibf(args, ibf_args);

    EXPECT_EQ(expected, medians);
}

TEST(ibf, genom_mean)
{
    arguments args{};
    ibf_arguments ibf_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};
    ibf_args.genome_file = std::string(DATA_DIR) + "/mini_genom.fasta";
    ibf_args.normalization_method = "mean";

    std::vector<uint32_t> expected{4};

    std::vector<uint32_t> means = ibf(args, ibf_args);

    EXPECT_EQ(expected, means);
}

TEST(search, small_example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{1};
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_DIR) + "/mini_gen.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_gene_not_found)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{0};
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_DIR) + "/mini_gen2.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, small_example_own_cutoffs)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    initialization_args(args);
    initialization_ibf_args(ibf_args);
    std::vector<uint32_t> expected{0};
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/mini_example.fasta"};
    ibf_args.expression_levels = {0};
    ibf_args.cutoffs = {2};

    ibf(args, ibf_args);

    search_args.search_file = std::string(DATA_DIR) + "/mini_gen3.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 1;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}

TEST(search, example)
{
    arguments args{};
    ibf_arguments ibf_args{};
    search_arguments search_args{};
    std::vector<uint32_t> expected{0,1};
    ibf_args.sequence_files = {std::string(DATA_DIR) + "/exp_01.fasta", std::string(DATA_DIR) + "/exp_02.fasta",
                               std::string(DATA_DIR) + "/exp_11.fasta", std::string(DATA_DIR) + "/exp_12.fasta"};
    ibf_args.samples = {2,2};
    ibf_args.expression_levels = {0.5};
    ibf_args.bin_size = {100000};
    ibf_args.path_out = std::string(DATA_DIR) + "/";
    args.compressed = false;
    ibf(args, ibf_args);

    // ./needle-search DATA_DIR"+"/gene.fasta -i DATA_DIR"+"/ -e 0.5 -c
    search_args.search_file = std::string(DATA_DIR) + "/gene.fasta";
    search_args.path_in = ibf_args.path_out;
    search_args.expression = 0.5;

    std::vector<uint32_t> results{search(args, search_args)};

    EXPECT_EQ(expected, results);
}