#include <gtest/gtest.h>
#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "ibf.h"
#include "shared.h"
#include "estimate.h"

#ifndef DATA_INPUT_DIR
#  define DATA_INPUT_DIR @DATA_INPUT_DIR@
#endif

using seqan3::operator""_shape;

void initialization_args(estimate_ibf_arguments & args)
{
    args.compressed = true;
    args.k = 4;
    args.shape = seqan3::ungapped{args.k};
    args.w_size = seqan3::window_size{4};
    args.s = seqan3::seed{0};
}

TEST(estimate, small_example)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.expression_thresholds = {1, 2, 4};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    estimate_args.path_in = ibf_args.path_out;
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_4");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_uncompressed)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.compressed = false;
    ibf_args.expression_thresholds = {1, 2, 4};
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    estimate_args.path_in = ibf_args.path_out;
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_4");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_gene_not_found)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.expression_thresholds = {2, 4};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen2.fasta";
    estimate_args.path_in = ibf_args.path_out;
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};
    std::vector<uint8_t> cutoffs{};

    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen2\t0\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_2");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_4");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_different_expressions_per_level)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test0_";
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Estimate_Test0_mini_example.minimiser"};
    ibf_args.expression_thresholds= {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"Test0_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"Test0_expression.out");
    std::string line;
    std::string expected{"gen1\t3\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test0_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Estimate_Test0_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test0_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"Estimate_Test0_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Estimate_Test0_IBF_Data");
    std::filesystem::remove(tmp_dir/"Test0_expression.out");
    std::filesystem::remove(tmp_dir/"Estimate_Test0_mini_example.minimiser");
}

TEST(estimate, small_example_different_expressions_per_level_normalization_1)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    estimate_args.normalization_method = 1;
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.number_expression_thresholds = 2;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Estimate_Test_mini_example.minimiser"};
    ibf_args.expression_thresholds= {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    std::string expected{"gen1\t1\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, small_example_different_expressions_per_level_normalization_1_uncompressed)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    estimate_args.normalization_method = 1;
    initialization_args(ibf_args);
    ibf_args.path_out = tmp_dir/"Estimate_Test2_";
    ibf_args.number_expression_thresholds = 2;
    ibf_args.compressed = false;
    std::vector<double> fpr = {0.05};
    std::vector<uint8_t> cutoffs{};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "mini_example.fasta"};

    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Estimate_Test2_mini_example.minimiser"};
    ibf_args.expression_thresholds= {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "mini_gen.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"Test2_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"Test2_expression.out");
    std::string line;
    std::string expected{"gen1\t1\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
            EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test2_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Estimate_Test2_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test2_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"Estimate_Test2_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Estimate_Test2_IBF_Data");
    std::filesystem::remove(tmp_dir/"Test2_expression.out");
    std::filesystem::remove(tmp_dir/"Estimate_Test2_mini_example.minimiser");
}

TEST(estimate, example)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    std::vector<double> fpr = {0.05};
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                                                         std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    minimiser_args.samples = {2, 2};
    ibf_args.expression_thresholds = {4, 32};
    ibf_args.path_out = tmp_dir/"Estimate_Test_Single_";
    ibf_args.compressed = false;
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);

    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "gene.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"Single_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"Single_expression.out");
    std::string line;
    std::string expected{"GeneA\t9\t32\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
             EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_Single_IBF_4");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Single_IBF_32");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Single_IBF_Data");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Single_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Single_expression.out");
}

TEST(estimate, example_multiple_threads)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                                                         std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    minimiser_args.samples = {2,2};
    ibf_args.expression_thresholds = {4, 32};
    std::vector<double> fpr = {0.05};
    ibf_args.path_out = tmp_dir/"Estimate_Test_Multiple_";
    ibf_args.compressed = false;
    std::vector<uint8_t> cutoffs{};
    ibf(sequence_files, ibf_args, minimiser_args, fpr, cutoffs);
    ibf_args.threads = 2;

    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "gene4.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"Multiple_expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"Multiple_expression.out");
    std::string line;
    std::string expected{"GeneA\t9\t32\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
             EXPECT_EQ(expected,line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_Multiple_IBF_32");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Multiple_IBF_4");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Multiple_IBF_Data");
    std::filesystem::remove(tmp_dir/"Estimate_Test_Multiple_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Multiple_expression.out");
}

TEST(estimate, example_different_expressions_per_level)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                                                         std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.samples = {2,2};
    ibf_args.number_expression_thresholds = 4;
    std::vector<double> fpr = {0.05};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.compressed = false;
    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Estimate_Test_exp_01.minimiser", tmp_dir/"Estimate_Test_exp_11.minimiser"};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "gene.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    // Count would expect 6 and 34
    std::string expected{"GeneA\t7\t26\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file, line) )
        {
             EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_2");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"expression.out");
}

TEST(estimate, example_different_expressions_per_level_multiple_threads)
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory
    estimate_ibf_arguments ibf_args{};
    minimiser_arguments minimiser_args{};
    estimate_arguments estimate_args{};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    std::vector<std::filesystem::path> sequence_files = {std::string(DATA_INPUT_DIR) + "exp_01.fasta", std::string(DATA_INPUT_DIR) + "exp_02.fasta",
                                                         std::string(DATA_INPUT_DIR) + "exp_11.fasta", std::string(DATA_INPUT_DIR) + "exp_12.fasta"};
    std::vector<uint8_t> cutoffs = {0, 0};
    minimiser_args.samples = {2,2};
    ibf_args.number_expression_thresholds = 4;
    std::vector<double> fpr = {0.05};
    ibf_args.path_out = tmp_dir/"Estimate_Test_";
    ibf_args.compressed = false;
    minimiser(sequence_files, ibf_args, minimiser_args, cutoffs);
    std::vector<std::filesystem::path> minimiser_files{tmp_dir/"Estimate_Test_exp_01.minimiser", tmp_dir/"Estimate_Test_exp_11.minimiser"};
    ibf_args.expression_thresholds= {};
    ibf(minimiser_files, ibf_args, fpr);

    ibf_args.threads = 2;
    ibf_args.expression_thresholds = {0, 1, 2};
    estimate_args.search_file = std::string(DATA_INPUT_DIR) + "gene4.fasta";
    estimate_args.path_in = ibf_args.path_out;
    ibf_args.path_out = tmp_dir/"expression.out";
    call_estimate(ibf_args, estimate_args);

    std::ifstream output_file(tmp_dir/"expression.out");
    std::string line;
    // Count would expect 6 and 34
    std::string expected{"GeneA\t7\t26\t"};
    if (output_file.is_open())
    {
        while ( std::getline (output_file,line) )
        {
             EXPECT_EQ(expected, line);
        }
        output_file.close();
    }
    std::filesystem::remove(tmp_dir/"Estimate_Test_exp_01.minimiser");
    std::filesystem::remove(tmp_dir/"Estimate_Test_exp_11.minimiser");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_0");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_1");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_2");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Level_3");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Levels.levels");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_FPRs.fprs");
    std::filesystem::remove(tmp_dir/"Estimate_Test_IBF_Data");
    std::filesystem::remove(tmp_dir/"expression.out");
}
