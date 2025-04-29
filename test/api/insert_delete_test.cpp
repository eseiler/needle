// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/test/expect_range_eq.hpp>

#include "../app_test.hpp"
#include "ibf.hpp"
#include "shared.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct delete_test : public app_test
{
    void initialization_args(configuration & args)
    {
        args.compressed = true;
        args.k = 4;
        args.shape = seqan3::ungapped{args.k};
        args.w_size = seqan3::window_size{4};
        args.s = seqan3::seed{0};
    }
};

struct insert_test : public delete_test
{
    // Reads the level file ibf creates
    template <typename float_or_int>
    void read_levels(std::vector<std::vector<float_or_int>> & expressions, std::filesystem::path filename)
    {
        ASSERT_TRUE(std::filesystem::exists(filename)) << filename;
        std::ifstream fin{filename};
        auto stream_view = seqan3::detail::istreambuf(fin);
        auto stream_it = std::ranges::begin(stream_view);
        size_t j{0};
        std::vector<float_or_int> empty_vector{};

        std::string buffer{};

        // Read line = expression levels
        do
        {
            if (j == expressions.size())
                expressions.push_back(empty_vector);
            std::ranges::copy(stream_view | seqan3::detail::take_until_or_throw(seqan3::is_char<' '>),
                              std::back_inserter(buffer));
            if constexpr (std::same_as<uint16_t, float_or_int>)
                expressions[j].push_back((uint16_t)std::stoi(buffer));
            else
                expressions[j].push_back((double)std::stod(buffer));
            buffer.clear();
            if (*stream_it != '/')
                ++stream_it;

            if (*stream_it == '\n')
            {
                ++stream_it;
                j++;
            }
        }
        while (*stream_it != '/');
        ++stream_it;

        fin.close();
    }
};

TEST_F(delete_test, no_given_thresholds)
{
    configuration config{};
    initialization_args(config);
    config.fpr = {0.05};
    config.cutoffs = {0, 0};
    config.compressed = false;
    config.number_expression_thresholds = 2;
    config.experiment_names = false;
    config.path_out = "IBF_delete_Exp_";
    config.sequence_files = {data("mini_example.fasta"), data("mini_example.fasta")};
    ibf(config);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{};
    load_ibf(ibf, "IBF_delete_Exp_IBF_Level_0");
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_0{seqan3::bin_count{2u},
                                                                              seqan3::bin_size{ibf.bin_size()},
                                                                              seqan3::hash_function_count{1u}};
    load_ibf(ibf, "IBF_delete_Exp_IBF_Level_1");
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_1{seqan3::bin_count{2u},
                                                                              seqan3::bin_size{ibf.bin_size()},
                                                                              seqan3::hash_function_count{1u}};

    config.path_in = config.path_out;
    config.delete_files = {0, 1};
    delete_bin(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_delete{};

    load_ibf(ibf_delete, "IBF_delete_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf_0 == ibf_delete));

    load_ibf(ibf_delete, "IBF_delete_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf_1 == ibf_delete));
}

TEST_F(insert_test, ibf)
{
    configuration config{};
    initialization_args(config);
    config.compressed = false;
    config.path_out = "IBF_True_Exp_";
    config.expression_thresholds = {1, 2};
    config.experiment_names = false;
    config.sequence_files = {data("mini_example.fasta"), data("mini_example.fasta")};
    config.fpr = {0.05};
    config.cutoffs = {0, 0};

    ibf(config);

    std::vector<uint8_t> cutoffs_insert{0};
    config.compressed = false;
    config.path_out = "IBF_True_Exp_";
    config.expression_thresholds = {1, 2};
    config.experiment_names = false;
    config.path_out = "IBF_Insert_Exp_";
    config.path_in = config.path_out;
    config.sequence_files = {data("mini_example.fasta")};
    ibf(config);

    insert(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_no_given_thresholds)
{
    configuration config{};
    initialization_args(config);
    config.compressed = false;
    config.path_out = "IBF_True_Exp_";
    config.number_expression_thresholds = 2;
    config.experiment_names = false;
    config.sequence_files = {data("mini_example.fasta"), data("mini_example.fasta")};
    config.fpr = {0.05};

    std::vector<uint16_t> expected{};
    config.cutoffs = {0, 0};

    std::vector<uint16_t> medians = ibf(config);

    std::vector<uint8_t> cutoffs_insert{0};
    config.compressed = false;
    config.number_expression_thresholds = 2;
    config.experiment_names = false;
    config.path_out = "IBF_Insert_Exp_";
    config.sequence_files = {data("mini_example.fasta")};
    std::vector<uint16_t> medians_insert = ibf(config);
    config.path_in = config.path_out;
    insert(config);
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_0");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_delete)
{
    configuration config{};
    initialization_args(config);
    config.compressed = false;
    config.path_out = "IBF_True_Exp_";
    config.expression_thresholds = {1, 2};
    config.experiment_names = false;
    config.sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta"), data("mini_example.fasta")};
    config.fpr = {0.05, 0.05};
    config.cutoffs = {0, 0, 0};

    ibf(config);

    std::vector<uint8_t> cutoffs_insert{0};
    config.compressed = false;
    config.path_out = "IBF_Insert_Exp_";
    config.expression_thresholds = {1, 2};
    config.experiment_names = false;
    config.sequence_files = {data("mini_example2.fasta")};
    std::vector<std::filesystem::path> sequence_files_test = {data("mini_example.fasta"),
                                                              data("mini_example2.fasta"),
                                                              data("mini_example.fasta")};

    ibf(config);
    config.path_in = config.path_out;
    config.delete_files = {1};
    delete_bin(config);
    insert(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_2");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibf_delete_no_given_threshold)
{
    configuration config{};
    initialization_args(config);
    config.compressed = false;
    config.path_out = "IBF_True_Exp_";
    config.number_expression_thresholds = 2;
    config.experiment_names = false;
    config.sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta"), data("mini_example.fasta")};
    config.fpr = {0.05, 0.05};
    config.cutoffs = {0, 0, 0};

    ibf(config);

    std::vector<uint8_t> cutoffs_insert{0};
    config.compressed = false;
    config.path_out = "IBF_Insert_Exp_";
    config.number_expression_thresholds = 2;
    config.experiment_names = false;
    config.sequence_files = {data("mini_example2.fasta")};
    config.sequence_files = {data("mini_example.fasta"), data("mini_example2.fasta"), data("mini_example.fasta")};

    ibf(config);
    config.path_in = config.path_out;
    config.delete_files = {1};
    delete_bin(config);
    insert(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_0");

    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBF_True_Exp_IBF_Level_1");
    load_ibf(ibf_insert, "IBF_Insert_Exp_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBF_True_Exp_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBF_Insert_Exp_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBF_True_Exp_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBF_Insert_Exp_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin)
{
    configuration config{};
    initialization_args(config);
    config.expression_thresholds = {1, 2};
    config.fpr = {0.05, 0.05};
    config.path_out = "IBFMIN_Test_Given_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser"), data("mini_example.minimiser")};
    ibf_min(config);

    config.expression_thresholds = {1, 2};
    config.path_out = "IBFMIN_Insert_Given_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser")};
    ibf(config);
    insert(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin_delete)
{
    configuration config{};
    initialization_args(config);
    config.expression_thresholds = {1, 2};
    config.fpr = {0.05, 0.05};
    config.path_out = "IBFMIN_Test_Given_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser"),
                             data("mini_example.minimiser"),
                             data("mini_example.minimiser")};
    ibf_min(config);

    config.expression_thresholds = {1, 2};
    config.path_out = "IBFMIN_Insert_Given_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser")};
    ibf_min(config);
    config.delete_files = {1};
    delete_bin(config);
    insert_min(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_1");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_2");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_2");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, ibfmin_no_given_thresholds)
{
    configuration config{};
    initialization_args(config);
    config.number_expression_thresholds = 2;
    config.fpr = {0.05, 0.05};
    config.path_out = "IBFMIN_Test_Given_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser"), data("mini_example.minimiser")};

    ibf_min(config);

    config.number_expression_thresholds = 2;
    config.path_out = "IBFMIN_Insert_Given_";
    config.compressed = false;
    config.fpr = {0.05};
    config.sequence_files = {data("mini_example.minimiser")};
    ibf_min(config);
    config.path_in = "IBFMIN_Insert_Given_";
    insert_min(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_Level_0");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_IBF_Level_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBFMIN_Test_Given_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBFMIN_Insert_Given_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}

TEST_F(insert_test, delete_ibfmin_no_given_thresholds)
{
    configuration config{};
    initialization_args(config);
    config.number_expression_thresholds = 2;
    config.fpr = {0.05, 0.05};
    config.path_out = "IBFMIN_Test_Given_Del_";
    config.compressed = false;
    config.sequence_files = {data("mini_example.minimiser"),
                             data("mini_example.minimiser"),
                             data("mini_example.minimiser")};

    ibf_min(config);

    config.number_expression_thresholds = 2;
    config.path_out = "IBFMIN_Insert_Given_Del_";
    config.compressed = false;
    config.fpr = {0.05};
    config.sequence_files = {data("mini_example.minimiser")};
    ibf_min(config);
    config.delete_files = {1};
    delete_bin(config);
    config.path_in = "IBFMIN_Insert_Given_Del_";
    insert_min(config);

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_insert;

    load_ibf(ibf, "IBFMIN_Test_Given_Del_IBF_Level_0");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_Del_IBF_Level_0");
    EXPECT_TRUE((ibf == ibf_insert));

    load_ibf(ibf, "IBFMIN_Test_Given_Del_IBF_Level_1");
    load_ibf(ibf_insert, "IBFMIN_Insert_Given_Del_IBF_Level_1");
    EXPECT_TRUE((ibf == ibf_insert));

    std::vector<std::vector<uint16_t>> expressions_ibf{};
    read_levels<uint16_t>(expressions_ibf, "IBFMIN_Test_Given_Del_IBF_Levels.levels");
    std::vector<std::vector<uint16_t>> expressions_insert{};
    read_levels<uint16_t>(expressions_insert, "IBFMIN_Insert_Given_Del_IBF_Levels.levels");
    EXPECT_EQ(expressions_ibf, expressions_insert);

    std::vector<std::vector<double>> fpr_ibf{};
    read_levels<double>(fpr_ibf, "IBFMIN_Test_Given_Del_IBF_FPRs.fprs");
    std::vector<std::vector<double>> fpr_insert{};
    read_levels<double>(fpr_insert, "IBFMIN_Insert_Given_Del_IBF_FPRs.fprs");
    EXPECT_EQ(fpr_ibf, fpr_insert);
}
