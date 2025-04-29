// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <filesystem>

#include <sharg/parser.hpp>

#include <seqan3/core/debug_stream.hpp>

#include "estimate.hpp"
#include "ibf.hpp"
#include "shared.hpp"

uint32_t helper_window{60};
uint64_t helper_shape{};
uint64_t helper_seed{0x8F'3F'73'B5'CF'1C'9A'DEULL};

struct positive_integer_validator
{
    using option_value_type = size_t;

    void operator()(option_value_type const & value) const
    {
        if (value < 1)
            throw sharg::validation_error{"The value must greater than 0."};
    }

    std::string get_help_page_message() const
    {
        return "Value must be greater than 0.";
    }
};

void add_parser_meta(sharg::parser & parser)
{
    parser.info.author = "Mitra Darvish";
    parser.info.citation = "Needle: a fast and space-efficient prefilter for estimating the quantification of very "
                           "large collections of expression experiments; Mitra Darvish, Enrico Seiler, Svenja "
                           "Mehringer, René Rahn, and Knut Reinert; Bioinformatics, Volume 38, Issue 17, 1 September "
                           "2022, Pages 4100-4108. doi: https://doi.org/10.1093/bioinformatics/btac492";
    parser.info.date = NEEDLE_DATE;
    // REUSE-IgnoreStart
    parser.info.long_copyright = "This application uses SPDX identifiers.\n"
                                 "SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin\n"
                                 "SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik\n"
                                 "SPDX-License-Identifier: BSD-3-Clause";
    // REUSE-IgnoreEnd
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.url = "https://github.com/seqan/needle";
    parser.info.version = NEEDLE_VERSION;
}

void initialise_configuration(sharg::parser & parser, configuration & config)
{
    parser.add_option(
        config.k,
        sharg::config{.short_id = 'k', .long_id = "kmer", .description = "Define k-mer size for the minimisers."});
    parser.add_option(
        helper_window,
        sharg::config{.short_id = 'w', .long_id = "window", .description = "Define window size for the minimisers."});
    parser.add_option(
        helper_shape,
        sharg::config{.short_id = '\0',
                      .long_id = "shape",
                      .description =
                          "Define a shape for the minimisers by the decimal of a bitvector, where 0 symbolizes a "
                          "position to be ignored, 1 a position considered.",
                      .default_message = "Ungapped shape of length k."});
    parser.add_option(
        helper_seed,
        sharg::config{.short_id = '\0', .long_id = "seed", .description = "Define seed for the minimisers."});
    parser.add_option(config.path_out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "Directory where output files should be saved.",
                                    .validator = sharg::output_directory_validator{}});
    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});
}

void initialise_arguments_ibf(sharg::parser & parser, configuration & config)
{
    parser.add_flag(config.compressed,
                    sharg::config{.short_id = 'c', .long_id = "compressed", .description = "Use compressed IBFS."});

    parser.add_option(config.fpr,
                      sharg::config{.short_id = 'f',
                                    .long_id = "fpr",
                                    .description =
                                        "List of bin false positive rate per expression level. If only one is given"
                                        ", then that fpr is used for all expression levels."});

    parser.add_option(config.expression_thresholds,
                      sharg::config{.short_id = 'e',
                                    .long_id = "expression_thresholds",
                                    .description = "Which expression thresholds should be used for"
                                                   " constructing the IBFs."});

    parser.add_option(config.number_expression_thresholds,
                      sharg::config{.short_id = 'l',
                                    .long_id = "number_expression_thresholds",
                                    .description = "Number of expression thresholds. "
                                                   "Can be set alternatively to expression_thresholds, then "
                                                   "the expression thresholds are determined automatically."});

    parser.add_option(config.num_hash,
                      sharg::config{.short_id = 'n',
                                    .long_id = "hash",
                                    .description = "Number of hash functions that should be used when constructing "
                                                   "one IBF."});
}

void parsing(sharg::parser & parser, configuration & config)
{
    helper_window = config.w_size.get();
    helper_seed = config.s.get();
    parser.parse();
    config.w_size = seqan3::window_size{helper_window};
    if (helper_shape == 0)
        config.shape = seqan3::ungapped{config.k};
    else
        config.shape = seqan3::bin_literal{helper_shape};
    config.s = seqan3::seed{adjust_seed(config.k, helper_seed)};
}

// Initialize arguments for ibf and minimiser
void initialise_arguments_minimiser(sharg::parser & parser, configuration & config)
{
    parser.add_option(config.include_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "include",
                                    .description = "Sequence file containing minimizers, only those "
                                                   "minimizers will be considered.",
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.exclude_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "exclude",
                                    .description = "Sequence file containing minimizers that should "
                                                   "not be stored.",
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.samples,
                      sharg::config{.short_id = '\0',
                                    .long_id = "samples",
                                    .description = "Define which samples belong together, sum has to be "
                                                   "equal to number of sequence files.",
                                    .default_message = "Every sequence file is one sample from one experiment."});

    parser.add_flag(config.paired,
                    sharg::config{.short_id = 'p',
                                  .long_id = "paired",
                                  .description = "Experiments are paired, i.e. two consecutive entries in the include "
                                                 "file list belong to the same experiment."});

    parser.add_option(
        config.cutoffs,
        sharg::config{.short_id = '\0',
                      .long_id = "cutoff",
                      .description = "Define for each sample, what number of found minimisers "
                                     "should be considered the result of a sequencing error "
                                     "and therefore be ignored.",
                      .default_message =
                          "Every sample has an automatically generated cutoff, which is based on the file size."});
}

void read_input_file_list(configuration & config)
{
    if (config.sequence_files[0].extension() != ".lst")
        return;

    std::ifstream fin{config.sequence_files[0]};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open file " + config.sequence_files[0].string() + " for reading."};

    config.sequence_files.clear();

    std::string line;
    while (std::getline(fin, line))
    {
        config.sequence_files.push_back(line);
    }
}

int run_needle_count(sharg::parser & parser)
{
    configuration config;
    initialise_configuration(parser, config);

    parser.info.short_description = "Get expression value depending on minimizers. This function is an alternative to "
                                    "pseudoaligners like kallisto. It estimates the expression value "
                                    "for all sequences in the genome file based on the exact minimiser occurrences of "
                                    "the given sequence files. Please run genome beforehand to create the genome file.";
    parser.add_positional_option(
        config.sequence_files,
        sharg::config{.description = "Sequence files.", .validator = sharg::input_file_validator{}});
    parser.add_option(config.include_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "include",
                                    .description = "A sequence file with transcripts",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});
    parser.add_option(config.genome_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "genome",
                                    .description = "A *.genome file created with the genome command.",
                                    .validator = sharg::input_file_validator{}});
    parser.add_flag(config.paired,
                    sharg::config{.short_id = 'p',
                                  .long_id = "paired",
                                  .description = "Experiments are paired, i.e. two consecutive entries in the include "
                                                 "file list belong to the same experiment."});

    parsing(parser, config);

    count(config);

    return 0;
}

int run_needle_count_genome(sharg::parser & parser)
{
    configuration config;
    initialise_configuration(parser, config);

    parser.info.short_description = "Creates the genome file necessary as an input to count.";
    parser.add_positional_option(config.genome_file,
                                 sharg::config{.description = "Please provide one sequence file with transcripts",
                                               .validator = sharg::input_file_validator{}});

    parser.add_option(config.exclude_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "exclude",
                                    .description = "Please provide one sequence file with minimizers to ignore.",
                                    .validator = sharg::input_file_validator{}});

    parsing(parser, config);

    count_genome(config);

    return 0;
}

int run_needle_estimate(sharg::parser & parser)
{
    configuration config{};
    parser.info.short_description = "Estimate expression value of transcript based on the Needle index.";

    config.path_out = "expressions.out";

    parser.add_positional_option(
        config.search_file,
        sharg::config{.description = "A sequence file", .validator = sharg::input_file_validator{}});

    parser.add_option(
        config.path_in,
        sharg::config{.short_id = 'i', .long_id = "in", .description = "Directory where input files can be found."});

    parser.add_option(config.path_out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "Directory where output files should be saved.",
                                    .validator = sharg::output_directory_validator{}});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_option(
        config.batch_size,
        sharg::config{
            .short_id = 'b',
            .long_id = "batch",
            .description =
                "Process at most this many queries at once. For example, if there are 100K queries to be estimated, "
                "and the batch size is 10K, there will be 10 batches. More batches, i.e. a smaller batch size, will "
                "reduce/cap the memory usage, but might result in slower performance due to repetitive index-IO.",
            .validator = positive_integer_validator{}});

    parser.add_flag(config.normalization_method,
                    sharg::config{.short_id = 'm',
                                  .long_id = "normalization-mode",
                                  .description =
                                      "Use normalization. Normalization is achieved by "
                                      "dividing the expression value with the expression threshold of the first "
                                      "ibf. Only make sense if every bin has its own expression "
                                      "thresholds (which is the case if expression thresholds "
                                      "were generated automatically)."});

    parsing(parser, config);

    call_estimate(config);

    return 0;
}

int run_needle_ibf(sharg::parser & parser)
{
    configuration config{};

    initialise_configuration(parser, config);
    initialise_arguments_ibf(parser, config);
    initialise_arguments_minimiser(parser, config);

    parser.info.short_description = "Constructs the Needle index.";
    parser.add_positional_option(config.sequence_files,
                                 sharg::config{.description =
                                                   "Please provide at least one sequence file OR provide one file "
                                                   "containing all sequence files with the extension '.lst'.",
                                               .validator = sharg::input_file_validator{}});

    parser.add_option(config.experiment_names,
                      sharg::config{.short_id = '\0',
                                    .long_id = "experiment-names",
                                    .description = "Use experiment names from the given txt file."});

    parser.add_option(config.expression_by_genome_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "levels-by-genome",
                                    .description = "Sequence file containing minimizers, only those minimizers will be "
                                                   "considered for determining the expression thresholds.",
                                    .validator = sharg::input_file_validator{}});

    parser.add_flag(config.ram_friendly,
                    sharg::config{.short_id = '\0',
                                  .long_id = "ram",
                                  .description = "When multithreading, prioritize lower RAM usage over speed."});

    parsing(parser, config);
    read_input_file_list(config);

    ibf(config);

    return 0;
}

int run_needle_insert(sharg::parser & parser)
{
    configuration config{};

    initialise_arguments_minimiser(parser, config);

    parser.info.short_description = "Inserts into a given uncompressed Needle index.";

    parser.add_flag(config.compressed,
                    sharg::config{.short_id = 'c', .long_id = "compressed", .description = "Use compressed IBFS."});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_option(
        config.path_in,
        sharg::config{.short_id = 'i', .long_id = "in", .description = "Directory where input files can be found."});

    parser.add_positional_option(config.sequence_files,
                                 sharg::config{.description =
                                                   "Please provide at least one sequence file OR provide one file "
                                                   "containing all sequence files with the extension '.lst'.",
                                               .validator = sharg::input_file_validator{}});

    parser.add_option(config.experiment_names,
                      sharg::config{.short_id = '\0',
                                    .long_id = "experiment-names",
                                    .description = "Use experiment names from the given txt file."});

    parser.add_option(config.expression_by_genome_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "levels-by-genome",
                                    .description = "Sequence file containing minimizers, only those minimizers will be "
                                                   "considered for determining the expression thresholds.",
                                    .validator = sharg::input_file_validator{}});

    parser.add_flag(config.ram_friendly,
                    sharg::config{.short_id = '\0',
                                  .long_id = "ram",
                                  .description = "When multithreading, prioritize lower RAM usage over speed."});

    parser.parse();
    if (std::filesystem::exists(config.path_in.string() + "IBF_Level_0"))
        config.samplewise = true;
    read_input_file_list(config);

    insert(config);

    return 0;
}

int run_needle_ibf_min(sharg::parser & parser)
{
    configuration config{};

    parser.info.short_description = "Constructs the Needle index from the minimiser files created by needle minimiser.";
    parser.add_positional_option(config.sequence_files,
                                 sharg::config{.description =
                                                   "Please provide at least one minimiser file OR provide one file "
                                                   "containing all minimiser files with the extension '.lst'.",
                                               .validator = sharg::input_file_validator{}});

    parser.add_option(config.path_out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "Directory where output files should be saved.",
                                    .validator = sharg::output_directory_validator{}});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_option(config.expression_by_genome_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "levels-by-genome",
                                    .description = "Sequence file containing minimizers, only those minimizers will be "
                                                   "considered for determining the expression thresholds.",
                                    .validator = sharg::input_file_validator{}});

    initialise_arguments_ibf(parser, config);

    parsing(parser, config);
    read_input_file_list(config);

    ibf_min(config);

    return 0;
}

int run_needle_insert_min(sharg::parser & parser)
{
    configuration config{};

    parser.info.short_description = "Constructs the Needle index from the minimiser files created by needle minimiser.";

    parser.add_positional_option(config.sequence_files,
                                 sharg::config{.description =
                                                   "Please provide at least one minimiser file OR provide one file "
                                                   "containing all minimiser files with the extension '.lst'.",
                                               .validator = sharg::input_file_validator{}});

    parser.add_option(config.path_out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "Directory where output files should be saved.",
                                    .validator = sharg::output_directory_validator{}});

    parser.add_flag(config.compressed,
                    sharg::config{.short_id = 'c', .long_id = "compressed", .description = "Use compressed IBFS."});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_option(
        config.path_in,
        sharg::config{.short_id = 'i', .long_id = "in", .description = "Directory where input files can be found."});

    parser.add_option(config.expression_by_genome_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "levels-by-genome",
                                    .description = "Sequence file containing minimizers, only those minimizers will be "
                                                   "considered for determining the expression thresholds.",
                                    .validator = sharg::input_file_validator{}});

    parser.parse();
    if (std::filesystem::exists(config.path_in.string() + "IBF_Level_0"))
        config.samplewise = true;
    read_input_file_list(config);

    insert_min(config);

    return 0;
}

int run_needle_delete_bin(sharg::parser & parser)
{
    configuration config{};

    parser.info.short_description = "Delete experiments specified by their position from the Needle index.";

    parser.add_positional_option(config.delete_files,
                                 sharg::config{.description = "Please provide at least one position to be deleted."});

    parser.add_option(config.path_out,
                      sharg::config{.short_id = 'o',
                                    .long_id = "out",
                                    .description = "Directory where output files should be saved.",
                                    .validator = sharg::output_directory_validator{}});

    parser.add_flag(config.compressed,
                    sharg::config{.short_id = 'c', .long_id = "compressed", .description = "Use compressed IBFS."});

    parser.add_option(config.threads,
                      sharg::config{.short_id = 't',
                                    .long_id = "threads",
                                    .description = "Number of threads to use.",
                                    .validator = positive_integer_validator{}});

    parser.add_option(
        config.path_in,
        sharg::config{.short_id = 'i', .long_id = "in", .description = "Directory where input files can be found."});

    parser.parse();

    if (std::filesystem::exists(config.path_in.string() + "IBF_Level_0"))
        config.samplewise = true;
    delete_bin(config);

    return 0;
}

int run_needle_minimiser(sharg::parser & parser)
{
    configuration config{};
    initialise_configuration(parser, config);
    initialise_arguments_minimiser(parser, config);

    parser.info.short_description = "Calculates minimiser for given experiments.";

    parser.add_positional_option(config.sequence_files,
                                 sharg::config{.description =
                                                   "Please provide at least one sequence file OR provide one file "
                                                   "containing all sequence files with the extension '.lst'.",
                                               .validator = sharg::input_file_validator{}});

    parser.add_flag(config.ram_friendly,
                    sharg::config{.short_id = '\0',
                                  .long_id = "ram",
                                  .description = "When multithreading, prioritize lower RAM usage over speed."});

    parsing(parser, config);
    read_input_file_list(config);

    minimiser(config);

    return 0;
}

int main(int argc, char const ** argv)
{
    try
    {
        sharg::parser needle_parser{
            "needle",
            argc,
            argv,
            sharg::update_notifications::on,
            {"count", "delete", "estimate", "genome", "ibf", "ibfmin", "insert", "insertmin", "minimiser"}};
        add_parser_meta(needle_parser);
        needle_parser.info.description.push_back(
            "Needle allows you to build an Interleaved Bloom Filter (IBF) with the "
            "command ibf or estimate the expression of transcripts with the command "
            "estimate.");

        needle_parser.parse();

        sharg::parser & sub_parser = needle_parser.get_sub_parser();
        add_parser_meta(sub_parser);
        if (sub_parser.info.app_name == std::string_view{"needle-count"})
            run_needle_count(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-delete"})
            run_needle_delete_bin(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-genome"})
            run_needle_count_genome(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-estimate"})
            run_needle_estimate(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-ibf"})
            run_needle_ibf(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-ibfmin"})
            run_needle_ibf_min(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-insert"})
            run_needle_insert(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-insertmin"})
            run_needle_insert_min(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"needle-minimiser"})
            run_needle_minimiser(sub_parser);
    }
    catch (std::exception const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        std::exit(-1);
    }

    return 0;
}
