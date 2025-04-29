// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "shared.hpp"

/*!\brief Get the concrete expression values (= median of all counts of one transcript) for given experiments.
*         This function can be used to estimate how good the median approach can be, if all count values are available.
* \param config               The minimiser arguments to use (seed, shape, window size).
*/
void count(configuration & config);

/*!\brief Creates a set of minimizers to ignore, which should be used as an input to count.
* \param config               The minimiser arguments to use (seed, shape, window size).
*/
void count_genome(configuration & config);

/*!\brief Reads a binary file that needle minimiser creates.
* \param filename           The filename of the binary file.
* \param hash_table         The hash table to store minimisers into.

*/
void read_binary(std::filesystem::path filename, robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table);

/*!\brief Reads the beginning of a binary file that needle minimiser creates.
* \param config               Min arguments.
*/
void read_binary_start(configuration & config);

/*! \brief Creates IBFs.
 * \param sequence_files  A vector of sequence file paths.
 * \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                        struct ibf_arguments.
 * \param minimiser_args  The minimiser specific arguments to use.
 * \param fpr             The average false positive rate that should be used.
 * \param cutoffs         List of cutoffs.
 * \param expression_by_genome_file File that contains the only minimisers that should be considered for the
 *                                  determination of the expression thresholds.
 * \param num_hash        The number of hash functions to use.
 *  \returns The expression thresholds per experiment.
 */
std::vector<uint16_t> ibf(configuration & config);

/*! \brief Creates IBFs based on the minimiser files
 * \param minimiser_files A vector of minimiser file paths.
 * \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
 *                        struct ibf_arguments.
 * \param fpr             The average false positive rate that should be used.
 * \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
 *                                  determination of the expression_thresholds.
 * \param num_hash        The number of hash functions to use.
 *  \returns The expression thresholds per experiment.
 */
std::vector<uint16_t> ibf_min(configuration & config);

/*! \brief Create minimiser and header files.
* \param sequence_files  A vector of sequence file paths.
* \param config            The minimiser arguments to use (seed, shape, window size).
* \param minimiser_args  The minimiser specific arguments to use.
* \param cutoffs         List of cutoffs.
*/
void minimiser(configuration & config);

/*! \brief Insert into IBFs.
* \param sequence_files  A vector of sequence file paths.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param minimiser_args  The minimiser specific arguments to use.
* \param cutoffs         List of cutoffs.
* \param expression_by_genome_file File that contains the only minimisers that should be considered for the
*                                  determination of the expression thresholds.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*  \returns The expression thresholds per experiment.
*/
std::vector<uint16_t> insert(configuration & config);

/*! \brief Insert into IBFs based on the minimiser files
* \param minimiser_files A vector of minimiser file paths.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param expression_by_genome_file File that contains the only minimisers that should be comnsidered for the
*                                  determination of the expression_thresholds.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*  \returns The expression thresholds per experiment.
*/
std::vector<uint16_t> insert_min(configuration & config);

/*! \brief Delete bins from ibfs
* \param delete_files    A vector of integers specifiying the bins to delete.
* \param ibf_args        The IBF specific arguments to use (bin size, number of hash functions, ...). See
*                        struct ibf_arguments.
* \param path_in         Input directory.
* \param samplewise      True, if expression levels were set beforehand.
*/
void delete_bin(configuration & config);
