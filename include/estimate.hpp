// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <filesystem>

#include "shared.hpp"

/*! \brief Function, which calls the estimate function.
*  \param args          The arguments estimate and ibf use.
*  \param estimate_args The estimate arguments.
*/
void call_estimate(configuration & config);
