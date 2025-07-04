# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# cmake-format: off

# hibf
set (NEEDLE_HIBF_VERSION 36c223527d7d4eb2bcdcd88d47c71984f50a00a3 CACHE STRING "")
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${NEEDLE_HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF"
)

# seqan3
set (NEEDLE_SEQAN3_VERSION 6dfa1b442d1fabd07024edcc37a29b61d5beae8f CACHE STRING "")
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${NEEDLE_SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# sharg
set (NEEDLE_SHARG_VERSION be113bcffe49c0d62cbd65a191820f05386aa8da CACHE STRING "")
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${NEEDLE_SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "BUILD_TESTING OFF"
)

# googletest
set (NEEDLE_GOOGLETEST_VERSION 1.17.0 CACHE STRING "")
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${NEEDLE_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# robin-hood-hashing
set (NEEDLE_ROBIN_HOOD_VERSION 7697343363af4cc3f42cab17be49e6af9ab181e2 CACHE STRING "")
CPMDeclarePackage (robin-hood
                   NAME robin-hood
                   GIT_TAG ${NEEDLE_ROBIN_HOOD_VERSION} # master
                   GITHUB_REPOSITORY martinus/robin-hood-hashing
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37 CACHE STRING "")
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)

# cmake-format: on
