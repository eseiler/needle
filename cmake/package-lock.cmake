# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# seqan3
set (NEEDLE_SEQAN3_VERSION d4a7c88fde0311e12e98e7822da772b99c887cb5)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${NEEDLE_SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# googletest
set (NEEDLE_GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME googletest
                   VERSION ${NEEDLE_GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                           "CMAKE_CXX_STANDARD 20"
)

# robin-hood-hashing
set (NEEDLE_ROBIN_HOOD_VERSION 7697343363af4cc3f42cab17be49e6af9ab181e2)
CPMDeclarePackage (robin-hood
                   NAME robin-hood
                   GIT_TAG ${NEEDLE_ROBIN_HOOD_VERSION}
                   GITHUB_REPOSITORY martinus/robin-hood-hashing
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION}
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)