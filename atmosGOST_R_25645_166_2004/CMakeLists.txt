cmake_minimum_required(VERSION 3.19)

project(atmosGOST_R_25645_166_2004)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    atmosGOST_R_25645_166_2004.cpp
    atmosGOST_R_25645_166_2004.h
    table_1_GOST_R_25645_166_2004.h
    table_2_GOST_R_25645_166_2004.h
    table_3_GOST_R_25645_166_2004.h
    tables_4_9_GOST_R_25645_166_2004.h
    tables_10_11_GOST_R_25645_166_2004.h
)

#add_library(atmosGOST_R_25645_166_2004_tmp STATIC ${SOURCES})
#target_include_directories(atmosGOST_R_25645_166_2004_tmp PUBLIC .)

add_library(atmosGOST_R_25645_166_2004
    INTERFACE
    ${SOURCES})
target_include_directories(atmosGOST_R_25645_166_2004
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(atmosGOST_R_25645_166_2004
    INTERFACE
    ${SOURCES}
)
