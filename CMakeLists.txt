cmake_minimum_required(VERSION 3.19)

project(GOST_R_25645_166_2004)

enable_testing()

add_subdirectory(atmosGOST_R_25645_166_2004)
add_subdirectory(test)
add_subdirectory(nrlmsise00_brodowski )
add_subdirectory(example_app)
