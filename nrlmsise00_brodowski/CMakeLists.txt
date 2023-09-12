cmake_minimum_required(VERSION 3.19)

project(nrlmsise00_brodowski)

set(NRLMSISE00_SOURCES
    nrlmsise-00/nrlmsise-00.c
    nrlmsise-00/nrlmsise-00.h
    nrlmsise-00/nrlmsise-00_data.c
    # nrlmsise-00/nrlmsise-00_test.c
)

if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/nrlmsise-00")
    message("*** cloning nrlmsise-00")
    include(FetchContent)
    FetchContent_Declare(
        nrlmsise00_brodowski
        GIT_REPOSITORY  https://git.linta.de/~brodo/nrlmsise-00.git
        GIT_TAG         master
        # FETCHCONTENT_BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nrlmsise-00
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    FetchContent_MakeAvailable(nrlmsise00_brodowski)
else()
    message("*** nrlmsise-00 already cloned")
endif()

add_library(nrlmsise00_brodowski
    INTERFACE
    ${NRLMSISE00_SOURCES}
    )

target_include_directories(nrlmsise00_brodowski
    INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/nrlmsise-00)

target_sources(nrlmsise00_brodowski
    INTERFACE
    ${NRLMSISE00_SOURCES}
)
