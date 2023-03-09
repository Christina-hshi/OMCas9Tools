file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

# # Add include directories.
# include_directories (${SEQAN_INCLUDE_DIRS})
# include_directories()

add_library(lib_core ${sources} ${headers})
target_include_directories(lib_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add executable and link against SeqAn dependencies.
add_executable(OMCas9Tools main.cpp)
target_link_libraries (OMCas9Tools -static lib_core boost_program_options)
