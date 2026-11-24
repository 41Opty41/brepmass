cmake_minimum_required(VERSION 3.16)
project(brepmass VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(brepmass_static STATIC
    src/mass.c
    src/tessellate.c
    src/io_mesh.c
)
target_include_directories(brepmass_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(brepmass SHARED
    src/mass.c
    src/tessellate.c
    src/io_mesh.c
)
target_compile_definitions(brepmass PRIVATE BM_BUILDING_DLL)
target_include_directories(brepmass PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(UNIX AND NOT APPLE)
    target_link_libraries(brepmass_static PUBLIC m)
    target_link_libraries(brepmass PUBLIC m)
elseif(APPLE)
    target_link_libraries(brepmass_static PUBLIC m)
    target_link_libraries(brepmass PUBLIC m)
endif()

if(MSVC)
    target_compile_options(brepmass_static PRIVATE /W4 /WX)
    target_compile_options(brepmass PRIVATE /W4 /WX)
else()
    target_compile_options(brepmass_static PRIVATE -Wall -Wextra -Werror)
    target_compile_options(brepmass PRIVATE -Wall -Wextra -Werror)
endif()

add_executable(brepmass_cli cli/main.c)
set_target_properties(brepmass_cli PROPERTIES OUTPUT_NAME brepmass)
target_link_libraries(brepmass_cli PRIVATE brepmass_static)

enable_testing()
add_executable(test_brepmass tests/test_brepmass.c)
target_link_libraries(test_brepmass PRIVATE brepmass_static)
add_test(NAME test_brepmass COMMAND test_brepmass)
