cmake_minimum_required(VERSION 3.16)
project(pushpull_pp_sdl2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED sdl2)

add_executable(pushpull_pp src/main.cpp)
target_include_directories(pushpull_pp PRIVATE ${SDL2_INCLUDE_DIRS})
target_compile_options(pushpull_pp PRIVATE ${SDL2_CFLAGS_OTHER} -Wall -Wextra -pedantic)
target_link_libraries(pushpull_pp PRIVATE ${SDL2_LIBRARIES})
