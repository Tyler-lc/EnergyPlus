include_directories("${CMAKE_SOURCE_DIR}/third_party/kiva/src")
include_directories("${CMAKE_BINARY_DIR}/third_party/kiva/src/libkiva")
include_directories( SYSTEM "${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/boost-1.61.0/")
include_directories( SYSTEM "${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/lis-1.5.66/include/")

configure_file("${CMAKE_SOURCE_DIR}/third_party/cmake/CMakeLists-kiva.txt" "${CMAKE_BINARY_DIR}/third_party/cmake/CMakeLists.txt")
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/third_party/cmake)
execute_process(COMMAND ${CMAKE_COMMAND} --build . WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/third_party/cmake)

add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/lis-1.5.66/")
add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/src/libkiva/")
if (BUILD_GROUND_PLOT)
  add_definitions("-DBOOST_ALL_NO_LIB")  
  include_directories( SYSTEM "${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/mathgl-2.3.5.1/include/")
  include_directories( SYSTEM "${CMAKE_BINARY_DIR}/third_party/kiva/vendor/mathgl-2.3.5.1/include/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/zlib-1.2.8/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/libpng-1.6.23/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/mathgl-2.3.5.1/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/boost-1.61.0/boost/filesystem/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/vendor/boost-1.61.0/boost/system/")
  add_subdirectory("${CMAKE_SOURCE_DIR}/third_party/kiva/src/libgroundplot/")
endif()
