find_package(Boost REQUIRED CONFIG)

# Clang 22 reports each __COUNTER__ in the Boost.Test macros as a C2y extension
# (-Wc2y-extensions), which -Werror turns into an error. The option is set on the
# source files because CMake puts source options after the -Wpedantic that
# cartogram_lib passes on, and only a later option overrides an earlier one.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wc2y-extensions HAVE_WARNING_C2Y_EXTENSIONS)

file(GLOB UNIT_SRCS CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp)

foreach(src IN LISTS UNIT_SRCS)
  get_filename_component(test_name ${src} NAME_WE)
  add_executable(${test_name} ${src})
  target_compile_definitions(${test_name}
    PRIVATE BOOST_TEST_MODULE=${test_name}
  )
  target_link_libraries(${test_name}
    PRIVATE
    boost::boost
    cartogram_lib
  )
  if(HAVE_WARNING_C2Y_EXTENSIONS)
    set_property(SOURCE ${src} APPEND PROPERTY COMPILE_OPTIONS
      -Wno-c2y-extensions
    )
  endif()
  add_test(NAME ${test_name} COMMAND ${test_name})
  set_tests_properties(${test_name} PROPERTIES LABELS "unit" TIMEOUT 1)
endforeach()
