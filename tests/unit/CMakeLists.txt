find_package(Boost REQUIRED CONFIG)

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
  add_test(NAME ${test_name} COMMAND ${test_name})
  set_tests_properties(${test_name} PROPERTIES LABELS "unit" TIMEOUT 1)
endforeach()
