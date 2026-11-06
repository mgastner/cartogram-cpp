find_package(Boost REQUIRED CONFIG)

# Slower numerical tests of the conformality core (e.g., optimizer recovery
# tests). Fast checks belong in tests/unit, whose tests time out after 1 s.
file(GLOB CONFORMALITY_TEST_SRCS CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp)

foreach(src IN LISTS CONFORMALITY_TEST_SRCS)
  get_filename_component(test_name ${src} NAME_WE)
  add_executable(${test_name} ${src})
  target_compile_definitions(${test_name}
    PRIVATE BOOST_TEST_MODULE=${test_name}
  )
  target_link_libraries(${test_name}
    PRIVATE
    boost::boost
    conformality_core
  )
  add_test(NAME ${test_name} COMMAND ${test_name})
  set_tests_properties(${test_name}
    PROPERTIES LABELS "conformality" TIMEOUT 300)
endforeach()
