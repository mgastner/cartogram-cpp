# How to Create and Run Tests

## Creating Tests

* Create a new file in the `tests` directory with the `.cpp` extension.
* Add the source file of the functionalities you want to test to the `CARTOGRAM_TEST_SOURCES_FROM_SRC` list in the `CMakeLists.txt` file in the `tests` directory.
    ```cmake
    set(CARTOGRAM_TEST_SOURCES_FROM_SRC
      "src/misc/string_to_decimal_converter.cpp"

      # Add additional test sources from src here if necessary
    )
    ```
  For example, `test_string_to_decimal_converter.cpp` tests the functionality of `string_to_decimal_converter.cpp` file in the `src/misc` directory; so, the source file is added to the list.

* Use Boost Test to write your tests. You may find the documentation [here](https://www.boost.org/doc/libs/1_84_0/libs/test/doc/html/index.html).

* Use the `test_string_to_decimal_converter.cpp` file as a template for your tests.

## Running Tests

* The tests are built with the project. To build the tests, run:
    ```bash
    make
    ```
    
* To run all the tests, run:
    ```bash
    ctest
    ```

* To get detailed information about the tests, run:
    ```bash
    ctest --verbose
    ```

* To see available tests, run:
    ```bash
    ctest --show-only
    ```

* To run a specific test, run:
    ```bash
    ctest -R <test_name>
    ```
    
    For example, to run the `test_string_to_decimal_converter` test, run:
    ```bash
    ctest -R test_string_to_decimal_converter
    ```
* To detail the output of a specific test, run:
    ```bash
    ctest -R <test_name> --verbose
    ```
    
    For example, to detail the output of the `test_string_to_decimal_converter` test, run:
    ```bash
    ctest -R test_string_to_decimal_converter --verbose
    ```
