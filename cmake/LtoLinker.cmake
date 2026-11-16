# Linker for full LTO with clang on Linux.
#
# cartogram_find_lto_linker(<options_var> <linker_var> [<compile_option>...])
#
# Looks for ld.lld in the directory of the compiler first and then on the
# search path, the order in which clang looks for it. The cache variable
# CARTOGRAM_LTO_LD_LLD overrides the search. A test program, compiled with the
# given compile options, is then linked with LTO and this ld.lld, because an
# lld older than clang may fail to read clang's bitcode.
#
# On success, <options_var> holds the link options, and <linker_var> the path
# of ld.lld. The option --ld-path makes clang run exactly this ld.lld, whatever
# -B, COMPILER_PATH, or the PATH say, and -fuse-ld=lld tells clang that the
# linker is lld. If no ld.lld is found, or if the test link fails, both
# variables are empty and a warning says why; the caller then builds without
# LTO.

include_guard(GLOBAL)
include(CheckCXXSourceCompiles)

function(cartogram_find_lto_linker options_var linker_var)
  set(${options_var} "" PARENT_SCOPE)
  set(${linker_var} "" PARENT_SCOPE)

  file(REAL_PATH "${CMAKE_CXX_COMPILER}" compiler_real_path)
  cmake_path(GET compiler_real_path PARENT_PATH compiler_dir)
  string(REGEX MATCH "^[0-9]+" clang_major "${CMAKE_CXX_COMPILER_VERSION}")

  # Only the name ld.lld: lld chooses its flavor by the name of the file, and
  # --ld-path passes that name on
  find_program(CARTOGRAM_LTO_LD_LLD NAMES ld.lld
    PATHS "${compiler_dir}" NO_DEFAULT_PATH)
  find_program(CARTOGRAM_LTO_LD_LLD NAMES ld.lld)

  if(NOT CARTOGRAM_LTO_LD_LLD)
    message(WARNING
      "ld.lld not found in ${compiler_dir} or on the PATH; disabling IPO. "
      "Install the lld of the same version as clang, for example: "
      "sudo apt install lld-${clang_major}")
    return()
  endif()

  set(link_options -fuse-ld=lld "--ld-path=${CARTOGRAM_LTO_LD_LLD}")

  list(JOIN ARGN " " CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LINK_OPTIONS -flto=full ${link_options})
  set(CMAKE_REQUIRED_QUIET TRUE)

  # Run the test link at every configuration, so that a change of the
  # compiler or of the flags cannot leave a stale result in the cache
  unset(CARTOGRAM_LTO_LINK_WORKS CACHE)
  check_cxx_source_compiles([[
    struct Shape {
      virtual ~Shape() = default;
      virtual int sides() const { return 0; }
    };
    struct Square final : Shape {
      int sides() const override { return 4; }
    };
    int main()
    {
      const Square square;
      const Shape &shape = square;
      return shape.sides() == 4 ? 0 : 1;
    }
  ]] CARTOGRAM_LTO_LINK_WORKS)

  if(NOT CARTOGRAM_LTO_LINK_WORKS)
    message(WARNING
      "${CARTOGRAM_LTO_LD_LLD} failed to link a test program with LTO; "
      "disabling IPO. Install the lld of the same version as clang, for "
      "example: sudo apt install lld-${clang_major}")
    return()
  endif()

  set(${options_var} ${link_options} PARENT_SCOPE)
  set(${linker_var} "${CARTOGRAM_LTO_LD_LLD}" PARENT_SCOPE)
endfunction()
