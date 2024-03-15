#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2014-2016, Steven Vancoillie                           *
#               2014-2022, Ignacio Fdez. Galván                        *
#               2020, Gerardo Raggi                                    *
#***********************************************************************

########################################################################
#                                                                      #
# Compiler-specific settings                                           #
#                                                                      #
########################################################################

# Set the default compiler flags to pick from
#================================================
# For new compilers, check for the proper compiler ID on the following web page:
# http://www.cmake.org/cmake/help/v2.8.12/cmake.html#variable:CMAKE_LANG_COMPILER_ID
# The regular compiler options have to be compatible with other flags (so they
# might be added or removed depending on certain options), while the flags for
# build targets are mutually exclusive (only one of them will be active).

# Ideally, build targets should accomplish the features in the table below.
# legend: opt = optimizations; bt = backtrace (call stack); trap = catch
# floating point exceptions; warn = produce warnings during compilation
#                    opt  bt   trap  warn
#                    ===  ===  ====  ====
# - Debug            NO   YES  NO    YES
# - Garble           YES  YES  YES   YES
# - RelWithDebInfo   YES  YES  NO    YES
# - Release          YES  NO   NO    NO
# - Fast             YES  NO   NO    NO
#
# Note: trapping floating point exceptions is not compatible with LAPACK
# (which produces and handles the exceptions), so a better alternative
# has to be found

#======================================#
# GNU Compiler Collection (GCC)
#======================================#
# defines __GNUC__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "GNU")
set (CXXFLAGS_GNU_BASIC "")
set (CXXFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELEASE "-O2")
set (CXXFLAGS_GNU_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "GNU")
set (CFLAGS_GNU_BASIC "-std=gnu99")
set (CFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_GNU_RELEASE "-O2")
set (CFLAGS_GNU_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "GNU")
set (FFLAGS_GNU_BASIC "-fno-aggressive-loop-optimizations") # Or Work is assumed to be 8 elements long
set (FFLAGS_GNU_PREPROCESS "-cpp")
set (FFLAGS_GNU_OPENMP "-fopenmp")
set (FFLAGS_GNU_NO_OPENMP "-fmax-stack-var-size=1048576")
set (FFLAGS_GNU_ILP64 "-fdefault-integer-8")
# build targets
set (FFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (FFLAGS_GNU_GARBLE "-O2 -g -Wall -finit-real=snan -finit-integer=730432726 -frounding-math -fsignaling-nans")
set (FFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (FFLAGS_GNU_RELEASE "-O2")
set (FFLAGS_GNU_FAST "-O3")
set (FFLAGS_GNU_BOUNDS "-fsanitize=address -fno-omit-frame-pointer")

set (MIN_GCC_VERSION "4.8")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS ${MIN_GCC_VERSION})
  message (FATAL_ERROR "At least version ${MIN_GCC_VERSION} is required for the GNU Fortran compiler.")
endif ()

# Add runtime checks
# (No boundary checks because of the Work arrays)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "6")
    set (FFLAGS_GNU_GARBLE "${FFLAGS_GNU_GARBLE} -fcheck=array-temps,do,mem,pointer,recursion")
  else ()
    set (FFLAGS_GNU_GARBLE "${FFLAGS_GNU_GARBLE} -fcheck=all,no-bounds")
  endif ()
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "8")
    set (FFLAGS_GNU_GARBLE "${FFLAGS_GNU_GARBLE} -finit-derived")
  endif ()
endif ()

set (CXXFLAGS_GNU_BIGOT "-Wall -Wextra -Wpedantic -Werror")
set (CFLAGS_GNU_BIGOT "-Wall -Wextra -Wpedantic -Werror")
set (FFLAGS_GNU_BIGOT "-Wall -Wrealloc-lhs -Werror")
set (GNU_TRAP "-ffpe-trap=invalid,zero,overflow")

#======================================#
# LLVM C frontend (Clang)
#======================================#
# defines __clang__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Clang")
set (CXXFLAGS_Clang_BASIC "-std=gnu99")
set (CXXFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELEASE "-O2")
set (CXXFLAGS_Clang_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Clang")
set (CFLAGS_Clang_BASIC "-std=gnu99")
set (CFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_Clang_RELEASE "-O2")
set (CFLAGS_Clang_FAST "-O2")

set (CXXFLAGS_Clang_BIGOT "-Wall -Werror")
set (CFLAGS_Clang_BIGOT "-Wall -Werror")
#set (Clang_TRAP "???")

#======================================#
# LLVM C frontend CMake 3.0 (AppleClang)
#======================================#
# defines __clang__ and __apple_build_version__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "AppleClang")
set (CXXFLAGS_AppleClang_BASIC "-std=gnu99")
set (CXXFLAGS_AppleClang_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_AppleClang_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_AppleClang_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_AppleClang_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_AppleClang_RELEASE "-O2")
set (CXXFLAGS_AppleClang_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "AppleClang")
set (CFLAGS_AppleClang_BASIC "-std=gnu99")
set (CFLAGS_AppleClang_OPENMP "-fopenmp")
# build targets
set (CFLAGS_AppleClang_DEBUG "-O0 -g -Wall")
set (CFLAGS_AppleClang_GARBLE "-O2 -g -Wall")
set (CFLAGS_AppleClang_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_AppleClang_RELEASE "-O2")
set (CFLAGS_AppleClang_FAST "-O2")

set (CXXFLAGS_AppleClang_BIGOT "-Wall -Werror")
set (CFLAGS_AppleClang_BIGOT "-Wall -Werror")
#set (AppleClang_TRAP "???")

#======================================#
# Intel "classic" C++/C/Fortran Compilers
#======================================#
# defines __INTEL_COMPILER
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Intel")
set (CXXFLAGS_Intel_BASIC "")
set (CXXFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CXXFLAGS_Intel_DEBUG "-debug -w3")
set (CXXFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELEASE "-O2")
set (CXXFLAGS_Intel_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Intel")
set (CFLAGS_Intel_BASIC "-std=gnu99 -gcc-sys")
set (CFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CFLAGS_Intel_DEBUG "-debug -w3")
set (CFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CFLAGS_Intel_RELEASE "-O2")
set (CFLAGS_Intel_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Intel")
set (FFLAGS_Intel_BASIC "")
set (FFLAGS_Intel_PREPROCESS "-fpp")
set (FFLAGS_Intel_OPENMP "-qopenmp")
set (FFLAGS_Intel_ILP64 "-i8 -heap-arrays")
# build targets
set (FFLAGS_Intel_DEBUG "-debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_GARBLE "-O2 -debug -traceback -warn all,nodeclarations -check all,nobounds,noarg_temp_created")
set (FFLAGS_Intel_RELWITHDEBINFO "-O2 -fno-alias -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_RELEASE "-O2 -fno-alias -traceback")
set (FFLAGS_Intel_FAST "-Ofast -fno-alias")

set (MIN_IFORT_VERSION "13")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS ${MIN_IFORT_VERSION})
  message (FATAL_ERROR "At least version ${MIN_IFORT_VERSION} is required for the Intel Fortran compiler.")
endif ()

# Intel versions prior to 15 used -openmp
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CXXFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_C_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (FFLAGS_Intel_OPENMP "-openmp")
endif ()

#fix for MAC OS X
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  set (CXXFLAGS_Intel_BASIC "-headerpad_max_install_names ${CXXFLAGS_Intel_BASIC}")
  set (CFLAGS_Intel_BASIC "-headerpad_max_install_names ${CFLAGS_Intel_BASIC}")
  set (FFLAGS_Intel_BASIC "-headerpad_max_install_names ${FFLAGS_Intel_BASIC}")
endif ()

#note that "declarations" means "implicit none", which we don't want (yet)
#          "externals" means "implicit none (external)", ditto (only 19.1 and later)
#v.13 puts the full path of included files as a comment and then complains about truncated source lines
set (CXXFLAGS_Intel_BIGOT "-w3 -Werror -diag-disable remark")
set (CFLAGS_Intel_BIGOT "-w3 -Werror -diag-disable remark")
set (warnings "all,nodeclarations")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14")
    set (warnings "${warnings},notruncated_source")
  elseif (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "19.1")
    set (warnings "${warnings},noexternals")
  endif ()
endif ()
set (FFLAGS_Intel_BIGOT "-warn ${warnings} -warn errors -diag-disable remark")
set (Intel_TRAP "-fpe0")

#======================================#
# IntelLLVM C++/C/Fortran Compilers
#======================================#
# defines __INTEL_LLVM_COMPILER
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "IntelLLVM")
set (CXXFLAGS_IntelLLVM_BASIC "")
set (CXXFLAGS_IntelLLVM_OPENMP "-qopenmp")
# build targets
set (CXXFLAGS_IntelLLVM_DEBUG "-debug -Wall")
set (CXXFLAGS_IntelLLVM_GARBLE "-O2 -debug -Wall")
set (CXXFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -debug -Wall")
set (CXXFLAGS_IntelLLVM_RELEASE "-O2")
set (CXXFLAGS_IntelLLVM_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "IntelLLVM")
set (CFLAGS_IntelLLVM_BASIC "-std=gnu99")
set (CFLAGS_IntelLLVM_OPENMP "-qopenmp")
# build targets
set (CFLAGS_IntelLLVM_DEBUG "-debug -Wall")
set (CFLAGS_IntelLLVM_GARBLE "-O2 -debug -Wall")
set (CFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -debug -Wall")
set (CFLAGS_IntelLLVM_RELEASE "-O2")
set (CFLAGS_IntelLLVM_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "IntelLLVM")
set (FFLAGS_IntelLLVM_BASIC "-standard-semantics -fltconsistency")
set (FFLAGS_IntelLLVM_PREPROCESS "-fpp")
set (FFLAGS_IntelLLVM_OPENMP "-qopenmp")
set (FFLAGS_IntelLLVM_ILP64 "-i8 -heap-arrays")
# build targets
set (FFLAGS_IntelLLVM_DEBUG "-debug -traceback -warn all,nodeclarations")
set (FFLAGS_IntelLLVM_GARBLE "-O2 -debug -traceback -warn all,nodeclarations -check all,nobounds,noarg_temp_created")
set (FFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -fno-alias -debug -traceback -warn all,nodeclarations")
set (FFLAGS_IntelLLVM_RELEASE "-O2 -fno-alias -traceback")
set (FFLAGS_IntelLLVM_FAST "-Ofast -fno-alias")

#fix for MAC OS X
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  set (CXXFLAGS_IntelLLVM_BASIC "-headerpad_max_install_names ${CXXFLAGS_Intel_BASIC}")
  set (CFLAGS_IntelLLVM_BASIC "-headerpad_max_install_names ${CFLAGS_Intel_BASIC}")
endif ()

#note that "declarations" means "implicit none", which we don't want (yet)
#          "externals" means "implicit none (external)", ditto
set (CXXFLAGS_IntelLLVM_BIGOT "-Wall -Werror -diag-disable remark")
set (CFLAGS_IntelLLVM_BIGOT "-Wall -Werror -diag-disable remark")
set (FFLAGS_IntelLLVM_BIGOT "-warn all,nodeclarations,noexternals -warn errors -diag-disable remark")
set (IntelLLVM_TRAP "-fpe0")

#======================================#
# Oracle Developer Studio (Sun Studio)
#======================================#
# defines __SUNPRO_F90
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "SunPro")
set (CFLAGS_SunPro_BASIC "-std=gnu99")
set (CFLAGS_SunPro_OPENMP "-xopenmp")
# build targets
set (CFLAGS_SunPro_DEBUG "-g -traceback -ftrap=%none")
set (CFLAGS_SunPro_GARBLE "-O -g -traceback -ftrap=%none -xcheck=init_local")
set (CFLAGS_SunPro_RELWITHDEBINFO "-O -g -traceback -ftrap=%none")
set (CFLAGS_SunPro_RELEASE "-O -ftrap=%none")
set (CFLAGS_SunPro_FAST "-O -ftrap=%none")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "SunPro")
set (FFLAGS_SunPro_BASIC "") # 32bit: -ftrap=common,no%division,%none
set (FFLAGS_SunPro_PREPROCESS "-fpp")
set (FFLAGS_SunPro_OPENMP "-xopenmp")
set (FFLAGS_SunPro_ILP64 "-dalign -xtypemap=real:64,double:64,integer:64")
# build targets
set (FFLAGS_SunPro_DEBUG "-g -traceback -ftrap=%none")
set (FFLAGS_SunPro_GARBLE "-O -g -traceback ftrap=%none -xcheck=init_local")
set (FFLAGS_SunPro_RELWITHDEBINFO "-O -g -traceback -ftrap=%none")
set (FFLAGS_SunPro_RELEASE "-O -ftrap=%none")
set (FFLAGS_SunPro_FAST "-fast -dbl_align_all=yes -ftrap=%none")

# older versions do not suport gnu99
if (CMAKE_C_COMPILER_ID STREQUAL "SunPro" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "5.15.0")
  set (CFLAGS_SunPro_BASIC "-std=c99")
endif ()

set (CFLAGS_SunPro_BIGOT "-errwarn=%all")
set (FFLAGS_SunPro_BIGOT "-errwarn=%all")
set (SunPro_TRAP "-ftrap=common")

#======================================#
# NVIDIA HPC Compilers, a.k.a. Portland Group (PGI)
#======================================#
# defines __PGI
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "PGI" "NVHPC")
set (CXXFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_PGI_OPENMP "-mp")
# build targets
set (CXXFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELEASE "-O2")
set (CXXFLAGS_PGI_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "PGI" "NVHPC")
set (CFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CFLAGS_PGI_OPENMP "-mp")
# build targets
set (CFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELEASE "-O2")
set (CFLAGS_PGI_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "PGI" "NVHPC")
set (FFLAGS_PGI_BASIC "-Mbackslash")
set (FFLAGS_PGI_PREPROCESS "-Mpreprocess")
set (FFLAGS_PGI_OPENMP "-mp")
set (FFLAGS_PGI_ILP64 "-i8")
# build targets
set (FFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELEASE "-O2")
set (FFLAGS_PGI_FAST "-fast")

# The option does not exist, or I cannot find it
#set (CXXFLAGS_PGI_BIGOT "???")
#set (CFLAGS_PGI_BIGOT "???")
#set (FFLAGS_PGI_BIGOT "???")
#set (PGI_TRAP "???")

#======================================#
# Numerical Algorithms Group (NAG)
#======================================#
# defines NAGFOR
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "NAG")
set (FFLAGS_NAG_BASIC "-kind=byte")
set (FFLAGS_NAG_PREPROCESS "-fpp")
set (FFLAGS_NAG_OPENMP "-openmp")
set (FFLAGS_NAG_ILP64 "-i8")
# build targets
set (FFLAGS_NAG_DEBUG "-O0 -g -gline -ieee=full")
set (FFLAGS_NAG_GARBLE "-O2 -g -gline -nan")
set (FFLAGS_NAG_RELWITHDEBINFO "-O2 -g -ieee=full")
set (FFLAGS_NAG_RELEASE "-O2 -ieee=full")
set (FFLAGS_NAG_FAST "-O4 -ieee=full")

# Add runtime checks
# (No boundary checks [-C=array -C=calls] because of the Work arrays. [-C=dangling] causes segfault)
set (FFLAGS_NAG_GARBLE "${FFLAGS_NAG_GARBLE} -C=alias -C=bits -C=do -C=intovf -C=present -C=pointer -C=recursion")

# Fix for NAG compiler (similar to GNU Fortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" AND CMAKE_C_COMPILER_ID STREQUAL "GNU" AND (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin"))
  set (FFLAGS_NAG_BASIC "${FFLAGS_NAG_BASIC} -Wc,-fno-aggressive-loop-optimizations")
endif ()

# The option does not exist, or I cannot find it
#set (FFLAGS_NAG_BIGOT "???")
set (NAG_TRAP "-ieee=stop") # not compatible with -ieee=full

########################################################################
#                                                                      #
# Pick up the specified compiler                                       #
#                                                                      #
########################################################################

# C++ compiler (optional)
foreach (CXX_COMPILER ${SUPPORTED_CXX_COMPILERS})
  if (CXX_COMPILER STREQUAL "${CMAKE_CXX_COMPILER_ID}")
    set (CXXCID ${CMAKE_CXX_COMPILER_ID})
  endif ()
endforeach ()
if (CXXCID STREQUAL "NVHPC")
  set (CXXCID "PGI")
endif ()

# C compiler
foreach (C_COMPILER ${SUPPORTED_C_COMPILERS})
  if (C_COMPILER STREQUAL "${CMAKE_C_COMPILER_ID}")
    set (CCID ${CMAKE_C_COMPILER_ID})
  endif ()
endforeach ()
if (CCID STREQUAL "NVHPC")
  set (CCID "PGI")
endif ()

# Fortran compiler
foreach (Fortran_COMPILER ${SUPPORTED_Fortran_COMPILERS})
  if (Fortran_COMPILER STREQUAL "${CMAKE_Fortran_COMPILER_ID}")
    set (FCID ${CMAKE_Fortran_COMPILER_ID})
  endif ()
endforeach ()
if (FCID STREQUAL "NVHPC")
  set (FCID "PGI")
endif ()

########################################################################
#                                                                      #
# Overwrite any CMake defaults (some are strange, so we set our own)   #
#                                                                      #
########################################################################

if (NOT CMAKE_CXX_FLAGS_USER)
  set (CMAKE_CXX_FLAGS_DEFAULT "$ENV{CXXFLAGS} ${CXXFLAGS_${CXXCID}_BASIC}"
       CACHE STRING "C++ compiler flags." FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG ${CXXFLAGS_${CXXCID}_DEBUG}
       CACHE STRING "C++ compiler flags used for debug builds." FORCE)
  set (CMAKE_CXX_FLAGS_GARBLE ${CXXFLAGS_${CXXCID}_GARBLE}
       CACHE STRING "C++ compiler flags used for garble builds." FORCE)
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CXXFLAGS_${CXXCID}_RELWITHDEBINFO}
       CACHE STRING "C++ compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_CXX_FLAGS_RELEASE ${CXXFLAGS_${CXXCID}_RELEASE}
       CACHE STRING "C++ compiler flags used for release builds." FORCE)
  set (CMAKE_CXX_FLAGS_FAST ${CXXFLAGS_${CXXCID}_FAST}
       CACHE STRING "C++ compiler flags used for fast release builds." FORCE)
  set (CMAKE_CXX_FLAGS_USER "TRUE"
       CACHE INTERNAL "set to FALSE to reset to default C++ flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced (CLEAR
    CMAKE_CXX_FLAGS_DEFAULT
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_GARBLE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_FAST
  )
else ()
  mark_as_advanced (FORCE
    CMAKE_CXX_FLAGS_DEFAULT
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_GARBLE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_FAST
  )
endif ()
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEFAULT}")

if (NOT CMAKE_C_FLAGS_USER)
  set (CMAKE_C_FLAGS_DEFAULT "$ENV{CFLAGS} ${CFLAGS_${CCID}_BASIC}"
       CACHE STRING "C compiler flags." FORCE)
  set (CMAKE_C_FLAGS_DEBUG ${CFLAGS_${CCID}_DEBUG}
       CACHE STRING "C compiler flags used for debug builds." FORCE)
  set (CMAKE_C_FLAGS_GARBLE ${CFLAGS_${CCID}_GARBLE}
       CACHE STRING "C compiler flags used for garble builds." FORCE)
  set (CMAKE_C_FLAGS_RELWITHDEBINFO ${CFLAGS_${CCID}_RELWITHDEBINFO}
       CACHE STRING "C compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_C_FLAGS_RELEASE ${CFLAGS_${CCID}_RELEASE}
       CACHE STRING "C compiler flags used for release builds." FORCE)
  set (CMAKE_C_FLAGS_FAST ${CFLAGS_${CCID}_FAST}
       CACHE STRING "C compiler flags used for fast release builds." FORCE)
  set (CMAKE_C_FLAGS_USER "TRUE"
       CACHE INTERNAL "set to FALSE to reset to default C flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced (CLEAR
    CMAKE_C_FLAGS_DEFAULT
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_GARBLE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_FAST
  )
else ()
  mark_as_advanced (FORCE
    CMAKE_C_FLAGS_DEFAULT
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_GARBLE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_FAST
  )
endif ()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEFAULT}")

if (NOT CMAKE_Fortran_FLAGS_USER)
  set (CMAKE_Fortran_FLAGS_DEFAULT "$ENV{FFLAGS} ${FFLAGS_${FCID}_BASIC}"
       CACHE STRING "Fortran compiler flags." FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG ${FFLAGS_${FCID}_DEBUG}
       CACHE STRING "Fortran compiler flags used for debug builds." FORCE)
  set (CMAKE_Fortran_FLAGS_GARBLE ${FFLAGS_${FCID}_GARBLE}
       CACHE STRING "Fortran compiler flags used for garble builds." FORCE)
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${FFLAGS_${FCID}_RELWITHDEBINFO}
       CACHE STRING "Fortran compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE ${FFLAGS_${FCID}_RELEASE}
       CACHE STRING "Fortran compiler flags used for release builds." FORCE)
  set (CMAKE_Fortran_FLAGS_FAST ${FFLAGS_${FCID}_FAST}
       CACHE STRING "Fortran compiler flags used for fast release builds." FORCE)
  set (CMAKE_Fortran_FLAGS_USER "TRUE"
       CACHE INTERNAL "set to FALSE to reset to default Fortran flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced (CLEAR
    CMAKE_Fortran_FLAGS_DEFAULT
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_GARBLE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_FAST
  )
else ()
  mark_as_advanced (FORCE
    CMAKE_Fortran_FLAGS_DEFAULT
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_GARBLE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_FAST
  )
endif ()
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEFAULT}")
