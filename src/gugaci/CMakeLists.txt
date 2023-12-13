#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

set (sources
     main.F90
     actextloop.F90
     actploop.F90
     actsegment.F90
     arrangeorb.F90
     calh0.F90
     cidiagelm.F90
     cidiagelm_g.F90
     cigrad.F90
     cipro.F90
     complete_ext_loop.F90
     davdiag.F90
     dblcloop.F90
     dblcloop_g.F90
     dblextloop.F90
     dblploop1.F90
     dblploop1_g.F90
     dblploop2.F90
     dblploop3.F90
     dblploop4.F90
     dblploop5.F90
     dblplpval.F90
     dbl_upwalk.F90
     ddsstt.F90
     denmat.F90
     drtrel.F90
     dvsdtd.F90
     extloop1.F90
     extploop.F90
     extploop_g.F90
     gpcsf.F90
     gugaci.F90
     gugaci_global.F90
     gugasyscom.F90
     h0diagonal.F90
     innloop.F90
     innloop_g.F90
     intsort.F90
     molinf.F90
     mrcimath.F90
     natureorb.F90
     paras_calculate.F90
     prod.F90
     prod_g.F90
     prod_util.F90
     readint.F90
     segdrt.F90
     segvalue.F90
     sysdep.F90
     tvsv.F90
     vddsdt.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/prog_template.cmake)

# Workaround for gcc spurious warning about large variables,
# versions are approximate
# (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98411)
# This could also be removed if large arrays are made allocatable
if (BIGOT AND NOT OPENMP AND
    CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "10.4")
  target_compile_options (${prog}_obj PRIVATE "-fmax-stack-var-size=8388608")
endif ()
