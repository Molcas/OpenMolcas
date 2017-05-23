************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine xbacktrace
CSVC: this routine tries to print a backtrace
      implicit none
C     use backtrace intrinsic (introduced since gfortran 4.8)
#if   defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8))
      call backtrace
C     use tracebackqq intrinsic for ifort
#elif defined(__INTEL_COMPILER)
      call tracebackqq
#endif
      end
