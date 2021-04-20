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
        subroutine MakeWwdHlp2 (Ww,W1,
     c                       dima,dimbe,dimga)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c       for the case a"=b" , beSGrp.ne.gaSGrp
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",be",ga" (I)
c
        implicit none
        integer dima,dimbe,dimga
        real*8 Ww(1:dima,1:dimbe,1:dimga)
        real*8 W1(1:dima,1:dimbe,1:dima,1:dimga)
c
c       help variables
        integer a,be,ga
c
c
CVpV 2014 Fix for Intel Compiler v14.*
#ifdef __INTEL_COMPILER
          do a=1,dima
          do ga=1,dimga
cDEC$ VECTOR UNALIGNED
          do be=1,dimbe
            Ww(a,be,ga)=W1(a,be,a,ga)
          end do
          end do
          end do
#else
          do ga=1,dimga
          do be=1,dimbe
          do a=1,dima
            Ww(a,be,ga)=W1(a,be,a,ga)
          end do
          end do
          end do
#endif
c
        return
        end
