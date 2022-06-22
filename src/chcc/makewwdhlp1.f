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
        subroutine MakeWwdHlp1 (Ww,W1,
     c                          dima,dimbe,dimbega)
c
c       this routine do:
c       Make  Ww(+)((aa)",(bega)") from W1(a",be",b",ga")
c       for the case beSGrp=gaSGrp
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",ga",bega" (I)
c
        implicit none
        integer dima,dimbe,dimbega
        real*8 Ww(1:dima,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dima,1:dimbe)
c
c       help variables
        integer a,be,ga,bega
c
c
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          do a=1,dima
            Ww(a,bega)=W1(a,be,a,ga)
          end do
          end do
          end do
c
        return
        end
