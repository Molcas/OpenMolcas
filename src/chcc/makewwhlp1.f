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
        subroutine MakeWwHlp1 (Ww,W1,
     c                       dima,dimb,dimab,dimbe,dimga,dimbega,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c       for the case a"=b" , be"=ga"
c       N.B. algoritmus nieje prilis vymakany
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       dimx   - dimension of a",b",ab",be",ga",bega" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimab,dimbe,dimga,dimbega,key
        real*8 Ww(1:dimab,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
c
c       help variables
        integer a,b,ab,be,ga,bega
c
c
        if (key.eq.1) then
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,bega)=W1(a,be,b,ga)+W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          bega=0
          do be=2,dimbe
          do ga=1,be-1
          bega=bega+1
          ab=0
          do a=2,dima
          do b=1,a-1
          ab=ab+1
            Ww(ab,bega)=W1(a,be,b,ga)-W1(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dimab*dimbega,dimab*dimbega,
c    c              Ww(1,1),0.5d0)
c
        return
        end
