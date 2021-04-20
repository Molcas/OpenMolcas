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
        subroutine MakeWwHlp3 (Ww,W1,W2,
     c                       dima,dimb,dimbe,dimga,dimbega,key)
c
c       this routine do:
c       Make  Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
c                                     and W2(b",be",a",ga")
c       for the case a".ne.b" , be"=ga"
c       N.B. algoritmus nieje prilis vymakany, prva cast sa da
c            urobit pomodou dcopy
c
c       parameter description
c       Ww     - array for Ww+(-) (O)
c       W1     - array for W1(a",be",b",ga") (I)
c       W2     - array for W1(b",be",a",ga") (I)
c       dimx   - dimension of a",b",be",ga",be>(>=)ga" (I)
c       key    - 1 - calc Ww+, 2 - calc Ww- (I)
c
        implicit none
        integer dima,dimb,dimbe,dimga,dimbega,key
        real*8 Ww(1:dima,1:dimb,1:dimbega)
        real*8 W1(1:dima,1:dimbe,1:dimb,1:dimga)
        real*8 W2(1:dimb,1:dimbe,1:dima,1:dimga)
c
c       help variables
        integer a,b,be,ga,bega
c
c
        if (key.eq.1) then
          bega=0
          do be=1,dimbe
          do ga=1,be
          bega=bega+1
          do b=1,dimb
          do a=1,dima
            Ww(a,b,bega)=W1(a,be,b,ga)+W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        else
          bega=0
          do be=2,dimbe
          do ga=1,be-1
          bega=bega+1
          do b=1,dimb
          do a=1,dima
            Ww(a,b,bega)=W1(a,be,b,ga)-W2(b,be,a,ga)
          end do
          end do
          end do
          end do
        end if
c
c        Cely clen ma Faktor 2, tu teda nevydelim 2
c        call mv0sv (dima*dimb*dimbega,dima*dimb*dimbega,
c    c              Ww(1,1,1),0.5d0)
c
        return
        end
