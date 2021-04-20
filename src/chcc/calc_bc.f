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
        subroutine Calc_Bc
c
c        this routine calc Bc
c        Bc(a,b,be,ga) =        (a,be|b,ga)
c                     - S(i)   (a,be,|b,i).t1(ga,i)
c                     - S(i)   (a,i,|b,ga).t1(be,i)
c                     + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
c
        implicit none
#include "chcc1.fh"
c
c        help var
        integer i,a,b,be,ga
c       integer j
        real*8 s
c
        do ga=1,nv
        do be=1,nv
        do b=1,nv
        do a=1,nv
c
c
c1            (a,be|b,ga)
          s=Q4(be,a,ga,b)
c
c2,3          - S(i)   (a,be,|b,i).t1(ga,i)
c         - S(i)   (a,i,|b,ga).t1(be,i)
          do i=1,no
           s=s-Q3(a,be,b,i)*T1c(ga,i)
           s=s-Q3(b,ga,a,i)*T1c(be,i)
          end do

c4        + S(i,j) (a,i|b,j).t1(be,i).t1(ga,j)
c          do j=1,no
c          do i=1,no
c          s=s+Q21(a,i,b,j)*T1c(be,i)*T1c(ga,j)
c          end do
c          end do
c
c
          Bc(a,b,be,ga)=s
c
        end do
        end do
        end do
        end do
c
c
        return
        end
