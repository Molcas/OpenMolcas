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
        subroutine MkT_exp (T2red,T2full,dimbe,no)
c
c        this routine produce expanded set of amplitudes from reduced set
c        (used in parallel case)
c       T2red(be'ga',u,v) -> T2full(be',ga',u,v)
c        for beGrp=gaGrp
c
        implicit none
        integer dimbe,no
        real*8 T2full(1:dimbe,1:dimbe,1:no,1:no)
        real*8 T2red(1:dimbe*(dimbe+1)/2,1:no,1:no)
c
c        help variables
        integer i,j,be,ga,bega
c
        do j=1,no
        do i=1,no
c
        bega=0
        do be=1,dimbe
        do ga=1,be
        bega=bega+1
c
cdir          T2red(bega,i,j)=T2Full(be,ga,i,j)
cinv
          T2Full(be,ga,i,j)=T2red(bega,i,j)
          T2Full(ga,be,j,i)=T2red(bega,i,j)
c
        end do
        end do
c
        end do
        end do
c
        return
        end
