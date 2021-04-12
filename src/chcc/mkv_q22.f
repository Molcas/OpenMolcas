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
        subroutine MkV_Q22 (W2,W1,dima)
c
c       this routine do
c       W1(j,u,i,a') = 2W2(i,u,j,a')-W2(j,u,i,a')
c
c       N.B. Kvajt odflaknute
c
        implicit none
#include "chcc1.fh"
        integer dima
        real*8 W1(1:no,1:no,1:no,1:dima)
        real*8 W2(1:no,1:no,1:no,1:dima)
c
c       help variables
        integer i,j,u,a
c
c
        do a=1,dima
        do i=1,no
        do u=1,no
        do j=1,no
          W1(j,u,i,a)=W2(j,u,i,a)-2.0d0*W2(i,u,j,a)
        end do
        end do
        end do
        end do
c
c
        return
        end
