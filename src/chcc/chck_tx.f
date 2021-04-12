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
        subroutine Chck_Tx (T)
c
c        check T(a,b,i,j)
c
        implicit none
#include "chcc1.fh"
c        real*8 T(1:nv,1:nv,1:no,1:no)
         real*8 T(1:nv,1:no,1:nv,1:no)
c        real*8 T(1:nv,1:nv,1:no)
c
        integer b,j,a,i,bad
        real*8 s
c
        bad=0
        do j=1,no
        do i=1,no
        do b=1,nv
        do a=1,nv
c
           s=T2c(a,b,i,j)
c
          if (abs(T(b,i,a,j)-s).gt.1.0d-10) then
            bad=bad+1
          end if
c
        end do
        end do
        end do
        end do
c
        write (6,*) ' Chck T2 :',bad
c
        return
        end
