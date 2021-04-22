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
        subroutine AdT_T17 (T1,Q,dimat,dimaq,no,adda,f)
c
c       this routine do:
c       T1(a,i) <<- f . Q(a',i)
c
        implicit none
        integer dimat,dimaq,no,adda
        real*8 T1(1:dimat,1:no)
        real*8 Q(1:dimaq,no)
        real*8 f
c
c       help variables
        integer i,a
c
        do i=1,no
          do a=1,dimaq
            t1(adda+a,i)=t1(adda+a,i)+f*q(a,i)
          end do
        end do
c
        return
        end
