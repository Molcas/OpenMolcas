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
        subroutine ExA_X4 (A,Ap,no)
c
c       this routine do:
c       Ap(i,u,v) <- A(ii,u,v)
c
        implicit none
        integer no
        real*8 Ap(1:no,1:no*no)
        real*8 A(1:no*(no+1)/2,1:no*no)
c
c       help variables
        integer uv,i,ii
c
        do uv=1,no*no
        do i=1,no
          ii=i*(i+1)/2
            Ap(i,uv)=A(ii,uv)
          end do
        end do
c
c
        return
        end
