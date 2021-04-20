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
        subroutine ExV_X43 (Vp,V,dimab,no)
c
c       this routine do:
c       Vp(a_b,ij) <- V(a_b,j,i) for i>=j
c
        implicit none
        integer dimab,no
        real*8 Vp(1:dimab,1:no*(no+1)/2)
        real*8 V(1:dimab,1:no,1:no)
c
c       help variables
        integer i,j,ij,ab
c
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
          do ab=1,dimab
            Vp(ab,ij)=V(ab,j,i)
          end do
        end do
        end do
c
c
        return
        end
