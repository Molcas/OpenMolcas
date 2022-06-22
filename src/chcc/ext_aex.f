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
        subroutine Ext_Aex (Aex,VV,no)
c
c        this routine def:
c        VV(i,u,v,j) <- Aex(ij,u,v)
c        for Aex(i,j,u,v) = Aex(j,i,v,u), Aex stored only for i>=j
c
        implicit none
        integer no
        real*8 VV(1:no,1:no,1:no,1:no)
        real*8 Aex(1:no*(no+1)/2,1:no,1:no)
c
c        help variables
        integer i,j,ij,u,v
c
c
        do u=1,no
        do v=1,no
          ij=0
          do i=1,no
          do j=1,i
          ij=ij+1
c
            VV(i,u,v,j)=Aex(ij,u,v)
            VV(j,v,u,i)=Aex(ij,u,v)
c
          end do
          end do
        end do
        end do
c
c
        return
        end
