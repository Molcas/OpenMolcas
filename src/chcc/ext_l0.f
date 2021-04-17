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
        subroutine  Ext_L0 (V1,V2,
     c                      no,dimij,dimc,nbs)
c
c       this routine do:
c       V2(i,j,m') <- V1(p,q,m')
c
        implicit none
        integer no,dimij,dimc,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:dimij,1:dimc)
c
c       help variables
        integer i,j,ij,m
c
        do m=1,dimc
        ij=0
        do i=1,no
        do j=1,i
        ij=ij+1
        V2(ij,m)=V1(i,j,m)
        end do
        end do
        end do
c
        return
        end
