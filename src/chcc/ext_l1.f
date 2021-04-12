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
        subroutine  Ext_L1 (V1,V2,
     c                      no,dima,dimc,adda,nbs)
c
c       this routine do:
c       V2(i,a',m') <- V1(p,q,m')
c
        implicit none
        integer no,dima,dimc,adda,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:no,1:dima,1:dimc)
c
c       help variables
        integer i,a,m
c
        do m=1,dimc
        do a=1,dima
        do i=1,no
        V2(i,a,m)=V1(i,adda+a,m)
        end do
        end do
        end do
c
        return
        end
