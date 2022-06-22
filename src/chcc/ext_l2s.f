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
        subroutine Ext_L2s (V1,V2,
     c                      dima,dimab,dimc,adda,addb,nbs)
c
c       this routine do:
c       V2(a'b',m') <- V1(p,q,m') for aGrp=bGrp
c
        implicit none
        integer dima,dimab,dimc,adda,addb,nbs
        real*8 V1(1:nbs,1:nbs,1:dimc)
        real*8 V2(1:dimab,1:dimc)
c
c       help variables
        integer a,b,ab,m
c
        do m=1,dimc
        ab=0
        do a=1,dima
        do b=1,a
        ab=ab+1
        V2(ab,m)=V1(adda+a,addb+b,m)
        end do
        end do
        end do
c
        return
        end
