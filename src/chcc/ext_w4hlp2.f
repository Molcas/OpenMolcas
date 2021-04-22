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
        subroutine Ext_W4hlp2 (V2,M1,
     c                       nc,dima,dimab,dimapp,dimbpp,addapp,addbpp)
c
c        this routine do:
c        Extract M1(m,a",b") <- V2(m,a'b')
c
        implicit none
        integer nc,dima,dimab,dimapp,dimbpp,addapp,addbpp
        real*8 V2(1:nc,1:dimab)
        real*8 M1(1:nc,1:dimapp,1:dimbpp)
c
c        help variables
c
        integer a,ab,app,bpp,m
c
        do app=1,dimapp
        a=addapp+app
        ab=a*(a-1)/2+addbpp
        do bpp=1,dimbpp
        ab=ab+1
c
          do m=1,nc
          M1(m,app,bpp)=V2(m,ab)
          end do
c
        end do
        end do
c
        return
c Avoid unused argument warnings
        if (.false.) Call Unused_integer(dima)
        end
