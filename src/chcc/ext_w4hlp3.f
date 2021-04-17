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
        subroutine Ext_W4hlp3 (V2,M1,
     c                     nc,dima,dimb,dimapp,dimbpp,addapp,addbpp)
c
c        this routine do:
c        Extract M1(m,a",b") <- V2(m,a',b')
c
        implicit none
        integer nc,dima,dimb,dimapp,dimbpp,addapp,addbpp
        real*8 V2(1:nc,1:dima,1:dimb)
        real*8 M1(1:nc,1:dimapp,1:dimbpp)
c
c        help variables
c
        integer a,b,app,bpp,m
c
        do app=1,dimapp
        a=addapp+app
        do bpp=1,dimbpp
        b=addbpp+bpp
c
          do m=1,nc
          M1(m,app,bpp)=V2(m,a,b)
          end do
c
        end do
        end do
c
        return
        end
