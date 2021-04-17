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
        subroutine ExH_X2 (Gvv,H,dima,dimbe,nv,adda,addbe)
c
c        this routine do:
c        H(a',be') <- Gvv(be,a)
c
        implicit none
        integer dima,dimbe,nv,adda,addbe
        real*8 H(1:dima,1:dimbe)
        real*8 Gvv(1:nv,1:nv)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
        bev=addbe+be
          do a=1,dima
            H(a,be)=Gvv(bev,adda+a)
          end do
        end do
c
        return
        end
