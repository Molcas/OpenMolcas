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
        subroutine ExH_T13 (V,Hvv,dimbe,addbe,nv)
c
c        this routine do:
c       Extract V1(a,be') <- Hvv(a,be)
c        Hvv(a,be') <- Fvv(a,be)
c
        implicit none
        integer dimbe,nv,addbe
        real*8 Hvv(1:nv,1:nv)
        real*8 V(1:nv,1:dimbe)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
          bev=be+addbe
          do a=1,nv
            V(a,be)=Hvv(a,bev)
          end do
        end do
c
c
        return
        end
