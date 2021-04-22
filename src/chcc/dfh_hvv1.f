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
        subroutine DfH_Hvv1 (Hvv,Fvv,nv,dimbe,addbe)
c
c        this routine do:
c        Hvv(a,be') <- Fvv(a,be)
c
        implicit none
        integer dimbe,nv,addbe
        real*8 Hvv(1:nv,1:dimbe)
        real*8 Fvv(1:nv,1:nv)
c
c        help variables
        integer a,be,bev
c
        do be=1,dimbe
          bev=be+addbe
          do a=1,nv
            Hvv(a,be)=Fvv(a,bev)
          end do
        end do
c
c
        return
        end
