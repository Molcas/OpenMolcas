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
        subroutine  Mk_T1t (T1,H,dimbepp,no,nv,addbepp)
c
c        this routine do:
c        H(i,be") <- T1o(be,i)
c
        implicit none
        integer dimbepp,no,nv,addbepp
        real*8 T1(1:nv,1:no)
        real*8 H(1:no,1:dimbepp)
c
c        help variables
        integer i,bepp
c
        do i=1,no
        do bepp=1,dimbepp
          H(i,bepp)=T1(addbepp+bepp,i)
        end do
        end do
c
        return
        end
