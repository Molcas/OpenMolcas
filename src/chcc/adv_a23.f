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
        subroutine AdV_A23 (V1,A,dimij,no)
c
c       this routine do:
c       A(ij,u,v) <<- V1(j,iu,v) + V1(i,jv,u)
c
c        Velmi odflaknute, da sa to urobit podstatne lepsie, ale
c        o4 proces osrat fok
c
        implicit none
        integer dimij,no
        real*8 A(1:dimij,1:no,1:no)
        real*8 V1(1:no,1:dimij,1:no)
c
c       help variables
        integer i,j,ij,u,v,iu,jv
c
        do v=1,no
        do u=1,no
c
          ij=0
c
          do i=1,no
          if (i.ge.u) then
          iu=i*(i-1)/2+u
          else
          iu=u*(u-1)/2+i
          end if
c
          do j=1,i
          ij=ij+1
          if (j.ge.v) then
          jv=j*(j-1)/2+v
          else
          jv=v*(v-1)/2+j
          end if
c
            A(ij,u,v)=A(ij,u,v)+V1(j,iu,v)+V1(i,jv,u)
c
          end do
          end do
c
        end do
        end do
c
        return
        end
