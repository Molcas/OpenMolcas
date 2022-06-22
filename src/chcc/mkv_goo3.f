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
        subroutine MkV_Goo3 (V,V2,dima,no)
c
c       this routine do:
c       Make AntiSymetric integrals
c        V2(a',j,i,u) <- 2 V(a,j|iu) - (a,i|ju)
c
        implicit none
        integer dima,no
        real*8 V2(1:dima,1:no,1:no,1:no)
        real*8 V(1:dima,1:no,1:no*(no+1)/2)
c
c       help variables
        integer a,i,j,u,iu,ju
c
c
        do u=1,no
c
          do i=1,no
            if (i.gt.u) then
            iu=(i-1)*i/2+u
            else
            iu=(u-1)*u/2+i
            end if
c
            do j=1,no
c
              if (j.gt.u) then
              ju=(j-1)*j/2+u
              else
              ju=(u-1)*u/2+j
              end if
c
              do a=1,dima
c
                V2(a,j,i,u)=2.0d0*V(a,j,iu)-V(a,i,ju)
c
              end do
            end do
          end do
        end do
c
c
        return
        end
