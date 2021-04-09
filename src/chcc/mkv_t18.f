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
        subroutine MkV_T18 (Va,V,dima,no)
c
c       this routine do:
c        Va(a,j,i,u) = - [2V(ai|ju)-V(aj|iu)]
c        from V(o_a,P,Q,a)
c
c        N.B. Kvajto odflaknute, ozaj ze hnus, treba sa zamysliet
c        ci je nutva permutacia s a-ckom na konci - bod QK2.2 @@
c
c
        implicit none
        integer dima,no
        real*8 Va(1:dima,1:no,1:no,1:no)
        real*8 V(1:no,1:no,1:no,1:dima)
c
c       help variables
        integer a,j,i,u
c
        do u=1,no
          do i=1,no
            do j=1,no
              do a=1,dima
                Va(a,j,i,u)=V(j,i,u,a)-2.0d0*V(i,j,u,a)
              end do
            end do
          end do
        end do
c
c
        return
        end
