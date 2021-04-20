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
        subroutine MkQ1 (V)
c
c        Q1(a,j,k,l) <- V(aj,kl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv,1:no,1:no*(no+1)/2)
c
c        help variables
        integer a,j,k,l,kl
c
        kl=0
        do k=1,no
        do l=1,k
        kl=kl+1
          do j=1,no
          do a=1,nv
            Q1(a,j,k,l)=V(a,j,kl)
            Q1(a,j,l,k)=V(a,j,kl)
          end do
          end do
        end do
        end do
c
        return
        end
