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
        subroutine MkQ22 (V)
c
c        Q22(a,b,k,l) <- V(ij,kl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:no*(no+1)/2)
c
c        help variables
        integer a,b,k,l,ab,kl
c
        kl=0
        do k=1,no
        do l=1,k
        kl=kl+1
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q22(a,b,k,l)=V(ab,kl)
            Q22(a,b,l,k)=V(ab,kl)
            Q22(b,a,k,l)=V(ab,kl)
            Q22(b,a,l,k)=V(ab,kl)
          end do
          end do
        end do
        end do
c
        return
        end
