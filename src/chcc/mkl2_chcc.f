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
        subroutine MkL2_chcc (V)
c
c        L2(m,a,b) <- V(m,ab)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nc,1:nv*(nv+1)/2)
c
c        help var
        integer a,b,ab,m
c
        ab=0
        do a=1,nv
        do b=1,a
        ab=ab+1
          do m=1,nc
             L2k(m,a,b)=V(m,ab)
             L2k(m,b,a)=V(m,ab)
          end do
        end do
        end do
c
        return
        end
