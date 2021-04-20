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
        subroutine MkQ3 (V)
c
c        Q1(a,b,c,l) <- V(ab,cl)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:no,1:nv)
c
c        help variables
        integer a,b,ab,c,l
c
        do l=1,no
        do c=1,nv
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q3(a,b,c,l)=V(ab,l,c)
            Q3(b,a,c,l)=V(ab,l,c)
          end do
          end do
        end do
        end do
c
        return
        end
