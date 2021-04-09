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
        subroutine MkQ4 (V)
c
c        Q4(a,b,c,d) <- V(ab,cd)
c
        implicit none
#include "chcc1.fh"
        real*8 V(1:nv*(nv+1)/2,1:nv*(nv+1)/2)
c
c        help variables
        integer a,b,c,d,ab,cd
c
        cd=0
        do c=1,nv
        do d=1,c
        cd=cd+1
          ab=0
          do a=1,nv
          do b=1,a
          ab=ab+1
            Q4(a,b,c,d)=V(ab,cd)
            Q4(a,b,d,c)=V(ab,cd)
            Q4(b,a,c,d)=V(ab,cd)
            Q4(b,a,d,c)=V(ab,cd)
          end do
          end do
        end do
        end do
c
        return
        end
