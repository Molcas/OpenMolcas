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
        subroutine MkT1T2
c
c        T1(a,i) = 0 (mozno neskor ine)
c        T2(a,b,i,j) = (ai|bj)/Dabij
c
        implicit none
#include "chcc1.fh"
c
c        help variables
        integer a,b,i,j
c
        do i=1,no
        do a=1,nv
          T1c(a,i)=0.0d0
        end do
        end do
c
        do j=1,no
        do i=1,no
          do b=1,nv
          do a=1,nv
            T2c(a,b,i,j)=Q21(a,i,b,j)/(OEo(i)+OEo(j)-OEv(a)-OEv(b))
          end do
          end do
        end do
        end do
c
        return
        end
