!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine ExtT1 (H,T1,dima,adda)
!
!       this routine do:
!       Extract H(a',i) <- T1(a,i) for given aGrp
!
!       parameter description:
!       H       - Output file (O)
!       T1      - T1 amplitudes (I)
!       dima    - dimension of given Group (I)
!       adda    - shift of a' in full a set (I)
!
!       N.B. Kvajt odflaknute
!
        implicit none
#include "chcc1.fh"
        integer dima,adda
        real*8 T1(1:nv,1:no)
        real*8 H(1:dima,1:no)
!
!       help variables
        integer a,i
!
        do i=1,no
          do a=1,dima
          H(a,i)=T1(adda+a,i)
          end do
        end do
!
!
        return
        end
