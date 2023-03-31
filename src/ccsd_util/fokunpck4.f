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
       subroutine fokunpck4 (fok,fii,dimfok,dimfi)
!
!     this routine distribute (Fok - dp) to Fii
!     fok    - Fok matrix (I)
!     fii    - Fii matrix (O)
!     dimfok - dimension for Fok matrix - norb (I)
!     dimfi  - dimension of occupied - no (I)
!
       integer dimfok,dimfi
!
       real*8 fok(1:dimfok,1:dimfok)
       real*8 fii(1:dimfi,1:dimfi)
!
!     help variables
!
       integer i,j
!
!1    distribute Fok to Fii
       do 400 j=1,dimfi
       do 401 i=1,dimfi
       fii(i,j)=fok(i,j)
 401    continue
 400    continue
!
       return
       end
