!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine Orthox(S,C,nOrb,nBas)
!***********************************************************************
!                                                                      *
!     purpose: Perform Gram-Schmidt orthogonalization scheme           *
!                                                                      *
!     input:                                                           *
!       S       : overlap in non-orthonormal basis (square storage)    *
!       C       : matrix transforming non-orthonormal basis            *
!      nBas,nOrb: dimensions                                           *
!                                                                      *
!     output:                                                          *
!       C       : matrix transforming orthonormal basis                *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nOrb, nBas
real(kind=wp) :: S(nOrb,nOrb), C(nBas,nOrb)
integer(kind=iwp) :: iBas, iOrb, jOrb, kOrb
real(kind=wp) :: A, F

do iOrb=1,nOrb
  F = Zero
  if (S(iOrb,iOrb) > Zero) F = One/sqrt(S(iOrb,iOrb))
  do iBas=1,nBas
    C(iBas,iOrb) = F*C(iBas,iOrb)
  end do
  do jOrb=1,nOrb
    S(iOrb,jOrb) = F*S(iOrb,jOrb)
    S(jOrb,iOrb) = F*S(jOrb,iOrb)
  end do
  S(iOrb,iOrb) = One
  do jOrb=iOrb+1,nOrb
    A = S(iOrb,jOrb)
    do iBas=1,nBas
      C(iBas,jOrb) = C(iBas,jOrb)-A*C(iBas,iOrb)
    end do
    do kOrb=1,nOrb
      S(jOrb,kOrb) = S(jOrb,kOrb)-A*S(iOrb,kOrb)
    end do
    do kOrb=1,nOrb
      S(kOrb,jOrb) = S(kOrb,jOrb)-A*S(kOrb,iOrb)
    end do
  end do
end do

end subroutine Orthox
