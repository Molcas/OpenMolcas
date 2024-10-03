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
integer(kind=iwp), intent(in) :: nOrb, nBas
real(kind=wp), intent(inout) :: S(nOrb,nOrb), C(nBas,nOrb)
integer(kind=iwp) :: iOrb, jOrb
real(kind=wp) :: A, F

do iOrb=1,nOrb
  F = Zero
  if (S(iOrb,iOrb) > Zero) F = One/sqrt(S(iOrb,iOrb))
  C(:,iOrb) = F*C(:,iOrb)
  S(iOrb,:) = F*S(iOrb,:)
  S(:,iOrb) = F*S(:,iOrb)
  S(iOrb,iOrb) = One
  do jOrb=iOrb+1,nOrb
    A = S(iOrb,jOrb)
    C(:,jOrb) = C(:,jOrb)-A*C(:,iOrb)
    S(jOrb,:) = S(jOrb,:)-A*S(iOrb,:)
    S(:,jOrb) = S(:,jOrb)-A*S(:,iOrb)
  end do
end do

end subroutine Orthox
