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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Tra2A(ij_pair,ij_Bas_pairs,kl_Orb_pairs,kSym,lSym,kBas,lBas,kAsh,lAsh,CMO_k,CMO_l,IJKL,IJKX,IJVX,VXIJ)
!***********************************************************************
!                                                                      *
!     run the first two quarter AO --> MO transformations with         *
!     both transformed indices being active, i.e.                      *
!     (ij!kl) --> (ij!vx)                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ij_pair, ij_Bas_pairs, kl_Orb_pairs, kSym, lSym, kBas, lBas, kAsh, lAsh
real(kind=wp), intent(in) :: CMO_k(kBas,kAsh), CMO_l(lBas,lAsh), IJKL(lBas*kBas)
real(kind=wp), intent(out) :: IJKX(kBas*lAsh), IJVX(lAsh*kAsh)
real(kind=wp), intent(inout) :: VXIJ(ij_bas_pairs,kl_Orb_pairs)
integer(kind=iwp) :: kl

! (ij|kl) -> (ij|kx)
call DGEMM_('T','N',kBas,lAsh,lBas,One,IJKL,lBas,CMO_l,lBas,Zero,IJKX,kBas)

! (ij|kx) -> (ij|vx)
if (kSym == lSym) then
  ! Triangular storage of target
  call DGEMM_Tri('T','N',kAsh,kAsh,kBas,One,IJKX,kBas,CMO_k,kBas,Zero,IJVX,kAsh)
else
  ! Rectangular storage of target
  call DGEMM_('T','N',lAsh,kAsh,kBas,One,IJKX,kBas,CMO_k,kBas,Zero,IJVX,lAsh)
end if

! Scatter result onto half-transformed integral list
do kl=1,kl_Orb_pairs
  VXIJ(ij_pair,kl) = IJVX(kl)
end do

return

end subroutine Tra2A
