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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine RotateNxN(CMO,kappa,nOrb2Loc,nBasis,kappa_cnt,xkappa_cnt,unitary_mat,rotated_CMO)
! this subroutine rotates the orbitals as CMO = exp(-kappa) * CMO
!
! for technical reasons, auxiliary matrices (kappa_cnt, xkappa_cnt, unitary_mat, rotated_CMO) are allocated outside of the loop
! that calls this routine. However it would work perfectly fine if these are allocated and deallocated within this routine to reduce
! the number of arguments (not recommended for speed)

use definitions, only: wp,iwp
use constants, only: Zero,One

implicit none

integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(inout) :: kappa(nOrb2Loc,nOrb2Loc),kappa_cnt(nOrb2Loc,nOrb2Loc),xkappa_cnt(nOrb2Loc,nOrb2Loc),&
                             unitary_mat(nOrb2Loc,nOrb2Loc), rotated_CMO(nBasis,nOrb2Loc)
integer(kind=iwp) :: i,k, iBas

! get exp(-kappa)
call expkap_localisation(kappa,nOrb2Loc,kappa_cnt,xkappa_cnt,unitary_mat)


! transform the orbitals
rotated_CMO(:,:) = Zero
call dgemm_('N','N',nBasis,nOrb2Loc,nOrb2Loc,One,CMO,nBasis,unitary_mat,nOrb2Loc,Zero,rotated_CMO,nBasis)

!reset CMO to be updated
CMO(:,:) = rotated_CMO(:,:)

end subroutine RotateNxN

