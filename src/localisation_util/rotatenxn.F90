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

subroutine RotateNxN(CMO,nOrb2Loc,nBasis,unitary_mat,rotated_CMO)
! this subroutine rotates the orbitals as rotated_CMO = CMO * exp(-kappa)

use definitions, only: wp,iwp
use constants, only: Zero,One

implicit none

integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(in) :: CMO(nBasis,nOrb2Loc)
real(kind=wp), intent(in) :: unitary_mat(nOrb2Loc,nOrb2Loc)
real(kind=wp), intent(out) :: rotated_CMO(nBasis,nOrb2Loc)
!Call recprt('CMO','(5F20.10) ',CMO,nBasis,nOrb2loc)
!Call recprt('Unitary_mat(nxn)','(5F20.10) ',unitary_mat,nOrb2Loc,nOrb2loc)


! transform the orbitals
rotated_CMO(:,:) = Zero
call dgemm_('N','N',nBasis,nOrb2Loc,nOrb2Loc,&
                    One,CMO,nBasis,&
                        unitary_mat,nOrb2Loc,&
                    Zero,rotated_CMO,nBasis)
!Call recprt('Rotated_CMO','(5F20.10) ',Rotated_CMO,nBasis,nOrb2loc)

end subroutine RotateNxN

