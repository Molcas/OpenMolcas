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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine GetDens_Localisation(Dens,CMO,nBas,nOcc)
! Author: T.B. Pedersen
!
! Purpose: compute density from CMOs as Dens = CMO * (CMO)^T

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBas, nOcc
real(kind=wp), intent(out) :: Dens(nBas,nBas)
real(kind=wp), intent(in) :: CMO(nBas,nOcc)
integer(kind=iwp) :: nTBs

nTBs = max(nBas,1)
call DGEMM_('N','T',nBas,nBas,nOcc,One,CMO,nTBs,CMO,nTBs,Zero,Dens,nTBs)

end subroutine GetDens_Localisation
