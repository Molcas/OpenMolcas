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

implicit none
real*8 Dens(*), CMO(*)
integer nBas, nOcc
integer nTBs

nTBs = max(nBas,1)
call DGEMM_('N','T',nBas,nBas,nOcc,1.0d0,CMO,nTBs,CMO,nTBs,0.0d0,Dens,nTBs)

end subroutine GetDens_Localisation
