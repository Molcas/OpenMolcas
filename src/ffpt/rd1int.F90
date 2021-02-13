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

subroutine Rd1Int_FFPT()
!***********************************************************************
!                                                                      *
!     Objective: Read the header of the one-electron integral file     *
!                Extract also symmetry and basis set information.      *
!                In addition read the overlap, the nuclear attraction  *
!                and kinetic integrals.                                *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nAtoms, nBas, nSym, Coor, Header
use stdalloc, only: mma_allocate

implicit none

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Read one-electron integral file header etc.                      *
!                                                                      *
!----------------------------------------------------------------------*

call Get_cArray('Seward Title',Header,144)
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Coor,3,nAtoms,label='Coor')
call Get_dArray('Unique Coordinates',Coor,3*nAtoms)

return

end subroutine Rd1Int_FFPT
