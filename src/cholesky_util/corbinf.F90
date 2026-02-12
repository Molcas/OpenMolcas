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
!----------------------------------------------------------------------*
! Global Variables describing chemical system                          *
!----------------------------------------------------------------------*
! nSym    - number of symmetries                                       *
! nOrb(i) - (i = 1, nSym), number of orbitals                          *
! nOcc(i) - (i = 1, nSym), number of occupied orbitals                 *
! nFro(i) - (i = 1, nSym), number of frozen orbitals                   *
! nDel(i) - (i = 1, nSym), number of orbitals deleted by linear deps.  *
! nExt(i) - (i = 1, nSym), number of virtual (external) orbitals       *
!----------------------------------------------------------------------*

module COrbInf

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: nDel(8), nExt(8), nFro(8), nOcc(8), nOrb(8), nSym

public :: nDel, nExt, nFro, nOcc, nOrb, nSym

end module COrbInf
