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

module RPA_globals

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nTitle, nSym, nFreeze(2), nBas(8), nOrb(8), nFro(8,2), nOcc(8,2), nVir(8,2), nDel(8,2), l_CMO, l_EMO, &
                     l_OccEn(2), l_VirEn(2), iPrint
real(kind=wp) :: NuclearRepulsionEnergy(1)
logical(kind=iwp) :: dRPA, SOSEX, doCD, doDF, doLDF, LumOrb
character(len=3) :: Reference
character(len=8) :: RPAModel
character(len=80) :: DFTFunctional
real(kind=wp), allocatable :: CMO(:,:), EMO(:,:), OccEn(:,:), VirEn(:,:)
integer(kind=iwp), parameter :: mTitle = 10
character(len=80) :: Title(mTitle)

public :: nTitle, nSym, nFreeze, nBas, nORb, nFro, nOcc, nVir, nDel, CMO, l_CMO, EMO, l_EMO, OccEn, l_OccEn, VirEn, l_VirEn, &
          iPrint, NuclearRepulsionEnergy, dRPA, SOSEX, doCD, doDF, doLDF, LumOrb, Reference, RPAModel, DFTFunctional, mTitle, Title

end module RPA_globals
