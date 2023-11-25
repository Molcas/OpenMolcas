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

subroutine RPA_SetInc()

use RPA_globals, only: DFTFunctional, doCD, doDF, dRPA, iPrint, l_CMO, l_EMO, l_OccEn, l_VirEn, LumOrb, mTitle, nBas, nDel, &
                       nFreeze, nFro, nOcc, nOrb, nSym, nTitle, NuclearRepulsionEnergy, nVir, Reference, RPAModel, SOSEX, Title
use Constants, only: Zero
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i

! rpa_config
Reference = 'Non'
RPAModel = 'None@Non'
DFTFunctional = 'Not defined     '
dRPA = .false.
SOSEX = .false.
doCD = .false.
doDF = .false.
LumOrb = .false.
iPrint = 0
! rpa_data
do i=1,mTitle
  write(Title(i),'(A)') repeat(' ',80)
end do
nTitle = 0
nSym = 0
nFreeze(:) = 0
nBas(:) = 0
nOrb(:) = 0
nFro(:,:) = 0
nDel(:,:) = 0
nOcc(:,:) = 0
nVir(:,:) = 0
l_CMO = 0
l_EMO = 0
l_OccEn(:) = 0
l_VirEn(:) = 0
NuclearRepulsionEnergy(1) = Zero

end subroutine RPA_SetInc
