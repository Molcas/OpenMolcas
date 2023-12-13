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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine Check_Amp2(nSym,nOcc,nVir,iSkip)

use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nOcc(nSym), nVir(nSym)
integer(kind=iwp), intent(out) :: iSkip
integer(kind=iwp) :: iSym, iSyma, iSymi, nT1amTot, nT1am(8)

iSkip = 0
nT1amTot = 0
do iSym=1,nSym
  nT1am(iSym) = 0
  do iSymi=1,nSym
    iSyma = Mul(iSymi,iSym)
    nT1am(iSym) = nT1am(iSym)+nVir(iSyma)*nOcc(iSymi)
  end do
  nT1amTot = nT1amTot+nT1am(iSym)
end do

if (nT1amTot > 0) iSkip = 1

return

end subroutine Check_Amp2
