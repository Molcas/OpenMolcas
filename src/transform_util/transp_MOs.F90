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

subroutine Transp_MOs(CMO1,CMO2,nSym,nFro,nIsh,nAsh,nSsh,nBas)
!        CMO1(alpha,p) -> CMO2(p,alpha)
!
!        At the same time, exclude Fro and Del orbitals

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO1(*)
real(kind=wp), intent(_OUT_) :: CMO2(*)
integer(kind=iwp), intent(in) :: nSym, nFro(nSym), nIsh(nSym), nAsh(nSym), nSsh(nSym), nBas(nSym)
integer(kind=iwp) :: i, iCount, iSym, jCount, kOff1, kOff2, lCount, nOrb

iCount = 0
lCount = 0
do iSym=1,nSym

  jCount = iCount+nBas(iSym)*nFro(iSym)

  nOrb = nIsh(iSym)+nAsh(iSym)+nSsh(iSym)

  do i=1,nOrb
    kOff1 = jCount+nBas(iSym)*(i-1)+1
    kOff2 = lCount+i
    call dCopy_(nBas(iSym),CMO1(kOff1),1,CMO2(kOff2),nOrb)
  end do

  lCount = lCount+nOrb*nBas(iSym)
  iCount = iCount+nBas(iSym)**2

end do

end subroutine Transp_MOs
