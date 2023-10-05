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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************

subroutine ChoMP2g_MOReOrd(CMO,COrb1,COrb2,iMoType1,iMOType2)
!
! Jonas Bostrom, Jan. 2010. (modified from ChoMP2_MOReOrd)
!
! Purpose: Make CMO:s of appropriate length, transpose COrb2,
!
!          CMO(alpha,p) -> COrb1(p,alpha)
!          CMO(alpha,q) -> COrb2(alpha,q)

use Cholesky, only: nBas, nSym
use ChoMP2, only: iAoMo, iMoAo, nMo
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(_OUT_) :: COrb1(*), COrb2(*)
integer(kind=iwp), intent(in) :: iMoType1, iMOType2
integer(kind=iwp) :: i, iCount, iSym, jCount, kOff1, kOff2, nOffOrb1(8), nOffOrb2(8), nOrb1(8), nOrb2(8)

do iSym=1,nSym
  nOffOrb1(iSym) = sum(nMo(iSym,1:iMOType1-1))
  nOffOrb2(iSym) = sum(nMo(iSym,1:iMOType2-1))
  nOrb1(iSym) = nMo(iSym,iMOType1)
  nOrb2(iSym) = nMo(iSym,iMOType2)
end do

iCount = 0
do iSym=1,nSym

  jCount = iCount+nOffOrb1(iSym)*nBas(iSym)

  do i=1,nOrb1(iSym)
    kOff1 = jCount+nBas(iSym)*(i-1)+1
    kOff2 = iMoAo(iSym,iSym,iMoType1)+i
    call dCopy_(nBas(iSym),CMO(kOff1),1,COrb1(kOff2),nOrb1(iSym))
  end do

  kOff1 = iCount+nOffOrb2(iSym)*nBas(iSym)
  kOff2 = iAoMo(iSym,iSym,iMoType2)
  COrb2(kOff2+1:kOff2+nBas(iSym)*nOrb2(iSym)) = CMO(kOff1+1:kOff1+nBas(iSym)*nOrb2(iSym))

  iCount = iCount+nBas(iSym)*nBas(iSym)

end do

end subroutine ChoMP2g_MOReOrd
