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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
subroutine CalcOrbOff()

use nq_Info, only: iOff_Ash, iOff_Bas, iOff_BasAct, mBas, mIrrep, mOrb, NASH, NASHT, nFro, nIsh, nOrbt, nPot1, OffBas, OffBas2, &
                   OffBasFro, OffOrb, OffOrb2, OffOrbTri
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iIrrep, jOffA_, jOffB_, nTri

NASHT = 0
jOffA_ = 0
jOffB_ = 0
nPot1 = 0
nTri = 0
nOrbt = 0
do iIrrep=0,mIrrep-1
  mOrb(iIrrep) = mBas(iIrrep)-nFro(iIrrep)
  nPot1 = nPot1+mOrb(iIrrep)**2
  nOrbt = nOrbt+mOrb(iIrrep)
  NASHT = NASHT+NASH(iIrrep)
  iOff_Ash(iIrrep) = jOffA_
  iOff_Bas(iIrrep) = jOffB_
  OffBasFro(iIrrep) = jOffB_+nFro(iIrrep)
  iOff_BasAct(iIrrep) = jOffB_+nIsh(iIrrep)+nFro(iIrrep)
  OffOrbTri(iIrrep) = nTri
  nTri = nTri+mOrb(iIrrep)*(mOrb(iIrrep)+1)/2
  jOffA_ = jOffA_+nAsh(iIrrep)
  jOffB_ = jOffB_+mBas(iIrrep)
end do

OffOrb(0) = 0
OffBas(0) = 1
OffBas2(0) = 1
OffOrb2(0) = 0
do IIrrep=1,mIrrep-1
  OffBas(iIrrep) = OffBas(iIrrep-1)+mBas(iIrrep-1)
  OffOrb(iIrrep) = OffOrb(iIrrep-1)+mOrb(iIrrep-1)
  OffBas2(iIrrep) = OffBas2(iIrrep-1)+mBas(iIrrep-1)**2
  OffOrb2(iIrrep) = OffOrb2(iIrrep-1)+mOrb(iIrrep-1)**2
end do

return

end subroutine CalcOrbOff
