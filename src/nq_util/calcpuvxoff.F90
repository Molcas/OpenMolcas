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

subroutine CalcPUVXOff()

use nq_Info, only: mIrrep, mOrb, nAsh, nPot2, nUVX, nUVXt, nVX, nVXt, OffPUVX, OffUVX, OffVX
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iIrrep, ijIrrep, IOff1, iOrb, jAct, jIrrep, kAct, kIrrep, klIrrep, lAct, lIrrep, nklAct

IOff1 = 0
do kIrrep=0,mIrrep-1
  kAct = nAsh(kIrrep)
  do lIrrep=0,kIrrep
    lAct = nAsh(lIrrep)
    nklAct = kAct*lAct
    if (kIrrep == lIrrep) nklAct = kAct*(kAct+1)/2
    OffVX(lIrrep,kIrrep) = IOff1
    nVX(lIrrep,kIrrep) = nklAct
    IOff1 = IOff1+nklAct
  end do
end do
nVXt = iOff1

IOff1 = 0
do jIrrep=0,mIrrep-1
  jAct = nAsh(jIrrep)
  do kIrrep=0,mIrrep-1
    kAct = nAsh(kIrrep)
    do lIrrep=0,kIrrep
      lAct = nAsh(lIrrep)
      nklAct = kAct*lAct
      if (kIrrep == lIrrep) nklAct = kAct*(kAct+1)/2
      OffUVX(lIrrep,kIrrep,jIrrep) = IOff1
      nUVX(lIrrep,kIrrep,jIrrep) = jAct*nklAct
      IOff1 = iOff1+jAct*nklAct
    end do
  end do
end do
nUVXt = IOff1

IOff1 = 0
do iIrrep=0,mIrrep-1
  OffPUVX(iIrrep) = IOff1
  iOrb = mOrb(iIrrep)
  do jIrrep=0,mIrrep-1
    jAct = nAsh(jIrrep)
    ijIrrep = 1+ieor(iIrrep,jIrrep)
    do kIrrep=0,mIrrep-1
      kAct = nAsh(kIrrep)
      do lIrrep=0,kIrrep
        lAct = nAsh(lIrrep)
        klIrrep = 1+ieor(kIrrep,lIrrep)
        if (ijIrrep == klIrrep) then
          iOff1 = iOff1+iOrb*nUVX(lIrrep,kIrrep,jIrrep)
        end if
      end do
    end do
  end do
end do
nPot2 = IOff1

!write(u6,*) 'OffPUVX new method',nPot2,MaxUVX
!write(u6,'(8(I5,2X))') (OffPUVX(iIrrep),iIrrep=0,mIrrep-1)

return

end subroutine CalcPUVXOff
