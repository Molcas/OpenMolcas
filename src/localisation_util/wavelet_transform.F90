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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine Wavelet_Transform(irc,CMO,nSym,nBas,nFro,nOrb2Loc,inv,Silent,xNrm)
! Author: F. Aquilante
!
! Purpose: wavelet transform of the MO basis (inv=0)
!          "       backtransform (inv=1)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: irc, nSym, nBas(nSym), nFro(nSym), nOrb2Loc(nSym), inv
real(kind=wp) :: CMO(*), xNrm
logical(kind=iwp) :: Silent
#include "WrkSpc.fh"
integer(kind=iwp) :: ipScr, iScr, iSym, kOff1, kOff2, kOffC, l_Scr, njOrb
character(len=*), parameter :: SecNam = 'Wavelet_Transform'
integer(kind=iwp), external :: Log2
real(kind=r8), external :: ddot_

irc = 0
xNrm = Zero
if (.not. Silent) then
  if (inv == 0) write(u6,'(/,1X,A)') 'Wavelet transform of the MOs'
  if (inv == 1) write(u6,'(/,1X,A)') 'Inverse wavelet transform of the MOs'
  write(u6,'(1X,A,8(1X,I6))') 'Frozen orbitals      :',(nFro(iSym),iSym=1,nSym)
  write(u6,'(1X,A,8(1X,I6))') 'Orbitals to transform:',(nOrb2Loc(iSym),iSym=1,nSym)
end if

if (inv == 1) go to 1000 ! Inverse wavelet transform

njOrb = Log2(nOrb2Loc(1))
l_Scr = nBas(1)*(2**njOrb-1)
do iSym=2,nSym
  njOrb = Log2(nOrb2Loc(iSym))
  l_Scr = max(l_Scr,nBas(iSym)*(2**njOrb-1))
end do
call GetMem('Scratch','Allo','Real',ipScr,l_Scr)
kOffC = 1
do iSym=1,nSym
  if (nOrb2Loc(iSym) > 0) then
    kOff1 = kOffC+nBas(iSym)*nFro(iSym)
    kOff2 = kOff1
    njOrb = Log2(nOrb2Loc(iSym))
    do while (njOrb >= 1)
      call FWT_Haar(nBas(iSym),njOrb,Work(ipScr),CMO(kOff2))
      njOrb = 2**njOrb
      kOff2 = kOff2+nBas(iSym)*njOrb
      njOrb = Log2(nOrb2Loc(iSym)-njOrb)
    end do
    xNrm = xNrm+dDot_(nBas(iSym)*nOrb2Loc(iSym),CMO(kOff1),1,CMO(kOff1),1)
    if (irc /= 0) then
      irc = 1
      xNrm = -huge(xNrm)
      return
    end if
  end if
  kOffC = kOffC+nBas(iSym)**2
end do
xNrm = sqrt(xNrm)
call GetMem('Scratch','Free','Real',ipScr,l_Scr)
return

1000 continue
njOrb = Log2(nOrb2Loc(1))
l_Scr = nBas(1)*2**njOrb
do iSym=2,nSym
  njOrb = Log2(nOrb2Loc(iSym))
  l_Scr = max(l_Scr,nBas(iSym)*2**njOrb)
end do
call GetMem('Scratch','Allo','Real',iScr,l_Scr)
kOffC = 1
do iSym=1,nSym
  if (nOrb2Loc(iSym) > 0) then
    kOff1 = kOffC+nBas(iSym)*nFro(iSym)
    kOff2 = kOff1
    njOrb = Log2(nOrb2Loc(iSym))
    do while (njOrb >= 1)
      call Inv_FWT_Haar(nBas(iSym),njOrb,Work(iScr),CMO(kOff2))
      njOrb = 2**njOrb
      kOff2 = kOff2+nBas(iSym)*njOrb
      njOrb = Log2(nOrb2Loc(iSym)-njOrb)
    end do
    xNrm = xNrm+dDot_(nBas(iSym)*nOrb2Loc(iSym),CMO(kOff1),1,CMO(kOff1),1)
    if (irc /= 0) then
      irc = 1
      xNrm = -huge(xNrm)
      return
    end if
  end if
  kOffC = kOffC+nBas(iSym)**2
end do
xNrm = sqrt(xNrm)
call GetMem('Scratch','Free','Real',iScr,l_Scr)

end subroutine Wavelet_Transform
