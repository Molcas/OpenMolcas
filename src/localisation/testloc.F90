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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine TestLoc(irc)
! Author: T.B. Pedersen
!
! Purpose: test localisation:
!
!          1) CC^T = XX^T ?
!          2) U=C^TSX , UU^T = diag ?
!          3) If C^TSC=1, X^TSX=1 ?
!
!          where C are the original MOs and X the localised ones.
!
!          Return codes: irc=0 (all OK), irc=1 (failure).

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp), intent(out) :: irc
#include "inflocal.fh"
#include "debug.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iComp, iOpt, ip0, ipDdff, ipDenC, ipDenX, ipOaux, ipOvlp, ipScr, ipUmat, iSyLbl, iSym, j, jrc, kC, kC1, &
                     kDC, kDd, kDX, kO, kOff, kSqr, kTri, kU, kX, kX1, lDen, lOaux, lOvlp, lScr, lUmat, mErr, nB2, nErr, nTO
real(kind=wp) :: Tol, Tst, xErr, xNrm
character(len=80) :: Txt
character(len=8) :: Label
logical(kind=iwp) :: Prnt
character(len=7), parameter :: SecNam = 'TestLoc'
integer(kind=iwp), external :: iPrintLevel
real(kind=r8), external :: ddot_

! Set return code.
! ----------------

irc = 0

! Set tolerance; allow more deviation for PAOs.
! ---------------------------------------------

if (LocPAO) then
  Tol = 1.0e-4_wp
else
  Tol = 1.0e-6_wp
end if

! Get the AO overlap matrix.
! --------------------------

lOvlp = nBas(1)**2
lOaux = 4+nBas(1)*(nBas(1)+1)/2
do iSym=2,nSym
  lOvlp = lOvlp+nBas(iSym)**2
  lOaux = lOaux+nBas(iSym)*(nBas(iSym)+1)/2
end do
call GetMem('TstOvlp','Allo','Real',ipOvlp,lOvlp)
call GetMem('TstOaux','Allo','Real',ipOaux,lOaux)
jrc = -1
iOpt = 2
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(jrc,iOpt,Label,iComp,Work(ipOaux),iSyLbl)
if (jrc /= 0) then
  write(Txt,'(A,I4)') 'RdOne returned',jrc
  call SysAbendMsg(SecNam,'I/O error!',Txt)
end if
Prnt = Debug .and. (iPrintLevel(-1) >= 5)
kTri = ipOaux
kSqr = ipOvlp
do iSym=1,nSym
  call Tri2Rec(Work(kTri),Work(kSqr),nBas(iSym),Prnt)
  kTri = kTri+nBas(iSym)*(nBas(iSym)+1)/2
  kSqr = kSqr+nBas(iSym)**2
end do
call GetMem('TstOaux','Free','Real',ipOaux,lOaux)

! Memory allocation.
! ------------------

lUmat = nOrb2Loc(1)**2
lScr = nBas(1)*nOrb2Loc(1)
do iSym=2,nSym
  lUmat = lUmat+nOrb2Loc(iSym)**2
  lScr = max(lScr,nBas(iSym)*nOrb2Loc(iSym))
end do
lDen = lOvlp
call GetMem('DenC','Allo','Real',ipDenC,lDen)
call GetMem('DenX','Allo','Real',ipDenX,lDen)
call GetMem('Ddff','Allo','Real',ipDdff,lDen)
call GetMem('Scratch','Allo','Real',ipScr,lScr)
call GetMem('Umat','Allo','Real',ipUmat,lUmat)

! Test 1) density.
! ----------------

nErr = 0

kDC = ipDenC
kDX = ipDenX
kDd = ipDdff
kC = ipMOrig
kX = ipCMO
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetDens_Localisation(Work(kDC),Work(kC1),nBas(iSym),nOrb2Loc(iSym))
  kX1 = kX+nBas(iSym)*nFro(iSym)
  call GetDens_Localisation(Work(kDX),Work(kX1),nBas(iSym),nOrb2Loc(iSym))
  nB2 = nBas(iSym)**2
  call dCopy_(nB2,Work(kDC),1,Work(kDd),1)
  call dAXPY_(nB2,-One,Work(kDX),1,Work(kDd),1)
  xNrm = sqrt(dDot_(nB2,Work(kDd),1,Work(kDd),1))
  if (xNrm > Tol) then
    write(u6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: ||CC^T - XX^T|| = ',xNrm,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kDC = kDC+nB2
  kDX = kDX+nB2
  kDd = kDd+nB2
  kC = kC+nB2
  kX = kX+nB2
end do
if (nErr /= 0) then
  irc = 1
  Go To 1 ! exit, after de-allocations
end if

! Test 2) U^TU = diag.
! --------------------

nErr = 0

kU = ipUmat
kO = ipOvlp
kC = ipMOrig
kX = ipCMO
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  kX1 = kX+nBas(iSym)*nFro(iSym)
  call GetUmat_Localisation(Work(kU),Work(kC1),Work(kO),Work(kX1),Work(ipScr),lScr,nBas(iSym),nOrb2Loc(iSym))
  nTO = max(nOrb2Loc(iSym),1)
  call DGEMM_('T','N',nOrb2Loc(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),One,Work(kU),nTO,Work(kU),nTO,Zero,Work(ipScr),nTO)
  xErr = -huge(xErr)
  ip0 = ipScr-1
  do j=1,nOrb2Loc(iSym)
    kOff = ip0+nOrb2Loc(iSym)*(j-1)
    do i=j+1,nOrb2Loc(iSym)
      Tst = abs(Work(kOff+i))
      xErr = max(xErr,Tst)
    end do
  end do
  if (xErr > Tol) then
    write(u6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: max. U^TU off-diag. = ',xErr,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  nB2 = nBas(iSym)**2
  kU = kU+nOrb2Loc(iSym)**2
  kO = kO+nB2
  kC = kC+nB2
  kX = kX+nB2
end do
if (nErr /= 0) then
  irc = 1
  Go To 1 ! exit, after de-allocations
end if

! Test 3) X^TSX=1.
! ----------------

nErr = 0

kU = ipUmat
kO = ipOvlp
kC = ipMOrig
kX = ipCMO
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetUmat_Localisation(Work(kU),Work(kC1),Work(kO),Work(kC1),Work(ipScr),lScr,nBas(iSym),nOrb2Loc(iSym))
  mErr = 0
  xErr = -huge(xErr)
  ip0 = kU-1
  do j=1,nOrb2Loc(iSym)
    kOff = ip0+nOrb2Loc(iSym)*(j-1)
    Tst = abs(Work(kOff+j)-One)
    if (Tst > Tol) then
      mErr = mErr+1
    end if
    do i=j+1,nOrb2Loc(iSym)
      Tst = abs(Work(kOff+i))
      xErr = max(xErr,Tst)
    end do
  end do
  if (xErr > Tol) then
    mErr = mErr+1
  end if
  if (mErr == 0) then
    kX1 = kX+nBas(iSym)*nFro(iSym)
    call GetUmat_Localisation(Work(kU),Work(kX1),Work(kO),Work(kX1),Work(ipScr),lScr,nBas(iSym),nOrb2Loc(iSym))
    xErr = -huge(xErr)
    ip0 = kU-1
    do j=1,nOrb2Loc(iSym)
      kOff = ip0+nOrb2Loc(iSym)*(j-1)
      Tst = abs(Work(kOff+j)-One)
      if (Tst > Tol) then
        mErr = mErr+1
      end if
      do i=j+1,nOrb2Loc(iSym)
        Tst = abs(Work(kOff+i))
        xErr = max(xErr,Tst)
      end do
    end do
    if (xErr > Tol) then
      write(u6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: max. X^TSX off-diag. = ',xErr,' (sym.',iSym,')'
      mErr = mErr+1
    end if
    if (mErr /= 0) nErr = nErr+1
  else
    write(u6,*) SecNam,': original orbitals not orthonormal!'
    write(u6,*) 'Orthonormality test is skipped (sym. ',iSym,')'
  end if
  nB2 = nBas(iSym)**2
  kU = kU+nOrb2Loc(iSym)**2
  kO = kO+nB2
  kC = kC+nB2
  kX = kX+nB2
end do
if (nErr /= 0) then
  irc = 1
  Go To 1 ! exit, after de-allocations
end if

! De-allocations and return.
! --------------------------

1 continue
call GetMem('Umat','Free','Real',ipUmat,lUmat)
call GetMem('Scratch','Free','Real',ipScr,lScr)
call GetMem('Ddff','Free','Real',ipDdff,lDen)
call GetMem('DenX','Free','Real',ipDenX,lDen)
call GetMem('DenC','Free','Real',ipDenC,lDen)
call GetMem('TstOvlp','Free','Real',ipOvlp,lOvlp)

end subroutine TestLoc
