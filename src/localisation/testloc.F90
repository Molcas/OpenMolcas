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

use Localisation_globals, only: CMO, LocPAO, MOrig, nBas, nFro, nOrb2Loc, nSym
use OneDat, only: sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: i, iComp, iOpt, ip0, iSyLbl, iSym, j, jrc, kC, kC1, kD, kOff, kSqr, kTri, kU, lOaux, lOvlp, lScr, lUmat, &
                     mErr, nB2, nErr, nTO
real(kind=wp) :: Tol, Tst, xErr, xNrm
character(len=80) :: Txt
character(len=8) :: Label
logical(kind=iwp) :: Prnt
real(kind=wp), allocatable :: DenC(:), DenX(:), Ddff(:), Oaux(:), Ovlp(:), Scr(:), Umat(:)
character(len=*), parameter :: SecNam = 'TestLoc'
integer(kind=iwp), external :: iPrintLevel
real(kind=wp), external :: ddot_
logical(kind=iwp), parameter :: debug = .false.

call Untested('TestLoc')

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
call mma_allocate(Ovlp,lOvlp,label='TstOvlp')
call mma_allocate(Oaux,lOaux,label='TstOaux')
jrc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(jrc,iOpt,Label,iComp,Oaux,iSyLbl)
if (jrc /= 0) then
  write(Txt,'(A,I4)') 'RdOne returned',jrc
  call SysAbendMsg(SecNam,'I/O error!',Txt)
end if
Prnt = Debug .and. (iPrintLevel(-1) >= 5)
kTri = 1
kSqr = 1
do iSym=1,nSym
  call Tri2Rec(Oaux(kTri),Ovlp(kSqr),nBas(iSym),Prnt)
  kTri = kTri+nBas(iSym)*(nBas(iSym)+1)/2
  kSqr = kSqr+nBas(iSym)**2
end do
call mma_deallocate(Oaux)

! Memory allocation.
! ------------------

lUmat = nOrb2Loc(1)**2
lScr = nBas(1)*nOrb2Loc(1)
do iSym=2,nSym
  lUmat = lUmat+nOrb2Loc(iSym)**2
  lScr = max(lScr,nBas(iSym)*nOrb2Loc(iSym))
end do
call mma_allocate(DenC,lOvlp,label='DenC')
call mma_allocate(DenX,lOvlp,label='DenX')
call mma_allocate(Ddff,lOvlp,label='Ddff')
call mma_allocate(Scr,lScr,label='Scratch')
call mma_allocate(Umat,lUmat,label='Umat')

! Test 1) density.
! ----------------

nErr = 0

kD = 1
kC = 1
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetDens_Localisation(DenC(kD),MOrig(kC1),nBas(iSym),nOrb2Loc(iSym))
  call GetDens_Localisation(DenX(kD),CMO(kC1),nBas(iSym),nOrb2Loc(iSym))
  nB2 = nBas(iSym)**2
  call dCopy_(nB2,DenC(kD),1,Ddff(kD),1)
  call dAXPY_(nB2,-One,DenX(kD),1,Ddff(kD),1)
  xNrm = sqrt(dDot_(nB2,Ddff(kD),1,Ddff(kD),1))
  if (xNrm > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: ||CC^T - XX^T|| = ',xNrm,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kD = kD+nB2
  kC = kC+nB2
end do
if (nErr /= 0) then
  call Error(1) ! exit, after de-allocations
  return
end if

! Test 2) U^TU = diag.
! --------------------

nErr = 0

kU = 1
kC = 1
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetUmat_Localisation(Umat(kU),MOrig(kC1),Ovlp(kC),CMO(kC1),Scr,nBas(iSym),nOrb2Loc(iSym))
  nTO = max(nOrb2Loc(iSym),1)
  call DGEMM_('T','N',nOrb2Loc(iSym),nOrb2Loc(iSym),nOrb2Loc(iSym),One,Umat(kU),nTO,Umat(kU),nTO,Zero,Scr,nTO)
  xErr = -huge(xErr)
  ip0 = 0
  do j=1,nOrb2Loc(iSym)
    kOff = ip0+nOrb2Loc(iSym)*(j-1)
    do i=j+1,nOrb2Loc(iSym)
      Tst = abs(Scr(kOff+i))
      xErr = max(xErr,Tst)
    end do
  end do
  if (xErr > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: max. U^TU off-diag. = ',xErr,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  nB2 = nBas(iSym)**2
  kU = kU+nOrb2Loc(iSym)**2
  kC = kC+nB2
end do
if (nErr /= 0) then
  call Error(1) ! exit, after de-allocations
  return
end if

! Test 3) X^TSX=1.
! ----------------

nErr = 0

kU = 1
kC = 1
do iSym=1,nSym
  kC1 = kC+nBas(iSym)*nFro(iSym)
  call GetUmat_Localisation(Umat(kU),MOrig(kC1),Ovlp(kC),MOrig(kC1),Scr,nBas(iSym),nOrb2Loc(iSym))
  mErr = 0
  xErr = -huge(xErr)
  ip0 = kU-1
  do j=1,nOrb2Loc(iSym)
    kOff = ip0+nOrb2Loc(iSym)*(j-1)
    Tst = abs(Umat(kOff+j)-One)
    if (Tst > Tol) then
      mErr = mErr+1
    end if
    do i=j+1,nOrb2Loc(iSym)
      Tst = abs(Umat(kOff+i))
      xErr = max(xErr,Tst)
    end do
  end do
  if (xErr > Tol) then
    mErr = mErr+1
  end if
  if (mErr == 0) then
    call GetUmat_Localisation(Umat(kU),CMO(kC1),Ovlp(kC),CMO(kC1),Scr,nBas(iSym),nOrb2Loc(iSym))
    xErr = -huge(xErr)
    ip0 = kU-1
    do j=1,nOrb2Loc(iSym)
      kOff = ip0+nOrb2Loc(iSym)*(j-1)
      Tst = abs(Umat(kOff+j)-One)
      if (Tst > Tol) then
        mErr = mErr+1
      end if
      do i=j+1,nOrb2Loc(iSym)
        Tst = abs(Umat(kOff+i))
        xErr = max(xErr,Tst)
      end do
    end do
    if (xErr > Tol) then
      write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: max. X^TSX off-diag. = ',xErr,' (sym.',iSym,')'
      mErr = mErr+1
    end if
    if (mErr /= 0) nErr = nErr+1
  else
    write(u6,*) SecNam,': original orbitals not orthonormal!'
    write(u6,*) 'Orthonormality test is skipped (sym. ',iSym,')'
  end if
  nB2 = nBas(iSym)**2
  kU = kU+nOrb2Loc(iSym)**2
  kC = kC+nB2
end do
if (nErr /= 0) then
  call Error(1) ! exit, after de-allocations
  return
end if

! De-allocations and return.
! --------------------------

call Error(0)

contains

subroutine Error(code)
  integer(kind=iwp), intent(in) :: code
  if (code /= 0) irc = code
  call mma_deallocate(Ovlp)
  call mma_deallocate(DenC)
  call mma_deallocate(DenX)
  call mma_deallocate(Ddff)
  call mma_deallocate(Scr)
  call mma_deallocate(Umat)
end subroutine Error

end subroutine TestLoc
