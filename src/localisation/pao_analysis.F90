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
! Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine PAO_Analysis(D,R,X)
! Thomas Bondo Pedersen, December 2005.
! - revised January 2006 (Thomas Bondo Pedersen).
!
! Purpose: test and analysis of Cholesky PAOs.

use Localisation_globals, only: AnaNrm, BName, MOrig, nAtoms, nBas, nFro, nOrb2Loc, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: D(*), R(*), X(*)
integer(kind=iwp) :: i, iGetVecs, iSym, kC, kCO, kCR, kOff, kOffX, kSX, kX, l_S, l_SX, l_Tst, lDAt, lRAt, lXAt, nB, nB2, nErr, nF, &
                     nFO, nO, nR, nRest, nRO
real(kind=wp) :: xB2, xFO, xNrm, xRMS, xRO
character(len=14) :: FilNam
logical(kind=iwp) :: Debug
integer(kind=iwp), allocatable :: nBas_per_Atom(:), nBas_Start(:)
real(kind=wp), allocatable :: DAt(:), EigI(:), EigR(:), RAt(:), S(:), SX(:), Tst(:), XAt(:)
real(kind=wp), parameter :: Tol = 1.0e-10_wp
character(len=*), parameter :: SecNam = 'PAO_Analysis'
real(kind=wp), external :: ddot_

call Untested('PAO_Analysis')

l_S = nBas(1)**2
do iSym=2,nSym
  l_S = l_S+nBas(iSym)**2
end do
call mma_allocate(S,l_S,label='S')
call GetOvlp_Localisation(S,'Sqr',nBas,nSym)

! Initialization.
! ---------------

Debug = .false.

! Tests:
! 0) RR^T = XX^T = D (check of decomposition)
! 1) Co^TSX = 0 (check that X is orthogonal to the orthogonal
!                complement space).
! 2) C^TSX is non-singular (check that X spans primary space)
! ===========================================================

write(u6,*) 'Testing PAOs...'
nErr = 0

! Test 0.
! -------

l_Tst = nBas(1)**2
do iSym=2,nSym
  l_Tst = max(l_Tst,nBas(iSym)**2)
end do
call mma_allocate(Tst,l_Tst,label='TstDen')

kOff = 1
do iSym=1,nSym
  nB = nBas(iSym)
  nB2 = nB**2
  nF = nFro(iSym)
  nO = nOrb2Loc(iSym)
  call GetDens_Localisation(Tst,R(kOff),nB,nB)
  call dAXPY_(nB2,-One,D(kOff),1,Tst,1)
  if (nB2 > 0) then
    xB2 = One/real(nB2,kind=wp)
    xRMS = sqrt(xB2*dDot_(nB2,Tst,1,Tst,1))
  else
    xRMS = Zero
  end if
  if (xRMS > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: RMS(D-RR^T) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kOffX = kOff+nB*nF
  call GetDens_Localisation(Tst,X(kOffX),nB,nO)
  call dAXPY_(nB2,-One,D(kOff),1,Tst,1)
  if (nB2 > 0) then
    xB2 = One/real(nB2,kind=wp)
    xRMS = sqrt(xB2*dDot_(nB2,Tst,1,Tst,1))
  else
    xRMS = Zero
  end if
  if (xRMS > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: RMS(D-XX^T) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kOff = kOff+nB2
end do

! Compute SX.
! -----------

l_SX = nBas(1)*nOrb2Loc(1)
do iSym=2,nSym
  l_SX = l_SX+nBas(iSym)*nOrb2Loc(iSym)
end do
call mma_allocate(SX,l_SX,label='SX')

kOff = 1
kSX = 1
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  kX = kOff+nBas(iSym)*nFro(iSym)
  call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nBas(iSym),One,S(kOff),nB,X(kX),nB,Zero,SX(kSX),nB)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
  kOff = kOff+nBas(iSym)*nBas(iSym)
end do

! Test 1.
! First frozen part of orthogonal complement, then the rest.
! ----------------------------------------------------------

kC = 1
kSX = 1
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  nF = max(nFro(iSym),1)
  nRest = nBas(iSym)-nFro(iSym)-nOrb2Loc(iSym)
  nR = max(nRest,1)
  call DGEMM_('T','N',nFro(iSym),nOrb2Loc(iSym),nBas(iSym),One,MOrig(kC),nB,SX(kSX),nB,Zero,Tst,nF)
  nFO = nFro(iSym)*nOrb2Loc(iSym)
  if (nFO > 0) then
    xFO = One/real(nFO,kind=wp)
    xRMS = sqrt(xFO*dDot_(nFO,Tst,1,Tst,1))
  else
    xRMS = Zero
  end if
  if (xRMS > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: RMS(Co^TSX [Frozen]) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kCR = kC+nBas(iSym)*(nFro(iSym)+nOrb2Loc(iSym))
  call DGEMM_('T','N',nRest,nOrb2Loc(iSym),nBas(iSym),One,MOrig(kCR),nB,SX(kSX),nB,Zero,Tst,nR)
  nRO = nRest*nOrb2Loc(iSym)
  if (nRO > 0) then
    xRO = One/real(nRO,kind=wp)
    xRMS = sqrt(xRO*dDot_(nRO,Tst,1,Tst,1))
  else
    xRMS = Zero
  end if
  if (xRMS > Tol) then
    write(u6,'(A,A,ES16.8,A,I2,A)') SecNam,': ERROR: RMS(Co^TSX [Rest]) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kC = kC+nBas(iSym)*nBas(iSym)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
end do

! Test 2.
! -------

kC = 1
kSX = 1
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  nO = max(nOrb2Loc(iSym),1)
  kCO = kC+nBas(iSym)*nFro(iSym)
  call DGEMM_('T','N',nOrb2Loc(iSym),nOrb2Loc(iSym),nBas(iSym),One,MOrig(kCO),nB,SX(kSX),nB,Zero,Tst,nO)
  call mma_allocate(EigR,nOrb2Loc(iSym),label='EigR')
  call mma_allocate(EigI,nOrb2Loc(iSym),label='EigI')
  iGetVecs = 0
  call Diag_Localisation(Tst,EigR,EigI,nOrb2Loc(iSym),iGetVecs)
  do i=1,nOrb2Loc(iSym)
    xNrm = sqrt(EigR(i)**2+EigI(i)**2)
    if (xNrm < Tol) then
      write(u6,'(A,A,I6,A,A,ES16.8,A,I2,A)') SecNam,': ERROR: ||eigenvalue',i,'|| of ','C^TSX = ',xNrm,' (sym.',iSym,')'
      nErr = nErr+1
    end if
  end do
  call mma_deallocate(EigR)
  call mma_deallocate(EigI)
  kC = kC+nBas(iSym)*nBas(iSym)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
end do

call mma_deallocate(Tst)
call mma_deallocate(SX)

! Overall check.
! --------------

if (nErr /= 0) then
  call SysAbendMsg(SecNam,'PAO localization failed!',' ')
else
  write(u6,*) '...OK!'
end if

! Analysis.
! This part does not work with symmetry.
! ======================================

if (nSym /= 1) then
  write(u6,*)
  write(u6,*) SecNam,': symmetry not implemented for analysis of PAOs. Sorry!'
  return
end if

! Allocate atom based density and CMO matrices.
! ---------------------------------------------

lDAt = nAtoms**2
lRAt = nAtoms*nBas(1)
lXAt = nAtoms*nOrb2Loc(1)
call mma_allocate(DAt,lDAt,label='DAt')
call mma_allocate(RAt,lRAt,label='RAt')
call mma_allocate(XAt,lXAt,label='XAt')

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtoms,label='nB_per_Atom')
call mma_allocate(nBas_Start,nAtoms,label='nB_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,BName,nBas(1),nAtoms,Debug)

! Generate bitmaps.
! -----------------

call GetAt_Localisation(D,nBas(1),nBas(1),DAt,nAtoms,2,nBas_per_Atom,nBas_Start,AnaNrm)
call GetAt_Localisation(R,nBas(1),nBas(1),RAt,nAtoms,1,nBas_per_Atom,nBas_Start,AnaNrm)
kOffX = nBas(1)*nFro(1)+1
call GetAt_Localisation(X(kOffX),nBas(1),nOrb2Loc(1),XAt,nAtoms,1,nBas_per_Atom,nBas_Start,AnaNrm)
write(FilNam,'(A)') 'PAO_Dnsty1.bmp'
call GenBMp_Loc(DAt,nAtoms,nAtoms,FilNam,'r')
write(FilNam,'(A)') 'PAO_LnDep1.bmp'
call GenBMp_Loc(RAt,nAtoms,nBas(1),FilNam,'r')
write(FilNam,'(A)') 'PAO_Chole1.bmp'
call GenBMp_Loc(XAt,nAtoms,nOrb2Loc(1),FilNam,'r')
write(u6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! Deallocations.
! --------------

call mma_deallocate(DAt)
call mma_deallocate(RAt)
call mma_deallocate(XAt)
call mma_deallocate(nBas_per_Atom)
call mma_deallocate(nBas_Start)

call mma_deallocate(S)

end subroutine PAO_Analysis
