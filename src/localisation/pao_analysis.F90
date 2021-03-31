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

implicit real*8(a-h,o-z)
real*8 D(*), R(*), X(*)
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

l_S = nBas(1)**2
do iSym=2,nSym
  l_S = l_S+nBas(iSym)**2
end do
call GetMem('S','Allo','Real',ip_S,l_S)
call GetOvlp_Localisation(Work(ip_S),'Sqr',nBas,nSym)
call PAO_Ana1(D,R,X,Work(ipMOrig),Work(ip_S),nBas,nFro,nOrb2Loc,nSym,Name,nAtoms,AnaNrm)
call GetMem('S','Free','Real',ip_S,l_S)

end subroutine PAO_Analysis

subroutine PAO_Ana1(D,R,X,C,S,nBas,nFro,nOrb2Loc,nSym,Nam,nAtoms,AnaNrm)

implicit real*8(a-h,o-z)
#include "Molcas.fh"
integer nBas(nSym), nFro(nSym), nOrb2Loc(nSym)
real*8 D(*), R(*), X(*), C(*), S(*)
character*(LENIN8) Nam(*)
character*3 AnaNrm
#include "WrkSpc.fh"

character*8 SecNam
parameter(SecNam='PAO_Ana1')

character*14 FilNam
logical Debug

external ddot_

parameter(Tol=1.0d-10)

! Initialization.
! ---------------

Debug = .false.

! Tests:
! 0) RR^T = XX^T = D (check of decomposition)
! 1) Co^TSX = 0 (check that X is orthogonal to the orthogonal
!                complement space).
! 2) C^TSX is non-singular (check that X spans primary space)
! ===========================================================

write(6,*) 'Testing PAOs...'
nErr = 0

! Test 0.
! -------

l_Tst = nBas(1)**2
do iSym=2,nSym
  l_Tst = max(l_Tst,nBas(iSym)**2)
end do
call GetMem('TstDen','Allo','Real',ip_Tst,l_Tst)

kOff = 1
do iSym=1,nSym
  nB = nBas(iSym)
  nB2 = nB**2
  nF = nFro(iSym)
  nO = nOrb2Loc(iSym)
  call GetDens_Localisation(Work(ip_Tst),R(kOff),nB,nB)
  call dAXPY_(nB2,-1.0d0,D(kOff),1,Work(ip_Tst),1)
  if (nB2 > 0) then
    xB2 = 1.0d0/dble(nB2)
    xRMS = sqrt(xB2*dDot_(nB2,Work(ip_Tst),1,Work(ip_Tst),1))
  else
    xRMS = 0.0d0
  end if
  if (xRMS > Tol) then
    write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: RMS(D-RR^T) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kOffX = kOff+nB*nF
  call GetDens_Localisation(Work(ip_Tst),X(kOffX),nB,nO)
  call dAXPY_(nB2,-1.0d0,D(kOff),1,Work(ip_Tst),1)
  if (nB2 > 0) then
    xB2 = 1.0d0/dble(nB2)
    xRMS = sqrt(xB2*dDot_(nB2,Work(ip_Tst),1,Work(ip_Tst),1))
  else
    xRMS = 0.0d0
  end if
  if (xRMS > Tol) then
    write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: RMS(D-XX^T) = ',xRMS,' (sym.',iSym,')'
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
call GetMem('SX','Allo','Real',ip_SX,l_SX)

kOff = 1
kSX = ip_SX
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  kX = kOff+nBas(iSym)*nFro(iSym)
  call DGEMM_('N','N',nBas(iSym),nOrb2Loc(iSym),nBas(iSym),1.0d0,S(kOff),nB,X(kX),nB,0.0d0,Work(kSX),nB)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
  kOff = kOff+nBas(iSym)*nBas(iSym)
end do

! Test 1.
! First frozen part of orthogonal complement, then the rest.
! ----------------------------------------------------------

kC = 1
kSX = ip_SX
kTst = ip_Tst
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  nF = max(nFro(iSym),1)
  nRest = nBas(iSym)-nFro(iSym)-nOrb2Loc(iSym)
  nR = max(nRest,1)
  call DGEMM_('T','N',nFro(iSym),nOrb2Loc(iSym),nBas(iSym),1.0d0,C(kC),nB,Work(kSX),nB,0.0d0,Work(kTst),nF)
  nFO = nFro(iSym)*nOrb2Loc(iSym)
  if (nFO > 0) then
    xFO = 1.0d0/dble(nFO)
    xRMS = sqrt(xFO*dDot_(nFO,Work(kTst),1,Work(kTst),1))
  else
    xRMS = 0.0d0
  end if
  if (xRMS > Tol) then
    write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: RMS(Co^TSX [Frozen]) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kCR = kC+nBas(iSym)*(nFro(iSym)+nOrb2Loc(iSym))
  call DGEMM_('T','N',nRest,nOrb2Loc(iSym),nBas(iSym),1.0d0,C(kCR),nB,Work(kSX),nB,0.0d0,Work(kTst),nR)
  nRO = nRest*nOrb2Loc(iSym)
  if (nRO > 0) then
    xRO = 1.0d0/dble(nRO)
    xRMS = sqrt(xRO*dDot_(nRO,Work(kTst),1,Work(kTst),1))
  else
    xRMS = 0.0d0
  end if
  if (xRMS > Tol) then
    write(6,'(A,A,D16.8,A,I2,A)') SecNam,': ERROR: RMS(Co^TSX [Rest]) = ',xRMS,' (sym.',iSym,')'
    nErr = nErr+1
  end if
  kC = kC+nBas(iSym)*nBas(iSym)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
end do

! Test 2.
! -------

kC = 1
kSX = ip_SX
kTst = ip_Tst
do iSym=1,nSym
  nB = max(nBas(iSym),1)
  nO = max(nOrb2Loc(iSym),1)
  kCO = kC+nBas(iSym)*nFro(iSym)
  call DGEMM_('T','N',nOrb2Loc(iSym),nOrb2Loc(iSym),nBas(iSym),1.0d0,C(kCO),nB,Work(kSX),nB,0.0d0,Work(kTst),nO)
  l_EigR = nOrb2Loc(iSym)
  l_EigI = nOrb2Loc(iSym)
  call GetMem('EigR','Allo','Real',ip_EigR,l_EigR)
  call GetMem('EigI','Allo','Real',ip_EigI,l_EigI)
  iGetVecs = 0
  call Diag_Localisation(Work(kTst),Work(ip_EigR),Work(ip_EigI),nOrb2Loc(iSym),iGetVecs)
  kOffR = ip_EigR-1
  kOffI = ip_EigI-1
  do i=1,nOrb2Loc(iSym)
    xNrm = sqrt(Work(kOffR+i)**2+Work(kOffI+i)**2)
    if (xNrm < Tol) then
      write(6,'(A,A,I6,A,A,D16.8,A,I2,A)') SecNam,': ERROR: ||eigenvalue',i,'|| of ','C^TSX = ',xNrm,' (sym.',iSym,')'
      nErr = nErr+1
    end if
  end do
  call GetMem('EigI','Free','Real',ip_EigI,l_EigI)
  call GetMem('EigR','Free','Real',ip_EigR,l_EigR)
  kC = kC+nBas(iSym)*nBas(iSym)
  kSX = kSX+nBas(iSym)*nOrb2Loc(iSym)
end do

call GetMem('SX','Free','Real',ip_SX,l_SX)
call GetMem('TstDen','Free','Real',ip_Tst,l_Tst)

! Overall check.
! --------------

if (nErr /= 0) then
  call SysAbendMsg(SecNam,'PAO localization failed!',' ')
else
  write(6,*) '...OK!'
end if

! Analysis.
! This part does not work with symmetry.
! ======================================

if (nSym /= 1) then
  write(6,*)
  write(6,*) SecNam,': symmetry not implemented for analysis of PAOs. Sorry!'
  return
end if

! Allocate atom based density and CMO matrices.
! ---------------------------------------------

lDAt = nAtoms**2
lRAt = nAtoms*nBas(1)
lXAt = nAtoms*nOrb2Loc(1)
call GetMem('DAt','Allo','Real',ipDAt,lDAt)
call GetMem('RAt','Allo','Real',ipRAt,lRAt)
call GetMem('XAt','Allo','Real',ipXAt,lXAt)

! Allocate and get index arrays for basis functions per atom.
! -----------------------------------------------------------

l_nBas_per_Atom = nAtoms
l_nBas_Start = nAtoms
call GetMem('nB_per_Atom','Allo','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('nB_Start','Allo','Inte',ip_nBas_Start,l_nBas_Start)
call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),Nam,nBas(1),nAtoms,Debug)

! Generate bitmaps.
! -----------------

call GetAt_Localisation(D,nBas(1),nBas(1),Work(ipDAt),nAtoms,2,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
call GetAt_Localisation(R,nBas(1),nBas(1),Work(ipRAt),nAtoms,1,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
kOffX = nBas(1)*nFro(1)+1
call GetAt_Localisation(X(kOffX),nBas(1),nOrb2Loc(1),Work(ipXAt),nAtoms,1,iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),AnaNrm)
write(FilNam,'(A)') 'PAO_Dnsty1.bmp'
call GenBMp_Loc(Work(ipDat),nAtoms,nAtoms,FilNam,'r')
write(FilNam,'(A)') 'PAO_LnDep1.bmp'
call GenBMp_Loc(Work(ipRat),nAtoms,nBas(1),FilNam,'r')
write(FilNam,'(A)') 'PAO_Chole1.bmp'
call GenBMp_Loc(Work(ipXat),nAtoms,nOrb2Loc(1),FilNam,'r')
write(6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

! Deallocations.
! --------------

call GetMem('nB_Start','Free','Inte',ip_nBas_Start,l_nBas_Start)
call GetMem('nB_per_Atom','Free','Inte',ip_nBas_per_Atom,l_nBas_per_Atom)
call GetMem('XAt','Free','Real',ipXAt,lXAt)
call GetMem('RAt','Free','Real',ipRAt,lRAt)
call GetMem('DAt','Free','Real',ipDAt,lDAt)

end subroutine PAO_Ana1
