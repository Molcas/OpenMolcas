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

subroutine ChoMP2g_density2(irc,EOcc,EVir,EFro,Wrk,lWrk)
!
! Jonas Bostrom, Mars 2010
!
! Purpose: Solve the CPHF-equations to obtain occ-vir
!          contributions to the 1-pdm.

use Cholesky, only: nSym, NumCho
use ChoMP2, only: iAdrOff, iFro, iOcc, iVir, kFLagr, kLagr, kPab, kPai, kPaK, kPij, kPik, kWab, kWai, kWaK, kWij, kWiK, kWJK, &
                  lFLagr, lUnit_F, LuVVec, LuWVec, MP2D, MP2W, nFro, nMoMo, nOcc, nOrb, nVir
use Constants, only: Zero, One, Two, Three, Four, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: EOcc(*), EVir(*), EFro(*)
integer(kind=iwp), intent(in) :: lWrk
real(kind=wp), intent(inout) :: Wrk(lWrk)
integer(kind=iwp) :: i, iA, iAdr, iB, iBat, iClos, iI, iIter, iJ, indx, index1, iOff, iOff1, iOff2, iOffL, iOpt, iOrb, iOrbI, &
                     iSeed, iSym, iSym1, iSymOffOV, iTypL, iVec, iVec1, iVecFF, iVecFO, iVecFV, iVecOF, iVecOO, iVecOV, iVecVF, &
                     iVecVO, iVecVV, j, jOrb, kCGVec(9), kDiag(2), kEndCGVec(9), kEndDiag, kEndFDiag, kEndLab, kEndLfa, kEndLia, &
                     kEndLij, kEndLiK, kEndLip, kEndLJK, kEndLKa, kEndLKi, kEndU, kEndVip, kLab, kLfa, kLia, kLij, kLiK, kLip, &
                     kLJK, kLKa, kLKi, kU, kVip, kWij2, lCGFVec, lCGOVec, lCGVec, lDiag, lFDiag, lLab, lLfa, lLia, lLij, lLiK, &
                     lLip, lLJK, lLKa, lLKi, lScr, lTot, lU, lVip, lWij2, maxvalue, nA, nBatL, nCGvec, nI, nIter, nOccAll(8), &
                     NumVec, nVec
real(kind=wp) :: Ea, Eb, Ei, Ej, Eps, Res, term
logical(kind=iwp) :: Done
character(len=5) :: Fname
character(len=*), parameter :: SecNam = 'ChoMP2g_Density2'
integer(kind=iwp), external :: IsFreeUnit

! Do not delete vectors after use.
! --------------------------------
iClos = 2
iTypL = 1
iVecFF = 1
iVecOV = 6
iVecVV = 9
iVecOO = 5
iVecFV = 3
iVecFO = 2
iVecOF = 4

! Setup the conjugate gradient procedure
! --------------------------------------
Eps = 1.0e-8_wp
Done = .false.
nIter = 100

! Allocate memory for Frozen Diagonal of A
! ----------------------------------------
lFDiag = nMoMo(1,iVecFV)
kDiag(1) = kFLagr(1)+lFLagr
kEndFDiag = kDiag(1)+lFDiag
Wrk(kDiag(1):kEndFDiag-1) = Zero

! Allocate memory for Occupied Diagonal of A
! ------------------------------------------
lDiag = nMoMo(1,iVecOV)
kDiag(2) = kEndFDiag
kEndDiag = kDiag(2)+lDiag
Wrk(kDiag(2):kEndDiag-1) = Zero

iSym = 1
nOccAll(iSym) = nFro(iSym)+nOcc(iSym)

! Open Cholesky vector files.
! ---------------------------
call ChoMP2_OpenF(1,1,iSym)

maxvalue = 200
nVec = min(NumCho(iSym),maxvalue)
nBatL = (NumCho(iSym)-1)/nVec+1
if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')

! Allocate memory for Lia-vector and LIa-vector
! ---------------------------------------------

lLfa = nMoMo(iSym,iVecFV)*nVec
kLfa = kEndDiag
kEndLfa = kLfa+lLfa

lLia = nMoMo(iSym,iVecOV)*nVec
kLia = kEndLfa
kEndLia = kLia+lLia

do iBat=1,nBatL
  if (iBat == nBatL) then
    NumVec = NumCho(iSym)-nVec*(nBatL-1)
  else
    NumVec = nVec
  end if
  iVec = nVec*(iBat-1)+1

  ! Read Lfa-vectors
  ! ----------------
  iOpt = 2
  lTot = nMoMo(iSym,iVecFV)*NumVec
  iAdr = nMoMo(iSym,iVecFV)*(iVec-1)+1+iAdrOff(iSym,iVecFV)
  call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLfa),lTot,iAdr)

  ! Read Lia-vectors
  ! ----------------
  iOpt = 2
  lTot = nMoMo(iSym,iVecOV)*NumVec
  iAdr = nMoMo(iSym,iVecOV)*(iVec-1)+1+iAdrOff(iSym,iVecOV)
  call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,iAdr)

  ! Construct Diagonal of A
  ! -----------------------
  do iVec1=1,NumVec
    do i=1,nMoMo(iSym,iVecFV)
      iOffL = (iVec1-1)*nMoMo(iSym,iVecFV)
      Wrk(kDiag(1)+i-1) = Wrk(kDiag(1)+i-1)+Three*Wrk(kLfa+i-1+iOffL)**2
    end do
  end do
  do iVec1=1,NumVec
    do i=1,nMoMo(iSym,iVecOV)
      iOffL = (iVec1-1)*nMoMo(iSym,iVecOV)
      Wrk(kDiag(2)+i-1) = Wrk(kDiag(2)+i-1)+Three*Wrk(kLia+i-1+iOffL)**2
    end do
  end do

  index1 = 0
  do iSym1=1,nSym
    nA = nVir(iSym1)
    nI = nFro(iSym1)
    do iI=1,nI
      Ei = EFro(iFro(iSym1)+iI)
      do iA=1,nA
        indx = iA-1+(iI-1)*nA+index1
        Ea = EVir(iVir(iSym1)+iA)
        Wrk(kDiag(1)+indx) = Wrk(kDiag(1)+indx)+Ea-Ei
      end do
    end do
    nI = nOcc(iSym1)
    do iI=1,nI
      Ei = EOcc(iOcc(iSym1)+iI)
      do iA=1,nA
        indx = iA-1+(iI-1)*nA+index1
        Ea = EVir(iVir(iSym1)+iA)
        Wrk(kDiag(2)+indx) = Wrk(kDiag(2)+indx)+Ea-Ei
      end do
    end do
    index1 = index1+nOcc(iSym1)*nVir(iSym1)
  end do
  do i=1,lFDiag
    Wrk(kDiag(1)+i-1) = 1/Wrk(kDiag(1)+i-1)
  end do
  do i=1,lDiag
    Wrk(kDiag(2)+i-1) = 1/Wrk(kDiag(2)+i-1)
  end do
end do

call ChoMP2_OpenF(iClos,1,iSym)

! Allocate vectors needed for the PCG
! -----------------------------------
lCGFVec = sum(nFro(1:nSym)*nVir(1:nSym))
lCGOVec = sum(nOcc(1:nSym)*nVir(1:nSym))
lCGVec = lCGOVec+lCGFVec
! Vector Legend (as they are named in Conj_Grad):
!                Z = 1, Ztemp = 2
!                R = 3, Rtemp = 4
!                P = 5, Ptemp = 6
!                X = 7, Xtemp = 8
!                AP = 9

kCGVec(1) = kEndDiag
kEndCGVec(1) = kCGVec(1)+lCGVec
Wrk(kCGVec(1):kEndCGVec(1)-1) = Zero
nCGvec = 9

do i=2,nCGVec
  kCGVec(i) = kEndCGVec(i-1)
  kEndCGVec(i) = kCGVec(i)+lCGVec
  Wrk(kCGVec(i):kEndCGVec(i)-1) = Zero
end do

! Calculate initial values for the CG-vectors needing that
! --------------_-----------------------------------------

do i=1,lCGFVec
  Wrk(kCGVec(1)+i-1) = Wrk(kFLagr(1)+i-1)*Wrk(kDiag(1)+i-1)
  Wrk(kCGVec(3)+i-1) = Wrk(kFLagr(1)+i-1)
  Wrk(kCGVec(5)+i-1) = Wrk(kCGVec(1)+i-1)
end do
iOff = lCGFVec
do i=1,lCGOVec
  Wrk(kCGVec(1)+i-1+iOff) = Wrk(kLagr(1)+i-1)*Wrk(kDiag(2)+i-1)
  Wrk(kCGVec(3)+i-1+iOff) = Wrk(kLagr(1)+i-1)
  Wrk(kCGVec(5)+i-1+iOff) = Wrk(kCGVec(1)+i-1+iOff)
end do

do iIter=1,nIter
  Wrk(kCGVec(9):kCGVec(9)+lCGVec-1) = Zero

  ! Calculate A*p
  ! -------------
  do iSym=1,nSym

    ! Open Cholesky vector files.
    ! ---------------------------
    call ChoMP2_OpenF(1,1,iSym)

    iSeed = 7
    LuVVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',2
    call DaName_MF_WA(LuVVec,Fname)

    iSeed = 7
    LuWVec = IsFreeUnit(iSeed)
    write(Fname,'(A4,I1)') 'TMPV',3
    call DaName_MF_WA(LuWVec,Fname)

    nVec = min(NumCho(iSym),maxvalue)
    if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')

    lScr = lWrk-kEndCGVec(9)
    iOff = lCGFVec

    ! Construct The Frozen part of A*P
    ! --------------------------------
    call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,'fvvf',iSym,nVec,Wrk(kCGVec(9)),lCGFVec,Wrk(kCGVec(5)),lCGFVec,One)
    call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,'ovvf',iSym,nVec,Wrk(kCGVec(9)),lCGFVec,Wrk(kCGVec(5)+iOff),lCGOVec,One)
    ! Construct The Occupied part of A*P
    ! ----------------------------------
    call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,'fvvo',iSym,nVec,Wrk(kCGVec(9)+iOff),lCGOVec,Wrk(kCGVec(5)),lCGFVec,One)
    call ChoMP2g_Constrap(irc,Wrk(kEndCGVec(9)),lScr,'ovvo',iSym,nVec,Wrk(kCGVec(9)+iOff),lCGOVec,Wrk(kCGVec(5)+iOff),lCGOVec,One)

    call ChoMP2_OpenF(iClos,1,iSym)
    call DaClos(LuVVec)
    call DaClos(LuWVec)

  end do

  index1 = 0
  do iSym1=1,nSym
    nA = nVir(iSym1)
    nI = nFro(iSym1)
    do iI=1,nI
      Ei = EFro(iFro(iSym1)+iI)
      do iA=1,nA
        indx = iA-1+(iI-1)*nA+index1
        Ea = EVir(iVir(iSym1)+iA)
        Wrk(kCGVec(9)+indx) = Wrk(kCGVec(9)+indx)+(Ea-Ei)*Wrk(kCGVec(5)+indx)
      end do
    end do
    nI = nOcc(iSym1)
    do iI=1,nI
      Ei = EOcc(iOcc(iSym1)+iI)
      do iA=1,nA
        indx = iA-1+(iI-1)*nA+index1
        Ea = EVir(iVir(iSym1)+iA)
        Wrk(kCGVec(9)+indx+iOff) = Wrk(kCGVec(9)+indx+iOff)+(Ea-Ei)*Wrk(kCGVec(5)+indx+iOff)
      end do
    end do
    index1 = index1+(nFro(iSym1)+nOcc(iSym1))*nVir(iSym1)
  end do

  call Conj_Grad(Done,lCGVec,Wrk(kDiag(1)),Wrk(kCGVec(7)),Wrk(kCGVec(8)),Wrk(kCGVec(3)),Wrk(kCGVec(4)),Wrk(kCGVec(5)), &
                 Wrk(kCGVec(6)),Wrk(kCGVec(1)),Wrk(kCGVec(2)),Wrk(kCGVec(9)),Eps,Res)
  if (Done) exit
end do
if (.not. Done) then
  write(u6,*) '***************WARNING************************'
  write(u6,*) ''
  write(u6,*) 'Too many iterations, this is what you get after:'
  write(u6,*) nIter,' Iterations'
  write(u6,*) 'The residual is',res,'and not',Eps
  write(u6,*) '**********************************************'
end if

! Construct MP2 Density contribution from parts of the matrix
! -----------------------------------------------------------
iSymOffOV = 0
do iSym=1,nSym
  do i=1,nOcc(iSym)
    iOrb = nFro(iSym)+i
    do j=1,nFro(iSym)
      jOrb = j
      MP2D(iSym)%A(jOrb,iOrb) = Wrk(kPiK(iSym)+i-1+nOcc(iSym)*(j-1))
      MP2D(iSym)%A(iOrb,jOrb) = MP2D(iSym)%A(jOrb,iOrb)
    end do
  end do

  do i=1,nOcc(iSym)
    iOrb = nFro(iSym)+i
    do j=1,nOcc(iSym)
      jOrb = nFro(iSym)+j
      MP2D(iSym)%A(jOrb,iOrb) = Wrk(kPij(iSym)+j-1+nOcc(iSym)*(i-1))
    end do
  end do

  do iI=1,nFro(iSym)
    do iA=1,nVir(iSym)
      Wrk(kPaK(iSym)+iA-1+nVir(iSym)*(iI-1)) = Wrk(kCGVec(7)+iA-1+nVir(iSym)*(iI-1)+iSymOffOV)
    end do
  end do
  do iI=1,nOcc(iSym)
    iOrbI = nFro(iSym)+iI
    do iA=1,nVir(iSym)
      Wrk(kPai(iSym)+iA-1+nVir(iSym)*(iI-1)) = Wrk(kCGVec(7)+iA-1+nVir(iSym)*(iOrbI-1)+iSymOffOV)
    end do
  end do

  do i=1,nFro(iSym)+nOcc(iSym)
    iOrb = i
    do j=1,nVir(iSym)
      jOrb = nFro(iSym)+nOcc(iSym)+j
      MP2D(iSym)%A(jOrb,iOrb) = Wrk(kCGVec(7)+j-1+(nVir(iSym))*(i-1)+iSymOffOV)
      MP2D(iSym)%A(iOrb,jOrb) = MP2D(iSym)%A(jOrb,iOrb)

    end do
  end do

  do i=1,nVir(iSym)
    iOrb = nFro(iSym)+nOcc(iSym)+i
    do j=1,nVir(iSym)
      jOrb = nFro(iSym)+nOcc(iSym)+j
      MP2D(iSym)%A(jOrb,iOrb) = Wrk(kPab(iSym)+j-1+nVir(iSym)*(i-1))
    end do
  end do
  iSymOffOV = iSymOffOV+(nOcc(iSym)+nFro(iSym))*nVir(iSym)
end do

! Add type (II) terms to W-matrix
! -------------------------------

do iSym=1,nsym
  do iI=1,nOcc(iSym)
    Ei = EOcc(iOcc(iSym)+iI)
    do iJ=1,nFro(iSym)
      Ej = EFro(iFro(iSym)+iJ)
      Wrk(kWiK(iSym)+iI-1+nOcc(iSym)*(iJ-1)) = Wrk(kWiK(iSym)+iI-1+nOcc(iSym)*(iJ-1))-Wrk(kPiK(iSym)+iI-1+nOcc(iSym)*(iJ-1))*(Ei+Ej)
    end do
  end do
  do iI=1,nOcc(iSym)
    Ei = EOcc(iOcc(iSym)+iI)
    do iJ=1,nOcc(iSym)
      Ej = EOcc(iOcc(iSym)+iJ)
      Wrk(kWij(iSym)+iJ-1+nOcc(iSym)*(iI-1)) = Wrk(kWij(iSym)+iJ-1+nOcc(iSym)*(iI-1))- &
                                               Half*Wrk(kPij(iSym)+iJ-1+nOcc(iSym)*(iI-1))*(Ei+Ej)
    end do
  end do
  do iA=1,nVir(iSym)
    Ea = EVir(iVir(iSym)+iA)
    do iB=1,nVir(iSym)
      Eb = EVir(iVir(iSym)+iB)
      Wrk(kWab(iSym)+iB-1+nVir(iSym)*(iA-1)) = Wrk(kWab(iSym)+iB-1+nVir(iSym)*(iA-1))- &
                                               Half*Wrk(kPab(iSym)+iB-1+nVir(iSym)*(iA-1))*(Ea+Eb)
    end do
  end do
  do iI=1,nFro(iSym)
    Ei = Efro(iFro(iSym)+iI)
    iOrbI = iI
    do iA=1,nVir(iSym)
      Wrk(kWak(iSym)+iA-1+nVir(iSym)*(iI-1)) = Wrk(kWak(iSym)+iA-1+nVir(iSym)*(iI-1))-Two*Wrk(kPaK(iSym)+iA-1+nVir(iSym)*(iI-1))*Ei
    end do
  end do
  do iI=1,nOcc(iSym)
    Ei = EOcc(iOcc(iSym)+iI)
    iOrbI = nFro(iSYm)+iI
    do iA=1,nVir(iSym)
      Wrk(kWai(iSym)+iA-1+nVir(iSym)*(iI-1)) = Wrk(kWai(iSym)+iA-1+nVir(iSym)*(iI-1))-Two*Wrk(kPai(iSym)+iA-1+nVir(iSym)*(iI-1))*Ei
    end do
  end do

end do

iSym = 1

! Add type (III) terms to W-matrix
! -------------------------------

! Open Cholesky vector files.
! ---------------------------
call ChoMP2_OpenF(1,1,iSym)

nVec = min(NumCho(iSym),maxvalue)
if (nVec < 1) call SysAbendMsg(SecNam,'Insufficient memory','[1]')
nBatL = (NumCho(iSym)-1)/nVec+1

! Allocate memory for U-vector
! ----------------------------

lU = nVec
kU = kEndDiag
kEndU = kU+lU

! Allocate memory for Lia-vector and LIa-vector
! ---------------------------------------------

lLfa = nMoMo(iSym,iVecFV)*nVec
kLfa = kEndU
kEndLfa = kLfa+lLfa

lLia = nMoMo(iSym,iVecOV)*nVec
kLia = kEndLfa
kEndLia = kLia+lLia

lLij = nMoMo(iSym,iVecOO)*nVec
kLij = kEndLia
kEndLij = kLij+lLij

lLiK = nMoMo(iSym,iVecOF)*nVec
kLiK = kEndLij
kEndLiK = kLiK+lLiK

lLKi = nMoMo(iSym,iVecFO)*nVec
kLKi = kEndLiK
kEndLKi = kLKi+lLKi

lLab = nMoMo(iSym,iVecVV)*nVec
kLab = kEndLKi
kEndLab = kLab+lLab

lLJK = nMoMo(iSym,iVecFF)*nVec
kLJK = kEndLab
kEndLJK = kLJK+lLJK

do iBat=1,nBatL
  if (iBat == nBatL) then
    NumVec = NumCho(iSym)-nVec*(nBatL-1)
  else
    NumVec = nVec
  end if
  iVec = nVec*(iBat-1)+1
  Wrk(kU:kU+lU-1) = Zero

  ! Read Lij^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecOO)*NumVec
    iAdr = 1+nMoMo(iSym,iVecOO)*(iVec-1)+iAdrOff(iSym,iVecOO)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,iAdr)
  end if

  ! Read LiK^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecOF)*NumVec
    iAdr = 1+nMoMo(iSym,iVecOF)*(iVec-1)+iAdrOff(iSym,iVecOF)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,iAdr)
  end if

  ! Read LKi^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecFO)*NumVec
    iAdr = 1+nMoMo(iSym,iVecFO)*(iVec-1)+iAdrOff(iSym,iVecFO)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKi),lTot,iAdr)
  end if

  ! Read LIK^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecFF)*NumVec
    iAdr = 1+nMoMo(iSym,iVecFF)*(iVec-1)+iAdrOff(iSym,iVecFF)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLJK),lTot,iAdr)
  end if

  ! Read Lia^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecOV)*NumVec
    iAdr = 1+nMoMo(iSym,iVecOV)*(iVec-1)+iAdrOff(iSym,iVecOV)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,iAdr)
  end if

  ! Read Lfa^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecFV)*NumVec
    iAdr = 1+nMoMo(iSym,iVecFV)*(iVec-1)+iAdrOff(iSym,iVecFV)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLfa),lTot,iAdr)
  end if

  ! Read Lab^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecVV)*NumVec
    iAdr = 1+nMoMo(iSym,iVecVV)*(iVec-1)+iAdrOff(iSym,iVecVV)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLab),lTot,iAdr)
  end if

  ! Construct U^J intermediate vectors
  ! ----------------------------------
  if (NumVec*nMoMo(iSym,iVecOO) /= 0) call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecOO),One,Wrk(kPij(iSym)),1,Wrk(kLij), &
                                                  nMoMo(iSym,iVecOO),One,Wrk(kU),1)

  if (NumVec*nMoMo(iSym,iVecFO) /= 0) call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecFO),Two,Wrk(kPiK(iSym)),1,Wrk(kLKi), &
                                                  nMoMo(iSym,iVecFO),One,Wrk(kU),1)

  if (NumVec*nMoMo(iSym,iVecOV) /= 0) call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecOV),Two,Wrk(kPai(iSym)),1,Wrk(kLia), &
                                                  nMoMo(iSym,iVecOV),One,Wrk(kU),1)

  if (NumVec*nMoMo(iSym,iVecFV) /= 0) call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecFV),Two,Wrk(kPaK(iSym)),1,Wrk(kLfa), &
                                                  nMoMo(iSym,iVecFV),One,Wrk(kU),1)

  if (NumVec*nMoMo(iSym,iVecVV) /= 0) call dGemm_('N','N',1,NumVec,nMoMo(iSym,iVecVV),One,Wrk(kPab(iSym)),1,Wrk(kLab), &
                                                  nMoMo(iSym,iVecVV),One,Wrk(kU),1)

  ! Construct contribution to Wij
  ! -----------------------------

  if (nMoMo(iSym,iVecOO) /= 0) call dGemm_('N','N',nMoMo(iSym,iVecOO),1,NumVec,-Two,Wrk(kLij),nMoMo(iSym,iVecOO),Wrk(kU),NumVec, &
                                           One,Wrk(kWij(iSym)),nMoMo(iSym,iVecOO))

  ! Construct Contribution to WiK
  ! -----------------------------

  if (nMoMo(iSym,iVecOF) /= 0) call dGemm_('N','N',nMoMo(iSym,iVecOF),1,NumVec,-Four,Wrk(kLKi),nMoMo(iSym,iVecOF),Wrk(kU),NumVec, &
                                           One,Wrk(kWiK(iSym)),nMoMo(iSym,iVecOF))

  ! Construct Contribution to WJK
  ! -----------------------------

  if (nMoMo(iSym,iVecFF) /= 0) call dGemm_('N','N',nMoMo(iSym,iVecFF),1,NumVec,-Two,Wrk(kLJK),nMoMo(iSym,iVecFF),Wrk(kU),NumVec, &
                                           One,Wrk(kWJK(iSym)),nMoMo(iSym,iVecFF))

end do

iVecFF = 1
iVecOF = 2
iVecVF = 3
iVecFO = 4
iVecOO = 5
iVecVO = 6
iVecFV = 7
iVecOV = 8
iVecVV = 9

! Allocate memory for Lia-vector and LIa-vector
! ---------------------------------------------

lLJK = nMoMo(iSym,iVecFF)*nVec
kLJK = kEndDiag
kEndLJK = kLJK+lLJK

lLKi = nMoMo(iSym,iVecOF)*nVec
kLKi = kEndLJK
kEndLKi = kLKi+lLKi

lLKa = nMoMo(iSym,iVecVF)*nVec
kLKa = kEndLKi
kEndLKa = kLKa+lLKa

lLiK = nMoMo(iSym,iVecFO)*nVec
kLiK = kEndLKa
kEndLiK = kLiK+lLiK

lLij = nMoMo(iSym,iVecOO)*nVec
kLij = kEndLiK
kEndLij = kLij+lLij

lLia = nMoMo(iSym,iVecVO)*nVec
kLia = kEndLij
kEndLia = kLia+lLia

lLip = nOccAll(iSym)*nOrb(iSym)*nVec
kLip = kEndLia
kEndLip = kLip+lLip

lVip = nOccAll(iSym)*nOrb(iSym)*nVec
kVip = kEndLip
kEndVip = kVip+lVip

lWij2 = nOccAll(iSym)*nOccAll(iSym)
kWij2 = kEndVip
Wrk(kWij2:kWij2+lWij2-1) = Zero

do iBat=1,nBatL
  if (iBat == nBatL) then
    NumVec = NumCho(iSym)-nVec*(nBatL-1)
  else
    NumVec = nVec
  end if
  iVec = nVec*(iBat-1)+1

  ! Read Lpq^J-vectors from disk
  ! ----------------------------

  if (NumVec > 0) then
    iOpt = 2
    lTot = nMoMo(iSym,iVecFF)*NumVec
    iAdr = 1+nMoMo(iSym,iVecFF)*(iVec-1)+iAdrOff(iSym,iVecFF)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLJK),lTot,iAdr)

    lTot = nMoMo(iSym,iVecOF)*NumVec
    iAdr = 1+nMoMo(iSym,iVecOF)*(iVec-1)+iAdrOff(iSym,iVecOF)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKi),lTot,iAdr)

    lTot = nMoMo(iSym,iVecVF)*NumVec
    iAdr = 1+nMoMo(iSym,iVecVF)*(iVec-1)+iAdrOff(iSym,iVecVF)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLKa),lTot,iAdr)

    lTot = nMoMo(iSym,iVecFO)*NumVec
    iAdr = 1+nMoMo(iSym,iVecFO)*(iVec-1)+iAdrOff(iSym,iVecFO)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLiK),lTot,iAdr)

    lTot = nMoMo(iSym,iVecOO)*NumVec
    iAdr = 1+nMoMo(iSym,iVecOO)*(iVec-1)+iAdrOff(iSym,iVecOO)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLij),lTot,iAdr)

    lTot = nMoMo(iSym,iVecVO)*NumVec
    iAdr = 1+nMoMo(iSym,iVecVO)*(iVec-1)+iAdrOff(iSym,iVecVO)
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Wrk(kLia),lTot,iAdr)
  end if

  ! Put the Lpq-vectors together to one large vector
  ! ------------------------------------------------

  iOff2 = kLip
  do iVec1=1,NumVec
    do i=1,nFro(iSym)
      iOff1 = kLJK+nMoMo(iSym,iVecFF)*(iVec1-1)+(i-1)*nFro(iSym)
      Wrk(iOff2:iOff2+nFro(iSym)-1) = Wrk(iOff1:iOff1+nFro(iSym)-1)
      iOff2 = iOff2+nFro(iSym)
      iOff1 = kLKi+nMoMo(iSym,iVecOF)*(iVec1-1)+(i-1)*nOcc(iSym)
      Wrk(iOff2:iOff2+nOcc(iSym)-1) = Wrk(iOff1:iOff1+nOcc(iSym)-1)
      iOff2 = iOff2+nOcc(iSym)
      iOff1 = kLKa+nMoMo(iSym,iVecVF)*(iVec1-1)+(i-1)*nVir(iSym)
      Wrk(iOff2:iOff2+nVir(iSym)-1) = Wrk(iOff1:iOff1+nVir(iSym)-1)
      iOff2 = iOff2+nVir(iSym)
    end do
    do i=1,nOcc(iSym)
      iOff1 = kLiK+nMoMo(iSym,iVecFO)*(iVec1-1)+(i-1)*nFro(iSym)
      Wrk(iOff2:iOff2+nFro(iSym)-1) = Wrk(iOff1:iOff1+nFro(iSym)-1)
      iOff2 = iOff2+nFro(iSym)
      iOff1 = kLij+nMoMo(iSym,iVecOO)*(iVec1-1)+(i-1)*nOcc(iSym)
      Wrk(iOff2:iOff2+nOcc(iSym)-1) = Wrk(iOff1:iOff1+nOcc(iSym)-1)
      iOff2 = iOff2+nOcc(iSym)
      iOff1 = kLia+nMoMo(iSym,iVecVO)*(iVec1-1)+(i-1)*nVir(iSym)
      Wrk(iOff2:iOff2+nVir(iSym)-1) = Wrk(iOff1:iOff1+nVir(iSym)-1)
      iOff2 = iOff2+nVir(iSym)
    end do
  end do

  do iVec1=1,NumVec
    iOff = nOrb(iSym)*nOccAll(iSym)*(iVec1-1)
    call dGemm_('N','N',nOrb(iSym),nOccAll(iSym),nOrb(iSym),One,MP2D(iSym)%A,nOrb(iSym),Wrk(kLip+iOff),nOrb(iSym),Zero, &
                Wrk(kVip+iOff),nOrb(iSym))
  end do

  ! Construct exchange contribution to Wij
  ! --------------------------------------

  do iVec1=1,NumVec
    iOff = nOrb(iSym)*nOccAll(iSym)*(iVec1-1)
    call dGemm_('T','N',nOccAll(iSym),nOccAll(iSym),nOrb(iSym),One,Wrk(kLip+iOff),nOrb(iSym),Wrk(kVip+iOff),nOrb(iSym),One, &
                Wrk(kWij2),nOccAll(iSym))
  end do
end do

call ChoMP2_OpenF(iClos,1,iSym)

! Add SCF-density to MP2-density contribution
! -------------------------------------------
do iSym1=1,nSym
  do i=1,nOcc(iSym1)+nFro(iSym1)
    MP2D(iSym1)%A(i,i) = MP2D(iSym1)%A(i,i)+Two
  end do
end do

! Construct Mp2 + SCF W-density
! -----------------------------

do i=1,nFro(iSym)
  iOrb = i
  do j=1,nFro(iSym)
    jOrb = j
    term = Zero
    if (i == j) term = Efro(iFro(iSym)+i)
    MP2W(iSym)%A(jOrb,iOrb) = -Wrk(kWJK(iSym)+j-1+nFro(iSym)*(i-1))+Two*term-Wrk(kWij2+j-1+nOccAll(iSym)*(i-1))
  end do
end do
do i=1,nOcc(iSym)
  iOrb = nFro(iSym)+i
  do j=1,nFro(iSym)
    jOrb = j
    MP2W(iSym)%A(jOrb,iOrb) = -Half*Wrk(kWiK(iSym)+i-1+nOcc(iSym)*(j-1))-Wrk(kWij2+jOrb-1+nOccAll(iSym)*(iOrb-1))
    MP2W(iSym)%A(iOrb,jOrb) = MP2W(iSym)%A(jOrb,iOrb)
  end do
end do
do i=1,nOcc(iSym)
  iOrb = nFro(iSym)+i
  do j=1,nOcc(iSym)
    jOrb = nFro(iSym)+j
    term = Zero
    if (i == j) term = EOcc(iOcc(iSym)+i)
    MP2W(iSym)%A(jOrb,iOrb) = -Wrk(kWij(iSym)+j-1+nOcc(iSym)*(i-1))+Two*term-Wrk(kWij2+jOrb-1+nOccAll(iSym)*(iOrb-1))
  end do
end do

do i=1,nFro(iSym)
  iOrb = i
  do j=1,nVir(iSym)
    jOrb = nFro(iSym)+nOcc(iSym)+j
    MP2W(iSym)%A(jOrb,iOrb) = -Half*Wrk(kWaK(iSym)+j-1+nVir(iSym)*(i-1))
    MP2W(iSym)%A(iOrb,jOrb) = MP2W(iSym)%A(jOrb,iOrb)
  end do
end do

do i=1,nOcc(iSym)
  iOrb = nFro(iSym)+i
  do j=1,nVir(iSym)
    jOrb = nFro(iSym)+nOcc(iSym)+j
    MP2W(iSym)%A(jOrb,iOrb) = -Half*Wrk(kWai(iSym)+j-1+nVir(iSym)*(i-1))
    MP2W(iSym)%A(iOrb,jOrb) = MP2W(iSym)%A(jOrb,iOrb)
  end do
end do

do i=1,nVir(iSym)
  iOrb = nFro(iSym)+nOcc(iSym)+i
  do j=1,nVir(iSym)
    jOrb = nFro(iSym)+nOcc(iSym)+j
    MP2W(iSym)%A(jOrb,iOrb) = -Wrk(kWab(iSym)+j-1+nVir(iSym)*(i-1))
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'Full One-electron Density Matrix'
do iSym=1,nSym
  call RecPrt('MP2D',' ',MP2D(iSym)%A,nOrb(iSym),nOrb(iSym))
end do
write(u6,*) 'Full One-electron energy weighted D Matrix'
do iSym=1,nSym
  call RecPrt('MP2W',' ',MP2W(iSym)%A,nOrb(iSym),nOrb(iSym))
end do
#endif

return

end subroutine ChoMP2g_density2
