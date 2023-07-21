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

subroutine ChoMP2g_ConstrAP(irc,Scr,lScr,type,iSym,nVec,Ap,lAp,Dens,lDens,factor)

use ChoMP2g

implicit real*8(a-h,o-z)
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "cholesky.fh"
#include "choorb.fh"
real*8 Scr(lScr), Ap(lAp), Dens(lDens)
character type(1:4)
logical DoX, DoY, DoZ
integer iOffp(8), iOffAp(8)
integer nOrb1(8,4), iOrbType(4)
character*8 ThisNm
character*16 SecNam
parameter(SecNam='ChoMP2g_ConstrAP',ThisNm='ConstrAP')
! Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1
!************************************

! Setup some lengths
! ------------------
nTypes = 4

iTypL = 1

! Setup lengths of needed cholesky/intermediate vectors
! -----------------------------------------------------

! Read the char-array Type = 'pqri'
! ---------------------------------
do i=1,nTypes
  if (type(i) == 'f') then
    iOrbType(i) = 1
    do iSym1=1,nSym
      nOrb1(iSym1,i) = nFro(iSym1)
    end do
  else if (type(i) == 'o') then
    iOrbType(i) = 2
    do iSym1=1,nSym
      nOrb1(iSym1,i) = nOcc(iSym1)
    end do
  else if (type(i) == 'v') then
    iOrbType(i) = 3
    do iSym1=1,nSym
      nOrb1(iSym1,i) = nVir(iSym1)
    end do
  else
    write(6,*) 'Forbidden Type pqrs in',SecNam
    irc = -1
    return
  end if
end do

! Set vector type index
! ---------------------
iPQ = 3*(iOrbType(1)-1)+iOrbType(2)
iIR = 3*(iOrbType(4)-1)+iOrbType(3)
iIQ = 3*(iOrbType(4)-1)+iOrbType(2)
iRP = 3*(iOrbType(3)-1)+iOrbType(1)
iIP = 3*(iOrbType(4)-1)+iOrbType(1)
iRQ = 3*(iOrbType(3)-1)+iOrbType(2)

! Set vector lengths
! ------------------
nPQ = nMoMo(iSym,iPQ)
nIR = nMoMo(iSym,iIR)
nIQ = nMoMo(iSym,iIQ)
nRP = nMoMo(iSym,iRP)
nIP = nMoMo(iSym,iIP)
nRQ = nMoMo(iSym,iRQ)

! Set address-offset for cholesky vector address
! ----------------------------------------------
iAdrLpq = iAdrOff(iSym,iPQ)
iAdrLir = iAdrOff(iSym,iIR)
iAdrLiq = iAdrOff(iSym,iIQ)
iAdrLrp = iAdrOff(iSym,iRP)
iAdrLip = iAdrOff(iSym,iIP)
iAdrLrq = iAdrOff(iSym,iRQ)

! Set

#ifdef _DEBUGPRINT_
write(6,*) 'iPQ',iPQ
write(6,*) 'iIQ',iIQ
write(6,*) 'iIR',iIR
write(6,*) 'iIP',iIP
write(6,*) 'iRQ',iRQ
write(6,*) 'iRP',iRP

write(6,*) 'nPQ',nPQ
write(6,*) 'nIQ',nIQ
write(6,*) 'nIR',nIR
write(6,*) 'nIP',nIP
write(6,*) 'nRQ',nRQ
write(6,*) 'nRP',nRP
#endif
! Decide what type of intermediate vectors is needed
! --------------------------------------------------
DoX = .true.
DoY = .true.
DoZ = .true.

if (iSym /= 1) DoX = .false.
if ((nPQ == 0) .or. (nIR == 0)) DoX = .false.
if ((nRP == 0) .or. (nIQ == 0)) DoY = .false.
if ((nIP == 0) .or. (nRQ == 0)) DoZ = .false.
if (type(1) == type(2)) DoZ = .false.

if (DoY .or. DoZ) then
  iOffp(1) = 0
  iOffap(1) = 0
  do iSym1=2,nSym
    iOffp(iSym1) = iOffp(iSym1-1)+nOrb1(iSym1-1,1)*nOrb1(iSym1-1,2)
    iOffap(iSym1) = iOffap(iSym1-1)+nOrb1(iSym1-1,3)*nOrb1(iSym1-1,4)
  end do
end if

!     Allocate memory for L-vectors
!     ------------------------------

!        Lri:
lLir = nIR*nVec
kLir = 1
kEndLir = kLir+lLir

!        Lpq:
lLpq = nPQ*nVec
kLpq = kEndLir
kEndLpq = kLpq+lLpq

!        Lrp:
lLrp = nRP*nVec
kLrp = kEndLpq
kEndLrp = kLrp+lLrp

! Lrq:
lLrq = nRQ*nVec
kLrq = kEndLrp
kEndLrq = kLrq+lLrq

! Lip:
lLip = nIP*nVec
kLip = kEndLrq
kEndLip = kLip+lLip

! Liq:
lLiq = nIQ*nVec
kLiq = kEndLip
kEndLiq = kLiq+lLiq

! Allocate memory for Intermediate Vectors
! ----------------------------------------
! X-vector
lX = NumCho(iSym)
kX = kEndLiq
kEndX = kX+lX
! Y-vector
lY = nVec*nIP
kY = kEndX
kEndY = kY+lY
! Z-vector
kZ = kEndY

nBatL = (NumCho(iSym)-1)/nVec+1

! Construction of Intermediate vectors
! ------------------------------------
do iBat=1,nBatL
  if (iBat == nBatL) then
    NumVec = NumCho(iSym)-nVec*(nBatL-1)
  else
    NumVec = nVec
  end if
  iVec = nVec*(iBat-1)+1

  ! Read Lpq-vectors
  ! ----------------
  if (DoX) then
    iOpt = 2
    lTot = nPQ*NumVec
    iAdr = nPQ*(iVec-1)+1+iAdrLpq
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLpq),lTot,iAdr)
  end if

  ! Read Liq-vectors
  ! ----------------
  if (DoY) then
    iOpt = 2
    lTot = nIQ*NumVec
    iAdr = nIQ*(iVec-1)+1+iAdrLiq
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLiq),lTot,iAdr)
  end if

  if (DoZ) then
    iOpt = 2
    lTot = nIP*NumVec
    iAdr = nIP*(iVec-1)+1+iAdrLip
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLip),lTot,iAdr)
  end if

  ! Construct X-vector
  ! ------------------
  if (DoX) then
    nRow = 1
    call dGemm_('T','N',nRow,NumVec,nPQ,1.0d0,Dens(1),nPQ,Scr(kLpq),nPQ,0.0d0,Scr(kX+iVec-1),nRow)
  end if

  ! Construct Y-vector
  ! ------------------
  if (DoY) then
    do iVec1=1,NumVec
      iOffL1 = 0
      iOffY1 = 0
      do iSymI=1,nSym
        iSymP = MulD2h(iSym,iSymI)
        iSymQ = iSymP
        nI = nOrb1(iSymI,4)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymQ,2)
        nR = nOrb1(iSymI,3)
        if (nQ*nP*nI*nR == 0) Go To 101
        iOffL = (iVec1-1)*nIQ+iOffL1
        iOffY = (iVec1-1)*nIP+iOffY1
        call dGemm_('T','N',nP,nI,nQ,1.0d0,Dens(1+iOffP(iSymP)),nQ,Scr(kLiq+iOffL),nQ,0.0d0,Scr(kY+iOffY),nP)
101     continue
        iOffL1 = iOffL1+nQ*nI
        iOffY1 = iOffY1+nP*nI
      end do
    end do
    if (nBatL /= 1) then
      iOpt = 1
      lTot = nIP*NumVec
      iAdr = nIP*(iVec-1)+1
      call dDaFile(LuVVec,iOpt,Scr(kY),lTot,iAdr)
    end if
  end if
  ! Construct Z-vectors
  ! -------------------
  if (DoZ) then
    do iVec1=1,NumVec
      iOffL1 = 0
      iOffZ1 = 0
      do iSymI=1,nSym
        iSymP = MulD2h(iSym,iSymI)
        iSymQ = iSymP
        nI = nOrb1(iSymI,4)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymQ,2)
        nR = nOrb1(iSymI,3)
        if (nQ*nP*nI*nR == 0) Go To 102
        iOffL = (iVec1-1)*nIP+iOffL1
        iOffZ = (iVec1-1)*nIQ+iOffZ1
        call dGemm_('N','N',nQ,nI,nP,1.0d0,Dens(1+iOffP(iSymQ)),nQ,Scr(kLip+iOffL),nP,0.0d0,Scr(kZ+iOffZ),nQ)
102     continue
        iOffL1 = iOffL1+nP*nI
        iOffZ1 = iOffZ1+nQ*nI
      end do
    end do
    if (nBatL /= 1) then
      iOpt = 1
      lTot = nIQ*NumVec
      iAdr = nIQ*(iVec-1)+1
      call dDaFile(LuWVec,iOpt,Scr(kZ),lTot,iAdr)
    end if

  end if
end do

! Contracting to produce the final contribution to A*p
! ----------------------------------------------------
do iBat=1,nBatL
  if (iBat == nBatL) then
    NumVec = NumCho(iSym)-nVec*(nBatL-1)
  else
    NumVec = nVec
  end if
  iVec = nVec*(iBat-1)+1

  ! Read Lir-vectors
  ! ----------------
  if (DoX) then
    iOpt = 2
    lTot = nIR*NumVec
    iAdr = nIR*(iVec-1)+1+iAdrLir
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLir),lTot,iAdr)
  end if

  ! Read Lrp-vectors
  ! ----------------
  if (DoY) then
    iOpt = 2
    lTot = nRP*NumVec
    iAdr = nRP*(iVec-1)+1+iAdrLrp
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLrp),lTot,iAdr)
  end if

  ! Read Lrq-vectors
  ! ----------------
  if (DoZ) then
    iOpt = 2
    lTot = nRQ*NumVec
    iAdr = nRQ*(iVec-1)+1+iAdrLrq
    call dDaFile(lUnit_F(iSym,iTypL),iOpt,Scr(kLrq),lTot,iAdr)
  end if

  if (DoX) then
    iCol = 1
    call dGemm_('N','N',nIR,iCol,NumVec,4.0d0*factor,Scr(kLir),nIR,Scr(kX+iVec-1),NumVec,1.0d0,Ap(1),nIR)
  end if

  if (DoY) then
    ! Read Y-vector (if they are stored on disk and not in memory)
    ! ------------------------------------------------------------
    if (nBatL /= 1) then
      iOpt = 2
      lTot = nIP*NumVec
      iAdr = nIP*(iVec-1)+1
      call dDaFile(LuVVec,iOpt,Scr(kY),lTot,iAdr)
    end if

    yfactor = 1.0d0*factor
    if (.not. DoZ) yfactor = yfactor*2.0d0
    do iVec1=1,NumVec
      iOffL1 = 0
      iOffY1 = 0
      do iSymI=1,nSym
        iSymP = MulD2h(iSym,iSymI)
        iSymR = iSymI
        nI = nOrb1(iSymI,4)
        nR = nOrb1(iSymR,3)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymP,2)

        if (nI*nR*nP*nQ == 0) Go To 201
        iOffL = (iVec1-1)*nRP+iOffL1
        iOffY = (iVec1-1)*nIP+iOffY1
        call dGemm_('T','N',nR,nI,nP,-1.0d0*yfactor,Scr(kLrp+iOffL),nP,Scr(kY+iOffY),nP,1.0d0,Ap(1+iOffAP(iSymI)),nR)
201     continue
        iOffL1 = iOffL1+nR*nP
        iOffY1 = iOffY1+nI*nP
      end do
    end do

  end if

  if (DoZ) then
    ! Read Z-vector (if they are stored on disk and not in memory)
    ! ------------------------------------------------------------
    if (nBatL /= 1) then
      iOpt = 2
      lTot = nIQ*NumVec
      iAdr = nIQ*(iVec-1)+1
      call dDaFile(LuWVec,iOpt,Scr(kZ),lTot,iAdr)
    end if

    do iVec1=1,NumVec
      iOffL1 = 0
      iOffZ1 = 0
      do iSymI=1,nSym
        iSymQ = MulD2h(iSym,iSymI)
        iSymR = iSymI
        nI = nOrb1(iSymI,4)
        nR = nOrb1(iSymR,3)
        nQ = nOrb1(iSymQ,2)
        nP = nOrb1(iSymQ,1)
        if (nI*nR*nQ*nP == 0) Go To 202
        iOffL = (iVec1-1)*nRQ+iOffL1
        iOffZ = (iVec1-1)*nIQ+iOffZ1
        call dGemm_('T','N',nR,nI,nQ,-1.0d0*factor,Scr(kLrq+iOffL),nQ,Scr(kZ+iOffZ),nQ,1.0d0,Ap(1+iOffAP(iSymI)),nR)
202     continue
        iOffL1 = iOffL1+nR*nQ
        iOffZ1 = iOffZ1+nI*nQ
      end do
    end do
  end if
end do

return

end subroutine ChoMP2g_ConstrAP
