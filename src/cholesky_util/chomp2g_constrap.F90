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

subroutine ChoMP2g_ConstrAP(irc,Scr,lScr,typ,iSym,nVec,Ap,lAp,Dens,lDens,factor)

use Symmetry_Info, only: Mul
use Cholesky, only: nSym, NumCho
use ChoMP2, only: iAdrOff, lUnit_F, LuVVec, LuWVec, nFro, nMoMo, nOcc, nVir
use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: lScr, iSym, nVec, lAp, lDens
real(kind=wp), intent(out) :: Scr(lScr)
real(kind=wp), intent(inout) :: Ap(lAp)
real(kind=wp), intent(in) :: Dens(lDens), factor
character, intent(in) :: typ(4)
integer(kind=iwp) :: i, iAdr, iAdrLip, iAdrLiq, iAdrLir, iAdrLpq, iAdrLrp, iAdrLrq, iBat, iCol, iIP, iIQ, iIR, iOffAp(8), iOffL, &
                     iOffL1, iOffp(8), iOffY, iOffY1, iOffZ, iOffZ1, iOpt, iOrbType(4), iPQ, iRP, iRQ, iSym1, iSymI, iSymP, iSymQ, &
                     iSymR, iTypL, iVec, iVec1, kEndLip, kEndLiq, kEndLir, kEndLpq, kEndLrp, kEndLrq, kEndX, kEndY, kLip, kLiq, &
                     kLir, kLpq, kLrp, kLrq, kX, kY, kZ, lLip, lLiq, lLir, lLpq, lLrp, lLrq, lTot, lX, lY, nBatL, nI, nIP, nIQ, &
                     nIR, nOrb1(8,4), nP, nPQ, nQ, nR, nRow, nRP, nRQ, nTypes, NumVec
real(kind=wp) :: yfactor
logical(kind=iwp) :: DoX, DoY, DoZ
character(len=*), parameter :: SecNam = 'ChoMP2g_ConstrAP'

! Setup some lengths
! ------------------
nTypes = 4

iTypL = 1

! Setup lengths of needed cholesky/intermediate vectors
! -----------------------------------------------------

! Read the char-array Type = 'pqri'
! ---------------------------------
do i=1,nTypes
  if (typ(i) == 'f') then
    iOrbType(i) = 1
    nOrb1(1:nSym,i) = nFro(1:nSym)
  else if (typ(i) == 'o') then
    iOrbType(i) = 2
    nOrb1(1:nSym,i) = nOcc(1:nSym)
  else if (typ(i) == 'v') then
    iOrbType(i) = 3
    nOrb1(1:nSym,i) = nVir(1:nSym)
  else
    write(u6,*) 'Forbidden Type pqrs in',SecNam
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
write(u6,*) 'iPQ',iPQ
write(u6,*) 'iIQ',iIQ
write(u6,*) 'iIR',iIR
write(u6,*) 'iIP',iIP
write(u6,*) 'iRQ',iRQ
write(u6,*) 'iRP',iRP

write(u6,*) 'nPQ',nPQ
write(u6,*) 'nIQ',nIQ
write(u6,*) 'nIR',nIR
write(u6,*) 'nIP',nIP
write(u6,*) 'nRQ',nRQ
write(u6,*) 'nRP',nRP
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
if (typ(1) == typ(2)) DoZ = .false.

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
    call dGemm_('T','N',nRow,NumVec,nPQ,One,Dens(1),nPQ,Scr(kLpq),nPQ,Zero,Scr(kX+iVec-1),nRow)
  end if

  ! Construct Y-vector
  ! ------------------
  if (DoY) then
    do iVec1=1,NumVec
      iOffL1 = 0
      iOffY1 = 0
      do iSymI=1,nSym
        iSymP = Mul(iSym,iSymI)
        iSymQ = iSymP
        nI = nOrb1(iSymI,4)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymQ,2)
        nR = nOrb1(iSymI,3)
        if (nQ*nP*nI*nR /= 0) then
          iOffL = (iVec1-1)*nIQ+iOffL1
          iOffY = (iVec1-1)*nIP+iOffY1
          call dGemm_('T','N',nP,nI,nQ,One,Dens(1+iOffP(iSymP)),nQ,Scr(kLiq+iOffL),nQ,Zero,Scr(kY+iOffY),nP)
        end if
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
        iSymP = Mul(iSym,iSymI)
        iSymQ = iSymP
        nI = nOrb1(iSymI,4)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymQ,2)
        nR = nOrb1(iSymI,3)
        if (nQ*nP*nI*nR /= 0) then
          iOffL = (iVec1-1)*nIP+iOffL1
          iOffZ = (iVec1-1)*nIQ+iOffZ1
          call dGemm_('N','N',nQ,nI,nP,One,Dens(1+iOffP(iSymQ)),nQ,Scr(kLip+iOffL),nP,Zero,Scr(kZ+iOffZ),nQ)
        end if
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
    call dGemm_('N','N',nIR,iCol,NumVec,Four*factor,Scr(kLir),nIR,Scr(kX+iVec-1),NumVec,One,Ap(1),nIR)
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

    yfactor = factor
    if (.not. DoZ) yfactor = yfactor*Two
    do iVec1=1,NumVec
      iOffL1 = 0
      iOffY1 = 0
      do iSymI=1,nSym
        iSymP = Mul(iSym,iSymI)
        iSymR = iSymI
        nI = nOrb1(iSymI,4)
        nR = nOrb1(iSymR,3)
        nP = nOrb1(iSymP,1)
        nQ = nOrb1(iSymP,2)

        if (nI*nR*nP*nQ /= 0) then
          iOffL = (iVec1-1)*nRP+iOffL1
          iOffY = (iVec1-1)*nIP+iOffY1
          call dGemm_('T','N',nR,nI,nP,-yfactor,Scr(kLrp+iOffL),nP,Scr(kY+iOffY),nP,One,Ap(1+iOffAP(iSymI)),nR)
        end if
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
        iSymQ = Mul(iSym,iSymI)
        iSymR = iSymI
        nI = nOrb1(iSymI,4)
        nR = nOrb1(iSymR,3)
        nQ = nOrb1(iSymQ,2)
        nP = nOrb1(iSymQ,1)
        if (nI*nR*nQ*nP /= 0) then
          iOffL = (iVec1-1)*nRQ+iOffL1
          iOffZ = (iVec1-1)*nIQ+iOffZ1
          call dGemm_('T','N',nR,nI,nQ,-factor,Scr(kLrq+iOffL),nQ,Scr(kZ+iOffZ),nQ,One,Ap(1+iOffAP(iSymI)),nR)
        end if
        iOffL1 = iOffL1+nR*nQ
        iOffZ1 = iOffZ1+nI*nQ
      end do
    end do
  end if
end do

return

end subroutine ChoMP2g_ConstrAP
