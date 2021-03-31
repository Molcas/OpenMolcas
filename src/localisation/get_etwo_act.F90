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

subroutine Get_Etwo_act(Dma,Dmb,nBDT,nBas,nSym,Etwo)

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp), intent(in) :: nBDT, nBas(8), nSym
real(kind=wp), intent(in) :: Dma(nBDT), Dmb(nBDT)
real(kind=wp), intent(out) :: Etwo
#include "WrkSpc.fh"
#include "choscf.fh"
#include "choscreen.fh"
#include "chotime.fh"
integer(kind=iwp) :: i, iOff, ipDai, ipDbi, ipDm(2), ipFCNO,ipFLT(2), ipKLT(2), ipPLT, ipPorb(2), ipV, irc, nBB, nIorb(8,2)
real(kind=wp) :: ChFracMem
!character(len=16) :: KSDFT
real(kind=r8), external :: ddot_
!real(kind=wp), external :: Get_ExFac

timings = .false.
Estimate = .false.
REORD = .false.

Update = .true.
DECO = .true.
dmpk = One
dFKmat = Zero
ALGO = 4
NSCREEN = 10

!nDMat = 2
nBB = 0
do i=1,nSym
  !nForb(i,1) = 0
  !nForb(i,2) = 0
  nBB = nBB+nBas(i)**2
end do
!call Get_cArray('DFT functional',KSDFT,16)
!ExFac = Get_ExFac(KSDFT)
!FactXI = ExFac
!FactXI = One  ! always HF energy
call GetMem('PLTc','Allo','Real',ipPLT,nBDT)
call dcopy_(nBDT,Dma,1,Work(ipPLT),1)
call daxpy_(nBDT,One,Dmb,1,Work(ipPLT),1)

call GetMem('ChMc','Allo','Real',ipPorb(1),2*nBB)
ipPorb(2) = ipPorb(1)+nBB
call GetMem('DSQc','Allo','Real',ipDm(1),2*nBB)
ipDm(2) = ipDm(1)+nBB
call UnFold(Dma,nBDT,Work(ipDm(1)),nBB,nSym,nBas)
call UnFold(Dmb,nBDT,Work(ipDm(2)),nBB,nSym,nBas)
iOff = 0
do i=1,nSym
  ipV = ipPorb(1)+iOff
  ipDai = ipDm(1)+iOff
  call CD_InCore(Work(ipDai),nBas(i),Work(ipV),nBas(i),nIorb(i,1),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
    call Abend()
  end if
  ipV = ipPorb(2)+iOff
  ipDbi = ipDm(2)+iOff
  call CD_InCore(Work(ipDbi),nBas(i),Work(ipV),nBas(i),nIorb(i,2),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Beta density. Sym= ',i,'   rc= ',irc
    call Abend()
  end if
  iOff = iOff+nBas(i)**2
end do

call GetMem('FCNO','Allo','Real',ipFCNO,2*nBDT)
call FZero(Work(ipFCNO),2*nBDT)
ipFLT(1) = ipFCNO
ipFLT(2) = ipFLT(1)+nBDT
call GetMem('KLTc','Allo','Real',ipKLT(1),2*nBDT)
call FZero(Work(ipKLT(1)),2*nBDT)
ipKLT(2) = ipKLT(1)+nBDT

call Cho_X_init(irc,ChFracMem)
if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_init.')
  call Abend()
end if

! BIGOT FIXME
call WarningMessage(2,'There is probably a bug here, ipPLT should have two elements.')
call Abend()
!call CHO_LK_SCF(irc,nDMat,ipFLT,ipKLT,nForb,nIorb,ipPorb,ipPLT,FactXI,nSCReen,dmpk,dFmat)
if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_LK_scf.')
  call Abend()
end if

call Cho_X_Final(irc)
if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_Final.')
  call Abend()
end if

Etwo = Half*(ddot_(nBDT,Dma,1,Work(ipFLT(1)),1)+ddot_(nBDT,Dmb,1,Work(ipFLT(2)),1))

call GetMem('KLTc','Free','Real',ipKLT(1),2*nBDT)
call GetMem('FCNO','Free','Real',ipFCNO,2*nBDT)
call GetMem('DSQc','Free','Real',ipDm(1),2*nBB)
call GetMem('ChMc','Free','Real',ipPorb(1),2*nBB)
call GetMem('PLTc','Free','Real',ipPLT,nBDT)

return

end subroutine Get_Etwo_act
