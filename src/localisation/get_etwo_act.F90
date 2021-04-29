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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6, r8
use Data_structures, only: DSBA_type, Allocate_DSBA, Deallocate_DSBA

implicit none
integer(kind=iwp), intent(in) :: nBDT, nBas(8), nSym
real(kind=wp), intent(in) :: Dma(nBDT), Dmb(nBDT)
real(kind=wp), intent(out) :: Etwo
#include "choscf.fh"
#include "choscreen.fh"
#include "chotime.fh"
integer(kind=iwp) :: i, iOff, irc, nBB, nForb(8,2), nIorb(8,2)
real(kind=wp) :: ChFracMem, dFmat, FactXI
!character(len=16) :: KSDFT
real(kind=wp), allocatable :: Dm1(:), Dm2(:)
integer(kind=iwp), external :: ip_of_Work
real(kind=r8), external :: ddot_ !, Get_ExFac
type (DSBA_type) :: FLT(2), KLT(2), POrb(2), PLT(2)

timings = .false.
Estimate = .false.
REORD = .false.

Update = .true.
DECO = .true.
dmpk = One
dFKmat = Zero
ALGO = 4
NSCREEN = 10

nForb(:,:) = 0
nBB = 0
do i=1,nSym
  nBB = nBB+nBas(i)**2
end do
!call Get_cArray('DFT functional',KSDFT,16)
!ExFac = Get_ExFac(KSDFT)
!FactXI = ExFac
FactXI = One  ! always HF energy
Call Allocate_DSBA(PLT(1),nBas,nBas,nSym,Case='TRI')
PLT(1)%A0(:) = Dma(:)+Dmb(:)

Call Allocate_DSBA(POrb(1),nBas,nBas,nSym)
Call Allocate_DSBA(POrb(2),nBas,nBas,nSym)
call mma_allocate(Dm1,nBB,label='Dm1')
call mma_allocate(Dm2,nBB,label='Dm2')
call UnFold(Dma,nBDT,Dm1,nBB,nSym,nBas)
call UnFold(Dmb,nBDT,Dm2,nBB,nSym,nBas)
iOff = 1
do i=1,nSym
  call CD_InCore(Dm1(iOff),nBas(i),Porb(1)%SB(i)%A2,nBas(i),nIorb(i,1),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
    call Abend()
  end if
  call CD_InCore(Dm2(iOff),nBas(i),Porb(2)%SB(i)%A2,nBas(i),nIorb(i,2),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Beta density. Sym= ',i,'   rc= ',irc
    call Abend()
  end if
  iOff = iOff+nBas(i)**2
end do

Call Allocate_DSBA(FLT(1),nBas,nBas,nSym,Case='TRI')
Call Allocate_DSBA(FLT(2),nBas,nBas,nSym,Case='TRI')
FLT(1)%A0(:)=Zero
FLT(2)%A0(:)=Zero

Call Allocate_DSBA(KLT(1),nBas,nBas,nSym,Case='TRI')
Call Allocate_DSBA(KLT(2),nBas,nBas,nSym,Case='TRI')
KLT(1)%A0(:)=Zero
KLT(2)%A0(:)=Zero

call Cho_X_init(irc,ChFracMem)
if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_init.')
  call Abend()
end if


call CHO_LK_SCF(irc,2,FLT,KLT,nForb,nIorb,Porb,PLT,FactXI,nSCReen,dmpk,dFmat)

if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_LK_scf.')
  call Abend()
end if

call Cho_X_Final(irc)
if (irc /= 0) then
  call WarningMessage(2,'Get_CNOs. Non-zero rc in Cho_X_Final.')
  call Abend()
end if

Etwo = Half*(ddot_(nBDT,Dma,1,FLT(1)%A0,1)+ddot_(nBDT,Dmb,1,FLT(2)%A0,1))

call deallocate_DSBA(PLT(1))
call deallocate_DSBA(Porb(2))
call deallocate_DSBA(Porb(1))
call mma_deallocate(Dm1)
call mma_deallocate(Dm2)
call deallocate_DSBA(FLT(2))
call deallocate_DSBA(FLT(1))
call deallocate_DSBA(KLT(2))
call deallocate_DSBA(KLT(1))

return

end subroutine Get_Etwo_act
