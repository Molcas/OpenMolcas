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

use Fock_util_global, only: Estimate, Update
use Cholesky, only: timings
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBDT, nBas(8), nSym
real(kind=wp), intent(in) :: Dma(nBDT), Dmb(nBDT)
real(kind=wp), intent(out) :: Etwo
integer(kind=iwp) :: i, iOff, irc, nBB, nForb(8,2), nIorb(8,2), NSCREEN
real(kind=wp) :: ChFracMem, dFmat, dmpk, FactXI
!character(len=80) :: KSDFT
real(kind=wp), allocatable :: Dm1(:), Dm2(:)
real(kind=wp), external :: ddot_ !, Get_ExFac
type (DSBA_Type) :: FLT(2), KLT(2), POrb(2), PLT(2)

timings = .false.
Estimate = .false.

Update = .true.
dmpk = One
NSCREEN = 10

nForb(:,:) = 0
nBB = 0
do i=1,nSym
  nBB = nBB+nBas(i)**2
end do
!call Get_cArray('DFT functional',KSDFT,80)
!ExFac = Get_ExFac(KSDFT)
!FactXI = ExFac
FactXI = One  ! always HF energy
call Allocate_DT(PLT(1),nBas,nBas,nSym,aCase='TRI')
PLT(1)%A0(:) = Dma(:)+Dmb(:)

call Allocate_DT(POrb(1),nBas,nBas,nSym)
call Allocate_DT(POrb(2),nBas,nBas,nSym)
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

call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI')
call Allocate_DT(FLT(2),nBas,nBas,nSym,aCase='TRI')
FLT(1)%A0(:) = Zero
FLT(2)%A0(:) = Zero

call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI')
call Allocate_DT(KLT(2),nBas,nBas,nSym,aCase='TRI')
KLT(1)%A0(:) = Zero
KLT(2)%A0(:) = Zero

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

call Deallocate_DT(PLT(1))
call Deallocate_DT(Porb(2))
call Deallocate_DT(Porb(1))
call mma_deallocate(Dm1)
call mma_deallocate(Dm2)
call Deallocate_DT(FLT(2))
call Deallocate_DT(FLT(1))
call Deallocate_DT(KLT(2))
call Deallocate_DT(KLT(1))

return

end subroutine Get_Etwo_act
