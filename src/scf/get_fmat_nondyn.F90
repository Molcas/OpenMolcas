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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************

subroutine Get_Fmat_nondyn(Dma,Dmb,nBDT,DFTX)

use Fock_util_global, only: Deco
use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use SpinAV, only: Do_SpinAV, DSC
use InfSCF, only: dmpk, E_nondyn, Erest_xc, KSDFT, nBas, nBB, nScreen, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBDT
real(kind=wp), intent(inout) :: Dma(nBDT), Dmb(nBDT)
logical(kind=iwp), intent(in) :: DFTX
integer(kind=iwp) :: i, iOff, ipDai, ipDbi, iRC, nDMat, nForb(8,2), nIorb(8,2)
real(kind=wp) :: dFMat, E2act(1), FactXI
type(DSBA_Type) :: FLT(2), KLT(2), PLT(2), POrb(2)
real(kind=wp), allocatable :: Dm(:,:)
real(kind=wp), external :: DDot_, Get_ExFac

nDMat = 2
nForb(1:nSym,:) = 0
if (DFTX) then
  FactXI = Get_ExFac(KSDFT)-One ! note this trick
else
  FactXI = One
end if

call Allocate_DT(PLT(1),nBas,nBas,nSym,aCase='TRI')
call Allocate_DT(PLT(2),nBas,nBas,nSym,aCase='TRI',Ref=PLT(1)%A0)
if (DFTX) then
  PLT(1)%A0(:) = Zero
else
  PLT(1)%A0(:) = Dma(:)+Dmb(:)
end if

call Allocate_DT(POrb(1),nBas,nBas,nSym)
call Allocate_DT(POrb(2),nBas,nBas,nSym)

call mma_allocate(Dm,nBB,2,Label='Dm')
call UnFold(Dma,nBDT,Dm(:,1),nBB,nSym,nBas)
call UnFold(Dmb,nBDT,Dm(:,2),nBB,nSym,nBas)

if (Do_SpinAV) then
  if (.not. DECO) then
    write(u6,*) ' Keywords NODE and SAVE are incompatible. '
    write(u6,*) ' NODE will be reset to default. '
    DECO = .true.
  end if
  Dm(:,1) = Dm(:,1)-DSc(:)
  Dm(:,2) = Dm(:,1)+DSc(:)
end if

iOff = 0
do i=1,nSym
  ipDai = 1+iOff
  call CD_InCore(Dm(ipDai,1),nBas(i),Porb(1)%SB(i)%A2,nBas(i),nIorb(i,1),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
    call RecPrt('Dm',' ',Dm(ipDai,1),nBas(i),nBas(i))
    call Abend()
  end if
  ipDbi = 1+iOff
  call CD_InCore(Dm(ipDbi,2),nBas(i),Porb(2)%SB(i)%A2,nBas(i),nIorb(i,2),1.0e-12_wp,irc)
  if (irc /= 0) then
    write(u6,*) ' Beta density. Sym= ',i,'   rc= ',irc
    call RecPrt('Dm',' ',Dm(ipDbi,1),nBas(i),nBas(i))
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

dFmat = Zero
call CHO_LK_SCF(irc,nDMat,FLT,KLT,nForb,nIorb,Porb,PLT,FactXI,nSCReen,dmpk,dFmat)
if (irc /= 0) then
  call WarningMessage(2,'Start6. Non-zero rc in Cho_LK_scf.')
  call Abend()
end if

if (Do_SpinAV) then
  call UnFold(Dma,nBDT,Dm(:,1),nBB,nSym,nBas)
  call UnFold(Dmb,nBDT,Dm(:,2),nBB,nSym,nBas)
  Dm(:,1) = Dm(:,1)-DSc(:)
  Dm(:,2) = Dm(:,2)+DSc(:)
  call Fold(nSym,nBas,Dm(:,1),Dma)
  call Fold(nSym,nBas,Dm(:,2),Dmb)
end if

E2act(1) = Half*(ddot_(nBDT,Dma,1,FLT(1)%A0,1)+ddot_(nBDT,Dmb,1,FLT(2)%A0,1))
call GADSum(E2act(1),1)

if (DFTX) then
  Erest_xc = Erest_xc-E2act(1)
else
  E_nondyn = E_nondyn-E2act(1)
end if

call Deallocate_DT(KLT(2))
call Deallocate_DT(KLT(1))
call Deallocate_DT(FLT(2))
call Deallocate_DT(FLT(1))
call mma_deallocate(Dm)
call Deallocate_DT(POrb(2))
call Deallocate_DT(POrb(1))
call Deallocate_DT(PLT(2))
call Deallocate_DT(PLT(1))

return

end subroutine Get_Fmat_nondyn
