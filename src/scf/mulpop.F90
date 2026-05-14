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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************

subroutine MulPop(CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
!***********************************************************************
!                                                                      *
! This routine is a wrapper for Charge (which prints Mulliken popu-    *
! lation analyzes) since charge needs square CMO matrices.             *
!                                                                      *
!***********************************************************************

use InfSCF, only: BName, nBas, nBB, nnB, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mBB, nD, mBT, mmB
real(kind=wp), intent(in) :: CMO(mBB,nD), Ovrlp(mBT), OccNo(mmB,nD)
real(kind=wp), allocatable :: Aux1(:), Aux2(:)

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (all(nBas(1:nSym) == nOrb(1:nSym))) then
  if (nD == 1) then
    call Charge(nSym,nBas,BName,CMO(:,1),OccNo(:,1),Ovrlp,2,.false.,.false.)
  else
    call Charge(nSym,nBas,BName,CMO(:,1),OccNo(:,1),Ovrlp,0,.false.,.false.)
    call Charge(nSym,nBas,BName,CMO(:,2),OccNo(:,2),Ovrlp,1,.false.,.false.)
  end if
else
  call mma_allocate(Aux1,nBB,Label='Aux1')
  call mma_allocate(Aux2,nnB,Label='Aux2')
  if (nD == 1) then
    Aux1(1:mBB) = CMO(:,1)
    Aux2(1:mmB) = OccNo(:,1)
    call PadCMO(Aux1,nSym,nBas,nOrb)
    call PadEor(Aux2,nSym,nBas,nOrb)
    call Charge(nSym,nBas,BName,Aux1,Aux2,Ovrlp,2,.false.,.false.)
  else
    Aux1(1:mBB) = CMO(:,2)
    Aux2(1:mmB) = OccNo(:,1)
    call PadCMO(Aux1,nSym,nBas,nOrb)
    call PadEor(Aux2,nSym,nBas,nOrb)
    call Charge(nSym,nBas,BName,Aux1,Aux2,Ovrlp,0,.false.,.false.)
    Aux1(1:mBB) = CMO(:,2)
    Aux2(1:mmB) = OccNo(:,2)
    call PadCMO(Aux1,nSym,nBas,nOrb)
    call PadEor(Aux2,nSym,nBas,nOrb)
    call Charge(nSym,nBas,BName,Aux1,Aux2,Ovrlp,1,.false.,.false.)
  end if
  call mma_deallocate(Aux1)
  call mma_deallocate(Aux2)
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine MulPop
