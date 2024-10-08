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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2016,2017, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine OptClc(CInter,nCI,nD,Ind,nInd)
!***********************************************************************
!                                                                      *
! purpose: calculate optimal density matrix and two-electron hamil-    *
!          tonial as well as optimal difference from interpolation     *
!          or extrapolation coefficients.                              *
!                                                                      *
! input:                                                               *
!   Dens    : density matrices (nDT,NumDT)                             *
!   TwoHam  : two-electron hamiltonian matrices (nDT,NumDT)            *
!   Vxc     : external potential       matrices (nDT,NumDT)            *
!   CInter  : interpolation coefficients (nCI)                         *
!                                                                      *
! output:                                                              *
!   Dens and TwoHam                                                    *
!                                                                      *
!***********************************************************************

use InfSCF, only: Dens, iDisk, kOptim, MapDns, nBT, nDens, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCI, nD, nInd, Ind(nInd)
real(kind=wp), intent(in) :: CInter(nCI,nD)
integer(kind=iwp) :: i, iD, iMap, Iter_D, MatNO
real(kind=wp) :: C
real(kind=wp), allocatable :: DnsTmp(:,:), TwoTmp(:,:), VxcTmp(:,:)

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!                                                                      *
! Allocate memory for matrices that contribute to the optimal one
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
if (kOptim == 1) return
call mma_allocate(DnsTmp,nBT,nD,Label='DnsTmp')
call mma_allocate(TwoTmp,nBT,nD,Label='TwoTmp')
call mma_allocate(VxcTmp,nBT,nD,Label='VxcTmp')
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Accumulate linear combinations in position iPsLst
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Start with the last iteration
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
iter_d = Ind(kOptim)

iMap = MapDns(iter_d)
if (iMap < 0) then
  call RWDTG(-iMap,DnsTmp,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
  call RWDTG(-iMap,TwoTmp,nBT*nD,'R','TWOHAM',iDisk,size(iDisk,1))
  call RWDTG(-iMap,VxcTmp,nBT*nD,'R','dVxcdR',iDisk,size(iDisk,1))
else
  DnsTmp(:,:) = Dens(:,:,iMap)
  TwoTmp(:,:) = TwoHam(:,:,iMap)
  VxcTmp(:,:) = Vxc(:,:,iMap)
end if

do iD=1,nD

  C = CInter(kOptim,iD)
  DnsTmp(:,iD) = C*DnsTmp(:,iD)
  TwoTmp(:,iD) = C*TwoTmp(:,iD)
  VxcTmp(:,iD) = C*VxcTmp(:,iD)

end do

Dens(:,:,nDens) = DnsTmp(:,:)
TwoHam(:,:,nDens) = TwoTmp(:,:)
Vxc(:,:,nDens) = VxcTmp(:,:)

do i=1,kOptim-1
  C = CInter(i,1)
  MatNo = Ind(i)

  iMap = MapDns(MatNo)
  if (iMap < 0) then
    call RWDTG(-iMap,DnsTmp,nBT*nD,'R','DENS  ',iDisk,size(iDisk,1))
    call RWDTG(-iMap,TwoTmp,nBT*nD,'R','TWOHAM',iDisk,size(iDisk,1))
    call RWDTG(-iMap,VxcTmp,nBT*nD,'R','dVxcdR',iDisk,size(iDisk,1))
  else
    DnsTmp(:,:) = Dens(:,:,iMap)
    TwoTmp(:,:) = TwoHam(:,:,iMap)
    VxcTmp(:,:) = Vxc(:,:,iMap)
  end if

  do iD=1,nD
    C = CInter(i,iD)
    Dens(:,iD,nDens) = Dens(:,iD,nDens)+C*DnsTmp(:,iD)
    TwoHam(:,iD,nDens) = TwoHam(:,iD,nDens)+C*TwoTmp(:,iD)
    Vxc(:,iD,nDens) = Vxc(:,iD,nDens)+C*VxcTmp(:,iD)
  end do

end do ! i

! Deallocate memory

call mma_deallocate(DnsTmp)
call mma_deallocate(TwoTmp)
call mma_deallocate(VxcTmp)

#ifdef _DEBUGPRINT_
call NrmClc(Dens(1,1,nDens),nBT*nD,'OptClc','D in iPsLst. ')
call NrmClc(TwoHam(1,1,nDens),nBT*nD,'OptClc','T in iPsLst. ')
call NrmClc(Vxc(1,1,nDens),nBT*nD,'OptClc','V in iPsLst. ')
#endif
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
return

end subroutine OptClc
