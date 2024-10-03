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
!               2003, Valera Veryazov                                  *
!               2022, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DMat(XCf,nXCf,nD)
!***********************************************************************
!                                                                      *
! Purpose: Compute aufbau density matrix                               *
!                                                                      *
!***********************************************************************

use Interfaces_SCF, only: dOne_SCF, MinDns
use InfSCF, only: CMO, DDnOFF, Dens, DNorm, iDisk, InVec, iPsLst, Iter, MapDns, MiniDN, nBas, nBT, nDens, nFrz, nIter, nIterP, &
                  nMem, nOrb, nSym, OccNo, TwoHam, Vxc
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nXCF, nD
real(kind=wp), intent(out) :: XCf(nXCf,nD)
integer(kind=iwp) :: iD, iFrom, iOnDsk, nCMO
logical(kind=iwp) :: alpha_density
real(kind=wp), allocatable :: Aux(:,:)
real(kind=wp), external :: DDot_

!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
nCMO = size(CMO,1)
!call Timing(Cpu1,Tim1,Tim2,Tim3)

! Form proper MapDns vector

if (MapDns(iter) == 0) then   ! Position not defined

  ! Update MapDns and eventually write earlier Dens, TwoHam,
  ! and Vxc matrices to disk.

  if (iter <= nMem) then    ! keep the array in memory

    MapDns(iter) = iter

  else

    iFrom = MapDns(iter-nMem) ! get index of array to dump

    MapDns(iter) = iFrom      ! assign to new array

    if (iter-nMem == 1) then
      MapDns(iter-nMem) = -1 !  Initiate
    else
      MapDns(iter-nMem) = MapDns(iter-nMem-1)-1
    end if

    ! Dump the vectors

    iOnDsk = -MapDns(iter-nMem)
    call RWDTG(iOnDsk,Dens(1,1,iFrom),nBT*nD,'W','DENS  ',iDisk,size(iDisk,1))
    call RWDTG(iOnDsk,TwoHam(1,1,iFrom),nBT*nD,'W','TWOHAM',iDisk,size(iDisk,1))
    call RWDTG(iOnDsk,Vxc(1,1,iFrom),nBT*nD,'W','dVxcdR',iDisk,size(iDisk,1))
  end if

end if

! Check if MapDns is correct

iPsLst = MapDns(iter)
if (iPsLst <= 0) then
  write(u6,*) 'DMat: iPsLst <= 0'
  write(u6,*) 'iPsLst=',iPsLst
  call Abend()
end if

! Form i-th density matrix in the position iPsLst

if ((InVec == 3) .and. (iter == 1)) then

  ! First density matrix is actually in the first position
  ! (read from RUNFILE) on the first iteration.

else

  ! Using the CMOs generate the new density in position iPsLst

  alpha_density = .true.
  do iD=1,nD
    call DOne_SCF(nSym,nBas,nOrb,nFrz,CMO(1,iD),nCMO,OccNo(1,iD),Dens(1,iD,iPsLst),alpha_density)
    alpha_density = .false.
  end do ! iD
end if

do iD=1,nD
  call ChkTrD(nSym,nBas,nOrb,OccNo(:,iD),size(OccNo,1),Dens(:,iD,iPsLst),size(Dens,1))
end do ! iD

#ifdef _DEBUGPRINT_
call NrmClc(Dens(1,1,iPsLst),nBT*nD,'DMat  ','D Iter    ')
#endif

! Put the actual densities on the run file.

! Note, this is from where the DFT part of the DrvXV code gets the
! total density when it computes the DFT contributions to the
! total energy and the Fock matrix.

call DensAB(nBT,iPsLst,nD,Dens)

! Form density difference (normal or minimized). Notice, that
! for minimized differences i-th density is copied to position
! nDens and kept there till the end of iteration.
! For normal differences position nDens is occupied by the density
! of the previous iteration, until
! (i)  OptClc is called (DIIS-only case), or
! (ii) after Diis (QNR case)
! after that, the actual density is stored at that memory location

if (MiniDn .and. (max(0,nIter(nIterP)-1) > 0)) then

  ! Minimized density option

  Dens(:,:,nDens) = Dens(:,:,iPsLst)
  if (iter > 1) call MinDns(Dens,nBT,nDens,XCf,nXCf,nD)

else if (.not. DDnOFF) then

  ! Do the density difference, D(iPsLst)-D(nDens)=D(k+1)-D(k)

  call mma_allocate(Aux,nBT,nD,Label='Aux')
  Aux(:,:) = Dens(:,:,iPsLst)
  Dens(:,:,iPsLst) = Dens(:,:,iPsLst)-Dens(:,:,nDens)
  Dens(:,:,nDens) = Aux(:,:)
  call mma_deallocate(Aux)

else

  Dens(:,:,nDens) = Dens(:,:,iPsLst)

end if

DNorm = real(nD,kind=wp)*DDot_(nBT*nD,Dens(1,1,iPsLst),1,Dens(1,1,iPsLst),1)

#ifdef _DEBUGPRINT_
write(u6,*) 'DNorm=',DNorm
call NrmClc(Dens(1,1,iPsLst),nBT*nD,'DMat  ','D iPsLst  ')
call NrmClc(Dens(1,1,nDens),nBT*nD,'DMat  ','D nDens   ')
#endif

!call Timing(Cpu2,Tim1,Tim2,Tim3)
!TimFld(4) = TimFld(4)+(Cpu2-Cpu1)
!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine DMat
