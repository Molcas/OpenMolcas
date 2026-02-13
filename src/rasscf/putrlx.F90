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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine PutRlx(D,DS,P,DAO,C)

use spin_correlation, only: tRootGrad
use rasscf_global, only: CBLBM, ENER, ExFac, iBLBM, IPCMROOT, iPr, iRLXRoot, iSymBB, ITER, lRoots, lSquare, NACPAR, NACPR2, &
                         NewFock, nFint, nRoots, NSXS, NTOT4, RlxGrd, iAdr15, ISTORP, JBLBM
use PrintLevel, only: DEBUG, USUAL
use output_ras, only: IPRLOC
use general_data, only: NTOT1, NTOT2, NSYM, JOBIPH, NBAS
use DWSol, only: DWSolv, DWSol_wgt, W_SOLV
use rctfld_module, only: lRF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit none
real*8 D(*), DS(*), P(*), DAO(*), C(*)
character(len=8) Fmt2
character(len=16), parameter :: ROUTINE = 'PUTRLX  '
real*8 rdum(1)
integer i, iFinal, iPrLev, istmp, itmp, jDisk, jtmp, kDisk, left, NFSize, NZ
real*8 rTmp, wgt
real*8, external :: DNRM2_
real*8, allocatable :: DA(:), DA_ave(:), DI(:), DSX(:), DS_ave(:), DX(:), F(:), B(:), Q(:), FA(:), FI(:), PUVX(:), TUVX(:), PX(:), &
                       WRK1(:), WRK2(:)
logical, external :: PCM_On

IPRLEV = IPRLOC(3)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering ',ROUTINE

! Read in corresponding density matrixes

! empty line before root resolved orbital gradient
NZ = NTOT2
kDisk = IADR15(3)  ! we need to save this address to reset later
jDisk = IADR15(3)
call mma_allocate(DA,NZ,Label='DA')
call mma_allocate(DI,NZ,Label='DI')
call mma_allocate(DSX,NZ,Label='DSX')
if (tRootGrad) then
  write(u6,*)
  do i=1,LRoots
    call ddafile(JOBIPH,2,D,NACPAR,jDisk)
    call ddafile(JOBIPH,2,DS,NACPAR,jDisk)
    call ddafile(JOBIPH,2,P,NACPR2,jDisk)
    call ddafile(JOBIPH,0,rdum,NACPR2,jDisk)

    ! Construc D-ACTIVE AND D-INACTIVE IN AO BASIS

    call Get_D1I_RASSCF(C,DI)

    call mma_allocate(DX,NACPAR,Label='DX')

    call DCOPY_(NACPAR,DS,1,DX,1)
    call DBLOCK(DX)
    call Get_D1A_RASSCF(C,DX,DSX)

    call DCOPY_(NACPAR,D,1,DX,1)
    call DBLOCK(DX)
    call Get_D1A_RASSCF(C,DX,DA)

    ! Construct the Fock matrix used for the connection term.

    NFSIZE = max(NACPAR,NTOT4)
    call mma_allocate(F,NFSIZE,Label='F')
    call mma_allocate(B,NZ,Label='B')
    call mma_allocate(Q,NZ,Label='Q')
    call mma_allocate(FA,NZ,Label='FA')
    call mma_allocate(FI,NZ,Label='FI')
    call mma_allocate(PUVX,NFINT,Label='PUVX')
    call mma_allocate(tuvx,nacpr2,Label='tuvx')
    IPR = 0
    if (IPRLEV == 3) IPR = 1
    if (IPRLEV == 4) IPR = 5
    if (IPRLEV == 5) IPR = 10
    PUVX(:) = Zero
    call TraCtl2(C,PUVX,TUVX,DI,FI,DA,FA,ipr,lsquare,ExFac)
    call SGFCIN(C,F,FI,DI,DA,DSX)
    call dcopy_(ntot4,[Zero],0,F,1)
    call dcopy_(ntot4,[Zero],0,B,1)

    ! Prevent FMAT from changing Active fock matrix

    iTmp = newfock
    newFock = -99999

    call Fmat(C,PUVX,DX,DA,FI,FA)
    newFock = itmp
    IFINAL = 1
    rTmp = CBLBM
    itmp = IBLBM
    jtmp = JBLBM
    istmp = ISYMBB

    if (ISTORP(NSYM+1) > 0) then
      call mma_allocate(PX,ISTORP(NSYM+1),Label='PX')
      call PMAT_RASSCF(P,PX)
    else
      call mma_allocate(PX,1,Label='PX')
    end if

    call FOCK(F,B,FI,FA,DX,PX,Q,PUVX,IFINAL,C)
    CBLBM = rtmp
    IBLBM = itmp
    JBLBM = jtmp
    ISYMBB = istmp
    call mma_deallocate(PX)

    write(u6,'(6x,a,i3,5x,f12.10)') 'Norm of electronic gradient for root ',i,DNRM2_(NSXS,B,1)

    call mma_deallocate(DX)
    call mma_deallocate(tuvx)
    call mma_deallocate(puvx)
    call mma_deallocate(FI)
    call mma_deallocate(FA)
    call mma_deallocate(Q)
    call mma_deallocate(B)
    call mma_deallocate(F)
  end do
end if

do i=1,iRlxRoot-1
  call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
  call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
  call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
  call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
end do
call DDaFile(JOBIPH,2,D,NACPAR,kDisk)
call DDaFile(JOBIPH,2,DS,NACPAR,kDisk)
call DDaFile(JOBIPH,2,P,NACPR2,kDisk)
call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)

! Construct D-ACTIVE AND D-INACTIVE IN AO BASIS

call Get_D1I_RASSCF(C,DI)

call mma_allocate(DX,NACPAR,Label='DX')

call DCOPY_(NACPAR,DS,1,DX,1)
call DBLOCK(DX)
call Get_D1A_RASSCF(C,DX,DSX)

call DCOPY_(NACPAR,D,1,DX,1)
call DBLOCK(DX)
call Get_D1A_RASSCF(C,DX,DA)

! Construct the Fock matrix used for the connection term.

NFSIZE = max(NACPAR,NTOT4)
call mma_allocate(F,NFSIZE,Label='F')
call mma_allocate(B,NZ,Label='B')
call mma_allocate(Q,NZ,Label='Q')
call mma_allocate(FA,NZ,Label='FA')
call mma_allocate(FI,NZ,Label='FI')
call mma_allocate(PUVX,NFINT,Label='PUVX')
call mma_allocate(tuvx,nacpr2,Label='tuvx')
IPR = 0
if (IPRLEV == 3) IPR = 1
if (IPRLEV == 4) IPR = 5
if (IPRLEV == 5) IPR = 10
PUVX(:) = Zero
call TraCtl2(C,PUVX,TUVX,DI,FI,DA,FA,ipr,lsquare,ExFac)

! DA constructed above is the density for the RLXROOT state.
! We should reconstruct the AO density matrix to prevent the
! mismatch of states for reaction field and geometry optimization
! It seems that IPCMROOT may not be defined sometimes
if (lRF .and. PCM_On() .and. ((IPCMROOT /= iRLXROOT) .or. (IPCMROOT <= 0) .or. (DWSolv%DWZeta /= Zero))) then
  !! Polarize PCM etc with weighted density
  call mma_allocate(DA_ave,max(NACPAR,NZ),Label='DA_ave')
  call mma_allocate(DS_ave,max(NACPAR,NZ),Label='DS_ave')
  call mma_allocate(WRK1,max(NACPAR,NZ),Label='WRK1')
  call mma_allocate(WRK2,max(NACPAR,NZ),Label='WRK2')

  DA_ave(:) = Zero
  DS_ave(:) = Zero

  call DWSol_wgt(2,ENER(:,ITER))
  kDisk = IADR15(3)
  do i=1,nRoots
    wgt = W_SOLV(i)
    call DDaFile(JOBIPH,2,WRK1,NACPAR,kDisk)
    call DDaFile(JOBIPH,2,WRK2,NACPAR,kDisk)
    call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
    call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
    if (wgt < 1.0e-09_wp) cycle
    call daxpy_(NACPAR,wgt,WRK1,1,DA_ave,1)
    call daxpy_(NACPAR,wgt,WRK2,1,DS_ave,1)
  end do

  ! Construc D-ACTIVE AND D-INACTIVE IN AO BASIS

  call DCOPY_(NACPAR,DS_ave,1,DX,1)
  call dcopy_(NZ,DSX,1,DS_ave,1)
  call DBLOCK(DX)
  call Get_D1A_RASSCF(C,DX,DSX)

  call DCOPY_(NACPAR,DA_ave,1,DX,1)
  call dcopy_(NZ,DA,1,DA_ave,1)
  call DBLOCK(DX)
  call Get_D1A_RASSCF(C,DX,DA)

  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)
  if ((IPRLEV >= USUAL) .and. (DWSolv%DWZeta > Zero)) then
    left = 6
    write(u6,*)
    write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
    write(u6,Fmt2//'A)') 'Dynamically weighted solvation has been employed'
    write(u6,Fmt2//'A,(T45,10F6.3))') 'Final weights for the reaction field:',(W_SOLV(i),i=1,nRoots)
  end if
end if

call SGFCIN(C,F,FI,DI,DA,DSX)
call dcopy_(ntot4,[Zero],0,F,1)
call dcopy_(ntot4,[Zero],0,B,1)

! Prevent FMAT from changing Active fock matrix

iTmp = newfock
newFock = -99999

if (lRF .and. PCM_On() .and. ((IPCMROOT /= iRLXROOT) .or. (IPCMROOT <= 0) .or. (DWSolv%DWZeta /= Zero))) then
  ! The rest of the RASSCF program uses state-specific density
  ! (note that, it is iRlxRoot!), so restore the one constructed
  ! above, before TraCtl2
  call dcopy_(NZ,DA_ave,1,DA,1)
  call dcopy_(NZ,DS_ave,1,DSX,1)
  call mma_deallocate(DA_ave)
  call mma_deallocate(DS_ave)

  call mma_allocate(WRK1,nTot1,Label='WRK1')
  call mma_allocate(WRK2,nTot1,Label='WRK2')
  call Fold(nSym,nBas,DI,WRK1)
  call Fold(nSym,nBas,DA,WRK2)
  call Daxpy_(nTot1,One,WRK1,1,WRK2,1)
  call Put_dArray('D1ao',WRK2,nTot1)
  call Fold(nSym,nBas,DSX,WRK1)
  call Put_dArray('D1sao',WRK1,nTot1)
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)
end if

call Fmat(C,PUVX,DX,DA,FI,FA)
newFock = itmp
IFINAL = 1
rTmp = CBLBM
itmp = IBLBM
jtmp = JBLBM
istmp = ISYMBB

if (ISTORP(NSYM+1) > 0) then
  call mma_allocate(PX,ISTORP(NSYM+1),Label='PX')
  call PMAT_RASSCF(P,PX)
else
  call mma_allocate(PX,1,Label='PX')
end if

call FOCK(F,B,FI,FA,DX,PX,Q,PUVX,IFINAL,C)
CBLBM = rtmp
IBLBM = itmp
JBLBM = jtmp
ISYMBB = istmp
call mma_deallocate(PX)

RlxGrd = DNRM2_(NSXS,B,1)

call mma_deallocate(DX)
call mma_deallocate(tuvx)
call mma_deallocate(puvx)
call mma_deallocate(FI)
call mma_deallocate(FA)
call mma_deallocate(Q)
call mma_deallocate(B)
call mma_deallocate(F)

! Add up one electron densities

call daxpy_(nZ,One,DA,1,DI,1)
call Fold(nSym,nBas,DI,DAO)

call mma_deallocate(DSX)
call mma_deallocate(DA)
call mma_deallocate(DI)

end subroutine PutRlx
