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

subroutine FockTwo_Drv_scf(nSym,nBas,nAux,Keep,DLT,DSQ,FLT,nFLT,ExFac,nBSQT,nBMX,nD,nOcc,lOcc,iDummy_run)

use OFembed, only: Do_OFemb, FMaux, OFE_first, Rep_EN
use InfSCF, only: ALGO
use Cholesky, only: timings
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nAux(nSym), Keep(nSym), nFLT, nBSQT, nBMX, nD, lOcc, nOcc(lOcc,nD), iDummy_run
real(kind=wp), intent(inout) :: DLT(nFLT,nD), FLT(nFLT,nD)
real(kind=wp), intent(in) :: DSQ(nBSQT,nD), ExFac
integer(kind=iwp) :: iRC, lBuf
real(kind=wp) :: TotCPU, TotCPU1, TotCPU2, TotWall, TotWall1, TotWall2
logical(kind=iwp) :: DoCholesky, GenInt
character(len=50) :: CFmt
real(kind=wp), allocatable :: FSQ(:,:), tFLT(:,:), W1(:), W2(:)

GenInt = .false.
DoCholesky = .false.
if (ALGO == 0) GenInt = .true. !use GenInt to regenerate integrals
call DecideOnCholesky(DoCholesky)

!write(u6,*) '*************************'
!write(u6,*) 'ONLY COULOMB CONTRIBUTION'
!write(u6,*) '*************************'
!ExFac = Zero
!write(u6,*) 'ExFac= ',ExFac

if (Do_OFemb) then ! Coul. potential from subsys B
  if (OFE_first) call mma_allocate(FMaux,nFlt,Label='FMaux')

  call Coul_DMB(OFE_first,nD,Rep_EN,FMaux,DLT(:,1),DLT(:,nD),nFlt)
  OFE_first = .false.
end if

call mma_allocate(FSQ,nBSQT,nD,Label='FSQ')
FSQ(:,:) = Zero

if ((.not. DoCholesky) .or. GenInt) call mma_Allocate(W2,NBMX*NBMX,Label='W2')

! nFlt is the total dimension of the LT fock matrix
call mma_allocate(tFLT,nFLT,nD,Label='tFLT')
tFLT(:,:) = Zero

call mma_maxDBLE(LBUF)

! Standard building of the Fock matrix from Two-el integrals
! or
! Building of the Fock matrix regenerating the integrals on the fly

call CWTIME(TotCPU1,TotWALL1)

if ((.not. DoCholesky) .or. GenInt) then

  if (DoCholesky) then
    ! save some space for GenInt
    LBUF = max(LBUF-LBUF/10,0)
    ! Make sure that the ri/ch vectors are in reordered mode
    call Cho_X_ReOVec(irc)
  end if
  call mma_allocate(W1,LBUF,Label='W1')

  if (LBUF < NBMX**2) then
    write(u6,*) 'FockTwo_Drv_SCF Error: Too little memory remains for the call to FOCKTWO_SCF.'
    write(u6,*) ' Largest allocatable array size LBUF=',LBUF
    write(u6,*) ' Max nr of bf in any symmetry,  NBMX=',NBMX
    write(u6,*) ' Required minimum size       NBMX**2=',NBMX**2
    write(u6,*) '    (All in Real words)'
    call ABEND()
  end if

  call FOCKTWO_scf(nSym,nBas,nAux,Keep,DLT,DSQ,tFLT,nFlt,FSQ,W1,size(W1),W2,size(W2),ExFac,nD,nBSQT)

end if

call CWTIME(TotCPU2,TotWALL2)
TOTCPU = TotCPU2-TotCPU1
TOTWALL = TotWALL2-TotWALL1

! Timings information for conventional or Cholesky with ALGO=0
if ((.not. DoCholesky) .or. GenInt) then
  if (timings) then
    CFmt = '(2x,A)'
    write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
    if (DoCholesky) then
      write(u6,CFmt) '---    Cholesky SCF - Integral regeneration   ---'
    else
      write(u6,CFmt) '-----------     Conventional SCF     ------------'
    end if
    write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
    write(u6,CFmt) 'Fock matrix construction        CPU       WALL   '
    write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
    write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
    write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
    write(u6,*)
  end if
end if

! Building of the Fock matrix directly from Cholesky vectors

if (DoCholesky .and. (.not. GenInt) .and. (iDummy_run == 1)) then
  write(u6,*) '*** Warning: missing feature in Cholesky code'
  write(u6,*) 'Use the results with extra care!'
end if

if (DoCholesky .and. (.not. GenInt) .and. (iDummy_run == 0)) &
  call CHOscf_drv(nBSQT,nD,nSym,nBas,DSQ(:,1),DLT(:,1),DSQ(:,nD),DLT(:,nD),tFLT(:,1),tFLT(:,nD),nFLT,ExFac,FSQ,nOcc(:,1),nOcc(:,nD))

FLT(:,:) = FLT(:,:)+tFLT(:,:)

call mma_deallocate(tFLT)

if (Do_OFemb) then ! add FM from subsystem B
  FLT(:,1) = FLT(:,1)+FMaux(:)
  if (nD == 2) FLT(:,2) = FLT(:,2)+FMaux(:)
end if

if ((.not. DoCholesky) .or. GenInt) then
  call mma_deallocate(W1)
  call mma_deallocate(W2)
end if

call mma_deallocate(FSQ)

return

end subroutine FockTwo_Drv_scf
