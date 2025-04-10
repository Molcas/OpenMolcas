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

subroutine ci_reconstruct(istate,nSDET,vector,indexSD)

use dmrginfo, only: LRRAS2, RGRAS2, nEle_RGLR, nDets_RGLR, MS2_RGLR, nStates_RGLR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: u6

implicit none
character(len=100), allocatable :: checkpoint(:)! for many states
integer :: istate, nsdet
integer i, j, idet, idx_det
integer norb, norbLR
integer neletol
integer nele_mod
integer nele_alpha, nele_beta
integer irrep_pre, iorbLR0, iorbLR
integer irrep_diff(8)
! need to be rewrittren using mma_allocate and mma_deallocate
! soon ...
integer, allocatable :: ele_orb_alpha(:)
integer, allocatable :: ele_orb_beta(:)
integer, allocatable :: pre_ele(:)
integer ndets_mclr
character(len=200) tmp_run
real*8 dtmp
logical IFFILE
integer rc
integer :: indexSD(nsdet) ! index
real*8 :: vector(nsdet)   ! determinants
integer :: nDets_Total, lcheckpoint
integer, external :: IsFreeUnit
type Slater_determinant
  integer :: itype = 1                ! excitation type
  integer :: inum = 1                 ! determinant number
  integer :: isign = 1                ! determinant phase
  integer, allocatable :: electron(:) ! determinant
  integer, allocatable :: ele_conf(:) ! electron configuration
  real*8, allocatable :: dV(:)        ! for many states
end type Slater_determinant
type(Slater_determinant), allocatable :: SD_DMRG(:)

! The total electrons
neletol = 0
neletol = nele_RGLR

! Check if there is single electron
nele_mod = mod(neletol,2)

! electrons in alpha or beta orbitals
if (nele_mod == 0) then
  ! If no single electron
  nele_alpha = neletol/2+ms2_RGLR/2
  nele_beta = neletol/2-ms2_RGLR/2
else
  nele_alpha = neletol/2+ms2_RGLR/2+1
  nele_beta = neletol/2-ms2_RGLR/2
end if

! Read in all the name of checkpoint file
lcheckpoint = 20
lcheckpoint = isFreeUnit(lcheckpoint)
call Molcas_Open(lcheckpoint,'dmrg_for_mclr.parameters')
call mma_allocate(checkpoint,nstates_RGLR,label='checkpoint')
checkpoint = ''
do i=1,6
  read(lcheckpoint,*)
end do
do i=1,nstates_RGLR
  read(lcheckpoint,*) checkpoint(i)
  read(lcheckpoint,*)
end do
close(lcheckpoint)

! reconstructing dets for current state
do i=istate,istate
  write(u6,*) trim(checkpoint(i))
end do

! preparing for point group symmetry
norb = sum(RGras2(1:8))
norbLR = sum(LRras2(1:8))
irrep_diff(:) = RGras2(1:8)-LRras2(1:8)

call mma_allocate(pre_ele,norbLR,label='pre_ele'); pre_ele = 0
iorbLR0 = 1
iorbLR = 0
irrep_pre = 0
do i=1,8
  iorbLR = iorbLR+LRras2(i)
  if (i /= 1) irrep_pre = irrep_pre+irrep_diff(i-1)
  do j=iorbLR0,iorbLR,1
    pre_ele(j) = irrep_pre
    write(u6,*) 'j,pre_ele(j)',j,pre_ele(j)
  end do
  iorbLR0 = iorbLR+1 ! At least work for C1
end do

write(u6,*) 'pre_ele, ndets_RGLR',pre_ele,ndets_RGLR
write(u6,*) 'nalpha,  nbeta     ',nele_alpha,nele_beta

! DETs read from mclr_dets.initial
open(UNIT=117,file='mclr_dets.initial',status='OLD')
allocate(SD_DMRG(ndets_RGLR))
do idet=1,ndets_RGLR
  call mma_allocate(SD_DMRG(idet)%electron,neletol)
  call mma_allocate(SD_DMRG(idet)%ele_conf,norb)
  call mma_allocate(SD_DMRG(idet)%dv,nstates_RGLR)
  SD_DMRG(idet)%electron = 0
  SD_DMRG(idet)%ele_conf = 0
  SD_DMRG(idet)%dv = Zero
  read(117,'(1X,I8,6X)',advance='no') SD_DMRG(idet)%ITYPE
  do i=1,neletol
    read(117,'(1X,I5)',advance='no') SD_DMRG(idet)%electron(i)
  end do
  read(117,'(5X,I3)',advance='no') SD_DMRG(idet)%isign
  read(117,'(11X,I20)') SD_DMRG(idet)%inum
end do
close(117)

write(u6,*) 'before get the executable file'

! get the executable file
call systemf('cp /home/eth/yma/Maquis_MPSLR/build/applications/srcas/srcas $PWD',rc)

write(u6,*) 'before get the executable file'

call mma_allocate(ele_orb_alpha,norb,label='ele_orb_alpha')
call mma_allocate(ele_orb_beta,norb,label='ele_orb_beta')
! All of the DETs into Maquis format
open(unit=118,file='dets.mclr')
do idet=1,ndets_RGLR
  !write(u6,*) 'idet',idet
  ele_orb_alpha = 0
  ele_orb_beta = 0
  do i=1,neletol
    if (SD_DMRG(idet)%electron(i) > 0) then ! With preconditioner
      j = abs(SD_DMRG(idet)%electron(i))
      ele_orb_alpha(j+pre_ele(j)) = 1
    else
      j = abs(SD_DMRG(idet)%electron(i))
      ele_orb_beta(j+pre_ele(j)) = 1
    end if
  end do
  ! The same style also used in Maquis input
  do i=1,norb
    if ((ele_orb_alpha(i) == 1) .and. (ele_orb_beta(i) == 1)) SD_DMRG(idet)%ele_conf(i) = 4
    if ((ele_orb_alpha(i) == 1) .and. (ele_orb_beta(i) == 0)) SD_DMRG(idet)%ele_conf(i) = 3
    if ((ele_orb_alpha(i) == 0) .and. (ele_orb_beta(i) == 1)) SD_DMRG(idet)%ele_conf(i) = 2
    if ((ele_orb_alpha(i) == 0) .and. (ele_orb_beta(i) == 0)) SD_DMRG(idet)%ele_conf(i) = 1
    !write(u6,*) 'SD_DMRG(idet)%ele_conf(i)',SD_DMRG(idet)%ele_conf(i)
  end do
  do i=1,norb
    write(118,'(I1)',advance='no') SD_DMRG(idet)%ele_conf(i)
    !write(u6,*)'SD_DMRG(idet)%ele_conf(',i,')',SD_DMRG(idet)%ele_conf(i)
  end do
  write(118,*)
end do
close(118)
call mma_deallocate(ele_orb_alpha)
call mma_deallocate(ele_orb_beta)

write(u6,*) 'After write dets.mclr file'

! Test for part of DETS
! The way of "open file" need to be written
!                                  Yingjin 2015.8.13
call systemf('wc -l dets.mclr > dets.mclr.info',rc)
open(unit=118,file='dets.mclr.info')
read(118,*) ndets_total
close(118)
! If too many determinants,
! try to use the single, double, triple gradually untill 9999 (as the maximum)
if (ndets_total > 9999) call systemf('head -9999 dets.mclr > ELE_CISR_FOR_MCLR',rc)

! Recover the determinants, off-diagional multiply 2
! ========= should be improved by Hash etc. ==========
do i=istate,istate
  !write(u6,*) 'SD_DMRG(idet)%dv(i)',i
  open(unit=118,file='GET_COEFF_IN_LIST')
  call f_inquire('ELE_CISR_FOR_MCLR',IFFILE)
  if (IFFILE) then
    tmp_run = './srcas '//trim(checkpoint(i))//' dets.mclr 1.0 1.0 0 ELE_CISR_FOR_MCLR > CIRE.scratch'
  else
    tmp_run = './srcas '//trim(checkpoint(i))//' dets.mclr 1.0 1.0 0 dets.mclr > CIRE.scratch'
  end if
  write(118,*) trim(tmp_run)
  close(118)
  call systemf('chmod +x GET_COEFF_IN_LIST',rc)
  call systemf('./GET_COEFF_IN_LIST',rc)
  ! read in the dets-coefficients
  open(unit=118,file='det_coeff.tmp')
  read(118,*) ndets_mclr
  do idet=1,ndets_mclr
    read(118,*) idx_det,SD_DMRG(idx_det)%dv(i)
    !write(u6,*) idx_det,SD_DMRG(idx_det)%dv(i)
  end do
  close(118)
  ! off-diagional multiply 2
  do idet=1,ndets_RGLR
    if (SD_DMRG(idet)%ITYPE == 1) then
      SD_DMRG(idet)%dv(i) = -SD_DMRG(idet)%dv(i)
    else
      dtmp = sqrt((SD_DMRG(idet)%dv(i)**2)*Two)
      SD_DMRG(idet)%dv(i) = -sign(dtmp,SD_DMRG(idet)%dv(i))
    end if
    !write(u6,*) 'i,idet,dv',i,idet,SD_DMRG(idet)%dv(i)
  end do
end do

dtmp = Zero
indexSD = 0
vector = Zero
do i=1,ndets_RGLR
  indexSD(i) = i !SD_DMRG(i)%inum*SD_DMRG(i)%isign
  vector(i) = SD_DMRG(i)%dv(istate)
  dtmp = dtmp+vector(i)**2
end do

write(u6,*) 'nele_alpha,nele_beta',nele_alpha,nele_beta
write(u6,*) 'Total CI weight is ',dtmp
write(u6,*) ' ================================================'
write(u6,*) '  IF for frequency, the CI weight must be'
write(u6,*) '      very close to 1 (i.e. 0.9999)'
write(u6,*) ' ------------------------------------------------'
write(u6,*) '  IF gradients in state-averaged case,'
write(u6,*) '      even very few is stil OK (e.g. 0.1)'
write(u6,*) '      however, better around 0.9'
write(u6,*) ' ================================================'
call xflush(u6)

!stop

call mma_deallocate(checkpoint)
call mma_deallocate(pre_ele)
do idet=1,size(SD_DMRG)
  call mma_deallocate(SD_DMRG(idet)%electron)
  call mma_deallocate(SD_DMRG(idet)%ele_conf)
  call mma_deallocate(SD_DMRG(idet)%dv)
end do
deallocate(SD_DMRG)

end subroutine ci_reconstruct
