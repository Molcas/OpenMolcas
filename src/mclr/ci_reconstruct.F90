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

use dmrginfo, only: LRRAS2, MS2_RGLR, nDets_RGLR, nEle_RGLR, nStates_RGLR, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: istate, nsdet
real(kind=wp), intent(out) :: vector(nsdet) ! determinants
integer(kind=iwp), intent(out) :: indexSD(nsdet) ! index
integer(kind=iwp) :: i, idet, idx_det, iorbLR, iorbLR0, irrep_diff(8), irrep_pre, j, lcheckpoint, ndets_mclr, nDets_Total, &
                     nele_alpha, nele_beta, nele_mod, neletol, norb, norbLR, rc
real(kind=wp) :: dtmp
logical(kind=iwp) :: IFFILE
character(len=200) :: tmp_run
integer(kind=iwp), allocatable :: ele_conf(:,:), ele_orb_alpha(:), ele_orb_beta(:), electron(:,:), inum(:), isgn(:), itype(:), &
                                  pre_ele(:)
real(kind=wp), allocatable :: dV(:,:)
character(len=100), allocatable :: checkpoint(:) ! for many states
integer(kind=iwp), external :: IsFreeUnit

! Slater determinant data
!  itype    : excitation type
!  inum     : determinant number
!  isgn     : determinant phase
!  electron : determinant
!  ele_conf : electron configuration
!  dV       : for many states

! The total electrons
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
norb = sum(RGras2(:))
norbLR = sum(LRras2(:))
irrep_diff(:) = RGras2(:)-LRras2(:)

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
call mma_allocate(itype,ndets_RGLR,Label='itype')
call mma_allocate(inum,ndets_RGLR,Label='inum')
call mma_allocate(isgn,ndets_RGLR,Label='isgn')
call mma_allocate(electron,neletol,ndets_RGLR,Label='electron')
call mma_allocate(ele_conf,norb,ndets_RGLR,Label='ele_conf')
call mma_allocate(dV,nstates_RGLR,ndets_RGLR,Label='dV')
itype(:) = 1
inum(:) = 1
isgn(:) = 1
electron(:,:) = 0
ele_conf(:,:) = 0
dV(:,:) = Zero
do idet=1,ndets_RGLR
  read(117,'(1X,I8,6X)',advance='no') itype(idet)
  do i=1,neletol
    read(117,'(1X,I5)',advance='no') electron(i,idet)
  end do
  read(117,'(5X,I3)',advance='no') isgn(idet)
  read(117,'(11X,I20)') inum(idet)
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
    if (electron(i,idet) > 0) then ! With preconditioner
      j = abs(electron(i,idet))
      ele_orb_alpha(j+pre_ele(j)) = 1
    else
      j = abs(electron(i,idet))
      ele_orb_beta(j+pre_ele(j)) = 1
    end if
  end do
  ! The same style also used in Maquis input
  do i=1,norb
    if ((ele_orb_alpha(i) == 1) .and. (ele_orb_beta(i) == 1)) ele_conf(i,idet) = 4
    if ((ele_orb_alpha(i) == 1) .and. (ele_orb_beta(i) == 0)) ele_conf(i,idet) = 3
    if ((ele_orb_alpha(i) == 0) .and. (ele_orb_beta(i) == 1)) ele_conf(i,idet) = 2
    if ((ele_orb_alpha(i) == 0) .and. (ele_orb_beta(i) == 0)) ele_conf(i,idet) = 1
    !write(u6,*) 'ele_conf(i,idet)',ele_conf(i,idet)
  end do
  do i=1,norb
    write(118,'(I1)',advance='no') ele_conf(i,idet)
    !write(u6,*) 'ele_conf(',i,',idet)',ele_conf(i,idet)
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
  !write(u6,*) 'dV(i,idet)',i
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
    read(118,*) idx_det,dV(i,idx_det)
    !write(u6,*) idx_det,dV(i,idx_det)
  end do
  close(118)
  ! off-diagional multiply 2
  do idet=1,ndets_RGLR
    if (itype(idet) == 1) then
      dV(i,idet) = -dV(i,idet)
    else
      dtmp = sqrt((dV(i,idet)**2)*Two)
      dV(i,idet) = -sign(dtmp,dV(i,idet))
    end if
    !write(u6,*) 'i,idet,dV',i,idet,dV(i,idet)
  end do
end do

indexSD(1:ndets_RGLR) = [(i,i=1,ndets_RGLR)] ! inum(i)*isgn(i)
indexSD(ndets_RGLR+1:) = 0
vector(1:ndets_RGLR) = dV(istate,:)
vector(ndets_RGLR+1:) = Zero
dtmp = sum(vector(1:ndets_RGLR)**2)

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
call mma_deallocate(dV)
call mma_deallocate(ele_conf)
call mma_deallocate(electron)
call mma_deallocate(isgn)
call mma_deallocate(inum)
call mma_deallocate(itype)

end subroutine ci_reconstruct
