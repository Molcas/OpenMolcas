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

subroutine setup_cho(nSym,nIsh,nAsh,nSsh,NumCho,mode)
! -------------------------
! This subroutine uses the input variables to compute
! Stuff(:)%Unit(:), Stuff(:)%ip(:), Stuff(:)%np(:), Stuff%sp(:)
! nisplit, nasplit,lsplit, each dimensioned (1:nsym),
! in module chocaspt2,
! and allocates fields %Unit(:), %np(:), and %sp(:)
! each with sizes lsplit(1:nsym), and
! %ip(:), with teh sizes nsym*lsplit(1:nsym)
! -------------------------

use stdalloc, only: mma_MaxDBLE
use ChoCASPT2, only: Stuff, lsplit, nisplit, nasplit, nksh, nkes, npsh, npes
use stdalloc, only: mma_allocate, mma_deallocate
use constants, only: Two
use definitions, only: iwp, wp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nIsh(nSym), nAsh(nSym), nSsh(nSym)
integer(kind=iwp), intent(in) :: NumCho(nSym)
character(len=*), intent(in) :: mode
integer(kind=iwp) iAorb(8), iIorb(8), iKorb(8)
character(len=4) modecopy
integer(kind=iwp), external :: cho_irange
integer(kind=iwp) I, J, IfTest, iK, iK1, iK2, ioff, iS, iSP, iSym, iw, jfrac, jIAc, jS, jSym, kEnd, kEndSym, kFrac, kS, kSta, &
                  kStaSym, lS, mDiff, Mem1, MemMx, mRHS, nA, nAO, nI, nIAc, nIO, nKsp, nkSum, nMin, nO, nOkrb, nOrb, nP, nPMax, &
                  nPOrb, npSum, nS
real(kind=wp) xmb, xMemMx
! Statement functions
integer(kind=iwp) MulD2h, nIc, nAc
! **********************************************************************
nIc(i,j) = Stuff(j)%sp(i)
!*****
nAc(i,j) = Stuff(j)%sp(nisplit(j)+i)
!*****
MulD2h(i,j) = ieor(i-1,j-1)+1
! **********************************************************************

#ifdef _DEBUGPRINT_
IFTEST = 1
#else
IFTEST = 0
#endif

modecopy = mode(1:4)
call UpCase(modecopy)
if (modecopy == 'FREE') then
  do jSym=1,nSym
    if (NumCho(jSym) > 0) then
      call mma_deallocate(Stuff(jSym)%Unit)
      call mma_deallocate(Stuff(jSym)%ip)
      call mma_deallocate(Stuff(jSym)%np)
      call mma_deallocate(Stuff(jSym)%sp)
    end if
  end do
  return
end if

do i=1,8
  lsplit(i) = 0
  nisplit(i) = 0
  nasplit(i) = 0
  nksh(i) = 0
  nkes(i) = 0
  npsh(i) = 0
  npes(i) = 0
end do

! PAM07: New arrays: 'k shells' = inactive+active orbitals
! 'p shells' = active+secondary orbitals
!  nksh(isym)= Nr of shells by symmetry, nkes=nr of k shells in Earlier Symmetries:
nksum = 0
npsum = 0
do isym=1,nsym
  ni = nish(isym)
  na = nash(isym)
  ns = nssh(isym)
  nksh(isym) = ni+na
  npsh(isym) = na+ns
  nkes(isym) = nksum
  npes(isym) = npsum
  nksum = nksum+nksh(isym)
  npsum = npsum+npsh(isym)
end do

! Local arrays:
! iIorb(iSym) = Nr of inactive orbitals in earlier symmetries.
! iAorb(iSym) = Nr of   active orbitals in earlier symmetries.
! iKorb(iSym) = iIorb(iSym) + iAorb(iSym)
iIorb(1) = 0
iAorb(1) = 0
iKorb(1) = 0
do iSym=2,nSym
  iIorb(iSym) = iIorb(iSym-1)+nIsh(iSym-1)
  iAorb(iSym) = iAorb(iSym-1)+nAsh(iSym-1)
  iKorb(iSym) = iIorb(iSym)+iAorb(iSym)
end do

! nIO   = total number of inactive orbitals.
! nAO   = total number of active orbitals.
! nOkrb = total number of inactive+active orbitals.
! nOrb  = total number of orbitals
nOkrb = 0
nIO = 0
nAO = 0
nO = 0
do iSym=1,nSym
  nOkrb = nOkrb+nIsh(iSym)+nAsh(iSym)
  nOrb = nIsh(iSym)+nAsh(iSym)+nSsh(iSym)
  nO = nO+nOrb*(nOrb+1)/2
  nIO = nIO+nIsh(iSym)
  nAO = nAO+nAsh(iSym)
end do

!xO = dble(nO)

! Task: set up the necessary administration of transformed Cholesky
! vectors. Each such vector is characterized by the composite symmetry
! label JSYM. Vectors with common JSYM are stored on the same unit.
! For a given JSYM, there are NumCho(jSym) vectors. These are separated
! into subsets -- batches -- There will be nISplit(jSym) batches with
! vectors with a fixed inactive orbital index (Inactive vectors), and
! nASplit(jSym) batches with vectors with a fixed active index
! (Active vectors). Each batch has vectors with fixed index within a
! range of size Stuff(jSym)%sp(i) where i is in 1..nISplit(jSym) or
!  in nISplit(jSym)+1..nISplit(jSym)+nASplit(jSym).
! The vectors themselves have a size of Stuff(jSym)%np(i)
! The vectors are generally referred to as e.g. p#k, where #k stands
! for the fixed inactive or active orbital, while p is all the orbitals
! with a symmetry such that the compound symmetry label is JSYM.
! There will be a (possibly large) number of such vectors, each indexed
! with 'J' generally. This is the summation index in the definition of
! the cholesky vectors.
do jSym=1,nSym
  if (NumCho(jSym) < 1) goto 99

  call mma_MaxDBLE(MemMx)
  xMemMx = dble(MemMx)
  ! xMemMx = largest allocatable field.

  jfrac = 1

10 continue
  Mem1 = 2*NumCho(jSym)/jfrac  ! hold 2 vectors in memory
  if ((Mem1 == 0) .and. (jfrac > 1)) then
    write(6,*) ' Setup_cho fails to set up the data structures'
    write(6,*) ' used for the Cholesky vectors.'
    write(6,*) ' Too little memory is available at this point.'
    write(6,*) ' Details:'
    xmb = xMemMx/1048576.0e0_wp
    write(6,'(1x,a,1x,f10.3)') ' Largest contiguous allocatable memory (MB):',xmb
    xmb = Two*dble(NUMCHO(JSYM))/1048576.0e0_wp
    write(6,'(1x,a,1x,f10.3)') '                        2*NumCho(jSym) (MB):',xmb
    write(6,*) ' Divided up on jFrac pieces. jFrac=',jFrac
    write(6,*) ' If this seems odd, please tell Molcas programmers.'
    write(6,*) ' Right now, the allocated memory is:'
    call Cho_x_Quit('setup_cho',': Sorry! Too little memory!!',' ')
  end if

  kfrac = 0
  kfrac = kfrac+1
  ! Subdivide the set of orbitals into kfrac pieces of size <= nKsp
  nKsp = nOkrb/kfrac
  ! If kfrac has grown too large, start again using kfrac=1 but larger jfrac
  if (nKsp == 0) then
    jfrac = jfrac+1
    goto 10
  end if

  ! Loop over pieces of the orbital set:
  nPmax = 0
  do iK1=1,nOkrb,nKsp
    iK2 = max(iK1-1+nKsp,nOkrb)
    ! For this piece, compute the number of (p,k) orbital pairs
    ! with p active or secondary, and with joint symmetry jSym:
    nP = 0
    do iK=iK1,iK2
      kS = cho_irange(iK,iKorb,nSym,.false.)
      jS = MulD2h(kS,jSym)
      nP = nP+nAsh(jS)+nSsh(jS)
    end do
    ! nPmax will be the largest number of such pairs in any piece.
    nPmax = max(nPmax,nP)
  end do

  ! PAM:The following looks strange, since the max is over the same set
  ! of numbers no matter what jSym is.
  mRHS = 0
  do iS=1,nSym
    jS = MulD2h(iS,jSym)
    mRHS = max(mRHS,max(nAsh(jS),nSsh(jS)))
  end do

  ! Conversion to real(kind=wp) to avoid integer overflow on 32-bit machines

  ! PAM:Why would this be 'mem for right-hand side?'
  !xRHS = dble(mRHS**2)           ! mem. for right hand side
  !xLpk = dble(Mem1*nPmax*nKsp)   ! store Cholesky MO vectors
  !xPIQK = dble((nPmax*nKsp)**2)   ! store integrals
  !xmNeed = xO+xPIQK+max(xLpk,Two*xRHS) ! Fmat+integrals+rhs

  ! This also looks strange -- nIAc=all the inact+act orbitals no matter what.?
  nIAc = nOkrb

  jIAc = nIAc
  lsplit(jSym) = 1
  do while (jIAc < nOkrb)
    lsplit(jSym) = lsplit(jSym)+1
    jIAc = jIAc+nIAc
  end do

  ! What does this mean? Simply that lsplit(jSym) will be 'the next'
  ! of something?
  lsplit(jSym) = lsplit(jSym)+1 ! separate out Ina from Act

  ! Need xmNeedNow units of memory
  !xmNeedNow = xmNeed+dble((3+nSym)*lsplit(jSym))

  ! Allocate arrays, in all (3+nSym)*lsplit(jSym) elements:
  call mma_allocate(Stuff(jsym)%sp,lsplit(jSym),Label='%sp')
  call mma_allocate(Stuff(jsym)%np,lsplit(jSym),Label='%np')
  call mma_allocate(Stuff(jsym)%ip,nSym*lsplit(jSym),Label='%ip')
  call mma_allocate(Stuff(jsym)%Unit,lsplit(jSym),Label='%Unit')
  do i=1,lsplit(jSym)
    Stuff(jSym)%sp(i) = 0
    Stuff(jSym)%np(i) = 0
    Stuff(jSym)%Unit(i) = 0
    do j=1,nsym
      Stuff(jsym)%ip(j+nsym*(i-1)) = 0
    end do
  end do

  nmin = min(nIO,nIAc)
  jS = 0
  if (nmin > 0) jS = nIO/nmin
  do i=0,jS-1
    Stuff(jSym)%sp(1+i) = nMin
  end do
  mDiff = jS*nmin-nIO
  if (mDiff > 0) then
    Stuff(jSym)%sp(1+jS) = mDiff
    jS = jS+1
  end if
  nisplit(jSym) = jS

  nmin = min(nAO,nIAc)
  jS = 0
  if (nmin > 0) jS = nAO/nmin
  do i=0,jS-1
    Stuff(jSym)%sp(1+nisplit(jSym)+i) = nmin
  end do
  mDiff = jS*nmin-nAO
  if (mDiff > 0) then
    Stuff(jSym)%sp(1+nisplit(jSym)+jS) = mDiff
    jS = jS+1
  end if
  nasplit(jSym) = jS

  ! Note: The inactive orbitals will be partitioned.
  ! The isp=1..nisplit(jSym) is counter of the partitions.
  ! ioff will be the number of inactive orbitals in earlier partitions.
  ioff = 0
  do isp=1,nisplit(jSym)
    ! ioff+1 is the first orbital of partition isp.
    lS = cho_irange(ioff+1,iIorb,nSym,.false.)
    ! lS is its symmetry.
    iK = nIc(isp,jSym)+ioff
    ! Note: nIc(isp,jSym)=Stuff(jSym)%sp(isp)
    ! iK is the last orbital of the partition, and kS its symmetry
    kS = cho_irange(iK,iIorb,nSym,.false.)

    nPorb = 0
    ! Loop over the symmetry range of inactive orbitals in the partition
    do iS=lS,kS
      ! jS is then the symmetry of its companion orbital in the pair.
      jS = MulD2h(iS,jSym)
      ! ipip(jSym) is pointer to an array dimensioned iP(nSym,nisplit(jSym))
      ! which is used for offsets. In some other array, after the position
      ! iP(jS,isp) follows space for nAsh(jS)+nSsh(jS) items.
      Stuff(jSym)%ip(jS+nsym*(isp-1)) = nPorb
      ! Update nPorb.
      nPorb = nPorb+nAsh(jS)+nSsh(jS)
    end do

    ! ipnp(jSym) is pointer to an array dimensioned nP(nisplit(jSym))
    ! It gives the total size for the items mentioned above.
    Stuff(jSym)%np(isp) = nPorb
    ioff = ioff+nIc(isp,jSym)
  end do

  ! Here follows similar arrays for the active orbitals.
  ! Addressing is the same, except we use iP(jS,nisplit(jSym)+isp)
  ! and nP(nisplit(jSym)+isp)
  ioff = 0
  do isp=1,nasplit(jSym)
    lS = cho_irange(ioff+1,iAorb,nSym,.false.)
    iw = nAc(isp,jSym)+ioff
    kS = cho_irange(iw,iAorb,nSym,.false.)
    nPorb = 0
    do iS=lS,kS
      jS = MulD2h(iS,jSym)
      Stuff(jSym)%ip(jS+nsym*(nisplit(jSym)+isp-1)) = nPorb
      nPorb = nPorb+nAsh(jS)+nSsh(jS)
    end do
    Stuff(jSym)%np(nisplit(jSym)+isp) = nPorb
    ioff = ioff+nAc(isp,jSym)
  end do

99 continue

end do

if (iftest /= 0) then
  write(u6,*)
  write(u6,*) ' setup_cho report:'
  write(u6,*)
  write(u6,'(1x,a,8i4)') ' Inactive :',(nIsh(isym),isym=1,nSym)
  write(u6,'(1x,a,8i4)') ' Active   :',(nAsh(isym),isym=1,nSym)
  write(u6,'(1x,a,8i4)') ' Secondary:',(nSsh(isym),isym=1,nSym)
  write(u6,*)
  write(u6,'(1x,a,8i4)') ' NumCho   :',(NumCho(isym),isym=1,nSym)
  write(u6,*)
  write(u6,*) ' Partition  Fixed orbitals       nPorb space'
  do jsym=1,nsym
    write(u6,*) ' Symm:',jSym
    kend = 0
    do isp=1,nisplit(jsym)
      ksta = kend+1
      kend = kend+Stuff(jSym)%sp(isp)
      kstasym = cho_irange(ksta,iIorb,nSym,.false.)
      kendsym = cho_irange(kend,iIorb,nSym,.false.)
      nPorb = Stuff(jSym)%np(isp)
      write(u6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)') isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
      write(u6,'(1x,a,8i8)') ' iP Offsets:',(Stuff(jSym)%ip(nSym*(isp-1)+MulD2h(jSym,iS)),iS=1,nSym)
      write(u6,*)
    end do
    kend = 0
    do isp=1,nasplit(jsym)
      ksta = kend+1
      kend = kend+Stuff(jSym)%sp(nisplit(jSym)+isp)
      kstasym = cho_irange(ksta,iAorb,nSym,.false.)
      kendsym = cho_irange(kend,iAorb,nSym,.false.)
      nPorb = Stuff(jSym)%np(nisplit(jSym)+isp)
      write(u6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)') isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
      write(u6,'(1x,a,8i8)') ' iP Offsets:', (Stuff(jSym)%ip(nSym*(nisplit(jSym)+isp-1)+MulD2h(jSym,iS)),iS=1,nSym)
      write(u6,*)
    end do
  end do
end if

end subroutine setup_cho
