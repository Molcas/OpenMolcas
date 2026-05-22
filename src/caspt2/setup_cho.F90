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
! Unt, nisplit, nasplit, lsplit
! each dimensioned (1:nsym),
! in module chocaspt2
! -------------------------

use Symmetry_Info, only: Mul
use Data_Structures, only: Alloc1DiArray_Type
use ChoCASPT2, only: nasplit, nisplit, nksh, npsh, Unt
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nIsh(nSym), nAsh(nSym), nSsh(nSym), NumCho(nSym)
character(len=4), intent(in) :: mode
integer(kind=iwp) :: iAorb(8), IfTest, iIorb(8), iK, iK1, iK2, iKorb(8), ioff, iS, iSP, iSym, iw, jfrac, jIAc, jS, jSym, kEnd, &
                     kEndSym, kFrac, kS, kSta, kStaSym, lS, lsplit(8), mDiff, Mem1, MemMx, mRHS, nAO, nIAc, nIO, nKsp, nMin, &
                     nOkrb, nP, nPMax, nPOrb
real(kind=wp) :: xmb, xMemMx
type(Alloc1DiArray_Type) :: ip(8), mp(8), sp(8)
integer(kind=iwp), external :: cho_irange

#ifdef _DEBUGPRINT_
IFTEST = 1
#else
IFTEST = 0
#endif

if (mode == 'FREE') then
  do jSym=1,nSym
    call mma_deallocate(Unt(jSym)%A,safe='*')
  end do
  return
end if

lsplit(:) = 0
nisplit(:) = 0
nasplit(:) = 0
nksh(:) = 0
npsh(:) = 0

! PAM07: New arrays: 'k shells' = inactive+active orbitals
! 'p shells' = active+secondary orbitals
!  nksh(isym)= Nr of shells by symmetry:
nksh(1:nsym) = nish(1:nsym)+nash(1:nsym)
npsh(1:nsym) = nash(1:nsym)+nssh(1:nsym)

! Local arrays:
! iIorb(iSym) = Nr of inactive orbitals in earlier symmetries.
! iAorb(iSym) = Nr of   active orbitals in earlier symmetries.
! iKorb(iSym) = iIorb(iSym) + iAorb(iSym)
iIorb(1) = 0
iAorb(1) = 0
do iSym=2,nSym
  iIorb(iSym) = iIorb(iSym-1)+nIsh(iSym-1)
  iAorb(iSym) = iAorb(iSym-1)+nAsh(iSym-1)
end do
iKorb(1:nSym) = iIorb(1:nSym)+iAorb(1:nSym)

! nIO   = total number of inactive orbitals.
! nAO   = total number of active orbitals.
! nOkrb = total number of inactive+active orbitals.
! nOrb  = total number of orbitals
nIO = sum(nIsh(1:nSym))
nAO = sum(nAsh(1:nSym))
nOkrb = nIO+nAO

!xO = real(nO,kind=wp)

! Task: set up the necessary administration of transformed Cholesky
! vectors. Each such vector is characterized by the composite symmetry
! label JSYM. Vectors with common JSYM are stored on the same unit.
! For a given JSYM, there are NumCho(jSym) vectors. These are separated
! into subsets -- batches -- There will be nISplit(jSym) batches with
! vectors with a fixed inactive orbital index (Inactive vectors), and
! nASplit(jSym) batches with vectors with a fixed active index
! (Active vectors). Each batch has vectors with fixed index within a
! range of size sp(jSym)%A(i) where i is in 1..nISplit(jSym) or
!  in nISplit(jSym)+1..nISplit(jSym)+nASplit(jSym).
! The vectors themselves have a size of mp(jSym)%A(i)
! The vectors are generally referred to as e.g. p#k, where #k stands
! for the fixed inactive or active orbital, while p is all the orbitals
! with a symmetry such that the compound symmetry label is JSYM.
! There will be a (possibly large) number of such vectors, each indexed
! with 'J' generally. This is the summation index in the definition of
! the cholesky vectors.
do jSym=1,nSym
  if (NumCho(jSym) < 1) cycle

  call mma_MaxDBLE(MemMx)
  xMemMx = real(MemMx,kind=wp)
  ! xMemMx = largest allocatable field.

  jfrac = 1

  do
    Mem1 = 2*NumCho(jSym)/jfrac  ! hold 2 vectors in memory
    if ((Mem1 == 0) .and. (jfrac > 1)) then
      write(u6,*) ' Setup_cho fails to set up the data structures'
      write(u6,*) ' used for the Cholesky vectors.'
      write(u6,*) ' Too little memory is available at this point.'
      write(u6,*) ' Details:'
      xmb = xMemMx/(Two**20)
      write(u6,'(1x,a,1x,f10.3)') ' Largest contiguous allocatable memory (MB):',xmb
      xmb = Two*real(NUMCHO(JSYM),kind=wp)/(Two**20)
      write(u6,'(1x,a,1x,f10.3)') '                        2*NumCho(jSym) (MB):',xmb
      write(u6,*) ' Divided up on jFrac pieces. jFrac=',jFrac
      write(u6,*) ' If this seems odd, please tell Molcas programmers.'
      write(u6,*) ' Right now, the allocated memory is:'
      call Cho_x_Quit('setup_cho',': Sorry! Too little memory!!',' ')
    end if

    kfrac = 1
    ! Subdivide the set of orbitals into kfrac pieces of size <= nKsp
    nKsp = nOkrb/kfrac
    ! If kfrac has grown too large, start again using kfrac=1 but larger jfrac
    if (nKsp /= 0) exit
    jfrac = jfrac+1
  end do

  ! Loop over pieces of the orbital set:
  nPmax = 0
  do iK1=1,nOkrb,nKsp
    iK2 = max(iK1-1+nKsp,nOkrb)
    ! For this piece, compute the number of (p,k) orbital pairs
    ! with p active or secondary, and with joint symmetry jSym:
    nP = 0
    do iK=iK1,iK2
      kS = cho_irange(iK,iKorb,nSym,.false.)
      jS = Mul(kS,jSym)
      nP = nP+nAsh(jS)+nSsh(jS)
    end do
    ! nPmax will be the largest number of such pairs in any piece.
    nPmax = max(nPmax,nP)
  end do

  ! PAM:The following looks strange, since the max is over the same set
  ! of numbers no matter what jSym is.
  mRHS = 0
  do iS=1,nSym
    jS = Mul(iS,jSym)
    mRHS = max(mRHS,max(nAsh(jS),nSsh(jS)))
  end do

  ! Conversion to real(kind=wp) to avoid integer overflow on 32-bit machines

  ! PAM:Why would this be 'mem for right-hand side?'
  !xRHS = real(mRHS**2,kind=wp)          ! mem. for right hand side
  !xLpk = real(Mem1*nPmax*nKsp,kind=wp)  ! store Cholesky MO vectors
  !xPIQK = real((nPmax*nKsp)**2,kind=wp) ! store integrals
  !xmNeed = xO+xPIQK+max(xLpk,Two*xRHS)  ! Fmat+integrals+rhs

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
  !xmNeedNow = xmNeed+real((3+nSym,kind=wp)*lsplit(jSym))

  ! Allocate arrays, in all (3+nSym)*lsplit(jSym) elements:
  call mma_allocate(sp(jsym)%A,lsplit(jSym),Label='sp')
  call mma_allocate(mp(jsym)%A,lsplit(jSym),Label='np')
  call mma_allocate(ip(jsym)%A,nSym*lsplit(jSym),Label='ip')
  call mma_allocate(Unt(jsym)%A,lsplit(jSym),Label='Unt')
  sp(jSym)%A(:) = 0
  mp(jSym)%A(:) = 0
  Unt(jSym)%A(:) = 0
  ip(jsym)%A(:) = 0

  nmin = min(nIO,nIAc)
  jS = 0
  if (nmin > 0) jS = nIO/nmin
  sp(jSym)%A(1:jS) = nMin
  mDiff = jS*nmin-nIO
  if (mDiff > 0) then
    sp(jSym)%A(jS+1) = mDiff
    jS = jS+1
  end if
  nisplit(jSym) = jS

  nmin = min(nAO,nIAc)
  jS = 0
  if (nmin > 0) jS = nAO/nmin
  sp(jSym)%A(nisplit(jSym)+1:nisplit(jSym)+jS) = nmin
  mDiff = jS*nmin-nAO
  if (mDiff > 0) then
    sp(jSym)%A(nisplit(jSym)+jS+1) = mDiff
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
    iK = sp(jSym)%A(isp)+ioff
    ! iK is the last orbital of the partition, and kS its symmetry
    kS = cho_irange(iK,iIorb,nSym,.false.)

    nPorb = 0
    ! Loop over the symmetry range of inactive orbitals in the partition
    do iS=lS,kS
      ! jS is then the symmetry of its companion orbital in the pair.
      jS = Mul(iS,jSym)
      ! ipip(jSym) is pointer to an array dimensioned iP(nSym,nisplit(jSym))
      ! which is used for offsets. In some other array, after the position
      ! iP(jS,isp) follows space for nAsh(jS)+nSsh(jS) items.
      ip(jSym)%A(jS+nsym*(isp-1)) = nPorb
      ! Update nPorb.
      nPorb = nPorb+nAsh(jS)+nSsh(jS)
    end do

    ! ipnp(jSym) is pointer to an array dimensioned nP(nisplit(jSym))
    ! It gives the total size for the items mentioned above.
    mp(jSym)%A(isp) = nPorb
    ioff = ioff+sp(jSym)%A(isp)
  end do

  ! Here follows similar arrays for the active orbitals.
  ! Addressing is the same, except we use iP(jS,nisplit(jSym)+isp)
  ! and nP(nisplit(jSym)+isp)
  ioff = 0
  do isp=1,nasplit(jSym)
    lS = cho_irange(ioff+1,iAorb,nSym,.false.)
    iw = sp(jSym)%A(nisplit(jSym)+isp)
    kS = cho_irange(iw,iAorb,nSym,.false.)
    nPorb = 0
    do iS=lS,kS
      jS = Mul(iS,jSym)
      ip(jSym)%A(jS+nsym*(nisplit(jSym)+isp-1)) = nPorb
      nPorb = nPorb+nAsh(jS)+nSsh(jS)
    end do
    mp(jSym)%A(nisplit(jSym)+isp) = nPorb
    ioff = ioff+sp(jSym)%A(nisplit(jSym)+isp)
  end do

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
      kend = kend+sp(jSym)%A(isp)
      kstasym = cho_irange(ksta,iIorb,nSym,.false.)
      kendsym = cho_irange(kend,iIorb,nSym,.false.)
      nPorb = mp(jSym)%A(isp)
      write(u6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)') isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
      write(u6,'(1x,a,8i8)') ' iP Offsets:',(ip(jSym)%A(nSym*(isp-1)+Mul(jSym,iS)),iS=1,nSym)
      write(u6,*)
    end do
    kend = 0
    do isp=1,nasplit(jsym)
      ksta = kend+1
      kend = kend+sp(jSym)%A(nisplit(jSym)+isp)
      kstasym = cho_irange(ksta,iAorb,nSym,.false.)
      kendsym = cho_irange(kend,iAorb,nSym,.false.)
      nPorb = mp(jSym)%A(nisplit(jSym)+isp)
      write(u6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)') isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
      write(u6,'(1x,a,8i8)') ' iP Offsets:', (ip(jSym)%A(nSym*(nisplit(jSym)+isp-1)+Mul(jSym,iS)),iS=1,nSym)
      write(u6,*)
    end do
  end do
end if

do jSym=1,nSym
  call mma_deallocate(sp(jSym)%A,safe='*')
  call mma_deallocate(mp(jSym)%A,safe='*')
  call mma_deallocate(ip(jSym)%A,safe='*')
end do

end subroutine setup_cho
