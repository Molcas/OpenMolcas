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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The description of a CI space common to all CI back ends (Lucia, SGUGA,
!> Faroald): the active orbitals in nPart partitions (CAS: 1, RAS: 3, GAS:
!> any), with the cumulative bounds MinEl(k) <= N(1..k) <= MaxEl(k) on the
!> number of electrons in the partitions 1..k, and configurations removed
!> on top of that (HEXS/DEXS): those with at least RemOcc(r) electrons in
!> partition RemPart(r). The levels (active orbitals of the CI) are ordered
!> by partition, within a partition by symmetry (LevPart, LevSym).
!>
!> The space is filled by the CI interface from the input of the calling
!> program and handed to the back ends at their initialization; this module
!> is stateless (no program data), so that both sides may use it.
module CI_Space_Def

use Molcas, only: MxGAS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
private

type CI_Space
  integer(kind=iwp) :: nSym = 0, nActEl = 0, iSpin = 0, nPart = 0, nLev = 0, nRem = 0
  ! orbitals per symmetry and partition
  integer(kind=iwp), allocatable :: nOrb(:,:)
  ! cumulative bounds on the electrons in the partitions 1..k
  integer(kind=iwp), allocatable :: MinEl(:), MaxEl(:)
  ! removed configurations: occupation of partition RemPart(r) >= RemOcc(r)
  integer(kind=iwp), allocatable :: RemPart(:), RemOcc(:)
  ! per level: partition and symmetry
  integer(kind=iwp), allocatable :: LevPart(:), LevSym(:)
end type CI_Space

! The removals as requested in the input (RASSCF keywords HEXS and DEXS),
! numbered by the input spaces (GAS spaces, or RAS1..RAS3 in the GAS form):
! iHEXS(1:nHEXS) lose their configurations of maximum occupation, iDEXS(1:nDEXS)
! also those of maximum occupation - 1 (CI_Space_Setup turns them into the
! removal rules of the CI space).
type CI_Space_Input
  integer(kind=iwp) :: nHEXS = 0, nDEXS = 0
  integer(kind=iwp) :: iHEXS(MxGAS) = 0, iDEXS(MxGAS) = 0
end type CI_Space_Input

public :: CI_Space, CI_Space_Input, Space_Allowed, Space_Free, Space_From_GAS, Space_Max_Occ, Space_Print, Space_Remove

contains

!> The space from the GAS description of the input: nGAS spaces with
!> nGSSH(g,s) orbitals of symmetry s, cumulative bounds IGSOcc(g,1:2) on
!> the electrons in the spaces 1..g. Spaces without orbitals are dropped
!> (their bounds are merged into the previous space), the bounds are
!> tightened to the physically possible ones, and the last partition holds
!> all nActEl electrons. Partitions with orbitals are kept as they are, so
!> the level order (by partition, then symmetry) is that of the input.
subroutine Space_From_GAS(Space,nSym,nActEl,iSpin,nGAS,nGSSH,IGSOcc)
  type(CI_Space), intent(inout) :: Space
  integer(kind=iwp), intent(in) :: nSym, nActEl, iSpin, nGAS, nGSSH(:,:), IGSOcc(:,:)
  integer(kind=iwp) :: g, k, s, iLev, nAbove, nBelow
  integer(kind=iwp), allocatable :: nO(:,:), Mn(:), Mx(:)

  call Space_Free(Space)
  Space%nSym = nSym
  Space%nActEl = nActEl
  Space%iSpin = iSpin
  call mma_allocate(nO,nSym,max(1,nGAS),Label='Space_nO')
  call mma_allocate(Mn,max(1,nGAS),Label='Space_Mn')
  call mma_allocate(Mx,max(1,nGAS),Label='Space_Mx')
  ! the partitions with orbitals; the bounds of an empty space restrict
  ! the same partial sum as those of the previous partition
  k = 0
  do g=1,nGAS
    if (sum(nGSSH(g,1:nSym)) == 0) then
      if (k == 0) then
        ! before the first partition the partial sum is zero
        if ((IGSOcc(g,1) > 0) .or. (IGSOcc(g,2) < 0)) call Inconsistent(g)
      else
        Mn(k) = max(Mn(k),IGSOcc(g,1))
        Mx(k) = min(Mx(k),IGSOcc(g,2))
      end if
      cycle
    end if
    k = k+1
    nO(1:nSym,k) = nGSSH(g,1:nSym)
    Mn(k) = IGSOcc(g,1)
    Mx(k) = IGSOcc(g,2)
  end do
  Space%nPart = k
  call mma_allocate(Space%nOrb,nSym,max(1,k),Label='Space%nOrb')
  call mma_allocate(Space%MinEl,max(1,k),Label='Space%MinEl')
  call mma_allocate(Space%MaxEl,max(1,k),Label='Space%MaxEl')
  Space%nOrb(:,1:k) = nO(:,1:k)
  ! the physically possible bounds: at most 2 per orbital of 1..k, at
  ! least nActEl minus 2 per orbital above k
  nBelow = 0
  do k=1,Space%nPart
    nBelow = nBelow+sum(Space%nOrb(:,k))
    nAbove = sum(Space%nOrb(:,k+1:Space%nPart))
    Space%MinEl(k) = max(Mn(k),0,nActEl-2*nAbove)
    Space%MaxEl(k) = min(Mx(k),nActEl,2*nBelow)
  end do
  if (Space%nPart > 0) then
    Space%MinEl(Space%nPart) = nActEl
    Space%MaxEl(Space%nPart) = nActEl
  end if
  do k=1,Space%nPart
    if (Space%MinEl(k) > Space%MaxEl(k)) call Inconsistent(k)
  end do
  ! the levels
  Space%nLev = sum(Space%nOrb(:,1:Space%nPart))
  call mma_allocate(Space%LevPart,max(1,Space%nLev),Label='Space%LevPart')
  call mma_allocate(Space%LevSym,max(1,Space%nLev),Label='Space%LevSym')
  iLev = 0
  do k=1,Space%nPart
    do s=1,nSym
      Space%LevPart(iLev+1:iLev+Space%nOrb(s,k)) = k
      Space%LevSym(iLev+1:iLev+Space%nOrb(s,k)) = s
      iLev = iLev+Space%nOrb(s,k)
    end do
  end do
  Space%nRem = 0
  call mma_deallocate(nO)
  call mma_deallocate(Mn)
  call mma_deallocate(Mx)

contains

  subroutine Inconsistent(j)
    integer(kind=iwp), intent(in) :: j
    write(u6,*) 'Space_From_GAS: inconsistent occupation bounds of space/partition',j
    call Abend()
  end subroutine Inconsistent

end subroutine Space_From_GAS

!> Remove the configurations with at least nOcc electrons in partition iPart.
subroutine Space_Remove(Space,iPart,nOcc)
  type(CI_Space), intent(inout) :: Space
  integer(kind=iwp), intent(in) :: iPart, nOcc
  integer(kind=iwp), allocatable :: P(:), O(:)
  integer(kind=iwp) :: n
  n = Space%nRem
  call mma_allocate(P,n+1,Label='Space_RemP')
  call mma_allocate(O,n+1,Label='Space_RemO')
  if (n > 0) then
    P(1:n) = Space%RemPart(1:n)
    O(1:n) = Space%RemOcc(1:n)
  end if
  P(n+1) = iPart
  O(n+1) = nOcc
  call mma_deallocate(Space%RemPart,safe='*')
  call mma_deallocate(Space%RemOcc,safe='*')
  call mma_allocate(Space%RemPart,n+1,Label='Space%RemPart')
  call mma_allocate(Space%RemOcc,n+1,Label='Space%RemOcc')
  Space%RemPart(:) = P(:)
  Space%RemOcc(:) = O(:)
  Space%nRem = n+1
  call mma_deallocate(P)
  call mma_deallocate(O)
end subroutine Space_Remove

!> Is a configuration with Occ(k) electrons in partition k in the space?
pure function Space_Allowed(Space,Occ) result(Allowed)
  type(CI_Space), intent(in) :: Space
  integer(kind=iwp), intent(in) :: Occ(:)
  logical(kind=iwp) :: Allowed
  integer(kind=iwp) :: k, n, r
  Allowed = .false.
  n = 0
  do k=1,Space%nPart
    n = n+Occ(k)
    if ((n < Space%MinEl(k)) .or. (n > Space%MaxEl(k))) return
  end do
  do r=1,Space%nRem
    if (Occ(Space%RemPart(r)) >= Space%RemOcc(r)) return
  end do
  Allowed = .true.
end function Space_Allowed

!> The maximum occupation of each partition, NMx(k) = maxa(k) + maxb(k):
!> the largest numbers of alpha and beta electrons (M_S = S) in partition
!> k over the pairs of alpha and beta partition occupations allowed by the
!> bounds (Lucia's rule for its determinant blocks; removals not applied).
subroutine Space_Max_Occ(Space,NMx)
  type(CI_Space), intent(in) :: Space
  integer(kind=iwp), intent(out) :: NMx(Space%nPart)
  integer(kind=iwp) :: na, nb, nP, iA, iB, nVA, nVB, k
  integer(kind=iwp), allocatable :: VA(:,:), VB(:,:), MxA(:), MxB(:), Lim(:)
  logical(kind=iwp) :: ok
  nP = Space%nPart
  na = (Space%nActEl+Space%iSpin-1)/2
  nb = Space%nActEl-na
  call mma_allocate(Lim,max(1,nP),Label='Space_Lim')
  do k=1,nP
    Lim(k) = sum(Space%nOrb(:,k))
  end do
  call Spin_Vectors(na,nVA,VA)
  call Spin_Vectors(nb,nVB,VB)
  call mma_allocate(MxA,max(1,nP),Label='Space_MxA')
  call mma_allocate(MxB,max(1,nP),Label='Space_MxB')
  MxA(:) = -1
  MxB(:) = -1
  do iA=1,nVA
    do iB=1,nVB
      ok = .true.
      do k=1,nP
        if ((sum(VA(1:k,iA)+VB(1:k,iB)) < Space%MinEl(k)) .or. (sum(VA(1:k,iA)+VB(1:k,iB)) > Space%MaxEl(k))) then
          ok = .false.
          exit
        end if
      end do
      if (.not. ok) cycle
      MxA(1:nP) = max(MxA(1:nP),VA(1:nP,iA))
      MxB(1:nP) = max(MxB(1:nP),VB(1:nP,iB))
    end do
  end do
  NMx(1:nP) = MxA(1:nP)+MxB(1:nP)
  call mma_deallocate(VA)
  call mma_deallocate(VB)
  call mma_deallocate(MxA)
  call mma_deallocate(MxB)
  call mma_deallocate(Lim)

contains

  ! all distributions V(1:nP,i) of n electrons of one spin over the
  ! partitions, at most Lim(k) in partition k
  subroutine Spin_Vectors(n,nV,V)
    integer(kind=iwp), intent(in) :: n
    integer(kind=iwp), intent(out) :: nV
    integer(kind=iwp), allocatable, intent(out) :: V(:,:)
    integer(kind=iwp) :: iPass, j, kk
    integer(kind=iwp), allocatable :: W(:)
    call mma_allocate(W,max(1,nP),Label='Space_W')
    do iPass=1,2
      nV = 0
      W(:) = 0
      do
        ! the current vector: the last partition takes the rest
        if (nP > 0) then
          W(nP) = n-sum(W(1:nP-1))
          if ((W(nP) >= 0) .and. (W(nP) <= Lim(nP))) then
            nV = nV+1
            if (iPass == 2) V(1:nP,nV) = W(1:nP)
          end if
        end if
        ! next vector over the partitions 1..nP-1 (odometer)
        kk = 0
        do j=1,nP-1
          if ((W(j) < Lim(j)) .and. (sum(W(1:nP-1)) < n)) then
            W(j) = W(j)+1
            kk = j
            exit
          end if
          W(j) = 0
        end do
        if (kk == 0) exit
      end do
      if (iPass == 1) call mma_allocate(V,max(1,nP),max(1,nV),Label='Space_V')
    end do
    call mma_deallocate(W)
  end subroutine Spin_Vectors

end subroutine Space_Max_Occ

subroutine Space_Print(Space)
  type(CI_Space), intent(in) :: Space
  integer(kind=iwp) :: k, r
  write(u6,'(1X,A,I0,A,I0,A,I0,A)') 'CI space: ',Space%nActEl,' electrons, multiplicity ',Space%iSpin,', ',Space%nPart, &
                                    ' partition(s)'
  do k=1,Space%nPart
    write(u6,'(3X,A,I0,A,8I4)') 'partition ',k,', orbitals per symmetry:',Space%nOrb(1:Space%nSym,k)
    write(u6,'(3X,A,I0,A,I0)') '  electrons in partitions 1..k: ',Space%MinEl(k),' - ',Space%MaxEl(k)
  end do
  do r=1,Space%nRem
    write(u6,'(3X,A,I0,A,I0,A)') 'removed: configurations with at least ',Space%RemOcc(r),' electrons in partition ', &
                                 Space%RemPart(r),' (HEXS/DEXS)'
  end do
end subroutine Space_Print

subroutine Space_Free(Space)
  type(CI_Space), intent(inout) :: Space
  call mma_deallocate(Space%nOrb,safe='*')
  call mma_deallocate(Space%MinEl,safe='*')
  call mma_deallocate(Space%MaxEl,safe='*')
  call mma_deallocate(Space%RemPart,safe='*')
  call mma_deallocate(Space%RemOcc,safe='*')
  call mma_deallocate(Space%LevPart,safe='*')
  call mma_deallocate(Space%LevSym,safe='*')
  Space%nPart = 0
  Space%nLev = 0
  Space%nRem = 0
end subroutine Space_Free

end module CI_Space_Def
