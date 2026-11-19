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
! Copyright (C) 2014, Steven Vancoillie                                *
!***********************************************************************

module citrans
! written by Steven Vancoillie, summer 2014.
!
! This module implements operations that are performed on a CI vector
! in order to transform it from one basis to another, in particular
! from CSFs (GUGA ordering) to determinants and back. In order to
! efficiently do this, all CSFs/determinants are grouped by
! configurations, i.e., a fixed set of doubly+singly occupied
! orbitals.
!
! The configurations are grouped by the number of double/single
! occupied orbitals.  Consider e.g. the case of 7in6, mult = 2.  We
! have n = 6 (#orbitals), N = 7 (#electrons), a = 4 (#alpha), b = 3
! (#beta). The number of doubly occupied orbitals (d) varies from 1 to
! 3, the singly occupied orbitals (s) from 5 to 1. The orbitals are
! represented here in descending order such that the lexicographical
! lowest rank has the lowest orbitals occupied.
!
!                                                  / n \
! I used the notation nCk for binomial coefficient |   |
!                                                  \ k /
!
!       do configs  so configs  determinants  CSFs
! d s   nCd         n-dCs       sCa-d         sCa-d - sCa-d+1
! 1 5   1   000002  5C5 xxxxx
!       2   000020  5C5 xxxxx
!       ...
!       6C1 200000  5C5 xxxxx
!
! 2 3   1   000022  1    0xxx
!                   2    x0xx
!                   ...
!                   4C3  xxx0
!       2   000202  1    0xxx
!                   ...
!       ...         ...
!       6C2 220000  1    0xxx
!                   ...
!                   4C3  xxx0
! 3 1   1   000222  1     00x
!                   ...
!                   3C1   x00
!       ...         ...
!       6C3 222000  1     00x
!                   ...
!                   3C1   x00
!
! So, a specific group d,s has a fixed number of doubly occupied
! strings and a fixed number of singly occupied strings per doubly
! occupied string. Each configuration is identified by its group,
! its doubly occupied rank, and its singly occupied rank. The ordering
! within a group is done by traversing the singly occupied strings
! first, i.e., the index of a configurations is given by the formula
! nsoc*(rankdo-1)+rankso, with nsoc the number of singly occupied
! strings per doubly occupied string in a group, i.e., n-dCs.

use CSF_Graph, only: Free_SG_Graph, MkSGUGA, SG_Graph, SG_Walk_Next, SG_Walks
use Molcas, only: MxGAS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
private

! global module variables

integer(kind=iwp) :: ndo_min, ndo_max
! number of doubly/singly occupied configurations per group
integer(kind=iwp), allocatable :: ndoc_group(:), nsoc_group(:)
! number of determinants/CSFs per group
integer(kind=iwp), allocatable :: ndet_group(:), ncsf_group(:)

type spintable
  integer(kind=iwp) :: ndet, ncsf
  real(kind=wp), allocatable :: coef(:,:)
end type spintable

type(spintable), allocatable :: spintabs(:)

! The split graph of the CSF space in the order of the CI vectors of the
! files (CSF_Graph), built by citrans itself, so that citrans_sort does not
! need SGUGA.
type(SG_Graph), target :: FileG
type(SG_Walks) :: FileW

! Symmetry: the CSF space holds the configurations of symmetry ci_lsym only
! (the symmetry of a configuration is the product of the symmetries of its
! singly occupied orbitals, orbsym). Only the kept configurations are
! stored: per group ido, k = 1..nkeep_group(ido), the doubly occupied
! string conf_doub, the singly occupied string conf_sing (over the orbitals
! not doubly occupied) and the key conf_key = nsoc*(rankdo-1)+rankso, at
! conf_off(ido)+k, sorted by the key (the order of the full CAS list).
! RAS/GAS: the orbitals are the levels of the split graph (partition by
! partition, within each by symmetry), orbras their partition; a
! configuration is kept if the partitions 1..k hold at least ci_nrasel(k)
! and at most ci_nrasmx(k) electrons (k < ci_nrsprt), as the graph of
! SGUGA. CAS: ci_nrsprt = 1.
integer(kind=iwp) :: ci_nsym = 1, ci_lsym = 1, ci_nrsprt = 1, part_pos(MxGAS) = 0, part_n(MxGAS) = 0
integer(kind=iwp), allocatable :: orbsym(:), orbras(:), ci_nrasel(:), ci_nrasmx(:), conf_off(:), conf_key(:), conf_doub(:), &
                                  conf_sing(:), nkeep_group(:)

public :: citrans_conf_dets, citrans_csf2ras, citrans_csf2sd, citrans_csf_rank, citrans_graph_free, citrans_graph_init, &
          citrans_ncsf, citrans_ras2csf, citrans_sd2csf, citrans_sort, citrans_sym_free, citrans_sym_init, ncsf_group, ndet_group, &
          ndo_max, ndo_min, ndoc_group, nsoc_group, spintable_create, spintabs, spintabs_allocate, spintabs_free

! Private extensions to mma interfaces

interface mma_allocate
  module procedure :: spt_mma_allo_1D, spt_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: spt_mma_free_1D
end interface

#include "compiler_features.h"

contains

!> The split graph of the CAS or RAS space (nel electrons in norb orbitals
!> = levels of the symmetries orbsym and partitions orbras, multiplicity
!> mult) in the file order, as SGUGA sets it up for RASSCF (SG_Init_Light:
!> IFRAS = 1 for three partitions, plus the number of symmetries with
!> active orbitals).
subroutine citrans_graph_init(nel,norb,mult)
  integer(kind=iwp), intent(in) :: nel, norb, mult
  integer(kind=iwp) :: isym, iprt
  call Free_SG_Graph(FileG,FileW)
  FileG%nSym = ci_nsym
  FileG%nActEl = nel
  FileG%iSpin = mult
  FileG%nLev = norb
  FileG%nRsPrt = ci_nrsprt
  FileG%nRas(:,:) = 0
  FileG%nRasEl(:) = 0
  FileG%nRasEl(1:ci_nrsprt) = ci_nrasel(1:ci_nrsprt)
  FileG%nRasMx(:) = -1
  FileG%nRasMx(1:ci_nrsprt) = ci_nrasmx(1:ci_nrsprt)
  FileG%IFRAS = 0
  if (ci_nrsprt == 3) FileG%IFRAS = 1
  do isym=1,ci_nsym
    do iprt=1,ci_nrsprt
      FileG%nRas(isym,iprt) = count((orbsym(1:norb) == isym) .and. (orbras(1:norb) == iprt))
    end do
    if (sum(FileG%nRas(isym,1:ci_nrsprt)) /= 0) FileG%IFRAS = FileG%IFRAS+1
  end do
  call mma_allocate(FileG%ISm,norb,Label='FileG%ISm')
  FileG%ISm(:) = orbsym(1:norb)
  call MkSGUGA(FileG,FileW)
end subroutine citrans_graph_init

!> The CSF space: nsym irreps, the symmetries osym of the norb active
!> orbitals, the target symmetry lsym, and the RAS/GAS restrictions (nrsprt
!> partitions, oras the partition of each orbital, nrasel and nrasmx the
!> minimum and maximum number of electrons in the partitions 1..k); builds
!> the configuration maps (after the group sizes of citrans_init). Without
!> the optional arguments C1 CAS.
subroutine citrans_sym_init(nel,norb,nsym,osym,lsym,nrsprt,oras,nrasel,nrasmx)
  integer(kind=iwp), intent(in) :: nel, norb
  integer(kind=iwp), intent(in), optional :: nsym, osym(norb), lsym, nrsprt, oras(norb), nrasel(:), nrasmx(:)
  integer(kind=iwp) :: ido, ntot

  call citrans_sym_free()
  ci_nsym = 1
  ci_lsym = 1
  ci_nrsprt = 1
  call mma_allocate(orbsym,max(1,norb),label='orbsym')
  call mma_allocate(orbras,max(1,norb),label='orbras')
  orbsym(:) = 1
  orbras(:) = 1
  if (present(nsym)) ci_nsym = nsym
  if (present(osym)) orbsym(1:norb) = osym(1:norb)
  if (present(lsym)) ci_lsym = lsym
  if (present(nrsprt)) ci_nrsprt = nrsprt
  if (ci_nrsprt > MxGAS) then
    write(u6,*) 'citrans: too many partitions',ci_nrsprt
    call Abend()
  end if
  call mma_allocate(ci_nrasel,ci_nrsprt,label='ci_nrasel')
  call mma_allocate(ci_nrasmx,ci_nrsprt,label='ci_nrasmx')
  ci_nrasel(:) = nel
  ci_nrasmx(:) = nel
  if (present(nrasel)) ci_nrasel(1:ci_nrsprt) = nrasel(1:ci_nrsprt)
  if (present(nrasmx)) ci_nrasmx(1:ci_nrsprt) = nrasmx(1:ci_nrsprt)
  if (present(oras)) orbras(1:norb) = oras(1:norb)

  call mma_allocate(conf_off,[ndo_min,ndo_max],label='conf_off')
  call mma_allocate(nkeep_group,[ndo_min,ndo_max],label='nkeep_group')
  call ras_partitions(norb)
  ! upper bound of the number of kept configurations (RAS, no symmetry)
  ntot = 0
  do ido=ndo_min,ndo_max
    ntot = ntot+gen_confs(nel,norb,ido,.false.,0)
  end do
  call mma_allocate(conf_key,max(1,ntot),label='conf_key')
  call mma_allocate(conf_doub,max(1,ntot),label='conf_doub')
  call mma_allocate(conf_sing,max(1,ntot),label='conf_sing')
  ntot = 0
  do ido=ndo_min,ndo_max
    conf_off(ido) = ntot
    nkeep_group(ido) = gen_confs(nel,norb,ido,.true.,ntot)
    call sort_confs(ntot,nkeep_group(ido))
    ntot = ntot+nkeep_group(ido)
  end do

end subroutine citrans_sym_init

subroutine citrans_sym_free()
  call mma_deallocate(orbsym,safe='*')
  call mma_deallocate(orbras,safe='*')
  call mma_deallocate(ci_nrasel,safe='*')
  call mma_deallocate(ci_nrasmx,safe='*')
  call mma_deallocate(conf_off,safe='*')
  call mma_deallocate(conf_key,safe='*')
  call mma_deallocate(conf_doub,safe='*')
  call mma_deallocate(conf_sing,safe='*')
  call mma_deallocate(nkeep_group,safe='*')
end subroutine citrans_sym_free

!> The partitions of the orbitals (contiguous, as the levels): first
!> orbital minus one, part_pos, and number of orbitals, part_n, of
!> partition k = 1..ci_nrsprt.
subroutine ras_partitions(norb)
  integer(kind=iwp), intent(in) :: norb
  integer(kind=iwp) :: iorb, k
  part_n(:) = 0
  do iorb=1,norb
    if ((orbras(iorb) < 1) .or. (orbras(iorb) > ci_nrsprt)) then
      write(u6,*) 'citrans: illegal RAS partition of orbital',iorb,orbras(iorb)
      call Abend()
    end if
    if (iorb > 1) then
      if (orbras(iorb) < orbras(iorb-1)) then
        write(u6,*) 'citrans: the RAS partitions are not contiguous'
        call Abend()
      end if
    end if
    part_n(orbras(iorb)) = part_n(orbras(iorb))+1
  end do
  part_pos(1) = 0
  do k=2,ci_nrsprt
    part_pos(k) = part_pos(k-1)+part_n(k-1)
  end do
end subroutine ras_partitions

!> The configurations of group ido (ido doubly, nel-2*ido singly occupied
!> orbitals) that fulfil the RAS/GAS restrictions, generated partition by
!> partition (never the full CAS list): the numbers d(k), sx(k) of doubly
!> and singly occupied orbitals of each partition run over all
!> distributions (odometer over the partitions 1..P-1, the last partition
!> takes the rest), and for each the patterns of the partitions are
!> combined (odometer over the patterns). Returns their number; with store,
!> those of symmetry ci_lsym are stored from conf_* (off+1) on, unsorted,
!> and their number is returned.
function gen_confs(nel,norb,ido,store,off) result(n)
  use second_quantization, only: binom_coef, lexrank
  integer(kind=iwp), intent(in) :: nel, norb, ido, off
  logical(kind=iwp), intent(in) :: store
  integer(kind=iwp) :: n
  integer(kind=iwp) :: iso, nsoc, np, k, ne, cnt, doub, socc, sing, ip(MxGAS), d(MxGAS), sx(MxGAS)
  integer(kind=iwp), allocatable :: lst(:,:,:), nl(:)
  logical(kind=iwp) :: ok, more

  n = 0
  iso = nel-2*ido
  nsoc = binom_coef(iso,norb-ido)
  np = ci_nrsprt
  d(:) = 0
  sx(:) = 0
  do
    ! the last partition takes the remaining orbitals
    d(np) = ido-sum(d(1:np-1))
    sx(np) = iso-sum(sx(1:np-1))
    ok = (d(np) >= 0) .and. (d(np) <= part_n(np)) .and. (sx(np) >= 0) .and. (sx(np) <= part_n(np)-d(np))
    ! RAS/GAS: ci_nrasel(k) to ci_nrasmx(k) electrons in the partitions 1..k
    ne = 0
    do k=1,np-1
      if (.not. ok) exit
      ne = ne+2*d(k)+sx(k)
      if ((ne < ci_nrasel(k)) .or. (ne > ci_nrasmx(k))) ok = .false.
    end do
    if (ok) then
      if (.not. store) then
        cnt = 1
        do k=1,np
          cnt = cnt*binom_coef(d(k),part_n(k))*binom_coef(sx(k),part_n(k)-d(k))
        end do
        n = n+cnt
      else
        ! the (doubly, singly) occupied patterns of each partition, all
        ! their combinations
        call part_lists(np,d,sx,lst,nl)
        if (all(nl(1:np) > 0)) then
          ip(1:np) = 1
          do
            doub = 0
            socc = 0
            do k=1,np
              doub = ior(doub,lst(1,ip(k),k))
              socc = ior(socc,lst(2,ip(k),k))
            end do
            sing = pext(socc,ibits(not(doub),0,norb))
            if (conf_sym(doub,sing,norb) == ci_lsym) then
              n = n+1
              conf_doub(off+n) = doub
              conf_sing(off+n) = sing
              conf_key(off+n) = nsoc*(lexrank(doub)-1)+lexrank(sing)
            end if
            more = .false.
            do k=1,np
              if (ip(k) < nl(k)) then
                ip(k) = ip(k)+1
                more = .true.
                exit
              end if
              ip(k) = 1
            end do
            if (.not. more) exit
          end do
        end if
        call mma_deallocate(lst)
        call mma_deallocate(nl)
      end if
    end if
    ! next distribution over the partitions 1..np-1
    more = .false.
    do k=1,np-1
      if ((sx(k) < part_n(k)-d(k)) .and. (sum(sx(1:np-1)) < iso)) then
        sx(k) = sx(k)+1
        more = .true.
        exit
      end if
      sx(k) = 0
      if ((d(k) < part_n(k)) .and. (sum(d(1:np-1)) < ido)) then
        d(k) = d(k)+1
        more = .true.
        exit
      end if
      d(k) = 0
    end do
    if (.not. more) exit
  end do

end function gen_confs

!> lst(1:2,i,k): doubly and singly occupied bits (whole orbital range) of
!> the i-th pattern of partition k, with d(k) doubly and sx(k) singly
!> occupied orbitals; nl(k) patterns.
subroutine part_lists(np,d,sx,lst,nl)
  use second_quantization, only: binom_coef, lex_init, lex_next
  integer(kind=iwp), intent(in) :: np, d(np), sx(np)
  integer(kind=iwp), allocatable, intent(out) :: lst(:,:,:), nl(:)
  integer(kind=iwp) :: k, nd, ns, id, is, pd, ps, i, mx
  call mma_allocate(nl,np,label='nl')
  mx = 1
  do k=1,np
    nl(k) = binom_coef(d(k),part_n(k))*binom_coef(sx(k),part_n(k)-d(k))
    mx = max(mx,nl(k))
  end do
  call mma_allocate(lst,[1,2],[1,mx],[1,np],label='lst')
  do k=1,np
    nd = binom_coef(d(k),part_n(k))
    ns = binom_coef(sx(k),part_n(k)-d(k))
    i = 0
    pd = lex_init(d(k),part_n(k))
    do id=1,nd
      ps = lex_init(sx(k),part_n(k)-d(k))
      do is=1,ns
        i = i+1
        lst(1,i,k) = ishft(pd,part_pos(k))
        lst(2,i,k) = ishft(pdep(ps,ibits(not(pd),0,max(part_n(k),1))),part_pos(k))
        if (is < ns) ps = lex_next(ps)
      end do
      if (id < nd) pd = lex_next(pd)
    end do
  end do
end subroutine part_lists

!> Sort the n configurations from conf_* (off+1) on by their key (heap sort).
subroutine sort_confs(off,n)
  integer(kind=iwp), intent(in) :: off, n
  integer(kind=iwp) :: i, iend
  do i=n/2,1,-1
    call sift(i,n)
  end do
  do iend=n,2,-1
    call swap(1,iend)
    call sift(1,iend-1)
  end do

contains

  subroutine sift(i0,m)
    integer(kind=iwp), intent(in) :: i0, m
    integer(kind=iwp) :: r, c
    r = i0
    do while (2*r <= m)
      c = 2*r
      if (c < m) then
        if (conf_key(off+c+1) > conf_key(off+c)) c = c+1
      end if
      if (conf_key(off+r) >= conf_key(off+c)) return
      call swap(r,c)
      r = c
    end do
  end subroutine sift

  subroutine swap(a,b)
    integer(kind=iwp), intent(in) :: a, b
    integer(kind=iwp) :: t
    t = conf_key(off+a)
    conf_key(off+a) = conf_key(off+b)
    conf_key(off+b) = t
    t = conf_doub(off+a)
    conf_doub(off+a) = conf_doub(off+b)
    conf_doub(off+b) = t
    t = conf_sing(off+a)
    conf_sing(off+a) = conf_sing(off+b)
    conf_sing(off+b) = t
  end subroutine swap

end subroutine sort_confs

!> The kept position (1..nkeep_group(ido)) of the configuration with key,
!> 0 if it is not kept (binary search).
function find_conf(ido,key) result(k)
  integer(kind=iwp), intent(in) :: ido, key
  integer(kind=iwp) :: k
  integer(kind=iwp) :: lo, hi, mid
  k = 0
  lo = 1
  hi = nkeep_group(ido)
  do while (lo <= hi)
    mid = (lo+hi)/2
    if (conf_key(conf_off(ido)+mid) == key) then
      k = mid
      return
    else if (conf_key(conf_off(ido)+mid) < key) then
      lo = mid+1
    else
      hi = mid-1
    end if
  end do
end function find_conf

!> The symmetry of a configuration (doubly occupied string doub, singly
!> occupied string sing within the other orbitals): the product of the
!> symmetries of the singly occupied orbitals.
function conf_sym(doub,sing,norb) result(isym)
  integer(kind=iwp) :: isym
  integer(kind=iwp), intent(in) :: doub, sing, norb
  integer(kind=iwp) :: iorb, socc
  socc = pdep(sing,ibits(not(doub),0,norb))
  isym = 0
  do iorb=1,norb
    if (btest(socc,iorb-1)) isym = ieor(isym,orbsym(iorb)-1)
  end do
  isym = isym+1
end function conf_sym

!> The number of CSFs of the CSF space (of symmetry ci_lsym).
function citrans_ncsf() result(n)
  integer(kind=iwp) :: n
  integer(kind=iwp) :: ido
  n = 0
  do ido=ndo_min,ndo_max
    n = n+nkeep_group(ido)*ncsf_group(ido)
  end do
end function citrans_ncsf

subroutine citrans_graph_free()
  call Free_SG_Graph(FileG,FileW)
end subroutine citrans_graph_free

subroutine citrans_sort(mode,ciold,cinew)
! sort CSFs such that they are grouped by configurations
! and fill in the transformation matrices to determinants

  use faroald, only: my_norb

  ! 'C' for configuration order, 'O' for original order
  character, intent(in) :: mode
  real(kind=wp), intent(in) :: ciold(:)
  real(kind=wp), intent(out) :: cinew(:)
  ! array with coupling coefficients
  !real(kind=wp), intent(out) :: coef(:)
  integer(kind=iwp) :: icsf, ido, idwn, ioff_csf, iphase, isyup, iup, mv, ncsf
  integer(kind=iwp), parameter :: maxorb = 32, maxdown = 16
  real(kind=wp) :: wtab(0:maxorb,maxdown)
  integer(kind=iwp), allocatable :: csf_offset(:), stepvector(:)

  ! Compute offsets for addressing into the reordering and coefficient arrays.
  call mma_allocate(csf_offset,[ndo_min,ndo_max],label='csf_offset')

# ifdef _DEBUGPRINT_
  write(u6,*) 'citrans_Sort'
  write(u6,'(5(1x,a4))') 'ido','ndoc','nsoc','ndet','ncsf'
# endif

  ncsf = 0
  do ido=ndo_min,ndo_max
#   ifdef _DEBUGPRINT_
    write(u6,'(5(1x,i4))') ido,ndoc_group(ido),nsoc_group(ido),ndet_group(ido),ncsf_group(ido)
#   endif
    ! use nconf/ncsf/ndet as offsets
    csf_offset(ido) = ncsf
    ncsf = ncsf+nkeep_group(ido)*ncsf_group(ido)
  end do

  ! The algorithm loops over CSFs (stepvectors) in an order determined
  ! by an external function, so that we do not need to know how that
  ! is done. Each time, we get a specific stepvector, e.g. 2u20du. We
  ! compute rankdo using the lexrank() function on the bitpattern 000101,
  ! and similarly rankso by applying the lexrank() function on the pattern
  ! 1101 (ignoring the doubly occupied orbitals), remembering that for
  ! bitpatterns we use the reverse orbital ordering. So our stepvector
  ! 2u20du (normal orbital ordering) has 000202 and xx0x as do and so
  ! strings, i.e. bitstrings 101 and 1101 respectively.
  !
  ! Once the configuration rank is found, we put in the coefficient of
  ! that CSF and when the loop is ended we have now a reordered set of
  ! CSF coefficients.

  ! set up table to compute CSF ranks
  call mkwtab(maxorb,maxdown,wtab)

  call mma_allocate(stepvector,my_norb,label='stepvector')

  ! initialize variables that track the stepvector
  mv = 1
  isyup = 1
  idwn = 1
  iup = 1
  do icsf=1,ncsf
    if (mv == 0) then
      write(u6,*) 'citrans_sort: fewer walks than CSFs:',icsf-1,ncsf
      call Abend()
    end if
    ! obtain the stepvector
    call SG_Walk_Next(FileG,FileW,ci_lsym,mv,isyup,idwn,iup,stepvector)

    ! position of the CSF in the configuration order, and its phase
    ioff_csf = csf_rank(stepvector,csf_offset,wtab,iphase)

    ! reorder the values
    if (mode == 'C') then
      cinew(ioff_csf+1) = iphase*ciold(icsf)
    else
      cinew(icsf) = iphase*ciold(ioff_csf+1)
    end if

  end do

  call mma_deallocate(csf_offset)
  call mma_deallocate(stepvector)

end subroutine citrans_sort

subroutine citrans_csf_rank(stepvector,k,iphase)
! Position k (1-based) of the CSF with the given GUGA step vector in the
! configuration order of citrans_sort, and its phase iphase:
! cinew(k) = iphase*ciold(icsf) for mode 'C'.

  use faroald, only: my_norb

  integer(kind=iwp), intent(in) :: stepvector(my_norb)
  integer(kind=iwp), intent(out) :: k, iphase
  integer(kind=iwp) :: ido, ncsf
  integer(kind=iwp), parameter :: maxorb = 32, maxdown = 16
  real(kind=wp) :: wtab(0:maxorb,maxdown)
  integer(kind=iwp), allocatable :: csf_offset(:)

  call mma_allocate(csf_offset,[ndo_min,ndo_max],label='csf_offset')
  ncsf = 0
  do ido=ndo_min,ndo_max
    csf_offset(ido) = ncsf
    ncsf = ncsf+nkeep_group(ido)*ncsf_group(ido)
  end do
  call mkwtab(maxorb,maxdown,wtab)
  k = csf_rank(stepvector,csf_offset,wtab,iphase)+1
  call mma_deallocate(csf_offset)

end subroutine citrans_csf_rank

function csf_rank(stepvector,csf_offset,wtab,iphase) result(ioff_csf)
! Offset (0-based) of the CSF with the given step vector in the
! configuration order, and its phase (see citrans_sort).

  use second_quantization, only: lexrank
  use faroald, only: my_norb

  integer(kind=iwp) :: ioff_csf
  integer(kind=iwp), intent(in) :: stepvector(my_norb), csf_offset(ndo_min:ndo_max)
  real(kind=wp), intent(in) :: wtab(0:,:)
  integer(kind=iwp), intent(out) :: iphase
  integer(kind=iwp) :: doub, ido, idown, iorb, iso, kconf, ndown, nsoc, rankdo, rankso, sing
  integer(kind=iwp) :: downvector(my_norb)

  ! determine configuration group and rank
  doub = 0
  sing = 0
  ido = 0
  iso = 0
  ndown = 0
  iphase = 1
  do iorb=1,my_norb
    select case (stepvector(iorb))
      case (1)
        sing = ibset(sing,iorb-ido-1)
        iso = iso+1
      case (2)
        sing = ibset(sing,iorb-ido-1)
        iso = iso+1
        ndown = ndown+1
        downvector(ndown) = iso
      case (3)
        doub = ibset(doub,iorb-1)
        ido = ido+1
        ! Each time we encounter a '2', we determine the phase change
        ! needed to put it in front of the string, i.e., (-1)^n, with n
        ! the number of singly occupied orbitals already encountered.
        if (mod(iso,2) == 1) iphase = -iphase
    end select
  end do

  ! configuration index
  rankdo = lexrank(doub)
  rankso = lexrank(sing)
  nsoc = nsoc_group(ido)
  kconf = find_conf(ido,nsoc*(rankdo-1)+rankso)
  if (kconf == 0) then
    write(u6,*) 'csf_rank: the CSF is not of the symmetry of the CSF space'
    call Abend()
  end if
  ioff_csf = csf_offset(ido)+(kconf-1)*ncsf_group(ido)

  ! CSF offset within this configuration, use per-ake's magical wtab
  ! to get the rank of a ud string within its set.
  do idown=1,ndown
    ioff_csf = ioff_csf+nint(wtab(downvector(idown)-2*idown,idown))
  end do

end function csf_rank

subroutine citrans_csf2sd(ci,det)

  use second_quantization, only: lex_init, lex_next, lexrank
  use faroald, only: my_nel, my_norb, ndeta, ndetb, nela

  real(kind=wp), intent(in) :: ci(:)
  real(kind=wp), intent(out) :: det(ndeta,ndetb)
  integer(kind=iwp) :: ido, iso, isoa, idet, kconf, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
                       doub, sing, alfa, beta, & ! occupation substrings
                       deta, detb, phase         ! determinant strings
  integer(kind=iwp), allocatable :: stepvector(:)
  real(kind=wp), allocatable :: tmp(:,:)

  ! Loop through the do,so configuration groups. For each group, load
  ! the spin table. Loop through the configurations and perform a matrix
  ! multiplication of the coefficient array with the spin table and get
  ! the determinant coefficients. Then, loop through each configuration
  ! and put the determinants in the correct place with the right phase.
  ! (Only the kept configurations of the CSF space, conf_doub/conf_sing;
  ! the determinants of the others are zero.)

  call mma_allocate(stepvector,my_norb,label='stepvector')
  det(:,:) = Zero

# ifdef _DEBUGPRINT_
  write(u6,*) 'citrans_csf2sd'
  write(u6,'(5(1x,a4))') 'ido','ndoc','nsoc','ndet','ncsf'
# endif

  ioff_csf1 = 1
  do ido=ndo_min,ndo_max
#   ifdef _DEBUGPRINT_
    write(u6,'(5(1x,i4))') ido,ndoc_group(ido),nsoc_group(ido),ndet_group(ido),ncsf_group(ido)
#   endif
    nconf = nkeep_group(ido)
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')
    tmp = Zero

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('N','N',ndet,nconf,ncsf,One,spintabs(ido)%coef,ndet,ci(ioff_csf1:ioff_csf2),ncsf,Zero,tmp,ndet)

    ! Store the determinant coefficients with the right phase factor in
    ! the correct place in the determinant matrix. The loops runs over
    ! the determinants in the order they are stored in the tmp array.

    iso = my_nel-2*ido
    isoa = nela-ido

    do kconf=1,nconf
      doub = conf_doub(conf_off(ido)+kconf)
      sing = conf_sing(conf_off(ido)+kconf)
      alfa = lex_init(isoa,iso)
      do idet=1,ndet
        beta = ibits(not(alfa),0,iso)
        phase = ds2ab(doub,sing,alfa,beta,deta,detb)
        det(lexrank(deta),lexrank(detb)) = phase*tmp(idet,kconf)
        alfa = lex_next(alfa)
      end do
    end do

    call mma_deallocate(tmp)
    ioff_csf1 = ioff_csf2+1
  end do

  call mma_deallocate(stepvector)

end subroutine citrans_csf2sd

subroutine citrans_sd2csf(det,ci)

  use second_quantization, only: lex_init, lex_next, lexrank
  use faroald, only: my_nel, my_norb, ndeta, ndetb, nela

  real(kind=wp), intent(in) :: det(ndeta,ndetb)
  real(kind=wp), intent(out) :: ci(:)
  integer(kind=iwp) :: ido, iso, isoa, idet, kconf, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
                       doub, sing, alfa, beta, & ! occupation substrings
                       deta, detb, phase         ! determinant strings
  integer(kind=iwp), allocatable :: stepvector(:)
  real(kind=wp), allocatable :: tmp(:,:)

  ! Loop through the do,so configuration groups. For each group, load
  ! the spin table. Loop through the configurations and perform a matrix
  ! multiplication of the coefficient array with the spin table and get
  ! the determinant coefficients. Then, loop through each configuration
  ! and put the determinants in the correct place with the right phase.

  call mma_allocate(stepvector,my_norb,label='stepvector')

  ioff_csf1 = 1
  do ido=ndo_min,ndo_max
    nconf = nkeep_group(ido)
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')

    ! Collect the determinant coefficients with the right phase factor
    ! from the correct place in the determinant matrix. The loops runs
    ! over the determinants in the order they should end up in tmp.

    iso = my_nel-2*ido
    isoa = nela-ido

    do kconf=1,nconf
      doub = conf_doub(conf_off(ido)+kconf)
      sing = conf_sing(conf_off(ido)+kconf)
      alfa = lex_init(isoa,iso)
      do idet=1,ndet
        beta = ibits(not(alfa),0,iso)
        phase = ds2ab(doub,sing,alfa,beta,deta,detb)
        tmp(idet,kconf) = phase*det(lexrank(deta),lexrank(detb))
        alfa = lex_next(alfa)
      end do
    end do

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('T','N',ncsf,nconf,ndet,One,spintabs(ido)%coef,ndet,tmp,ndet,Zero,ci(ioff_csf1:ioff_csf2),ncsf)

    call mma_deallocate(tmp)
    ioff_csf1 = ioff_csf2+1
  end do

  call mma_deallocate(stepvector)

end subroutine citrans_sd2csf

!> As citrans_csf2sd and citrans_sd2csf, with the CI vector of the RAS
!> string classes (faroald_ras) instead of the full determinant matrix.
subroutine citrans_csf2ras(ci,vec)

  use second_quantization, only: lex_init, lex_next
  use faroald, only: my_nel, my_norb, nela
  use faroald_ras, only: n_ras, ras_addr

  real(kind=wp), intent(in) :: ci(:)
  real(kind=wp), intent(out) :: vec(n_ras)
  integer(kind=iwp) :: ido, iso, isoa, idet, kconf, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
                       doub, sing, alfa, beta, & ! occupation substrings
                       deta, detb, phase         ! determinant strings
  integer(kind=iwp), allocatable :: stepvector(:)
  real(kind=wp), allocatable :: tmp(:,:)

  ! Loop through the do,so configuration groups. For each group, load
  ! the spin table. Loop through the configurations and perform a matrix
  ! multiplication of the coefficient array with the spin table and get
  ! the determinant coefficients. Then, loop through each configuration
  ! and put the determinants in the correct place with the right phase.
  ! (Only the kept configurations of the CSF space, conf_doub/conf_sing;
  ! the determinants of the others are zero.)

  call mma_allocate(stepvector,my_norb,label='stepvector')
  vec(:) = Zero

# ifdef _DEBUGPRINT_
  write(u6,*) 'citrans_csf2ras'
  write(u6,'(5(1x,a4))') 'ido','ndoc','nsoc','ndet','ncsf'
# endif

  ioff_csf1 = 1
  do ido=ndo_min,ndo_max
#   ifdef _DEBUGPRINT_
    write(u6,'(5(1x,i4))') ido,ndoc_group(ido),nsoc_group(ido),ndet_group(ido),ncsf_group(ido)
#   endif
    nconf = nkeep_group(ido)
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')
    tmp = Zero

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('N','N',ndet,nconf,ncsf,One,spintabs(ido)%coef,ndet,ci(ioff_csf1:ioff_csf2),ncsf,Zero,tmp,ndet)

    ! Store the determinant coefficients with the right phase factor in
    ! the correct place in the determinant matrix. The loops runs over
    ! the determinants in the order they are stored in the tmp array.

    iso = my_nel-2*ido
    isoa = nela-ido

    do kconf=1,nconf
      doub = conf_doub(conf_off(ido)+kconf)
      sing = conf_sing(conf_off(ido)+kconf)
      alfa = lex_init(isoa,iso)
      do idet=1,ndet
        beta = ibits(not(alfa),0,iso)
        phase = ds2ab(doub,sing,alfa,beta,deta,detb)
        vec(ras_addr(deta,detb)) = phase*tmp(idet,kconf)
        alfa = lex_next(alfa)
      end do
    end do

    call mma_deallocate(tmp)
    ioff_csf1 = ioff_csf2+1
  end do

  call mma_deallocate(stepvector)

end subroutine citrans_csf2ras

subroutine citrans_ras2csf(vec,ci)

  use second_quantization, only: lex_init, lex_next
  use faroald, only: my_nel, my_norb, nela
  use faroald_ras, only: n_ras, ras_addr

  real(kind=wp), intent(in) :: vec(n_ras)
  real(kind=wp), intent(out) :: ci(:)
  integer(kind=iwp) :: ido, iso, isoa, idet, kconf, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
                       doub, sing, alfa, beta, & ! occupation substrings
                       deta, detb, phase         ! determinant strings
  integer(kind=iwp), allocatable :: stepvector(:)
  real(kind=wp), allocatable :: tmp(:,:)

  ! Loop through the do,so configuration groups. For each group, load
  ! the spin table. Loop through the configurations and perform a matrix
  ! multiplication of the coefficient array with the spin table and get
  ! the determinant coefficients. Then, loop through each configuration
  ! and put the determinants in the correct place with the right phase.

  call mma_allocate(stepvector,my_norb,label='stepvector')

  ioff_csf1 = 1
  do ido=ndo_min,ndo_max
    nconf = nkeep_group(ido)
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')

    ! Collect the determinant coefficients with the right phase factor
    ! from the correct place in the determinant matrix. The loops runs
    ! over the determinants in the order they should end up in tmp.

    iso = my_nel-2*ido
    isoa = nela-ido

    do kconf=1,nconf
      doub = conf_doub(conf_off(ido)+kconf)
      sing = conf_sing(conf_off(ido)+kconf)
      alfa = lex_init(isoa,iso)
      do idet=1,ndet
        beta = ibits(not(alfa),0,iso)
        phase = ds2ab(doub,sing,alfa,beta,deta,detb)
        tmp(idet,kconf) = phase*vec(ras_addr(deta,detb))
        alfa = lex_next(alfa)
      end do
    end do

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('T','N',ncsf,nconf,ndet,One,spintabs(ido)%coef,ndet,tmp,ndet,Zero,ci(ioff_csf1:ioff_csf2),ncsf)

    call mma_deallocate(tmp)
    ioff_csf1 = ioff_csf2+1
  end do

  call mma_deallocate(stepvector)

end subroutine citrans_ras2csf

subroutine citrans_conf_dets(kcsf,k1,ncsf,ndet,deta,detb,coef)
! Determinant expansion of the configuration holding the CSF kcsf, kcsf
! being an index of the configuration order (citrans_sort mode 'C'):
! CSF k1-1+j (j = 1..ncsf) of this configuration is
! sum_i coef(i,j) |deta(i) detb(i)>, with deta(i)/detb(i) the alpha/beta
! strings of the determinants (as in citrans_csf2sd, phase included).
! deta, detb and coef are allocated here, to be freed by the caller.

  use second_quantization, only: lex_init, lex_next
  use faroald, only: my_nel, nela

  integer(kind=iwp), intent(in) :: kcsf
  integer(kind=iwp), intent(out) :: k1, ncsf, ndet
  integer(kind=iwp), allocatable, intent(out) :: deta(:), detb(:)
  real(kind=wp), allocatable, intent(out) :: coef(:,:)
  integer(kind=iwp) :: ido, iso, isoa, ioff, idet, doub, sing, alfa, beta, phase, kconf

  ! locate the configuration group and the configuration within it (among
  ! those of the symmetry of the CSF space)
  ioff = 0
  do ido=ndo_min,ndo_max
    if (kcsf <= ioff+nkeep_group(ido)*ncsf_group(ido)) exit
    ioff = ioff+nkeep_group(ido)*ncsf_group(ido)
  end do
  if ((kcsf < 1) .or. (ido > ndo_max)) then
    write(u6,*) 'citrans_conf_dets: CSF index out of range:',kcsf
    call Abend()
  end if
  ncsf = ncsf_group(ido)
  ndet = ndet_group(ido)
  kconf = (kcsf-ioff-1)/ncsf+1
  k1 = ioff+(kconf-1)*ncsf+1
  doub = conf_doub(conf_off(ido)+kconf)
  sing = conf_sing(conf_off(ido)+kconf)

  iso = my_nel-2*ido
  isoa = nela-ido

  call mma_allocate(deta,ndet,label='deta')
  call mma_allocate(detb,ndet,label='detb')
  call mma_allocate(coef,ndet,ncsf,label='coef')
  alfa = lex_init(isoa,iso)
  do idet=1,ndet
    beta = ibits(not(alfa),0,iso)
    phase = ds2ab(doub,sing,alfa,beta,deta(idet),detb(idet))
    coef(idet,:) = phase*spintabs(ido)%coef(idet,:)
    alfa = lex_next(alfa)
  end do

end subroutine citrans_conf_dets

subroutine spintabs_allocate()

# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A) _BIND_C_
    import :: wp
    real(kind=wp), allocatable :: A(:,:)
    end subroutine c_null_alloc
    end interface
    integer(kind=iwp) :: i
# endif

    call mma_allocate(spintabs,[ndo_min,ndo_max],label='spintabs')

# ifdef _GARBLE_
    ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
    do i=ndo_min,ndo_max
      call c_null_alloc(spintabs(i)%coef)
    end do
# endif

# include "macros.fh"
    unused_proc(mma_allocate(spintabs,0))

    end subroutine spintabs_allocate

    subroutine spintabs_free()

      integer(kind=iwp) :: i

      If (allocated(Spintabs)) Then
         do i=lbound(spintabs,1),ubound(spintabs,1)
           call mma_deallocate(spintabs(i)%coef,safe='*')
         end do
         call mma_deallocate(spintabs,safe='*')
      End If

    end subroutine spintabs_free

    subroutine spintable_create(nso,ndown,spintab)
! given a number of singly occupied orbitals and down couplings,
! construct the table of CSFs and determinant transformation matrices.
!
! example: given 3 nso and 1 ndown, the matrix constructed is:
!     aab aba baa
! udu ..
! uud
!
! for now, just print the matrix for testing, later put it in some
! memory location to be used by the conversion routine

      integer(kind=iwp), intent(in) :: nso, ndown
      type(spintable), intent(inout) :: spintab
      integer(kind=iwp) :: idown, icsf, ncsf, ndet
      integer(kind=iwp), allocatable :: down_orb(:), udvec(:)

      ndet = spintab%ndet
      ncsf = spintab%ncsf

      call mma_allocate(down_orb,ndown+1,label='down_orb')
      call mma_allocate(udvec,nso,label='udvec')

      call mma_allocate(spintab%coef,ndet,ncsf,label='spintab%coef')

      ! The CSFs here are generated in order of ascending rank that matches
      ! wtab, the ranking table used to sort the CSFs. The most alternating
      ! udud...ud string comes first, and the u..ud...d string last.

      call csf_init(nso,ndown,down_orb)
      do icsf=1,ncsf
        ! create udvec of singly occupied orbitals
        udvec = 1
        do idown=1,ndown
          udvec(down_orb(idown)) = 2
        end do
        ! expand into determinants, store coefficients
        call ud2det(udvec,spintab%coef(:,icsf))
        call csf_next(ndown,down_orb)
      end do

      call mma_deallocate(down_orb)
      call mma_deallocate(udvec)

    end subroutine spintable_create

    subroutine ud2det(udvec,coef)
! A stepvector in this case is represented by a series of integers 1,2
! for each singly-occupied orbital and corresponds to the steps u, d.
! As such it is actually a limited stepvector, since doubly occupied
! or empty orbitals are omitted. The reason is that we only generate
! the stepvectors for a general configuration case, and then modify
! the phase factors later when doing the actual transformation.
!
! The determinants are traversed in lexicographic order, with all alpha
! orbitals the lowest orbitals.

      use second_quantization, only: binom_coef, lex_init, lex_next

      integer(kind=iwp), intent(in) :: udvec(:)
      real(kind=wp), intent(out) :: coef(:)
      ! the coupling coefficient
      integer(kind=iwp) :: phase, & ! phase factor
                           iso, nso, nsoa, nsob, ialfa, ibeta, ia, ib, idet, ndet, deta, ilev, nlev
      real(kind=wp) :: nom, den ! fraction holding the coefficient
      ! lex keeps track of alpha orbitals
      integer(kind=iwp), parameter :: maxorb = 64

      ! find number of singly occupied orbitals
      nsoa = 0
      nsob = 0
      nlev = size(udvec)
      do ilev=1,nlev
        select case (udvec(ilev))
          case (1)
            nsoa = nsoa+1
          case (2)
            nsob = nsob+1
        end select
      end do
      nso = nsoa+nsob

      ndet = binom_coef(nsoa,nso)
      ! loop over possible determinants
      deta = lex_init(nsoa,nso)
      do idet=1,ndet
        nom = One
        den = One
        phase = 1
        ialfa = 0
        ibeta = 0
        iso = 0
        ia = 0
        ib = 0
        do ilev=1,nlev
          select case (udvec(ilev))
            case (1)
              ib = ib+1
              if (btest(deta,iso)) then
                nom = nom*(ia+ib-ibeta)
                ialfa = ialfa+1
              else
                nom = nom*(ia+ib-ialfa)
                ibeta = ibeta+1
              end if
              iso = iso+1
              den = den*ib
            case (2)
              ia = ia+1
              ib = ib-1
              if (btest(deta,iso)) then
                nom = nom*(ibeta-ia+1)
                if (mod(ib,2) == 0) phase = -phase
                ialfa = ialfa+1
              else
                nom = nom*(ialfa-ia+1)
                if (mod(ib,2) /= 0) phase = -phase
                ibeta = ibeta+1
              end if
              iso = iso+1
              den = den*(ib+2)
            case default
              write(u6,'(1x,a)') 'ud2det: udvec element /= 1 or 2, fatal...'
              call AbEnd()
          end select
        end do

        ! Compute the determinant coefficient and position and add the value
        ! into the determinant array for this configuration. Not that the
        ! determinant coefficients are stored in lexicographic bit-order.
        coef(idet) = phase*sqrt(nom/den)

        deta = lex_next(deta)
      end do

    end subroutine ud2det

    function ds2ab(doub,sing,alfa,beta,deta,detb) result(phase)
! convert a determinant characterized by a doubly occupied, singly
! occupied, and alpha/beta substrings to an alpha and beta string.
!
! A doubly occupied string identifies the doubly occupied orbitals
! with a one-bit. A singly occupied string identifies the singly
! occupied orbitals with a one-bit within the non-doubly occupied
! orbitals. Finally, an alpha substring identifies the alpha orbitals
! within the singly occupied orbital space.
!
! example: 6in6, singlet, determinant 2a20ba would be represented as:
!     <-    <-   <-                  <-      <-
! 000101, 1101, 101 (dsa) and as 100111, 010101 (ab)
! as usual counting orbitals from the right with bits (lower-most bit)

      use faroald, only: my_norb

      integer(kind=iwp) :: phase
      integer(kind=iwp), intent(in) :: doub, sing, alfa, beta
      integer(kind=iwp), intent(out) :: deta, detb
      integer(kind=iwp) :: mask, not_doub, pos
      logical(kind=iwp) :: switch

      ! First, we have to get the a/b singly occupied orbitals. This can be
      ! easily done using successive bit scattering operations on the alpha
      ! and beta substrings, which apparently are present in e.g. Intel's
      ! Haswell architecture, called pext/pdep. Here, I implemented the bit
      ! scattering operation as a function pdep: res = pdep(val,mask), which
      ! scatters the bits of 'val' according to 'mask' into 'res'.

      not_doub = ibits(not(doub),0,my_norb)
      deta = pdep(pdep(alfa,sing),not_doub)
      detb = pdep(pdep(beta,sing),not_doub)

      ! Add together with the doubly occupied orbitals.
      deta = ior(deta,doub)
      detb = ior(detb,doub)

      ! Finally, determine the phase of the determinant. The coupling
      ! coefficients were constructed for an orbital ordering where alpha
      ! and beta alternate, e.g. 2a20ab -> 11'233'56'. However, we order alpha
      ! first then beta, so we have 1235|1'3'6'. In order to account for the
      ! change in phase, we should build our determinant by creating the
      ! electrons of the alternating-order determinant in reverse, i.e., 6'
      ! first, then 5, then 3', and so on. If we do this for the split-order
      ! determinant, each time we add a beta electron we should count the
      ! alpha electrons already added and change the phase if that number is
      ! odd.

      mask = 0
      switch = .false.
      pos = 0
      do while (ishft(deta,-pos) /= 0)
        if (switch) mask = ibset(mask,pos)
        if (iand(ishft(detb,-pos),1) == 1) switch = .not. switch
        pos = pos+1
      end do
      phase = 1-2*poppar(iand(deta,mask))

    end function ds2ab

    function pdep(val,mask) result(res)

      integer(kind=iwp) :: res
      integer(kind=iwp), intent(in) :: val, mask
      integer(kind=iwp) :: mask_bit, pos, tmp_mask, tmp_val, val_bit

      tmp_val = val
      tmp_mask = mask

      ! nothing set by default
      res = 0

      pos = 0
      do while (tmp_mask /= 0)
        mask_bit = iand(tmp_mask,1)
        if (mask_bit == 1) then
          val_bit = iand(tmp_val,1)
          res = ior(res,ishft(val_bit,pos))
          tmp_val = ishft(tmp_val,-1)
        end if
        tmp_mask = ishft(tmp_mask,-1)
        pos = pos+1
      end do

    end function pdep

    ! the inverse of pdep: gathers the bits of 'val' at the positions set in
    ! 'mask' into the lowest bits of 'res'
    function pext(val,mask) result(res)

      integer(kind=iwp) :: res
      integer(kind=iwp), intent(in) :: val, mask
      integer(kind=iwp) :: pos, k, tmp_mask

      res = 0
      tmp_mask = mask
      pos = 0
      k = 0
      do while (tmp_mask /= 0)
        if (iand(tmp_mask,1) == 1) then
          if (btest(val,pos)) res = ibset(res,k)
          k = k+1
        end if
        tmp_mask = ishft(tmp_mask,-1)
        pos = pos+1
      end do

    end function pext

! a CSF is a CSF consisting of singly occupied orbitals and is given
! by its down_orb string, that is the orbitals which are down coupled.

    subroutine csf_init(nso,ndown,down_orb)

      integer(kind=iwp), intent(in) :: nso, ndown
      integer(kind=iwp), intent(out) :: down_orb(ndown+1)
      integer(kind=iwp) :: i

      do i=1,ndown
        down_orb(i) = 2*i
      end do
      down_orb(ndown+1) = nso+1

    end subroutine csf_init

    subroutine csf_next(ndown,down_orb)

      integer(kind=iwp), intent(in) :: ndown
      integer(kind=iwp), intent(inout) :: down_orb(ndown+1)
      integer(kind=iwp) :: i, j

      do i=1,ndown
        if (down_orb(i) < down_orb(i+1)-1) then
          down_orb(i) = down_orb(i)+1
          do j=1,i-1
            down_orb(j) = 2*j
          end do
          return
        end if
      end do

    end subroutine csf_next

    subroutine mkwtab(mxn1,mxn2,wtab)

      use second_quantization, only: binom_coef

      integer(kind=iwp), intent(in) :: mxn1, mxn2
      real(kind=wp), intent(out) :: wtab(0:mxn1,mxn2)
      integer(kind=iwp) :: n1, n2

      do n1=0,mxn1
        do n2=1,mxn2
          wtab(n1,n2) = real(binom_coef(n1+n2,n1+2*n2),kind=wp)*real(n1,kind=wp)/real(n1+2*n2,kind=wp)
        end do
      end do

    end subroutine mkwtab

! Extensions to mma_interfaces, using preprocessor templates
! (see mma_util/stdalloc.F90)

! Define spt_mma_allo_1D, spt_mma_allo_1D_lim, spt_mma_free_1D
#define _TYPE_ type(spintable)
#  define _SUBR_NAME_ spt_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'spt_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

  end module citrans
