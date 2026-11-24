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

module faroald_ras
! String-driven RAS/GAS for Faroald (Olsen et al., J. Chem. Phys. 89, 2185
! (1988)): the alpha and beta strings are grouped in occupation classes
! (n_1,..,n_np), the numbers of electrons in the partitions (RAS1/2/3, or
! the GAS spaces); the CI vector holds only the blocks (alpha class, beta
! class) of the space, where the partitions 1..k hold ras_mn(k) to
! ras_mx(k) electrons. No string or determinant outside the space is ever
! addressed.
!
! The orbitals are those of Faroald, the levels of the split graph,
! partition by partition (contiguous bit ranges ras_pos(k)+1..
! ras_pos(k)+ras_o(k)). Within a class a string is addressed by the
! lexical ranks of its substrings, the first partition fastest. The
! element of determinant
! |alpha beta> is that of Faroald's determinant matrix psi(ia,ib), same
! phase convention; block (ca,cb) is stored column-major, na(ca) x nb(cb).
!
! H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx (E_tu E_vx - delta_uv E_tx),
! g_tuvx = (tu|vx), is applied as
!   same spin:     sum_tu h_tu a+_t a_u
!                  + sum_{t>v,u>x} (g_tuvx-g_txvu) a+_t a+_v a_x a_u
!   opposite spin: sum_tuvx g_tuvx a+_ta a_ua a+_vb a_xb
! so that no intermediate string outside the RAS space is needed.
!
! A CAS space is the special case of one partition (one class per spin).
!
! Singlets (spin_flip, nela = nelb, M_S = S = 0): the vectors are
! symmetric under the exchange of the alpha and beta strings,
! C(Ia,Ib) = C(Ib,Ia) (Faroald's phase convention, as in sigma_update).
! Then the alpha and beta strings are the same, block (ga,gb) is the
! transpose of block (gb,ga), and sigma = X + X^T with X = the beta-beta
! part plus half of the alpha-beta part (pairs tu >= vx of
! g_tuvx a+_ta a_ua a+_vb a_xb, tu = vx with weight 1/2); the alpha-alpha
! density equals the beta-beta one and the spin density is zero.

use second_quantization, only: ann, ann2, binom_coef, cre, cre2, ex1, fase, lex_init, lex_next, lexrank
use stdalloc, only: mma_allocate, mma_deallocate
use Molcas, only: MxGAS
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
private

type ras_spin
  integer(kind=iwp) :: nel = 0, nocc = 0, ncls = 0, nstr = 0, nex = 0
  ! occupation classes: electrons in the partitions, numbers of substrings,
  ! number of strings and offset in the rank order (0-based); occ_skey
  ! the packed occupations (occ_key) sorted, occ_sidx their classes
  integer(kind=iwp), allocatable :: occ_n(:,:), occ_cnt(:,:), occ_nstr(:), occ_off(:), occ_skey(:), occ_sidx(:)
  ! the groups (the blocks of the CI layout): occupation class and string
  ! symmetry; per group: electrons in the partitions (cls_n), symmetry (cls_sym),
  ! number of strings, offset of its strings in the global string list
  ! (0-based). The strings are numbered group by group.
  integer(kind=iwp), allocatable :: cls_n(:,:), cls_sym(:), cls_nstr(:), cls_off(:)
  ! perm(k): the global string of rank address k (occ_off+ranks)
  integer(kind=iwp), allocatable :: perm(:)
  ! per global string: bit string, group, index within the group
  integer(kind=iwp), allocatable :: str(:), str_cls(:), str_idx(:)
  ! single excitations J = a+_p a_q I (sign sgn), I = 1..nstr, entries
  ! ex_off(I)+1..ex_off(I+1); only J of this spin's classes
  integer(kind=iwp), allocatable :: ex_off(:), ex_p(:), ex_q(:), ex_sgn(:), ex_j(:)
  ! same-spin double excitations J = a+_t a+_v a_x a_u I (t > v, u > x),
  ! entries dex_off(I)+1..dex_off(I+1): J, sign, (t,u,v,x) packed (pack4);
  ! has_dex = .false. if they did not fit in memory (generated on the fly)
  logical(kind=iwp) :: has_dex = .false.
  integer(kind=iwp) :: ndex = 0
  integer(kind=iwp), allocatable :: dex_off(:), dex_j(:), dex_sgn(:), dex_tuvx(:)
end type ras_spin

! partitions: number, orbitals, first orbital - 1, cumulative bounds on
! the electrons in the partitions 1..k; ras_rdx the radices of occ_key
integer(kind=iwp) :: ras_norb = 0, ras_np = 0, ras_o(MxGAS) = 0, ras_pos(MxGAS) = 0, ras_mn(MxGAS) = 0, ras_mx(MxGAS) = 0, &
                     ras_rdx(MxGAS) = 0, n_ras = 0, ras_nsym = 1, ras_lsym = 1, ras_osym(30) = 1
type(ras_spin) :: sa, sb
! block offsets (0-based) of (alpha group, beta group), -1 outside the
! RAS space or of another symmetry; blk_offT: the same blocks transposed
! (beta x alpha)
integer(kind=iwp), allocatable :: blk_off(:,:), blk_offT(:,:)
logical(kind=iwp) :: ras_active = .false., ras_singlet = .false.
! the beta single excitations Jb -> Ib = a+_v a_x Jb in batches, one per
! beta group of Jb and (v,x): batch b holds the entries
! bat_off(b)+1..bat_off(b+1) (Jb and Ib as indices within their groups,
! sign); bat_cj/bat_ci the groups of Jb/Ib; mxbat the largest batch
integer(kind=iwp) :: nbat = 0, mxbat = 0
! (N-1)-string resolution of the alpha-beta sigma (Olsen; Lucia's rsbb2bn):
! <Ia|a+_t a_u|Ja> = sum_K <Ia|a+_t|K> <K|a_u|Ja> over the alpha strings K
! with one electron less that occur (nres_k), grouped in blocks of equal
! occupation of the partitions and symmetry (res_nblk; block B: K
! res_boff(B)+1..res_boff(B+1), symmetry res_bsym(B), occupations res_bocc(:,B)).
! res_cre(t,K) = +-(global alpha string a+_t K), 0 if t is occupied in K or
! a+_t K is not an alpha string. The orbitals sorted by symmetry: res_ord
! (sorted -> orbital), res_soff/res_ns (offset and number per symmetry).
! Coarse blocks (res_ncb): the fine blocks of one symmetry together.
! res_okb(B,cb) / res_cokb: the K of (fine/coarse) block B couple with the beta
! strings of group cb. res_mode: 0 the direct loops (opposite_spin), 1 the
! resolution with the fine blocks, 2 with the coarse blocks, -1 not yet
! chosen: the candidates res_try are timed in the first sigma calls and the
! fastest is used from then on (see res_choose).
integer(kind=iwp) :: nres_k = 0, res_nblk = 0, res_ncb = 0, res_mode = 0, res_ntry = 0, res_itry = 0, res_try(3) = 0
real(kind=wp) :: res_tim(0:2) = Zero
integer(kind=iwp), allocatable :: res_cre(:,:), res_boff(:), res_bsym(:), res_bocc(:,:), res_ord(:), res_soff(:), res_ns(:), &
                                  res_cboff(:), res_cbsym(:)
logical(kind=iwp), allocatable :: res_okb(:,:), res_cokb(:,:)
! the work arrays of the resolution, kept between the calls (no new pages
! per sigma); used with the current block size as leading dimension
real(kind=wp), allocatable :: res_Xw(:), res_Sw(:)
! the cost model of res_choose, in units of one element of a vector update
! of the direct loops: per update (Cost_Op), per gathered or scattered
! element (Cost_GS), per DGEMM call (Cost_Call) and multiply-add (Cost_Flop);
! fitted to the timings of CAS/RAS/GAS cases with and without symmetry
! (gfortran, MKL, one thread; right within a factor 3)
real(kind=wp), parameter :: Cost_Op = 56.0_wp, Cost_GS = 6.0_wp, Cost_Call = 78.0_wp, Cost_Flop = 0.35_wp
integer(kind=iwp), allocatable :: bat_cj(:), bat_ci(:), bat_v(:), bat_x(:), bat_off(:), be_j(:), be_i(:), be_sgn(:)

public :: faroald_ras_free, faroald_ras_init, n_ras, ras_active, ras_addr, ras_count_dets, ras_from_full, ras_sigma, &
          ras_to_full, ras_traci, ras_trans_one_pdm, ras_trans_two_pdm

contains

!> Set up the string groups, the single-excitation lists and the CI layout
!> for nela/nelb electrons in norb orbitals, o(k) of them in partition k
!> (np partitions, RAS: 3), the partitions 1..k holding mn(k) to mx(k)
!> electrons; the orbitals have the symmetries osym (nsym irreps), the CI
!> space the symmetry lsym (determinants sym(alpha) x sym(beta) = lsym
!> only). spin_flip: all vectors are singlets (nela = nelb), symmetric
!> under the exchange of the alpha and beta strings.
!> ab_timed: the path of the alpha-beta sigma is the fastest of the first
!> sigma calls, else (reproducible, for the verification) by the cost
!> estimate (see res_choose).
subroutine faroald_ras_init(nela,nelb,norb,np,o,mn,mx,nsym,osym,lsym,spin_flip,ab_timed)
  integer(kind=iwp), intent(in) :: nela, nelb, norb, np, o(np), mn(np), mx(np), nsym, osym(norb), lsym
  logical(kind=iwp), intent(in) :: spin_flip, ab_timed
  integer(kind=iwp) :: ca, cb, n, k

  call faroald_ras_free()
  if (sum(o) /= norb) then
    write(u6,*) 'faroald_ras_init: the orbitals of the partitions /= norb',o,norb
    call Abend()
  end if
  if (norb > 30) then
    write(u6,*) 'faroald_ras_init: at most 30 orbitals',norb
    call Abend()
  end if
  if ((np < 1) .or. (np > MxGAS)) then
    write(u6,*) 'faroald_ras_init: illegal number of partitions',np
    call Abend()
  end if
  ras_norb = norb
  ras_np = np
  ras_o(:) = 0
  ras_o(1:np) = o(:)
  ras_mn(1:np) = mn(:)
  ras_mx(1:np) = mx(:)
  ras_pos(1) = 0
  ras_rdx(1) = 1
  do k=2,np
    ras_pos(k) = ras_pos(k-1)+o(k-1)
    ras_rdx(k) = ras_rdx(k-1)*(o(k-1)+1)
  end do
  ras_nsym = nsym
  ras_lsym = lsym
  ras_osym(:) = 1
  ras_osym(1:norb) = osym(1:norb)
  if (spin_flip .and. (nela /= nelb)) then
    write(u6,*) 'faroald_ras_init: spin-flip symmetry needs nela = nelb',nela,nelb
    call Abend()
  end if
  ras_singlet = spin_flip

  call mk_classes(sa,nela,nelb)
  call mk_classes(sb,nelb,nela)
  call mk_strings(sa)
  call mk_strings(sb)
  call mk_singles(sa)
  call mk_singles(sb)
  call mk_doubles(sa)
  call mk_doubles(sb)
  call mk_batches()
  call mk_resolution()

  call mma_allocate(blk_off,sa%ncls,sb%ncls,label='blk_off')
  call mma_allocate(blk_offT,sb%ncls,sa%ncls,label='blk_offT')
  blk_off(:,:) = -1
  blk_offT(:,:) = -1
  n = 0
  do cb=1,sb%ncls
    do ca=1,sa%ncls
      if (allowed(sa%cls_n(:,ca),sb%cls_n(:,cb)) .and. (ieor(sa%cls_sym(ca)-1,sb%cls_sym(cb)-1)+1 == ras_lsym)) then
        blk_off(ca,cb) = n
        n = n+sa%cls_nstr(ca)*sb%cls_nstr(cb)
      end if
    end do
  end do
  n_ras = n
  n = 0
  do ca=1,sa%ncls
    do cb=1,sb%ncls
      if (blk_off(ca,cb) >= 0) then
        blk_offT(cb,ca) = n
        n = n+sa%cls_nstr(ca)*sb%cls_nstr(cb)
      end if
    end do
  end do
  call res_choose(ab_timed)
  ras_active = .true.

end subroutine faroald_ras_init

subroutine faroald_ras_free()
  call free_spin(sa)
  call free_spin(sb)
  call mma_deallocate(bat_cj,safe='*')
  call mma_deallocate(bat_ci,safe='*')
  call mma_deallocate(bat_v,safe='*')
  call mma_deallocate(bat_x,safe='*')
  call mma_deallocate(bat_off,safe='*')
  call mma_deallocate(be_j,safe='*')
  call mma_deallocate(be_i,safe='*')
  call mma_deallocate(be_sgn,safe='*')
  nbat = 0
  mxbat = 0
  call mma_deallocate(res_cre,safe='*')
  call mma_deallocate(res_boff,safe='*')
  call mma_deallocate(res_bsym,safe='*')
  call mma_deallocate(res_bocc,safe='*')
  call mma_deallocate(res_cboff,safe='*')
  call mma_deallocate(res_cbsym,safe='*')
  call mma_deallocate(res_okb,safe='*')
  call mma_deallocate(res_cokb,safe='*')
  call mma_deallocate(res_Xw,safe='*')
  call mma_deallocate(res_Sw,safe='*')
  res_nblk = 0
  res_ncb = 0
  res_mode = 0
  res_ntry = 0
  res_itry = 0
  call mma_deallocate(res_ord,safe='*')
  call mma_deallocate(res_soff,safe='*')
  call mma_deallocate(res_ns,safe='*')
  nres_k = 0
  call mma_deallocate(blk_off,safe='*')
  call mma_deallocate(blk_offT,safe='*')
  n_ras = 0
  ras_active = .false.
  ras_singlet = .false.
end subroutine faroald_ras_free

subroutine free_spin(s)
  type(ras_spin), intent(inout) :: s
  call mma_deallocate(s%occ_n,safe='*')
  call mma_deallocate(s%occ_cnt,safe='*')
  call mma_deallocate(s%occ_nstr,safe='*')
  call mma_deallocate(s%occ_off,safe='*')
  call mma_deallocate(s%occ_skey,safe='*')
  call mma_deallocate(s%occ_sidx,safe='*')
  call mma_deallocate(s%cls_n,safe='*')
  call mma_deallocate(s%cls_sym,safe='*')
  call mma_deallocate(s%cls_nstr,safe='*')
  call mma_deallocate(s%cls_off,safe='*')
  call mma_deallocate(s%perm,safe='*')
  call mma_deallocate(s%str,safe='*')
  call mma_deallocate(s%str_cls,safe='*')
  call mma_deallocate(s%str_idx,safe='*')
  call mma_deallocate(s%ex_off,safe='*')
  call mma_deallocate(s%ex_p,safe='*')
  call mma_deallocate(s%ex_q,safe='*')
  call mma_deallocate(s%ex_sgn,safe='*')
  call mma_deallocate(s%ex_j,safe='*')
  call mma_deallocate(s%dex_off,safe='*')
  call mma_deallocate(s%dex_j,safe='*')
  call mma_deallocate(s%dex_sgn,safe='*')
  call mma_deallocate(s%dex_tuvx,safe='*')
  s%has_dex = .false.
  s%ndex = 0
  s%nel = 0
  s%nocc = 0
  s%ncls = 0
  s%nstr = 0
  s%nex = 0
end subroutine free_spin

!> A determinant of the classes with electrons n_a, n_b in the partitions
!> is in the space: the partitions 1..k hold ras_mn(k) to ras_mx(k)
!> electrons.
pure function allowed(n_a,n_b) result(ok)
  integer(kind=iwp), intent(in) :: n_a(:), n_b(:)
  logical(kind=iwp) :: ok
  integer(kind=iwp) :: k, ne
  ok = .false.
  ne = 0
  do k=1,ras_np
    ne = ne+n_a(k)+n_b(k)
    if ((ne < ras_mn(k)) .or. (ne > ras_mx(k))) return
  end do
  ok = .true.
end function allowed

!> The packed occupation of a class (mixed radix ras_rdx).
pure function occ_key(n) result(key)
  integer(kind=iwp), intent(in) :: n(:)
  integer(kind=iwp) :: key
  key = sum(n(1:ras_np)*ras_rdx(1:ras_np))
end function occ_key

!> All distributions V(1:np,i) of nel electrons of one spin over the np
!> partitions with capacities o (at most o(k) in partition k), nv of them.
subroutine occ_vectors(np,o,nel,nv,V)
  integer(kind=iwp), intent(in) :: np, o(np), nel
  integer(kind=iwp), intent(out) :: nv
  integer(kind=iwp), allocatable, intent(out) :: V(:,:)
  integer(kind=iwp) :: pass, j, w(MxGAS)
  logical(kind=iwp) :: more
  do pass=1,2
    nv = 0
    w(:) = 0
    do
      w(np) = nel-sum(w(1:np-1))
      if ((w(np) >= 0) .and. (w(np) <= o(np))) then
        nv = nv+1
        if (pass == 2) V(1:np,nv) = w(1:np)
      end if
      more = .false.
      do j=1,np-1
        if ((w(j) < o(j)) .and. (sum(w(1:np-1)) < nel)) then
          w(j) = w(j)+1
          more = .true.
          exit
        end if
        w(j) = 0
      end do
      if (.not. more) exit
    end do
    if (pass == 1) call mma_allocate(V,np,max(1,nv),label='ras_occv')
  end do
end subroutine occ_vectors

!> The occupation classes of a spin with nel electrons (the other spin has
!> nelo): all occupation vectors that occur in at least one determinant of
!> the space; occ_skey/occ_sidx the sorted packed occupations.
subroutine mk_classes(s,nel,nelo)
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp), intent(in) :: nel, nelo
  integer(kind=iwp) :: np, nv, nw, iv, iw, k, kk, j, t
  integer(kind=iwp), allocatable :: V(:,:), W(:,:)
  logical(kind=iwp), allocatable :: keep(:)

  s%nel = nel
  np = ras_np
  call occ_vectors(np,ras_o(1:np),nel,nv,V)
  call occ_vectors(np,ras_o(1:np),nelo,nw,W)
  call mma_allocate(keep,max(1,nv),label='ras_keep')
  keep(:) = .false.
  do iv=1,nv
    do iw=1,nw
      if (allowed(V(:,iv),W(:,iw))) then
        keep(iv) = .true.
        exit
      end if
    end do
  end do
  k = count(keep(1:nv))
  s%nocc = k
  call mma_allocate(s%occ_n,np,max(1,k),label='occ_n')
  call mma_allocate(s%occ_cnt,np,max(1,k),label='occ_cnt')
  call mma_allocate(s%occ_nstr,max(1,k),label='occ_nstr')
  call mma_allocate(s%occ_off,max(1,k),label='occ_off')
  call mma_allocate(s%occ_skey,max(1,k),label='occ_skey')
  call mma_allocate(s%occ_sidx,max(1,k),label='occ_sidx')
  k = 0
  do iv=1,nv
    if (.not. keep(iv)) cycle
    k = k+1
    s%occ_n(:,k) = V(:,iv)
    do j=1,np
      s%occ_cnt(j,k) = binom_coef(V(j,iv),ras_o(j))
    end do
    s%occ_nstr(k) = product(s%occ_cnt(:,k))
    s%occ_skey(k) = occ_key(V(:,iv))
    s%occ_sidx(k) = k
  end do
  ! sort the keys (insertion sort, few classes)
  do k=2,s%nocc
    kk = k
    do while (kk > 1)
      if (s%occ_skey(kk-1) <= s%occ_skey(kk)) exit
      t = s%occ_skey(kk)
      s%occ_skey(kk) = s%occ_skey(kk-1)
      s%occ_skey(kk-1) = t
      t = s%occ_sidx(kk)
      s%occ_sidx(kk) = s%occ_sidx(kk-1)
      s%occ_sidx(kk-1) = t
      kk = kk-1
    end do
  end do
  s%nstr = 0
  do k=1,s%nocc
    s%occ_off(k) = s%nstr
    s%nstr = s%nstr+s%occ_nstr(k)
  end do
  call mma_deallocate(keep)
  call mma_deallocate(V)
  call mma_deallocate(W)
end subroutine mk_classes

!> The class with the packed occupation key (0 if none).
function occ_find(s,key) result(k)
  type(ras_spin), intent(in) :: s
  integer(kind=iwp), intent(in) :: key
  integer(kind=iwp) :: k
  integer(kind=iwp) :: iLo, iHi, iMid
  k = 0
  iLo = 1
  iHi = s%nocc
  do while (iLo <= iHi)
    iMid = (iLo+iHi)/2
    if (s%occ_skey(iMid) == key) then
      k = s%occ_sidx(iMid)
      return
    else if (s%occ_skey(iMid) < key) then
      iLo = iMid+1
    else
      iHi = iMid-1
    end if
  end do
end function occ_find

!> The strings of all occupation classes, grouped by class and symmetry:
!> the groups (cls_*), the bit strings in group order and the map perm
!> from the rank address to the group order.
subroutine mk_strings(s)
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp) :: k, g, j, isym, ng, c, gnew, ip(MxGAS), pt(MxGAS)
  integer(kind=iwp), allocatable :: bits(:), sym(:), cnt(:,:), pos(:,:)
  logical(kind=iwp) :: more

  ! the strings in rank order and their symmetries
  call mma_allocate(bits,max(1,s%nstr),label='ras_bits')
  call mma_allocate(sym,max(1,s%nstr),label='ras_sym')
  call mma_allocate(cnt,ras_nsym,max(1,s%nocc),label='ras_cnt')
  cnt(:,:) = 0
  g = 0
  do k=1,s%nocc
    ! odometer over the substrings of the partitions, the first fastest
    do j=1,ras_np
      ip(j) = 1
      pt(j) = lex_init(s%occ_n(j,k),ras_o(j))
    end do
    do
      g = g+1
      bits(g) = 0
      do j=1,ras_np
        bits(g) = ior(bits(g),ishft(pt(j),ras_pos(j)))
      end do
      sym(g) = str_sym(bits(g))
      cnt(sym(g),k) = cnt(sym(g),k)+1
      more = .false.
      do j=1,ras_np
        if (ip(j) < s%occ_cnt(j,k)) then
          ip(j) = ip(j)+1
          pt(j) = lex_next(pt(j))
          more = .true.
          exit
        end if
        ip(j) = 1
        pt(j) = lex_init(s%occ_n(j,k),ras_o(j))
      end do
      if (.not. more) exit
    end do
  end do

  ! the groups: (class, symmetry) with strings
  ng = count(cnt(:,:) > 0)
  s%ncls = ng
  call mma_allocate(s%cls_n,ras_np,max(1,ng),label='cls_n')
  call mma_allocate(s%cls_sym,max(1,ng),label='cls_sym')
  call mma_allocate(s%cls_nstr,max(1,ng),label='cls_nstr')
  call mma_allocate(s%cls_off,max(1,ng),label='cls_off')
  call mma_allocate(pos,ras_nsym,max(1,s%nocc),label='ras_pos')
  pos(:,:) = 0
  c = 0
  g = 0
  do k=1,s%nocc
    do isym=1,ras_nsym
      if (cnt(isym,k) == 0) cycle
      c = c+1
      s%cls_n(:,c) = s%occ_n(:,k)
      s%cls_sym(c) = isym
      s%cls_nstr(c) = cnt(isym,k)
      s%cls_off(c) = g
      pos(isym,k) = c
      g = g+cnt(isym,k)
    end do
  end do

  ! the strings in group order
  call mma_allocate(s%str,max(1,s%nstr),label='ras_str')
  call mma_allocate(s%str_cls,max(1,s%nstr),label='ras_str_cls')
  call mma_allocate(s%str_idx,max(1,s%nstr),label='ras_str_idx')
  call mma_allocate(s%perm,max(1,s%nstr),label='ras_perm')
  cnt(:,:) = 0
  do k=1,s%nocc
    do g=s%occ_off(k)+1,s%occ_off(k)+s%occ_nstr(k)
      c = pos(sym(g),k)
      cnt(sym(g),k) = cnt(sym(g),k)+1
      gnew = s%cls_off(c)+cnt(sym(g),k)
      s%perm(g) = gnew
      s%str(gnew) = bits(g)
      s%str_cls(gnew) = c
      s%str_idx(gnew) = cnt(sym(g),k)
    end do
  end do
  do g=1,s%nstr
    if (str_index(s,s%str(g)) /= g) then
      write(u6,*) 'faroald_ras: inconsistent string address',g,str_index(s,s%str(g))
      call Abend()
    end if
  end do
  call mma_deallocate(bits)
  call mma_deallocate(sym)
  call mma_deallocate(cnt)
  call mma_deallocate(pos)
end subroutine mk_strings

!> The symmetry of a string (product of the symmetries of its orbitals).
function str_sym(c) result(isym)
  integer(kind=iwp), intent(in) :: c
  integer(kind=iwp) :: isym
  integer(kind=iwp) :: iorb
  isym = 0
  do iorb=1,ras_norb
    if (btest(c,iorb-1)) isym = ieor(isym,ras_osym(iorb)-1)
  end do
  isym = isym+1
end function str_sym

!> Bits pos+1..pos+len of c (0 for len = 0).
pure function sub_bits(c,pos,len) result(d)
  integer(kind=iwp), intent(in) :: c, pos, len
  integer(kind=iwp) :: d
  if (len == 0) then
    d = 0
  else
    d = ibits(c,pos,len)
  end if
end function sub_bits

!> Global address of the string with bits c (sign bit ignored) among the
!> strings of spin s, 0 if its class is not one of s.
function str_index(s,c) result(g)
  type(ras_spin), intent(in) :: s
  integer(kind=iwp), intent(in) :: c
  integer(kind=iwp) :: g
  integer(kind=iwp) :: b(MxGAS), n(MxGAS), j, k, idx, mult

  g = 0
  do j=1,ras_np
    b(j) = sub_bits(c,ras_pos(j),ras_o(j))
    n(j) = popcnt(b(j))
  end do
  if (sum(n(1:ras_np)) /= s%nel) return
  k = occ_find(s,occ_key(n(1:ras_np)))
  if (k == 0) return
  ! the rank address, the first partition fastest
  idx = 0
  mult = 1
  do j=1,ras_np
    idx = idx+(lexrank(b(j))-1)*mult
    mult = mult*s%occ_cnt(j,k)
  end do
  g = s%perm(s%occ_off(k)+idx+1)
end function str_index

!> Single excitations a+_p a_q of each string within the strings of s.
subroutine mk_singles(s)
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp) :: g, p, q, d, j, n, pass

  call mma_allocate(s%ex_off,s%nstr+1,label='ras_ex_off')
  do pass=1,2
    n = 0
    do g=1,s%nstr
      s%ex_off(g) = n
      do q=1,ras_norb
        if (.not. btest(s%str(g),q-1)) cycle
        do p=1,ras_norb
          d = ex1(p,q,s%str(g))
          if (d == -1) cycle
          j = str_index(s,ibits(d,0,31))
          if (j == 0) cycle
          n = n+1
          if (pass == 2) then
            s%ex_p(n) = p
            s%ex_q(n) = q
            s%ex_sgn(n) = fase(d)
            s%ex_j(n) = j
          end if
        end do
      end do
    end do
    s%ex_off(s%nstr+1) = n
    if (pass == 1) then
      s%nex = n
      call mma_allocate(s%ex_p,max(1,n),label='ras_ex_p')
      call mma_allocate(s%ex_q,max(1,n),label='ras_ex_q')
      call mma_allocate(s%ex_sgn,max(1,n),label='ras_ex_sgn')
      call mma_allocate(s%ex_j,max(1,n),label='ras_ex_j')
    end if
  end do
end subroutine mk_singles

!> The batches of the beta single excitations (bat_*, be_*): per beta
!> group and pair (v,x), in the order of the pair key (x-1)*norb+v.
subroutine mk_batches()
  integer(kind=iwp) :: c, gb, e, k, n, b, pass
  integer(kind=iwp), allocatable :: cnt(:), pos(:)

  call mma_allocate(cnt,ras_norb**2,label='ras_bcnt')
  call mma_allocate(pos,ras_norb**2,label='ras_bpos')
  do pass=1,2
    b = 0
    n = 0
    do c=1,sb%ncls
      cnt(:) = 0
      do gb=sb%cls_off(c)+1,sb%cls_off(c)+sb%cls_nstr(c)
        do e=sb%ex_off(gb)+1,sb%ex_off(gb+1)
          k = (sb%ex_q(e)-1)*ras_norb+sb%ex_p(e)
          cnt(k) = cnt(k)+1
        end do
      end do
      do k=1,ras_norb**2
        if (cnt(k) == 0) cycle
        b = b+1
        if (pass == 2) then
          bat_cj(b) = c
          bat_ci(b) = 0
          bat_v(b) = mod(k-1,ras_norb)+1
          bat_x(b) = (k-1)/ras_norb+1
          bat_off(b) = n
          mxbat = max(mxbat,cnt(k))
        end if
        pos(k) = b
        n = n+cnt(k)
      end do
      if (pass == 1) cycle
      do gb=sb%cls_off(c)+1,sb%cls_off(c)+sb%cls_nstr(c)
        do e=sb%ex_off(gb)+1,sb%ex_off(gb+1)
          k = pos((sb%ex_q(e)-1)*ras_norb+sb%ex_p(e))
          bat_off(k) = bat_off(k)+1
          be_j(bat_off(k)) = sb%str_idx(gb)
          be_i(bat_off(k)) = sb%str_idx(sb%ex_j(e))
          be_sgn(bat_off(k)) = sb%ex_sgn(e)
          if (bat_ci(k) == 0) bat_ci(k) = sb%str_cls(sb%ex_j(e))
          if (bat_ci(k) /= sb%str_cls(sb%ex_j(e))) then
            write(u6,*) 'faroald_ras: beta excitation batch with two target groups',k
            call Abend()
          end if
        end do
      end do
    end do
    if (pass == 1) then
      nbat = b
      mxbat = 0
      call mma_allocate(bat_cj,max(1,b),label='ras_bat_cj')
      call mma_allocate(bat_ci,max(1,b),label='ras_bat_ci')
      call mma_allocate(bat_v,max(1,b),label='ras_bat_v')
      call mma_allocate(bat_x,max(1,b),label='ras_bat_x')
      call mma_allocate(bat_off,b+1,label='ras_bat_off')
      call mma_allocate(be_j,max(1,n),label='ras_be_j')
      call mma_allocate(be_i,max(1,n),label='ras_be_i')
      call mma_allocate(be_sgn,max(1,n),label='ras_be_sgn')
    end if
  end do
  ! bat_off(b) was advanced to the end of batch b: shift back
  do b=nbat,1,-1
    bat_off(b+1) = bat_off(b)
  end do
  bat_off(1) = 0
  call mma_deallocate(cnt)
  call mma_deallocate(pos)
end subroutine mk_batches

!> Number of beta excitations per chunk of a batch for the gathered
!> work arrays (kc x nstr alpha, two of them, at most a quarter of the
!> available memory).
!> The (N-1)-string resolution tables of the alpha strings: the K strings
!> that occur (a_u Ja for some alpha string Ja), grouped by symmetry, and
!> res_cre. Not set up (nres_k = 0) if the full K space is too large for
!> the rank table.
subroutine mk_resolution()
  use sorting, only: argsort
  use sorting_funcs, only: leq_i
  integer(kind=iwp), parameter :: MaxRank = 50000000
  integer(kind=iwp) :: ga, u, t, d, k, nfull, nk, sym, kc, kbits, i, j, nb, km, cb
  integer(kind=iwp), allocatable :: na(:)
  integer(kind=iwp), allocatable :: r2c(:), kb(:), ksym(:), kocc(:,:), key(:), ordk(:)
  nres_k = 0
  if (sa%nel < 1) return
  nfull = binom_coef(sa%nel-1,ras_norb)
  if (nfull > MaxRank) return
  ! the K that occur, in order of first appearance
  call mma_allocate(r2c,nfull,label='res_r2c')
  call mma_allocate(kb,nfull,label='res_kb')
  r2c(:) = 0
  nk = 0
  do ga=1,sa%nstr
    do u=1,ras_norb
      d = ann(u,sa%str(ga))
      if (d == -1) cycle
      kbits = ibclr(d,31)
      k = lexrank(kbits)
      if (r2c(k) /= 0) cycle
      nk = nk+1
      r2c(k) = nk
      kb(nk) = kbits
    end do
  end do
  call mma_deallocate(r2c)
  ! occupations and symmetry of each K, and the order by block (occupation, symmetry)
  call mma_allocate(ksym,nk,label='res_ksym')
  call mma_allocate(kocc,ras_np,nk,label='res_kocc')
  call mma_allocate(key,nk,label='res_key')
  do kc=1,nk
    sym = 1
    do u=1,ras_norb
      if (btest(kb(kc),u-1)) sym = ieor(sym-1,ras_osym(u)-1)+1
    end do
    ksym(kc) = sym
    do j=1,ras_np
      kocc(j,kc) = popcnt(sub_bits(kb(kc),ras_pos(j),ras_o(j)))
    end do
  end do
  km = 0
  do kc=1,nk
    km = max(km,occ_key(kocc(:,kc)))
  end do
  do kc=1,nk
    key(kc) = (ksym(kc)-1)*(km+1)+occ_key(kocc(:,kc))
  end do
  call mma_allocate(ordk,nk,label='res_ordk')
  ordk(:) = argsort(key,leq_i)
  nb = 0
  do i=1,nk
    if (i == 1) then
      nb = 1
    else if (key(ordk(i)) /= key(ordk(i-1))) then
      nb = nb+1
    end if
  end do
  call mma_allocate(res_boff,nb+1,label='res_boff')
  call mma_allocate(res_bsym,nb,label='res_bsym')
  call mma_allocate(res_bocc,ras_np,nb,label='res_bocc')
  nb = 0
  do i=1,nk
    if ((i == 1) .or. (key(ordk(max(1,i-1))) /= key(ordk(i)))) then
      nb = nb+1
      res_boff(nb) = i-1
      res_bsym(nb) = ksym(ordk(i))
      res_bocc(:,nb) = kocc(:,ordk(i))
    end if
  end do
  res_boff(nb+1) = nk
  res_nblk = nb
  ! the coarse blocks: the fine blocks of one symmetry
  call mma_allocate(res_cboff,ras_nsym+1,label='res_cboff')
  call mma_allocate(res_cbsym,ras_nsym,label='res_cbsym')
  res_ncb = 0
  do i=1,nb
    if ((i == 1) .or. (res_bsym(max(1,i-1)) /= res_bsym(i))) then
      res_ncb = res_ncb+1
      res_cboff(res_ncb) = res_boff(i)
      res_cbsym(res_ncb) = res_bsym(i)
    end if
  end do
  res_cboff(res_ncb+1) = nk
  ! okb: the K of block B couple with beta group cb if an alpha class
  ! a+_u K (u in some partition) is allowed with cb
  call mma_allocate(res_okb,nb,sb%ncls,label='res_okb')
  call mma_allocate(res_cokb,res_ncb,sb%ncls,label='res_cokb')
  call mma_allocate(na,ras_np,label='res_na')
  res_okb(:,:) = .false.
  do i=1,nb
    do j=1,ras_np
      if (res_bocc(j,i) >= ras_o(j)) cycle
      na(:) = res_bocc(:,i)
      na(j) = na(j)+1
      if (occ_find(sa,occ_key(na)) == 0) cycle
      do cb=1,sb%ncls
        if (allowed(na,sb%cls_n(:,cb))) res_okb(i,cb) = .true.
      end do
    end do
  end do
  call mma_deallocate(na)
  res_cokb(:,:) = .false.
  j = 0
  do i=1,nb
    if ((i == 1) .or. (res_bsym(max(1,i-1)) /= res_bsym(i))) j = j+1
    res_cokb(j,:) = res_cokb(j,:) .or. res_okb(i,:)
  end do
  ! res_cre in the symmetry order of K
  call mma_allocate(res_cre,ras_norb,max(1,nk),label='res_cre')
  res_cre(:,:) = 0
  do i=1,nk
    kbits = kb(ordk(i))
    do t=1,ras_norb
      d = cre(t,kbits)
      if (d == -1) cycle
      res_cre(t,i) = fase(d)*str_index(sa,ibclr(d,31))
    end do
  end do
  call mma_deallocate(ordk)
  call mma_deallocate(key)
  call mma_deallocate(kocc)
  call mma_deallocate(ksym)
  call mma_deallocate(kb)
  ! the orbitals sorted by symmetry
  call mma_allocate(res_ord,ras_norb,label='res_ord')
  call mma_allocate(res_soff,ras_nsym,label='res_soff')
  call mma_allocate(res_ns,ras_nsym,label='res_ns')
  i = 0
  do sym=1,ras_nsym
    res_soff(sym) = i
    do u=1,ras_norb
      if (ras_osym(u) /= sym) cycle
      i = i+1
      res_ord(i) = u
    end do
    res_ns(sym) = i-res_soff(sym)
  end do
  nres_k = nk
end subroutine mk_resolution

!> The path of the alpha-beta sigma (res_mode): the direct loops, or the
!> resolution with fine or coarse blocks. timed: the candidates (within a
!> factor 8 of the cheapest estimate) are timed in the first sigma calls and
!> the fastest is used (different paths differ by round-off, so the results
!> are not reproducible to the last digit); else
!> the cheapest by an estimate of the cost (dir_counts, res_counts; ties in
!> the order direct, fine, coarse), which depends on the CI space only.
subroutine res_choose(timed)
  logical(kind=iwp), intent(in) :: timed
  real(kind=wp) :: c_dir, c_fine, c_coarse, cmin, d_el, d_op, d_gs, f_gs, f_call, f_fl, c_gs, c_call, c_fl
  logical(kind=iwp) :: coarse
  character(len=*), parameter :: Path(0:2) = ['direct loops               ','resolution, fine blocks    ', &
                                              'resolution, coarse blocks  ']
  res_mode = 0
  res_ntry = 0
  res_itry = 0
  res_tim(:) = Zero
  if (nres_k == 0) return
  call dir_counts(d_el,d_op,d_gs)
  call res_counts(res_nblk,res_boff,res_bsym,res_okb,f_gs,f_call,f_fl)
  call res_counts(res_ncb,res_cboff,res_cbsym,res_cokb,c_gs,c_call,c_fl)
  c_dir = d_el+Cost_Op*d_op+Cost_GS*d_gs
  c_fine = Cost_GS*f_gs+Cost_Call*f_call+Cost_Flop*f_fl
  c_coarse = Cost_GS*c_gs+Cost_Call*c_call+Cost_Flop*c_fl
  ! (CAS: the coarse blocks are the fine ones)
  coarse = res_ncb /= res_nblk
  if (.not. coarse) c_coarse = c_fine
  cmin = min(c_dir,c_fine,c_coarse)
  if (timed) then
    ! the candidates within a factor 8 of the cheapest estimate
    res_ntry = 0
    if (c_dir <= 8.0_wp*cmin) then
      res_ntry = res_ntry+1
      res_try(res_ntry) = 0
    end if
    if (c_fine <= 8.0_wp*cmin) then
      res_ntry = res_ntry+1
      res_try(res_ntry) = 1
    end if
    if (coarse .and. (c_coarse <= 8.0_wp*cmin)) then
      res_ntry = res_ntry+1
      res_try(res_ntry) = 2
    end if
    res_mode = merge(-1,res_try(1),res_ntry > 1)
  else if (c_dir <= cmin) then
    res_mode = 0
  else if (c_fine <= cmin) then
    res_mode = 1
  else
    res_mode = 2
  end if
  write(u6,'(1X,A,3ES11.3)') 'Faroald alpha-beta sigma, cost estimate (direct, resolution fine, coarse):',c_dir,c_fine,c_coarse
  if (res_mode >= 0) then
    write(u6,'(1X,2A)') 'Faroald alpha-beta sigma uses: ',trim(Path(res_mode))
  else
    write(u6,'(1X,A)') 'Faroald alpha-beta sigma: path timed in the first sigma calls'
  end if
end subroutine res_choose

!> The work of the direct loops (opposite_spin_dir): el the elements of the
!> vector updates, op their number (one per alpha excitation whose integral
!> is allowed by symmetry, for singlets only ka >= kb), gs the gathered and
!> scattered elements
subroutine dir_counts(el,op,gs)
  real(kind=wp), intent(out) :: el, op, gs
  integer(kind=iwp) :: n2, b, ca, ga, ea, ka, kb, nk, kc, nch, sk
  real(kind=wp) :: na
  real(kind=wp), allocatable :: h(:,:), m(:,:)
  integer(kind=iwp), allocatable :: psym(:)
  n2 = ras_norb**2
  call mma_allocate(psym,n2,label='res_psym')
  do kb=1,n2
    psym(kb) = ieor(ras_osym(mod(kb-1,ras_norb)+1)-1,ras_osym((kb-1)/ras_norb+1)-1)
  end do
  ! h(ka,ca): the alpha excitations of pair ka of the strings of group ca;
  ! m(kb,ca): those that couple with the beta pair kb
  call mma_allocate(h,n2,sa%ncls,label='res_h')
  call mma_allocate(m,n2,sa%ncls,label='res_m')
  h(:,:) = Zero
  do ca=1,sa%ncls
    do ga=sa%cls_off(ca)+1,sa%cls_off(ca)+sa%cls_nstr(ca)
      do ea=sa%ex_off(ga)+1,sa%ex_off(ga+1)
        ka = (sa%ex_q(ea)-1)*ras_norb+sa%ex_p(ea)
        h(ka,ca) = h(ka,ca)+One
      end do
    end do
  end do
  m(:,:) = Zero
  do ca=1,sa%ncls
    do kb=1,n2
      sk = psym(kb)
      do ka=merge(kb,1,ras_singlet),n2
        if (psym(ka) == sk) m(kb,ca) = m(kb,ca)+h(ka,ca)
      end do
    end do
  end do
  kc = chunk_size()
  el = Zero
  op = Zero
  gs = Zero
  do b=1,nbat
    kb = (bat_x(b)-1)*ras_norb+bat_v(b)
    nk = bat_off(b+1)-bat_off(b)
    nch = (nk+kc-1)/kc
    na = Zero
    do ca=1,sa%ncls
      if (blk_off(ca,bat_cj(b)) >= 0) then
        na = na+m(kb,ca)
        gs = gs+real(nk,kind=wp)*real(sa%cls_nstr(ca),kind=wp)
      end if
      if (blk_off(ca,bat_ci(b)) >= 0) gs = gs+real(2*nk,kind=wp)*real(sa%cls_nstr(ca),kind=wp)
    end do
    el = el+real(nk,kind=wp)*na
    op = op+real(nch,kind=wp)*na
  end do
  call mma_deallocate(m)
  call mma_deallocate(h)
  call mma_deallocate(psym)
end subroutine dir_counts

!> The work of the resolution with the blocks boff/bsym/okb
!> (opposite_spin_res): gs the gathered and scattered elements, nc the
!> DGEMM calls, fl their multiply-adds
subroutine res_counts(nb,boff,bsym,okb,gs,nc,fl)
  integer(kind=iwp), intent(in) :: nb, boff(nb+1), bsym(nb)
  logical(kind=iwp), intent(in) :: okb(nb,sb%ncls)
  real(kind=wp), intent(out) :: gs, nc, fl
  integer(kind=iwp) :: iblk, nk, jb, ib, e, cbj, usym, nu, tsym, nt
  gs = Zero
  nc = Zero
  fl = Zero
  do iblk=1,nb
    nk = boff(iblk+1)-boff(iblk)
    do jb=1,sb%nstr
      cbj = sb%str_cls(jb)
      if (.not. okb(iblk,cbj)) cycle
      usym = ieor(ieor(ras_lsym-1,sb%cls_sym(cbj)-1),bsym(iblk)-1)+1
      nu = res_ns(usym)
      if (nu == 0) cycle
      ! gather and scatter (the scatter of Ib = jb's share: n per K)
      gs = gs+real(nk,kind=wp)*real(nu+2*ras_norb,kind=wp)
      do e=sb%ex_off(jb)+1,sb%ex_off(jb+1)
        ib = sb%ex_j(e)
        if (.not. okb(iblk,sb%str_cls(ib))) cycle
        tsym = ieor(ieor(usym-1,ras_osym(sb%ex_p(e))-1),ras_osym(sb%ex_q(e))-1)+1
        nt = res_ns(tsym)
        if (nt == 0) cycle
        nc = nc+One
        fl = fl+real(nk,kind=wp)*real(nu,kind=wp)*real(nt,kind=wp)
      end do
    end do
  end do
end subroutine res_counts

function chunk_size() result(kc)
  use stdalloc, only: mma_maxDBLE
  integer(kind=iwp) :: kc
  integer(kind=iwp) :: mxmem
  call mma_maxDBLE(mxmem)
  kc = max(1,min(mxbat,512,mxmem/(8*max(1,sa%nstr))))
end function chunk_size

!> The same-spin double excitations of each string within the strings of s
!> (lists, if they fit in a quarter of the available memory).
subroutine mk_doubles(s)
  use stdalloc, only: mma_maxDBLE
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp) :: g, t, u, v, x, d, dd, j, n, pass, mxmem

  call mma_allocate(s%dex_off,s%nstr+1,label='ras_dex_off')
  do pass=1,2
    n = 0
    do g=1,s%nstr
      s%dex_off(g) = n
      do u=2,ras_norb
        if (.not. btest(s%str(g),u-1)) cycle
        do x=1,u-1
          d = ann2(x,u,s%str(g))
          if (d == -1) cycle
          do t=2,ras_norb
            if (btest(d,t-1)) cycle
            do v=1,t-1
              dd = cre2(t,v,d)
              if (dd == -1) cycle
              j = str_index(s,ibits(dd,0,31))
              if (j == 0) cycle
              n = n+1
              if (pass == 2) then
                s%dex_j(n) = j
                s%dex_sgn(n) = fase(dd)
                s%dex_tuvx(n) = pack4(t,u,v,x)
              end if
            end do
          end do
        end do
      end do
    end do
    s%dex_off(s%nstr+1) = n
    if (pass == 1) then
      call mma_maxDBLE(mxmem)
      if (3*n > mxmem/4) then
        call mma_deallocate(s%dex_off)
        s%has_dex = .false.
        return
      end if
      s%ndex = n
      call mma_allocate(s%dex_j,max(1,n),label='ras_dex_j')
      call mma_allocate(s%dex_sgn,max(1,n),label='ras_dex_sgn')
      call mma_allocate(s%dex_tuvx,max(1,n),label='ras_dex_tuvx')
    end if
  end do
  s%has_dex = .true.
end subroutine mk_doubles

!> (t,u,v,x) packed into one integer, and back.
pure function pack4(t,u,v,x) result(k)
  integer(kind=iwp), intent(in) :: t, u, v, x
  integer(kind=iwp) :: k
  k = (t-1)+ras_norb*((u-1)+ras_norb*((v-1)+ras_norb*(x-1)))
end function pack4

pure subroutine unpack4(k,t,u,v,x)
  integer(kind=iwp), intent(in) :: k
  integer(kind=iwp), intent(out) :: t, u, v, x
  integer(kind=iwp) :: r
  r = k
  t = mod(r,ras_norb)+1
  r = r/ras_norb
  u = mod(r,ras_norb)+1
  r = r/ras_norb
  v = mod(r,ras_norb)+1
  x = r/ras_norb+1
end subroutine unpack4

!> Position (1-based) of the determinant |deta detb> in the RAS vector, 0
!> if it is not in the RAS space.
function ras_addr(deta,detb) result(k)
  integer(kind=iwp), intent(in) :: deta, detb
  integer(kind=iwp) :: k
  integer(kind=iwp) :: ga, gb, ca, cb
  k = 0
  ga = str_index(sa,deta)
  if (ga == 0) return
  gb = str_index(sb,detb)
  if (gb == 0) return
  ca = sa%str_cls(ga)
  cb = sb%str_cls(gb)
  if (blk_off(ca,cb) < 0) return
  k = blk_off(ca,cb)+sa%str_idx(ga)+(sb%str_idx(gb)-1)*sa%cls_nstr(ca)
end function ras_addr

!> The RAS vector from Faroald's full determinant matrix psi(ndeta,ndetb)
!> (lexical ranks of the whole strings); for verification.
subroutine ras_from_full(psi,vec)
  real(kind=wp), intent(in) :: psi(:,:)
  real(kind=wp), intent(out) :: vec(n_ras)
  integer(kind=iwp) :: ca, cb, ia, ib, k
  do cb=1,sb%ncls
    do ca=1,sa%ncls
      if (blk_off(ca,cb) < 0) cycle
      k = blk_off(ca,cb)
      do ib=1,sb%cls_nstr(cb)
        do ia=1,sa%cls_nstr(ca)
          k = k+1
          vec(k) = psi(lexrank(sa%str(sa%cls_off(ca)+ia)),lexrank(sb%str(sb%cls_off(cb)+ib)))
        end do
      end do
    end do
  end do
end subroutine ras_from_full

!> Faroald's full determinant matrix from the RAS vector (zero outside the
!> RAS space); for verification.
subroutine ras_to_full(vec,psi)
  real(kind=wp), intent(in) :: vec(n_ras)
  real(kind=wp), intent(out) :: psi(:,:)
  integer(kind=iwp) :: ca, cb, ia, ib, k
  psi(:,:) = Zero
  do cb=1,sb%ncls
    do ca=1,sa%ncls
      if (blk_off(ca,cb) < 0) cycle
      k = blk_off(ca,cb)
      do ib=1,sb%cls_nstr(cb)
        do ia=1,sa%cls_nstr(ca)
          k = k+1
          psi(lexrank(sa%str(sa%cls_off(ca)+ia)),lexrank(sb%str(sb%cls_off(cb)+ib))) = vec(k)
        end do
      end do
    end do
  end do
end subroutine ras_to_full

!> sgm = H vec in the RAS space (h, g over the orbitals of Faroald, as for
!> sigma_update; no core energy).
subroutine ras_sigma(h,g,vec,sgm)
  real(kind=wp), intent(in) :: h(ras_norb,ras_norb), g(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  real(kind=wp), intent(out) :: sgm(n_ras)
  real(kind=wp), allocatable :: vT(:), sT(:)

  sgm(:) = Zero
  ! beta-beta: the columns of the blocks
  call same_spin(sa,sb,blk_off,h,g,vec,sgm)
  if (ras_singlet) then
    ! half of alpha-beta, then sgm := sgm + sgm^T
    call opposite_spin(g,vec,sgm)
    call sym_add(sgm)
    return
  end if
  ! alpha-alpha: the columns of the transposed blocks
  call mma_allocate(vT,max(1,n_ras),label='ras_vT')
  call mma_allocate(sT,max(1,n_ras),label='ras_sT')
  call transpose_blocks(vec,vT,.true.)
  sT(:) = Zero
  call same_spin(sb,sa,blk_offT,h,g,vT,sT)
  call transpose_blocks(sT,vT,.false.)
  sgm(:) = sgm(:)+vT(1:n_ras)
  call mma_deallocate(vT)
  call mma_deallocate(sT)
  ! alpha-beta
  call opposite_spin(g,vec,sgm)

end subroutine ras_sigma

!> Singlets: v := v + v^T (v^T: alpha and beta strings exchanged).
subroutine sym_add(v)
  real(kind=wp), intent(inout) :: v(n_ras)
  integer(kind=iwp) :: ca, cb, na, nb, ia, ib, o, oT
  real(kind=wp) :: w
  do cb=1,sb%ncls
    do ca=cb,sa%ncls
      if (blk_off(ca,cb) < 0) cycle
      na = sa%cls_nstr(ca)
      nb = sb%cls_nstr(cb)
      o = blk_off(ca,cb)
      oT = blk_off(cb,ca)
      if (ca == cb) then
        do ib=1,nb
          v(o+ib+(ib-1)*na) = Two*v(o+ib+(ib-1)*na)
          do ia=ib+1,na
            w = v(o+ia+(ib-1)*na)+v(o+ib+(ia-1)*na)
            v(o+ia+(ib-1)*na) = w
            v(o+ib+(ia-1)*na) = w
          end do
        end do
      else
        do ib=1,nb
          do ia=1,na
            w = v(o+ia+(ib-1)*na)+v(oT+ib+(ia-1)*nb)
            v(o+ia+(ib-1)*na) = w
            v(oT+ib+(ia-1)*nb) = w
          end do
        end do
      end if
    end do
  end do
end subroutine sym_add

!> forward: v (blocks na x nb) -> w (blocks nb x na); backward the inverse.
subroutine transpose_blocks(v,w,forward)
  real(kind=wp), intent(in) :: v(:)
  real(kind=wp), intent(out) :: w(:)
  logical(kind=iwp), intent(in) :: forward
  integer(kind=iwp) :: ca, cb, na, nb, ia, ib, o, oT
  do cb=1,sb%ncls
    do ca=1,sa%ncls
      if (blk_off(ca,cb) < 0) cycle
      na = sa%cls_nstr(ca)
      nb = sb%cls_nstr(cb)
      o = blk_off(ca,cb)
      oT = blk_offT(cb,ca)
      do ib=1,nb
        do ia=1,na
          if (forward) then
            w(oT+ib+(ia-1)*nb) = v(o+ia+(ib-1)*na)
          else
            w(o+ia+(ib-1)*na) = v(oT+ib+(ia-1)*nb)
          end if
        end do
      end do
    end do
  end do
end subroutine transpose_blocks

!> The same-spin part of sigma for the column spin c (row spin r): for each
!> column string I the coefficients f(J) = <I|H_c|J> (one-electron and
!> normal-ordered two-electron same-spin terms), then
!> sgm(:,I) += f(J) vec(:,J) over the row classes of both blocks. off are
!> the block offsets (row class, column class).
!> comp(cj,ci): some row class cr is allowed (off >= 0) with both column
!> classes ci and cj, i.e. an excitation between strings of ci and cj of
!> the column spin can contribute; the rest is skipped by the same-spin
!> loops (RAS: most pairs of occupation classes have no common partner).
subroutine col_pairs(r,c,off,comp)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  logical(kind=iwp), allocatable, intent(out) :: comp(:,:)
  integer(kind=iwp) :: ci, cj, cr
  call mma_allocate(comp,max(1,c%ncls),max(1,c%ncls),label='ras_comp')
  comp(:,:) = .false.
  do ci=1,c%ncls
    do cj=1,c%ncls
      do cr=1,r%ncls
        if ((off(cr,ci) >= 0) .and. (off(cr,cj) >= 0)) then
          comp(cj,ci) = .true.
          exit
        end if
      end do
    end do
  end do
end subroutine col_pairs

subroutine same_spin(r,c,off,h,g,vec,sgm)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: h(ras_norb,ras_norb), g(ras_norb,ras_norb,ras_norb,ras_norb), vec(:)
  real(kind=wp), intent(inout) :: sgm(:)
  integer(kind=iwp) :: gi, gj, ci, cj, cr, e, nl, il, nr, oi, oj, t, u, v, x, d, dd, ii, ij
  integer(kind=iwp), allocatable :: lst(:)
  real(kind=wp), allocatable :: f(:)
  logical(kind=iwp), allocatable :: used(:), comp(:,:)

  call col_pairs(r,c,off,comp)
  call mma_allocate(f,max(1,c%nstr),label='ras_f')
  call mma_allocate(lst,max(1,c%nstr),label='ras_lst')
  call mma_allocate(used,max(1,c%nstr),label='ras_used')
  f(:) = Zero
  used(:) = .false.
  do gi=1,c%nstr
    ci = c%str_cls(gi)
    nl = 0
    ! (no row class allowed with the class of I: no contribution)
    if (.not. any(comp(:,ci))) cycle
    ! one-electron: J = a+_p a_q I
    do e=c%ex_off(gi)+1,c%ex_off(gi+1)
      gj = c%ex_j(e)
      if (.not. comp(c%str_cls(gj),ci)) cycle
      call add(gj,c%ex_sgn(e)*h(c%ex_p(e),c%ex_q(e)))
    end do
    ! two-electron, same spin: J = a+_t a+_v a_x a_u I, t > v, u > x
    if (c%has_dex) then
      do e=c%dex_off(gi)+1,c%dex_off(gi+1)
        gj = c%dex_j(e)
        if (.not. comp(c%str_cls(gj),ci)) cycle
        call unpack4(c%dex_tuvx(e),t,u,v,x)
        call add(gj,c%dex_sgn(e)*(g(t,u,v,x)-g(t,x,v,u)))
      end do
    else
      do u=2,ras_norb
        if (.not. btest(c%str(gi),u-1)) cycle
        do x=1,u-1
          d = ann2(x,u,c%str(gi))
          if (d == -1) cycle
          do t=2,ras_norb
            if (btest(d,t-1)) cycle
            do v=1,t-1
              dd = cre2(t,v,d)
              if (dd == -1) cycle
              gj = str_index(c,ibits(dd,0,31))
              if (gj == 0) cycle
              if (.not. comp(c%str_cls(gj),ci)) cycle
              call add(gj,fase(dd)*(g(t,u,v,x)-g(t,x,v,u)))
            end do
          end do
        end do
      end do
    end if
    ! sigma update
    ii = c%str_idx(gi)
    do il=1,nl
      gj = lst(il)
      if (f(gj) /= Zero) then
        cj = c%str_cls(gj)
        ij = c%str_idx(gj)
        do cr=1,r%ncls
          if ((off(cr,ci) < 0) .or. (off(cr,cj) < 0)) cycle
          nr = r%cls_nstr(cr)
          oi = off(cr,ci)+(ii-1)*nr
          oj = off(cr,cj)+(ij-1)*nr
          sgm(oi+1:oi+nr) = sgm(oi+1:oi+nr)+f(gj)*vec(oj+1:oj+nr)
        end do
      end if
      f(gj) = Zero
      used(gj) = .false.
    end do
  end do
  call mma_deallocate(f)
  call mma_deallocate(lst)
  call mma_deallocate(used)
  call mma_deallocate(comp)

contains

  subroutine add(k,val)
    integer(kind=iwp), intent(in) :: k
    real(kind=wp), intent(in) :: val
    if (.not. used(k)) then
      used(k) = .true.
      nl = nl+1
      lst(nl) = k
    end if
    f(k) = f(k)+val
  end subroutine add

end subroutine same_spin

!> The alpha-beta part: sgm(Ia,Ib) += sum_tuvx g_tuvx <Ia|e^a_tu|Ja>
!> <Ib|e^b_vx|Jb> vec(Ja,Jb), over the determinants of the RAS space only.
!> Per batch of beta excitations Jb_k -> Ib_k (one group, one (v,x)):
!> Cp(k,Ja) = sgn_k vec(Ja,Jb_k) over the alpha strings allowed with Jb,
!> Vw(:,Ia) += sgn g_tuvx Cp(:,Ja) over the alpha excitations Ja -> Ia
!> (vector loops over k), sgm(Ia,Ib_k) += Vw(k,Ia). Singlets: only the
!> pairs tu >= vx (pair key (u-1)*norb+t, the order of the excitation
!> lists), tu = vx with weight 1/2.
!> The alpha-beta part of sigma (as opposite_spin) through the (N-1)-string
!> resolution, in Lucia's loop order: for each batch of K strings (of one
!> symmetry), gather X(K,u) = <K|a_u|Ja> C(Ja,Jb) once per beta string Jb,
!> S(K,t,Ib) += sgn(Ib<-Jb) X G(vx) for every beta excitation Ib = a+_v a_x Jb
!> (one DGEMM over the orbitals u and t of the allowed symmetries), then
!> sigma(a+_t K,Ib) += <Ia|a+_t|K> S(K,t,Ib). done = .false. if not even one
!> K per batch fits in memory (then opposite_spin does it directly).
subroutine opposite_spin_res(mode,g,vec,sgm,done)
  use stdalloc, only: mma_maxDBLE
  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: g(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  real(kind=wp), intent(inout) :: sgm(n_ras)
  logical(kind=iwp), intent(out) :: done
  integer(kind=iwp) :: n, ul, tl, u, t, v, x, ka, kb, mxmem, maxkb, nsb
  real(kind=wp) :: w
  real(kind=wp), allocatable :: Gp(:,:,:,:)

  done = .false.
  n = ras_norb
  nsb = sb%nstr
  call mma_maxDBLE(mxmem)
  maxkb = int((real(mxmem,kind=wp)/3.0_wp-real(n,kind=wp)**4)/(real(n,kind=wp)*real(nsb+1,kind=wp)),kind=iwp)
  if (maxkb < 1) return
  done = .true.
  ! Gp(u',t',v,x): the weight of a+_t a_u a+_v a_x (u', t' in the symmetry order)
  call mma_allocate(Gp,n,n,n,n,label='res_Gp')
  Gp(:,:,:,:) = Zero
  do x=1,n
    do v=1,n
      kb = (x-1)*n+v
      do tl=1,n
        t = res_ord(tl)
        do ul=1,n
          u = res_ord(ul)
          ka = (u-1)*n+t
          if (ras_singlet .and. (ka < kb)) cycle
          w = g(t,u,v,x)
          if (ras_singlet .and. (ka == kb)) w = Half*w
          Gp(ul,tl,v,x) = w
        end do
      end do
    end do
  end do
  if (mode == 1) then
    maxkb = min(maxkb,maxval(res_boff(2:res_nblk+1)-res_boff(1:res_nblk)))
  else
    maxkb = min(maxkb,maxval(res_cboff(2:res_ncb+1)-res_cboff(1:res_ncb)))
  end if
  if (allocated(res_Sw)) then
    if (size(res_Sw) < maxkb*n*nsb) call mma_deallocate(res_Sw)
  end if
  if (.not. allocated(res_Sw)) call mma_allocate(res_Sw,maxkb*n*nsb,label='res_Sw')
  if (allocated(res_Xw)) then
    if (size(res_Xw) < maxkb*n) call mma_deallocate(res_Xw)
  end if
  if (.not. allocated(res_Xw)) call mma_allocate(res_Xw,maxkb*n,label='res_Xw')
  if (mode == 1) then
    call res_blocks(res_nblk,res_boff,res_bsym,res_okb,Gp,size(res_Xw),res_Xw,size(res_Sw),res_Sw,maxkb,vec,sgm)
  else
    call res_blocks(res_ncb,res_cboff,res_cbsym,res_cokb,Gp,size(res_Xw),res_Xw,size(res_Sw),res_Sw,maxkb,vec,sgm)
  end if
  call mma_deallocate(Gp)
end subroutine opposite_spin_res

!> The blocks of K of opposite_spin_res (fine or coarse, see res_mode).
subroutine res_blocks(nb,boff,bsym,okb,Gp,nX,Xs,nS,S,maxkb,vec,sgm)
  integer(kind=iwp), intent(in) :: nb, boff(nb+1), bsym(nb), nX, nS, maxkb
  logical(kind=iwp), intent(in) :: okb(nb,sb%ncls)
  real(kind=wp), intent(in) :: Gp(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  ! X(kk,ul) at kk+(ul-1)*nkb, S(kk,t',ib) at kk+(t'-1)*nkb+(ib-1)*nkb*norb
  real(kind=wp), intent(inout) :: Xs(nX), S(nS), sgm(n_ras)
  integer(kind=iwp) :: n, nsb, iblk, ks, k0, k1, nkb, ib, jb, cbj, cbi, sja, usym, nu, u0, ul, u, kk, kc, r, ia, ca, o, e, v, x, &
                       tsym, nt, t0, tl, t
  logical(kind=iwp) :: any_x
  n = ras_norb
  nsb = sb%nstr
  do iblk=1,nb
    ks = bsym(iblk)
    do k0=boff(iblk)+1,boff(iblk+1),maxkb
      k1 = min(k0+maxkb-1,boff(iblk+1))
      nkb = k1-k0+1
      do ib=1,nsb
        if (okb(iblk,sb%str_cls(ib))) S((ib-1)*nkb*n+1:ib*nkb*n) = Zero
      end do
      do jb=1,nsb
        cbj = sb%str_cls(jb)
        if (.not. okb(iblk,cbj)) cycle
        sja = ieor(ras_lsym-1,sb%cls_sym(cbj)-1)+1
        usym = ieor(sja-1,ks-1)+1
        nu = res_ns(usym)
        if (nu == 0) cycle
        u0 = res_soff(usym)
        ! X(K,u) = <K|a_u|Ja> C(Ja,Jb), Ja = a+_u K
        any_x = .false.
        do ul=1,nu
          u = res_ord(u0+ul)
          do kk=1,nkb
            kc = k0+kk-1
            r = res_cre(u,kc)
            Xs(kk+(ul-1)*nkb) = Zero
            if (r == 0) cycle
            ia = abs(r)
            ca = sa%str_cls(ia)
            if (blk_off(ca,cbj) < 0) cycle
            o = blk_off(ca,cbj)+sa%str_idx(ia)+(sb%str_idx(jb)-1)*sa%cls_nstr(ca)
            if (r > 0) then
              Xs(kk+(ul-1)*nkb) = vec(o)
            else
              Xs(kk+(ul-1)*nkb) = -vec(o)
            end if
            any_x = .true.
          end do
        end do
        if (.not. any_x) cycle
        ! S(:,t,Ib) += sgn X G(u,t;vx) for each beta excitation Ib = a+_v a_x Jb
        do e=sb%ex_off(jb)+1,sb%ex_off(jb+1)
          ib = sb%ex_j(e)
          v = sb%ex_p(e)
          x = sb%ex_q(e)
          tsym = ieor(ieor(usym-1,ras_osym(v)-1),ras_osym(x)-1)+1
          nt = res_ns(tsym)
          if (nt == 0) cycle
          t0 = res_soff(tsym)
          if (.not. okb(iblk,sb%str_cls(ib))) cycle
          call DGEMM_('N','N',nkb,nt,nu,real(sb%ex_sgn(e),kind=wp),Xs,nkb,Gp(u0+1,t0+1,v,x),n,One,S(1+t0*nkb+(ib-1)*nkb*n),nkb)
        end do
      end do
      ! sigma(Ia,Ib) += <Ia|a+_t|K> S(K,t,Ib)
      do ib=1,nsb
        cbi = sb%str_cls(ib)
        if (.not. okb(iblk,cbi)) cycle
        do tl=1,n
          t = res_ord(tl)
          do kk=1,nkb
            r = res_cre(t,k0+kk-1)
            if (r == 0) cycle
            ia = abs(r)
            ca = sa%str_cls(ia)
            if (blk_off(ca,cbi) < 0) cycle
            o = blk_off(ca,cbi)+sa%str_idx(ia)+(sb%str_idx(ib)-1)*sa%cls_nstr(ca)
            if (r > 0) then
              sgm(o) = sgm(o)+S(kk+(tl-1)*nkb+(ib-1)*nkb*n)
            else
              sgm(o) = sgm(o)-S(kk+(tl-1)*nkb+(ib-1)*nkb*n)
            end if
          end do
        end do
      end do
    end do
  end do
end subroutine res_blocks


subroutine opposite_spin(g,vec,sgm)
  real(kind=wp), intent(in) :: g(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  real(kind=wp), intent(inout) :: sgm(n_ras)
  integer(kind=iwp) :: m, k
  real(kind=wp) :: c0, c1
  logical(kind=iwp) :: done

  ! the method: chosen (res_mode >= 0), or the next candidate to be timed
  ! (CPU time of the intrinsic clock: CWTime has a resolution of 10 ms,
  ! too coarse for small CI spaces)
  if (res_mode < 0) then
    res_itry = res_itry+1
    m = res_try(res_itry)
    call cpu_time(c0)
  else
    m = res_mode
  end if
  done = .false.
  if (m > 0) call opposite_spin_res(m,g,vec,sgm,done)
  if (.not. done) call opposite_spin_dir(g,vec,sgm)
  if (res_mode < 0) then
    call cpu_time(c1)
    res_tim(m) = c1-c0
    if (res_itry == res_ntry) then
      ! every process computes the whole sigma: all must take the same path
      ! (different paths differ by round-off, and the processes would then
      ! hold different CI vectors), so the choice uses the sum of the timings
      call gadgop(res_tim,3_iwp,'+')
      res_mode = res_try(1)
      do k=2,res_ntry
        if (res_tim(res_try(k)) < res_tim(res_mode)) res_mode = res_try(k)
      end do
      write(u6,'(1X,A,3ES11.3)') 'Faroald alpha-beta sigma, s per call (direct, resolution fine, coarse; 0 = not tried):', &
                                res_tim(0:2)
      write(u6,'(1X,A,I0,A)') 'Faroald alpha-beta sigma uses path ',res_mode,' (0 direct, 1 fine, 2 coarse)'
    end if
  end if

end subroutine opposite_spin

!> The alpha-beta part of sigma by the direct loops over the beta and alpha
!> excitations (batched over beta excitations).
subroutine opposite_spin_dir(g,vec,sgm)
  real(kind=wp), intent(in) :: g(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  real(kind=wp), intent(inout) :: sgm(n_ras)
  integer(kind=iwp) :: b, cbj, cbi, v, x, kb, kc, k0, nk, k, caj, cai, naj, nai, o, ga, ea, ia, ka
  real(kind=wp) :: w
  real(kind=wp), allocatable :: Cp(:,:), Vw(:,:)

  kc = chunk_size()
  call mma_allocate(Cp,kc,max(1,sa%nstr),label='ras_Cp')
  call mma_allocate(Vw,kc,max(1,sa%nstr),label='ras_Vw')
  do b=1,nbat
    cbj = bat_cj(b)
    cbi = bat_ci(b)
    v = bat_v(b)
    x = bat_x(b)
    kb = (x-1)*ras_norb+v
    do k0=bat_off(b),bat_off(b+1)-1,kc
      nk = min(kc,bat_off(b+1)-k0)
      ! gather
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        naj = sa%cls_nstr(caj)
        do k=1,nk
          o = blk_off(caj,cbj)+(be_j(k0+k)-1)*naj
          Cp(k,sa%cls_off(caj)+1:sa%cls_off(caj)+naj) = be_sgn(k0+k)*vec(o+1:o+naj)
        end do
      end do
      do cai=1,sa%ncls
        if (blk_off(cai,cbi) < 0) cycle
        Vw(1:nk,sa%cls_off(cai)+1:sa%cls_off(cai)+sa%cls_nstr(cai)) = Zero
      end do
      ! alpha excitations
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        do ga=sa%cls_off(caj)+1,sa%cls_off(caj)+sa%cls_nstr(caj)
          do ea=sa%ex_off(ga+1),sa%ex_off(ga)+1,-1
            ka = (sa%ex_q(ea)-1)*ras_norb+sa%ex_p(ea)
            if (ras_singlet .and. (ka < kb)) exit
            w = g(sa%ex_p(ea),sa%ex_q(ea),v,x)
            if (w == Zero) cycle
            ia = sa%ex_j(ea)
            if (blk_off(sa%str_cls(ia),cbi) < 0) cycle
            w = sa%ex_sgn(ea)*w
            if (ras_singlet .and. (ka == kb)) w = Half*w
            Vw(1:nk,ia) = Vw(1:nk,ia)+w*Cp(1:nk,ga)
          end do
        end do
      end do
      ! scatter
      do cai=1,sa%ncls
        if (blk_off(cai,cbi) < 0) cycle
        nai = sa%cls_nstr(cai)
        do k=1,nk
          o = blk_off(cai,cbi)+(be_i(k0+k)-1)*nai
          sgm(o+1:o+nai) = sgm(o+1:o+nai)+Vw(k,sa%cls_off(cai)+1:sa%cls_off(cai)+nai)
        end do
      end do
    end do
  end do
  call mma_deallocate(Cp)
  call mma_deallocate(Vw)

end subroutine opposite_spin_dir

!> Transition one-particle densities between RAS vectors:
!> da(p,q) = <bra|a+_pa a_qa|ket>, db(p,q) = <bra|a+_pb a_qb|ket>,
!> d = da+db, sd = da-db (as Faroald's transition_one_pdm).
subroutine ras_trans_one_pdm(bra,ket,d,sd)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(out) :: d(ras_norb,ras_norb), sd(ras_norb,ras_norb)
  real(kind=wp), allocatable :: da(:,:), db(:,:), bT(:), kT(:)

  call mma_allocate(da,ras_norb,ras_norb,label='ras_da')
  call mma_allocate(db,ras_norb,ras_norb,label='ras_db')
  call mma_allocate(bT,max(1,n_ras),label='ras_bT')
  call mma_allocate(kT,max(1,n_ras),label='ras_kT')
  da(:,:) = Zero
  db(:,:) = Zero
  call col_one_pdm(sa,sb,blk_off,bra,ket,db)
  if (ras_singlet) then
    ! alpha = beta
    d(:,:) = Two*db(:,:)
    sd(:,:) = Zero
    call mma_deallocate(da)
    call mma_deallocate(db)
    call mma_deallocate(bT)
    call mma_deallocate(kT)
    return
  end if
  call transpose_blocks(bra,bT,.true.)
  call transpose_blocks(ket,kT,.true.)
  call col_one_pdm(sb,sa,blk_offT,bT,kT,da)
  d(:,:) = da(:,:)+db(:,:)
  sd(:,:) = da(:,:)-db(:,:)
  call mma_deallocate(da)
  call mma_deallocate(db)
  call mma_deallocate(bT)
  call mma_deallocate(kT)
end subroutine ras_trans_one_pdm

!> dm(p,q) += <bra|a+_p a_q|ket> for the column spin c (row spin r).
subroutine col_one_pdm(r,c,off,bra,ket,dm)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: bra(:), ket(:)
  real(kind=wp), intent(inout) :: dm(ras_norb,ras_norb)
  integer(kind=iwp) :: gi, gj, e, cr, nr, oi, oj
  do gi=1,c%nstr
    do e=c%ex_off(gi)+1,c%ex_off(gi+1)
      gj = c%ex_j(e)
      do cr=1,r%ncls
        if ((off(cr,c%str_cls(gi)) < 0) .or. (off(cr,c%str_cls(gj)) < 0)) cycle
        nr = r%cls_nstr(cr)
        oi = off(cr,c%str_cls(gi))+(c%str_idx(gi)-1)*nr
        oj = off(cr,c%str_cls(gj))+(c%str_idx(gj)-1)*nr
        dm(c%ex_p(e),c%ex_q(e)) = dm(c%ex_p(e),c%ex_q(e))+c%ex_sgn(e)*dot_product(bra(oj+1:oj+nr),ket(oi+1:oi+nr))
      end do
    end do
  end do
end subroutine col_one_pdm

!> Spin-free transition two-particle density between RAS vectors,
!> p2(t,u,v,x) = <bra|E_tu E_vx - delta_uv E_tx|ket>
!>             = sum_{s,s'} <bra|a+_ts a+_vs' a_xs' a_us|ket>
!> (as Faroald's transition_two_pdm).
subroutine ras_trans_two_pdm(bra,ket,p2)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(out) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  real(kind=wp), allocatable :: bT(:), kT(:)

  p2(:,:,:,:) = Zero
  ! same spin
  call col_two_pdm(sa,sb,blk_off,bra,ket,p2)
  if (ras_singlet) then
    ! alpha-alpha = beta-beta; alpha-beta only for the pairs tu >= vx
    p2(:,:,:,:) = Two*p2(:,:,:,:)
    call ab_two_pdm_any(bra,ket,p2)
    return
  end if
  call mma_allocate(bT,max(1,n_ras),label='ras_bT')
  call mma_allocate(kT,max(1,n_ras),label='ras_kT')
  call transpose_blocks(bra,bT,.true.)
  call transpose_blocks(ket,kT,.true.)
  call col_two_pdm(sb,sa,blk_offT,bT,kT,p2)
  call mma_deallocate(bT)
  call mma_deallocate(kT)
  ! opposite spin
  call ab_two_pdm_any(bra,ket,p2)
end subroutine ras_trans_two_pdm

!> The alpha-beta part of the 2-RDM: through the (N-1)-string resolution
!> with the blocks the sigma uses (res_mode 1/2), else (direct loops chosen
!> or not enough memory) by ab_two_pdm.
subroutine ab_two_pdm_any(bra,ket,p2)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  logical(kind=iwp) :: done
  done = .false.
  if (res_mode > 0) call ab_two_pdm_res(res_mode,bra,ket,p2,done)
  if (.not. done) call ab_two_pdm(bra,ket,p2)
end subroutine ab_two_pdm_any

!> The alpha-beta part of the 2-RDM (as ab_two_pdm) through the
!> (N-1)-string resolution: per batch of K strings, X(K,u,Jb) =
!> <K|a_u|Ja> ket(Ja,Jb) and Y(K,t,Ib) = <K|a_t|Ia> bra(Ia,Ib) for all beta
!> strings, M(t,u,vx) += sgn Y(:,:,Ib)^T X(:,:,Jb) for every beta excitation
!> Ib = a+_v a_x Jb (one DGEMM over the orbitals of the allowed symmetries),
!> and p2(t,u,v,x) += M(t,u,v,x) + M(v,x,t,u) (for singlets M is symmetric
!> under the exchange of the pairs, which ab_two_pdm uses instead).
!> done = .false. if not even one K per batch fits.
subroutine ab_two_pdm_res(mode,bra,ket,p2,done)
  use stdalloc, only: mma_maxDBLE
  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  logical(kind=iwp), intent(out) :: done
  integer(kind=iwp) :: n, nsb, mxmem, maxkb, t, u, v, x, tl, ul
  real(kind=wp), allocatable :: M(:,:,:,:)
  done = .false.
  n = ras_norb
  nsb = sb%nstr
  call mma_maxDBLE(mxmem)
  if (allocated(res_Xw)) mxmem = mxmem+size(res_Xw)
  if (allocated(res_Sw)) mxmem = mxmem+size(res_Sw)
  maxkb = int((real(mxmem,kind=wp)/3.0_wp-real(n,kind=wp)**4)/(Two*real(n,kind=wp)*real(nsb+1,kind=wp)),kind=iwp)
  if (maxkb < 1) return
  done = .true.
  if (mode == 1) then
    maxkb = min(maxkb,maxval(res_boff(2:res_nblk+1)-res_boff(1:res_nblk)))
  else
    maxkb = min(maxkb,maxval(res_cboff(2:res_ncb+1)-res_cboff(1:res_ncb)))
  end if
  if (allocated(res_Sw)) then
    if (size(res_Sw) < maxkb*n*nsb) call mma_deallocate(res_Sw)
  end if
  if (.not. allocated(res_Sw)) call mma_allocate(res_Sw,maxkb*n*nsb,label='res_Sw')
  if (allocated(res_Xw)) then
    if (size(res_Xw) < maxkb*n*nsb) call mma_deallocate(res_Xw)
  end if
  if (.not. allocated(res_Xw)) call mma_allocate(res_Xw,maxkb*n*nsb,label='res_Xw')
  ! M(t',u',v,x), t' and u' in the symmetry order of the orbitals
  call mma_allocate(M,n,n,n,n,label='res_M')
  M(:,:,:,:) = Zero
  if (mode == 1) then
    call ab_res_blocks(res_nblk,res_boff,res_bsym,res_okb,size(res_Xw),res_Xw,size(res_Sw),res_Sw,maxkb,bra,ket,M)
  else
    call ab_res_blocks(res_ncb,res_cboff,res_cbsym,res_cokb,size(res_Xw),res_Xw,size(res_Sw),res_Sw,maxkb,bra,ket,M)
  end if
  do x=1,n
    do v=1,n
      do ul=1,n
        u = res_ord(ul)
        do tl=1,n
          t = res_ord(tl)
          p2(t,u,v,x) = p2(t,u,v,x)+M(tl,ul,v,x)
          p2(v,x,t,u) = p2(v,x,t,u)+M(tl,ul,v,x)
        end do
      end do
    end do
  end do
  call mma_deallocate(M)
end subroutine ab_two_pdm_res

!> The blocks of K of ab_two_pdm_res (fine or coarse).
subroutine ab_res_blocks(nb,boff,bsym,okb,nX,Xk,nY,Yb,maxkb,bra,ket,M)
  integer(kind=iwp), intent(in) :: nb, boff(nb+1), bsym(nb), nX, nY, maxkb
  logical(kind=iwp), intent(in) :: okb(nb,sb%ncls)
  ! Xk(kk,u',jb) at kk+(u'-1)*nkb+(jb-1)*nkb*norb (ket), Yb(kk,t',ib) likewise (bra)
  real(kind=wp), intent(inout) :: Xk(nX), Yb(nY), M(ras_norb,ras_norb,ras_norb,ras_norb)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  integer(kind=iwp) :: n, nsb, iblk, ks, k0, k1, nkb, jb, ib, cb, sa_sym, usym, tsym, nu, nt, u0, t0, ul, u, kk, r, ia, ca, &
                       o, e, v, x
  n = ras_norb
  nsb = sb%nstr
  do iblk=1,nb
    ks = bsym(iblk)
    do k0=boff(iblk)+1,boff(iblk+1),maxkb
      k1 = min(k0+maxkb-1,boff(iblk+1))
      nkb = k1-k0+1
      ! the gathers of ket (Xk) and bra (Yb) for the beta strings that couple
      do jb=1,nsb
        cb = sb%str_cls(jb)
        if (.not. okb(iblk,cb)) cycle
        sa_sym = ieor(ras_lsym-1,sb%cls_sym(cb)-1)+1
        usym = ieor(sa_sym-1,ks-1)+1
        nu = res_ns(usym)
        if (nu == 0) cycle
        u0 = res_soff(usym)
        do ul=1,nu
          u = res_ord(u0+ul)
          do kk=1,nkb
            o = kk+(u0+ul-1)*nkb+(jb-1)*nkb*n
            Xk(o) = Zero
            Yb(o) = Zero
            r = res_cre(u,k0+kk-1)
            if (r == 0) cycle
            ia = abs(r)
            ca = sa%str_cls(ia)
            if (blk_off(ca,cb) < 0) cycle
            if (r > 0) then
              Xk(o) = ket(blk_off(ca,cb)+sa%str_idx(ia)+(sb%str_idx(jb)-1)*sa%cls_nstr(ca))
              Yb(o) = bra(blk_off(ca,cb)+sa%str_idx(ia)+(sb%str_idx(jb)-1)*sa%cls_nstr(ca))
            else
              Xk(o) = -ket(blk_off(ca,cb)+sa%str_idx(ia)+(sb%str_idx(jb)-1)*sa%cls_nstr(ca))
              Yb(o) = -bra(blk_off(ca,cb)+sa%str_idx(ia)+(sb%str_idx(jb)-1)*sa%cls_nstr(ca))
            end if
          end do
        end do
      end do
      ! M(t',u',v,x) += sgn Yb(:,t',Ib)^T Xk(:,u',Jb)
      do jb=1,nsb
        cb = sb%str_cls(jb)
        if (.not. okb(iblk,cb)) cycle
        sa_sym = ieor(ras_lsym-1,sb%cls_sym(cb)-1)+1
        usym = ieor(sa_sym-1,ks-1)+1
        nu = res_ns(usym)
        if (nu == 0) cycle
        u0 = res_soff(usym)
        do e=sb%ex_off(jb)+1,sb%ex_off(jb+1)
          ib = sb%ex_j(e)
          if (.not. okb(iblk,sb%str_cls(ib))) cycle
          v = sb%ex_p(e)
          x = sb%ex_q(e)
          tsym = ieor(ieor(usym-1,ras_osym(v)-1),ras_osym(x)-1)+1
          nt = res_ns(tsym)
          if (nt == 0) cycle
          t0 = res_soff(tsym)
          call DGEMM_('T','N',nt,nu,nkb,real(sb%ex_sgn(e),kind=wp),Yb(1+t0*nkb+(ib-1)*nkb*n),nkb,Xk(1+u0*nkb+(jb-1)*nkb*n), &
                      nkb,One,M(t0+1,u0+1,v,x),n)
        end do
      end do
    end do
  end do
end subroutine ab_res_blocks

!> The alpha-beta part of ras_trans_two_pdm, <bra|a+_ta a_ua a+_vb a_xb|ket>
!> added to (t,u,v,x) and (v,x,t,u), batched as opposite_spin:
!> Kp(k,Ja) = sgn_k ket(Ja,Jb_k), Bp(k,Ia) = bra(Ia,Ib_k), val = sgn
!> Bp(:,Ia).Kp(:,Ja). Singlets: P_ab(tuvx) = P_ab(vxtu), only the pairs
!> tu >= vx, 2 val to (t,u,v,x) and (v,x,t,u) (once for tu = vx).
subroutine ab_two_pdm(bra,ket,p2)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  integer(kind=iwp) :: b, cbj, cbi, v, x, kb, kc, k0, nk, k, caj, cai, naj, nai, o, ga, ea, ia, ka, t, u
  real(kind=wp) :: val
  real(kind=wp), allocatable :: Kp(:,:), Bp(:,:)

  kc = chunk_size()
  call mma_allocate(Kp,kc,max(1,sa%nstr),label='ras_Kp')
  call mma_allocate(Bp,kc,max(1,sa%nstr),label='ras_Bp')
  do b=1,nbat
    cbj = bat_cj(b)
    cbi = bat_ci(b)
    v = bat_v(b)
    x = bat_x(b)
    kb = (x-1)*ras_norb+v
    do k0=bat_off(b),bat_off(b+1)-1,kc
      nk = min(kc,bat_off(b+1)-k0)
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        naj = sa%cls_nstr(caj)
        do k=1,nk
          o = blk_off(caj,cbj)+(be_j(k0+k)-1)*naj
          Kp(k,sa%cls_off(caj)+1:sa%cls_off(caj)+naj) = be_sgn(k0+k)*ket(o+1:o+naj)
        end do
      end do
      do cai=1,sa%ncls
        if (blk_off(cai,cbi) < 0) cycle
        nai = sa%cls_nstr(cai)
        do k=1,nk
          o = blk_off(cai,cbi)+(be_i(k0+k)-1)*nai
          Bp(k,sa%cls_off(cai)+1:sa%cls_off(cai)+nai) = bra(o+1:o+nai)
        end do
      end do
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        do ga=sa%cls_off(caj)+1,sa%cls_off(caj)+sa%cls_nstr(caj)
          do ea=sa%ex_off(ga+1),sa%ex_off(ga)+1,-1
            t = sa%ex_p(ea)
            u = sa%ex_q(ea)
            ka = (u-1)*ras_norb+t
            if (ras_singlet .and. (ka < kb)) exit
            ia = sa%ex_j(ea)
            if (blk_off(sa%str_cls(ia),cbi) < 0) cycle
            val = sa%ex_sgn(ea)*dot_product(Bp(1:nk,ia),Kp(1:nk,ga))
            if (ras_singlet) then
              val = Two*val
              p2(t,u,v,x) = p2(t,u,v,x)+val
              if (ka /= kb) p2(v,x,t,u) = p2(v,x,t,u)+val
            else
              p2(t,u,v,x) = p2(t,u,v,x)+val
              p2(v,x,t,u) = p2(v,x,t,u)+val
            end if
          end do
        end do
      end do
    end do
  end do
  call mma_deallocate(Kp)
  call mma_deallocate(Bp)

end subroutine ab_two_pdm

!> p2 += <bra|a+_t a+_v a_x a_u|ket> (all index orders) for the column
!> spin c (row spin r); J = a+_t a+_v a_x a_u I with t > v, u > x.
subroutine col_two_pdm(r,c,off,bra,ket,p2)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: bra(:), ket(:)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  integer(kind=iwp) :: gi, e, t, u, v, x, d, dd, gj
  logical(kind=iwp), allocatable :: comp(:,:)
  call col_pairs(r,c,off,comp)
  do gi=1,c%nstr
    if (.not. any(comp(:,c%str_cls(gi)))) cycle
    if (c%has_dex) then
      do e=c%dex_off(gi)+1,c%dex_off(gi+1)
        gj = c%dex_j(e)
        if (.not. comp(c%str_cls(gj),c%str_cls(gi))) cycle
        call unpack4(c%dex_tuvx(e),t,u,v,x)
        call put(gi,gj,c%dex_sgn(e),t,u,v,x)
      end do
    else
      do u=2,ras_norb
        if (.not. btest(c%str(gi),u-1)) cycle
        do x=1,u-1
          d = ann2(x,u,c%str(gi))
          if (d == -1) cycle
          do t=2,ras_norb
            if (btest(d,t-1)) cycle
            do v=1,t-1
              dd = cre2(t,v,d)
              if (dd == -1) cycle
              gj = str_index(c,ibits(dd,0,31))
              if (gj == 0) cycle
              if (.not. comp(c%str_cls(gj),c%str_cls(gi))) cycle
              call put(gi,gj,fase(dd),t,u,v,x)
            end do
          end do
        end do
      end do
    end if
  end do
  call mma_deallocate(comp)

contains

  ! the contribution of J = a+_t a+_v a_x a_u I (sign sg) in all index orders
  subroutine put(gi,gj,sg,t,u,v,x)
    integer(kind=iwp), intent(in) :: gi, gj, sg, t, u, v, x
    integer(kind=iwp) :: cr, nr, oi, oj
    real(kind=wp) :: val
    val = Zero
    do cr=1,r%ncls
      if ((off(cr,c%str_cls(gi)) < 0) .or. (off(cr,c%str_cls(gj)) < 0)) cycle
      nr = r%cls_nstr(cr)
      oi = off(cr,c%str_cls(gi))+(c%str_idx(gi)-1)*nr
      oj = off(cr,c%str_cls(gj))+(c%str_idx(gj)-1)*nr
      val = val+dot_product(bra(oj+1:oj+nr),ket(oi+1:oi+nr))
    end do
    if (val == Zero) return
    val = sg*val
    p2(t,u,v,x) = p2(t,u,v,x)+val
    p2(v,x,t,u) = p2(v,x,t,u)+val
    p2(t,x,v,u) = p2(t,x,v,u)-val
    p2(v,u,t,x) = p2(v,u,t,x)-val
  end subroutine put

end subroutine col_two_pdm

!> out = sum_pq h1(p,q) E_pq vin (both spins).
subroutine ras_apply_1e(h1,vin,vout)
  real(kind=wp), intent(in) :: h1(ras_norb,ras_norb), vin(n_ras)
  real(kind=wp), intent(out) :: vout(n_ras)
  real(kind=wp), allocatable :: iT(:), oT(:)
  vout(:) = Zero
  call col_apply_1e(sa,sb,blk_off,h1,vin,vout)
  call mma_allocate(iT,max(1,n_ras),label='ras_iT')
  call mma_allocate(oT,max(1,n_ras),label='ras_oT')
  call transpose_blocks(vin,iT,.true.)
  oT(:) = Zero
  call col_apply_1e(sb,sa,blk_offT,h1,iT,oT)
  call transpose_blocks(oT,iT,.false.)
  vout(:) = vout(:)+iT(1:n_ras)
  call mma_deallocate(iT)
  call mma_deallocate(oT)
end subroutine ras_apply_1e

!> vout(:,J) += h1(p,q) <J|a+_p a_q|I> vin(:,I) for the column spin c.
subroutine col_apply_1e(r,c,off,h1,vin,vout)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: h1(ras_norb,ras_norb), vin(:)
  real(kind=wp), intent(inout) :: vout(:)
  integer(kind=iwp) :: gi, gj, e, cr, nr, oi, oj
  real(kind=wp) :: f
  do gi=1,c%nstr
    do e=c%ex_off(gi)+1,c%ex_off(gi+1)
      f = c%ex_sgn(e)*h1(c%ex_p(e),c%ex_q(e))
      if (f == Zero) cycle
      gj = c%ex_j(e)
      do cr=1,r%ncls
        if ((off(cr,c%str_cls(gi)) < 0) .or. (off(cr,c%str_cls(gj)) < 0)) cycle
        nr = r%cls_nstr(cr)
        oi = off(cr,c%str_cls(gi))+(c%str_idx(gi)-1)*nr
        oj = off(cr,c%str_cls(gj))+(c%str_idx(gj)-1)*nr
        vout(oj+1:oj+nr) = vout(oj+1:oj+nr)+f*vin(oi+1:oi+nr)
      end do
    end do
  end do
end subroutine col_apply_1e

!> The CI transformation of traci_faroald (Malmqvist, one orbital at a
!> time) in the RAS space; T (norb x norb over the orbitals of Faroald)
!> must not mix the RAS partitions, so that every step stays in the RAS
!> space: elements between partitions up to ThrMix (rounding noise) are set
!> to zero in a local copy Tc, larger ones are an error.
subroutine ras_traci(T,nT,vec)
  integer(kind=iwp), intent(in) :: nT
  real(kind=wp), intent(in) :: T(nT)
  real(kind=wp), intent(inout) :: vec(n_ras)
  real(kind=wp), allocatable :: Tc(:), h1(:,:), v0(:), v1(:), v2(:)
  integer(kind=iwp) :: k, p, prt(30), ca, cb, ia, ib, nk, o
  real(kind=wp) :: tkk, fac(0:2)
  real(kind=wp), parameter :: ThrMix = 1.0e-10_wp

  do k=1,ras_np
    prt(ras_pos(k)+1:ras_pos(k)+ras_o(k)) = k
  end do
  call mma_allocate(Tc,nT,label='ras_Tc')
  Tc(:) = T(:)
  do k=1,ras_norb
    do p=1,ras_norb
      if (prt(p) == prt(k)) cycle
      if (abs(Tc((k-1)*ras_norb+p)) > ThrMix) then
        write(u6,*) 'ras_traci: the transformation mixes RAS partitions',p,k,Tc((k-1)*ras_norb+p)
        call Abend()
      end if
      Tc((k-1)*ras_norb+p) = Zero
    end do
  end do

  call mma_allocate(h1,ras_norb,ras_norb,label='ras_h1')
  call mma_allocate(v0,max(1,n_ras),label='ras_v0')
  call mma_allocate(v1,max(1,n_ras),label='ras_v1')
  call mma_allocate(v2,max(1,n_ras),label='ras_v2')
  do k=1,ras_norb
    tkk = Tc((k-1)*ras_norb+k)
    h1(:,:) = Zero
    do p=1,ras_norb
      if (p /= k) h1(p,k) = Tc((k-1)*ras_norb+p)/tkk
    end do
    fac(0) = One
    fac(1) = tkk
    fac(2) = tkk*tkk
    ! v0 = tkk^n_k vec
    do cb=1,sb%ncls
      do ca=1,sa%ncls
        if (blk_off(ca,cb) < 0) cycle
        o = blk_off(ca,cb)
        do ib=1,sb%cls_nstr(cb)
          do ia=1,sa%cls_nstr(ca)
            o = o+1
            nk = 0
            if (btest(sa%str(sa%cls_off(ca)+ia),k-1)) nk = nk+1
            if (btest(sb%str(sb%cls_off(cb)+ib),k-1)) nk = nk+1
            v0(o) = fac(nk)*vec(o)
          end do
        end do
      end do
    end do
    call ras_apply_1e(h1,v0,v1)
    call ras_apply_1e(h1,v1,v2)
    vec(:) = v0(1:n_ras)+v1(1:n_ras)+Half*v2(1:n_ras)
  end do
  call mma_deallocate(Tc)
  call mma_deallocate(h1)
  call mma_deallocate(v0)
  call mma_deallocate(v1)
  call mma_deallocate(v2)
end subroutine ras_traci

!> The number of determinants of the RAS/GAS space (all symmetries),
!> without setting anything up: nela/nelb electrons, o(k) orbitals in
!> partition k (np partitions), mn(k) to mx(k) electrons in the partitions
!> 1..k; in floating point (no overflow), for the memory estimate of RASSCF.
function ras_count_dets(nela,nelb,np,o,mn,mx) result(n)
  integer(kind=iwp), intent(in) :: nela, nelb, np, o(np), mn(np), mx(np)
  real(kind=wp) :: n
  integer(kind=iwp) :: nva, nvb, ia, ib, k, ne
  integer(kind=iwp), allocatable :: VA(:,:), VB(:,:)
  real(kind=wp) :: fa
  logical(kind=iwp) :: ok
  n = Zero
  call occ_vectors(np,o,nela,nva,VA)
  call occ_vectors(np,o,nelb,nvb,VB)
  do ia=1,nva
    fa = One
    do k=1,np
      fa = fa*bin_r(VA(k,ia),o(k))
    end do
    do ib=1,nvb
      ok = .true.
      ne = 0
      do k=1,np
        ne = ne+VA(k,ia)+VB(k,ib)
        if ((ne < mn(k)) .or. (ne > mx(k))) ok = .false.
      end do
      if (.not. ok) cycle
      do k=1,np
        fa = fa*bin_r(VB(k,ib),o(k))
      end do
      n = n+fa
      fa = One
      do k=1,np
        fa = fa*bin_r(VA(k,ia),o(k))
      end do
    end do
  end do
  call mma_deallocate(VA)
  call mma_deallocate(VB)

contains

  pure function bin_r(k,m) result(b)
    integer(kind=iwp), intent(in) :: k, m
    real(kind=wp) :: b
    integer(kind=iwp) :: i
    b = Zero
    if ((k < 0) .or. (k > m)) return
    b = One
    do i=1,min(k,m-k)
      b = b*real(m-min(k,m-k)+i,kind=wp)/real(i,kind=wp)
    end do
  end function bin_r

end function ras_count_dets

end module faroald_ras
