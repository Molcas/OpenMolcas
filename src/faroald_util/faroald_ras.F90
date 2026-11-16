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
! String-driven RAS for Faroald (Olsen et al., J. Chem. Phys. 89, 2185
! (1988)): the alpha and beta strings are grouped in occupation classes
! (n1,n2,n3), the numbers of electrons in RAS1, RAS2 and RAS3; the CI
! vector holds only the blocks (alpha class, beta class) of the RAS space,
! n1a+n1b >= minr1 and n3a+n3b <= maxr3. No string or determinant outside
! the RAS space is ever addressed.
!
! The orbitals are those of Faroald, the levels of the split graph: RAS1
! = 1..o1, RAS2 = o1+1..o1+o2, RAS3 = o1+o2+1..norb (contiguous bit
! ranges). Within a class a string is addressed by the lexical ranks of
! its three substrings, r1 fastest. The element of determinant
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
! A CAS space is the special case o1 = o3 = 0 (one class per spin).
!
! Singlets (spin_flip, nela = nelb, M_S = S = 0): the vectors are
! symmetric under the exchange of the alpha and beta strings,
! C(Ia,Ib) = C(Ib,Ia) (Faroald's phase convention, as in sigma_update).
! Then the alpha and beta strings are the same, block (ga,gb) is the
! transpose of block (gb,ga), and sigma = X + X^T with X = the beta-beta
! part plus half of the alpha-beta part (pairs tu >= vx of
! g_tuvx a+_ta a_ua a+_vb a_xb, tu = vx with weight 1/2); the alpha-alpha
! density equals the beta-beta one and the spin density is zero.

use second_quantization, only: ann2, binom_coef, cre2, ex1, fase, lex_init, lex_next, lexrank
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
private

type ras_spin
  integer(kind=iwp) :: nel = 0, nocc = 0, ncls = 0, nstr = 0, nex = 0
  ! occupation classes: electrons in RAS1/2/3, number of substrings,
  ! number of strings and offset in the rank order (0-based); occ_of(n1,n3)
  ! the class of (n1,n3), 0 if not one of this spin
  integer(kind=iwp), allocatable :: occ_n(:,:), occ_cnt(:,:), occ_nstr(:), occ_off(:), occ_of(:,:)
  ! the groups (the blocks of the CI layout): occupation class and string
  ! symmetry; per group: electrons in RAS1/2/3 (cls_n), symmetry (cls_sym),
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

integer(kind=iwp) :: ras_norb = 0, ras_o(3) = 0, ras_pos(3) = 0, ras_minr1 = 0, ras_maxr3 = 0, n_ras = 0, ras_nsym = 1, &
                     ras_lsym = 1, ras_osym(30) = 1
type(ras_spin) :: sa, sb
! block offsets (0-based) of (alpha group, beta group), -1 outside the
! RAS space or of another symmetry; blk_offT: the same blocks transposed
! (beta x alpha)
integer(kind=iwp), allocatable :: blk_off(:,:), blk_offT(:,:)
logical(kind=iwp) :: ras_active = .false., ras_singlet = .false.

public :: faroald_ras_free, faroald_ras_init, n_ras, ras_active, ras_addr, ras_count_dets, ras_from_full, ras_sigma, &
          ras_to_full, ras_traci, ras_trans_one_pdm, ras_trans_two_pdm

contains

!> Set up the string groups, the single-excitation lists and the CI layout
!> for nela/nelb electrons in norb orbitals, o(1:3) of them in RAS1/2/3,
!> with at least minr1 electrons in RAS1 and at most maxr3 in RAS3; the
!> orbitals have the symmetries osym (nsym irreps), the CI space the
!> symmetry lsym (determinants sym(alpha) x sym(beta) = lsym only).
!> spin_flip: all vectors are singlets (nela = nelb), symmetric under the
!> exchange of the alpha and beta strings.
subroutine faroald_ras_init(nela,nelb,norb,o,minr1,maxr3,nsym,osym,lsym,spin_flip)
  integer(kind=iwp), intent(in) :: nela, nelb, norb, o(3), minr1, maxr3, nsym, osym(norb), lsym
  logical(kind=iwp), intent(in) :: spin_flip
  integer(kind=iwp) :: ca, cb, n

  call faroald_ras_free()
  if (sum(o) /= norb) then
    write(u6,*) 'faroald_ras_init: RAS1+RAS2+RAS3 /= norb',o,norb
    call Abend()
  end if
  if (norb > 30) then
    write(u6,*) 'faroald_ras_init: at most 30 orbitals',norb
    call Abend()
  end if
  ras_norb = norb
  ras_o(:) = o(:)
  ras_pos(1) = 0
  ras_pos(2) = o(1)
  ras_pos(3) = o(1)+o(2)
  ras_minr1 = minr1
  ras_maxr3 = maxr3
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
  ras_active = .true.

end subroutine faroald_ras_init

subroutine faroald_ras_free()
  call free_spin(sa)
  call free_spin(sb)
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
  call mma_deallocate(s%occ_of,safe='*')
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

!> A determinant of the classes with electrons n_a(1:3), n_b(1:3) is in the
!> RAS space.
pure function allowed(n_a,n_b) result(ok)
  integer(kind=iwp), intent(in) :: n_a(3), n_b(3)
  logical(kind=iwp) :: ok
  ok = (n_a(1)+n_b(1) >= ras_minr1) .and. (n_a(3)+n_b(3) <= ras_maxr3)
end function allowed

!> The occupation classes of a spin with nel electrons (the other spin has
!> nelo): all (n1,n2,n3) that occur in at least one determinant of the RAS
!> space.
subroutine mk_classes(s,nel,nelo)
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp), intent(in) :: nel, nelo
  integer(kind=iwp) :: n1, n3, m1, m3, k, pass
  integer(kind=iwp) :: n(3), m(3)
  logical(kind=iwp) :: ok

  s%nel = nel
  call mma_allocate(s%occ_of,[0,ras_o(1)],[0,ras_o(3)],label='occ_of')
  s%occ_of(:,:) = 0
  do pass=1,2
    k = 0
    do n3=0,min(ras_o(3),nel)
      do n1=0,min(ras_o(1),nel-n3)
        n = [n1,nel-n1-n3,n3]
        if (n(2) > ras_o(2)) cycle
        ok = .false.
        do m3=0,min(ras_o(3),nelo)
          do m1=0,min(ras_o(1),nelo-m3)
            m = [m1,nelo-m1-m3,m3]
            if (m(2) > ras_o(2)) cycle
            if (allowed(n,m)) ok = .true.
          end do
        end do
        if (.not. ok) cycle
        k = k+1
        if (pass == 2) then
          s%occ_n(:,k) = n(:)
          s%occ_cnt(1,k) = binom_coef(n(1),ras_o(1))
          s%occ_cnt(2,k) = binom_coef(n(2),ras_o(2))
          s%occ_cnt(3,k) = binom_coef(n(3),ras_o(3))
          s%occ_nstr(k) = product(s%occ_cnt(:,k))
          s%occ_of(n1,n3) = k
        end if
      end do
    end do
    if (pass == 1) then
      s%nocc = k
      call mma_allocate(s%occ_n,3,max(1,k),label='occ_n')
      call mma_allocate(s%occ_cnt,3,max(1,k),label='occ_cnt')
      call mma_allocate(s%occ_nstr,max(1,k),label='occ_nstr')
      call mma_allocate(s%occ_off,max(1,k),label='occ_off')
    end if
  end do
  s%nstr = 0
  do k=1,s%nocc
    s%occ_off(k) = s%nstr
    s%nstr = s%nstr+s%occ_nstr(k)
  end do
end subroutine mk_classes

!> The strings of all occupation classes, grouped by class and symmetry:
!> the groups (cls_*), the bit strings in group order and the map perm
!> from the rank address to the group order.
subroutine mk_strings(s)
  type(ras_spin), intent(inout) :: s
  integer(kind=iwp) :: k, g, i1, i2, i3, p1, p2, p3, isym, ng, c, gnew
  integer(kind=iwp), allocatable :: bits(:), sym(:), cnt(:,:), pos(:,:)

  ! the strings in rank order and their symmetries
  call mma_allocate(bits,max(1,s%nstr),label='ras_bits')
  call mma_allocate(sym,max(1,s%nstr),label='ras_sym')
  call mma_allocate(cnt,ras_nsym,max(1,s%nocc),label='ras_cnt')
  cnt(:,:) = 0
  g = 0
  do k=1,s%nocc
    p3 = lex_init(s%occ_n(3,k),ras_o(3))
    do i3=1,s%occ_cnt(3,k)
      p2 = lex_init(s%occ_n(2,k),ras_o(2))
      do i2=1,s%occ_cnt(2,k)
        p1 = lex_init(s%occ_n(1,k),ras_o(1))
        do i1=1,s%occ_cnt(1,k)
          g = g+1
          bits(g) = ior(ior(p1,ishft(p2,ras_pos(2))),ishft(p3,ras_pos(3)))
          sym(g) = str_sym(bits(g))
          cnt(sym(g),k) = cnt(sym(g),k)+1
          if (i1 < s%occ_cnt(1,k)) p1 = lex_next(p1)
        end do
        if (i2 < s%occ_cnt(2,k)) p2 = lex_next(p2)
      end do
      if (i3 < s%occ_cnt(3,k)) p3 = lex_next(p3)
    end do
  end do

  ! the groups: (class, symmetry) with strings
  ng = count(cnt(:,:) > 0)
  s%ncls = ng
  call mma_allocate(s%cls_n,3,max(1,ng),label='cls_n')
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
  integer(kind=iwp) :: b1, b2, b3, k

  g = 0
  b1 = sub_bits(c,ras_pos(1),ras_o(1))
  b2 = sub_bits(c,ras_pos(2),ras_o(2))
  b3 = sub_bits(c,ras_pos(3),ras_o(3))
  if (popcnt(b1)+popcnt(b2)+popcnt(b3) /= s%nel) return
  k = s%occ_of(popcnt(b1),popcnt(b3))
  if (k == 0) return
  g = s%perm(s%occ_off(k)+lexrank(b1)+s%occ_cnt(1,k)*((lexrank(b2)-1)+s%occ_cnt(2,k)*(lexrank(b3)-1)))
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
subroutine same_spin(r,c,off,h,g,vec,sgm)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: h(ras_norb,ras_norb), g(ras_norb,ras_norb,ras_norb,ras_norb), vec(:)
  real(kind=wp), intent(inout) :: sgm(:)
  integer(kind=iwp) :: gi, gj, ci, cj, cr, e, nl, il, nr, oi, oj, t, u, v, x, d, dd, ii, ij
  integer(kind=iwp), allocatable :: lst(:)
  real(kind=wp), allocatable :: f(:)
  logical(kind=iwp), allocatable :: used(:)

  call mma_allocate(f,max(1,c%nstr),label='ras_f')
  call mma_allocate(lst,max(1,c%nstr),label='ras_lst')
  call mma_allocate(used,max(1,c%nstr),label='ras_used')
  f(:) = Zero
  used(:) = .false.
  do gi=1,c%nstr
    ci = c%str_cls(gi)
    nl = 0
    ! one-electron: J = a+_p a_q I
    do e=c%ex_off(gi)+1,c%ex_off(gi+1)
      gj = c%ex_j(e)
      call add(gj,c%ex_sgn(e)*h(c%ex_p(e),c%ex_q(e)))
    end do
    ! two-electron, same spin: J = a+_t a+_v a_x a_u I, t > v, u > x
    if (c%has_dex) then
      do e=c%dex_off(gi)+1,c%dex_off(gi+1)
        call unpack4(c%dex_tuvx(e),t,u,v,x)
        call add(c%dex_j(e),c%dex_sgn(e)*(g(t,u,v,x)-g(t,x,v,u)))
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
!> <Ib|e^b_vx|Jb> vec(Ja,Jb), over the determinants (Ja,Jb) of the RAS
!> space only: for each beta string Jb and its single excitations
!> Jb -> Ib (v,x), the alpha strings Ja of the allowed blocks and their
!> single excitations Ja -> Ia (t,u). Singlets: only the pairs tu >= vx
!> (pair key (u-1)*norb+t, the order of the excitation lists), tu = vx
!> with weight 1/2.
subroutine opposite_spin(g,vec,sgm)
  real(kind=wp), intent(in) :: g(ras_norb,ras_norb,ras_norb,ras_norb), vec(n_ras)
  real(kind=wp), intent(inout) :: sgm(n_ras)
  integer(kind=iwp) :: gb, cbj, cbi, eb, ib, v, x, caj, cai, ja, ga, ea, ia, oj, oi, sgb, kb, ea1
  real(kind=wp) :: c, w

  kb = 0
  do gb=1,sb%nstr
    cbj = sb%str_cls(gb)
    do eb=sb%ex_off(gb)+1,sb%ex_off(gb+1)
      ib = sb%ex_j(eb)
      cbi = sb%str_cls(ib)
      v = sb%ex_p(eb)
      x = sb%ex_q(eb)
      sgb = sb%ex_sgn(eb)
      if (ras_singlet) kb = (x-1)*ras_norb+v
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        oj = blk_off(caj,cbj)+(sb%str_idx(gb)-1)*sa%cls_nstr(caj)
        do ja=1,sa%cls_nstr(caj)
          c = vec(oj+ja)
          if (c == Zero) cycle
          c = sgb*c
          ga = sa%cls_off(caj)+ja
          ea1 = sa%ex_off(ga)+1
          if (ras_singlet) then
            ! the list is sorted by the pair key: skip the pairs tu < vx
            do while (ea1 <= sa%ex_off(ga+1))
              if ((sa%ex_q(ea1)-1)*ras_norb+sa%ex_p(ea1) >= kb) exit
              ea1 = ea1+1
            end do
          end if
          do ea=ea1,sa%ex_off(ga+1)
            ia = sa%ex_j(ea)
            cai = sa%str_cls(ia)
            if (blk_off(cai,cbi) < 0) cycle
            oi = blk_off(cai,cbi)+sa%str_idx(ia)+(sb%str_idx(ib)-1)*sa%cls_nstr(cai)
            w = sa%ex_sgn(ea)*g(sa%ex_p(ea),sa%ex_q(ea),v,x)*c
            if (ras_singlet .and. (ea == ea1) .and. (sa%ex_p(ea) == v) .and. (sa%ex_q(ea) == x)) w = Half*w
            sgm(oi) = sgm(oi)+w
          end do
        end do
      end do
    end do
  end do

end subroutine opposite_spin

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
  integer(kind=iwp) :: ea, eb, ja, ia, jb, ib, caj, cai, cbj, cbi, t, u, v, x
  real(kind=wp) :: c, val

  p2(:,:,:,:) = Zero
  ! same spin
  call col_two_pdm(sa,sb,blk_off,bra,ket,p2)
  if (ras_singlet) then
    ! alpha-alpha = beta-beta; alpha-beta only for the pairs tu >= vx
    p2(:,:,:,:) = Two*p2(:,:,:,:)
    call singlet_ab_two_pdm(bra,ket,p2)
    return
  end if
  call mma_allocate(bT,max(1,n_ras),label='ras_bT')
  call mma_allocate(kT,max(1,n_ras),label='ras_kT')
  call transpose_blocks(bra,bT,.true.)
  call transpose_blocks(ket,kT,.true.)
  call col_two_pdm(sb,sa,blk_offT,bT,kT,p2)
  call mma_deallocate(bT)
  call mma_deallocate(kT)
  ! opposite spin: <a+_ta a_ua a+_vb a_xb>, also as (v,x,t,u); over the
  ! ket determinants (Ja,Jb) of the RAS space
  do jb=1,sb%nstr
    cbj = sb%str_cls(jb)
    do eb=sb%ex_off(jb)+1,sb%ex_off(jb+1)
      ib = sb%ex_j(eb)
      cbi = sb%str_cls(ib)
      v = sb%ex_p(eb)
      x = sb%ex_q(eb)
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        do ja=sa%cls_off(caj)+1,sa%cls_off(caj)+sa%cls_nstr(caj)
          c = ket(blk_off(caj,cbj)+sa%str_idx(ja)+(sb%str_idx(jb)-1)*sa%cls_nstr(caj))
          if (c == Zero) cycle
          c = sb%ex_sgn(eb)*c
          do ea=sa%ex_off(ja)+1,sa%ex_off(ja+1)
            ia = sa%ex_j(ea)
            cai = sa%str_cls(ia)
            if (blk_off(cai,cbi) < 0) cycle
            t = sa%ex_p(ea)
            u = sa%ex_q(ea)
            val = sa%ex_sgn(ea)*c*bra(blk_off(cai,cbi)+sa%str_idx(ia)+(sb%str_idx(ib)-1)*sa%cls_nstr(cai))
            p2(t,u,v,x) = p2(t,u,v,x)+val
            p2(v,x,t,u) = p2(v,x,t,u)+val
          end do
        end do
      end do
    end do
  end do
end subroutine ras_trans_two_pdm

!> Singlets: the alpha-beta part of ras_trans_two_pdm from the pairs
!> tu >= vx: P_ab(tuvx) = P_ab(vxtu), so 2 val goes to (t,u,v,x) and
!> (v,x,t,u) (once for tu = vx).
subroutine singlet_ab_two_pdm(bra,ket,p2)
  real(kind=wp), intent(in) :: bra(n_ras), ket(n_ras)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  integer(kind=iwp) :: ea, eb, ja, ia, jb, ib, caj, cai, cbj, cbi, t, u, v, x, kb
  real(kind=wp) :: c, val

  do jb=1,sb%nstr
    cbj = sb%str_cls(jb)
    do eb=sb%ex_off(jb)+1,sb%ex_off(jb+1)
      ib = sb%ex_j(eb)
      cbi = sb%str_cls(ib)
      v = sb%ex_p(eb)
      x = sb%ex_q(eb)
      kb = (x-1)*ras_norb+v
      do caj=1,sa%ncls
        if (blk_off(caj,cbj) < 0) cycle
        do ja=sa%cls_off(caj)+1,sa%cls_off(caj)+sa%cls_nstr(caj)
          c = ket(blk_off(caj,cbj)+sa%str_idx(ja)+(sb%str_idx(jb)-1)*sa%cls_nstr(caj))
          if (c == Zero) cycle
          c = Two*sb%ex_sgn(eb)*c
          do ea=sa%ex_off(ja+1),sa%ex_off(ja)+1,-1
            t = sa%ex_p(ea)
            u = sa%ex_q(ea)
            if ((u-1)*ras_norb+t < kb) exit
            ia = sa%ex_j(ea)
            cai = sa%str_cls(ia)
            if (blk_off(cai,cbi) < 0) cycle
            val = sa%ex_sgn(ea)*c*bra(blk_off(cai,cbi)+sa%str_idx(ia)+(sb%str_idx(ib)-1)*sa%cls_nstr(cai))
            p2(t,u,v,x) = p2(t,u,v,x)+val
            if ((t /= v) .or. (u /= x)) p2(v,x,t,u) = p2(v,x,t,u)+val
          end do
        end do
      end do
    end do
  end do
end subroutine singlet_ab_two_pdm

!> p2 += <bra|a+_t a+_v a_x a_u|ket> (all index orders) for the column
!> spin c (row spin r); J = a+_t a+_v a_x a_u I with t > v, u > x.
subroutine col_two_pdm(r,c,off,bra,ket,p2)
  type(ras_spin), intent(in) :: r, c
  integer(kind=iwp), intent(in) :: off(:,:)
  real(kind=wp), intent(in) :: bra(:), ket(:)
  real(kind=wp), intent(inout) :: p2(ras_norb,ras_norb,ras_norb,ras_norb)
  integer(kind=iwp) :: gi, e, t, u, v, x, d, dd, gj
  do gi=1,c%nstr
    if (c%has_dex) then
      do e=c%dex_off(gi)+1,c%dex_off(gi+1)
        call unpack4(c%dex_tuvx(e),t,u,v,x)
        call put(gi,c%dex_j(e),c%dex_sgn(e),t,u,v,x)
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
              call put(gi,gj,fase(dd),t,u,v,x)
            end do
          end do
        end do
      end do
    end if
  end do

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

  do p=1,ras_norb
    if (p <= ras_pos(2)) then
      prt(p) = 1
    else if (p <= ras_pos(3)) then
      prt(p) = 2
    else
      prt(p) = 3
    end if
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

!> The number of determinants of the RAS space (all symmetries), without
!> setting anything up; in floating point (no overflow), for the memory
!> estimate of RASSCF.
function ras_count_dets(nela,nelb,o,minr1,maxr3) result(n)
  integer(kind=iwp), intent(in) :: nela, nelb, o(3), minr1, maxr3
  real(kind=wp) :: n
  integer(kind=iwp) :: a1, a3, b1, b3
  n = Zero
  do a3=0,min(o(3),nela)
    do a1=0,min(o(1),nela-a3)
      if (nela-a1-a3 > o(2)) cycle
      do b3=0,min(o(3),nelb)
        do b1=0,min(o(1),nelb-b3)
          if (nelb-b1-b3 > o(2)) cycle
          if ((a1+b1 < minr1) .or. (a3+b3 > maxr3)) cycle
          n = n+bin_r(a1,o(1))*bin_r(nela-a1-a3,o(2))*bin_r(a3,o(3))*bin_r(b1,o(1))*bin_r(nelb-b1-b3,o(2))*bin_r(b3,o(3))
        end do
      end do
    end do
  end do

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
