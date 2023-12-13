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

public :: citrans_csf2sd, citrans_sd2csf, citrans_sort, ncsf_group, ndet_group, ndo_max, ndo_min, ndoc_group, nsoc_group, &
          spintable_create, spintabs, spintabs_allocate, spintabs_free

! Private extensions to mma interfaces

interface cptr2loff
  module procedure spt_cptr2loff
end interface
interface mma_allocate
  module procedure spt_mma_allo_1D, spt_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure spt_mma_free_1D
end interface

contains

subroutine citrans_sort(mode,ciold,cinew)
! sort CSFs such that they are grouped by configurations
! and fill in the transformation matrices to determinants

  use second_quantization, only: lexrank
  use faroald, only: my_norb

  ! 'C' for configuration order, 'O' for original order
  character, intent(in) :: mode
  real(kind=wp), intent(in) :: ciold(:)
  real(kind=wp), intent(out) :: cinew(:)
  ! array with coupling coefficients
  !real(kind=wp), intent(out) :: coef(:)
  integer(kind=iwp) :: doub, icsf, ido, idown, idwn, ioff_csf, iorb, iphase, iso, iup, mv, ncsf, ndoc, ndown, nsoc, rankdo, &
                       rankso, sing
  integer(kind=iwp), parameter :: maxorb = 32, maxdown = 16
  real(kind=wp) :: wtab(0:maxorb,maxdown)
  integer(kind=iwp), allocatable :: csf_offset(:), downvector(:), stepvector(:)

  ! Compute offsets for addressing into the reordering and coefficient arrays.
  call mma_allocate(csf_offset,[ndo_min,ndo_max],label='csf_offset')

# ifdef _DEBUGPRINT_
  write(u6,'(5(1x,a4))') 'ido','ndoc','nsoc','ndet','ncsf'
# endif

  ncsf = 0
  do ido=ndo_min,ndo_max
#   ifdef _DEBUGPRINT_
    write(u6,'(5(1x,i4))') ido,ndoc_group(ido),nsoc_group(ido),ndet_group(ido),ncsf_group(ido)
#   endif
    ndoc = ndoc_group(ido)
    nsoc = nsoc_group(ido)
    ! use nconf/ncsf/ndet as offsets
    csf_offset(ido) = ncsf
    ncsf = ncsf+ndoc*nsoc*ncsf_group(ido)
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
  call mma_allocate(downvector,my_norb,label='downvector')
  ! initialize variables that track the stepvector
  mv = 1
  idwn = 1
  iup = 1
  do icsf=1,ncsf
    ! obtain the stepvector
    call stepvector_next(mv,idwn,iup,stepvector)

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
    ndoc = ndoc_group(ido)
    nsoc = nsoc_group(ido)
    ioff_csf = csf_offset(ido)+(nsoc*(rankdo-1)+rankso-1)*ncsf_group(ido)

    ! CSF offset within this configuration, use per-ake's magical wtab
    ! to get the rank of a ud string within its set.
    do idown=1,ndown
      ioff_csf = ioff_csf+nint(wtab(downvector(idown)-2*idown,idown))
    end do

    ! reorder the values
    if (mode == 'C') then
      cinew(ioff_csf+1) = iphase*ciold(icsf)
    else
      cinew(icsf) = iphase*ciold(ioff_csf+1)
    end if

  end do

  call mma_deallocate(csf_offset)
  call mma_deallocate(stepvector)
  call mma_deallocate(downvector)

end subroutine citrans_sort

subroutine citrans_csf2sd(ci,det)

  use second_quantization, only: lex_init, lex_next, lexrank
  use faroald, only: my_nel, my_norb, ndeta, ndetb, nela

  real(kind=wp), intent(in) :: ci(:)
  real(kind=wp), intent(out) :: det(ndeta,ndetb)
  integer(kind=iwp) :: ido, iso, isoa, idoc, isoc, iconf, idet, ndoc, nsoc, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
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
    ndoc = ndoc_group(ido)
    nsoc = nsoc_group(ido)
    nconf = ndoc*nsoc
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('N','N',ndet,nconf,ncsf,One,spintabs(ido)%coef,ndet,ci(ioff_csf1:ioff_csf2),ncsf,Zero,tmp,ndet)

    ! Store the determinant coefficients with the right phase factor in
    ! the correct place in the determinant matrix. The loops runs over
    ! the determinants in the order they are stored in the tmp array.

    iso = my_nel-2*ido
    isoa = nela-ido

    iconf = 0
    doub = lex_init(ido,my_norb)
    do idoc=1,ndoc
      sing = lex_init(iso,my_norb-ido)
      do isoc=1,nsoc
        iconf = iconf+1
        alfa = lex_init(isoa,iso)
        do idet=1,ndet
          beta = ibits(not(alfa),0,iso)
          phase = ds2ab(doub,sing,alfa,beta,deta,detb)
          det(lexrank(deta),lexrank(detb)) = phase*tmp(idet,iconf)
          alfa = lex_next(alfa)
        end do
        sing = lex_next(sing)
      end do
      doub = lex_next(doub)
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
  integer(kind=iwp) :: ido, iso, isoa, idoc, isoc, iconf, idet, ndoc, nsoc, ndet, ncsf, nconf, ioff_csf1, ioff_csf2, &
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
    ndoc = ndoc_group(ido)
    nsoc = nsoc_group(ido)
    nconf = ndoc*nsoc
    ndet = ndet_group(ido)
    ncsf = ncsf_group(ido)
    ioff_csf2 = ioff_csf1+ncsf*nconf-1

    call mma_allocate(tmp,ndet,nconf,label='tmp')

    ! Collect the determinant coefficients with the right phase factor
    ! from the correct place in the determinant matrix. The loops runs
    ! over the determinants in the order they should end up in tmp.

    iso = my_nel-2*ido
    isoa = nela-ido

    iconf = 0
    doub = lex_init(ido,my_norb)
    do idoc=1,ndoc
      sing = lex_init(iso,my_norb-ido)
      do isoc=1,nsoc
        iconf = iconf+1
        alfa = lex_init(isoa,iso)
        do idet=1,ndet
          beta = ibits(not(alfa),0,iso)
          phase = ds2ab(doub,sing,alfa,beta,deta,detb)
          tmp(idet,iconf) = phase*det(lexrank(deta),lexrank(detb))
          alfa = lex_next(alfa)
        end do
        sing = lex_next(sing)
      end do
      doub = lex_next(doub)
    end do

    ! Compute the determinant coefficients from the CSF coefficients.
    call dgemm_('T','N',ncsf,nconf,ndet,One,spintabs(ido)%coef,ndet,tmp,ndet,Zero,ci(ioff_csf1:ioff_csf2),ncsf)

    call mma_deallocate(tmp)
    ioff_csf1 = ioff_csf2+1
  end do

  call mma_deallocate(stepvector)

end subroutine citrans_sd2csf

subroutine spintabs_allocate()

# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
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

  call mma_deallocate(spintabs)

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
! (see src/mma_util/stdalloc.f)

! Define spt_cptr2loff, spt_mma_allo_1D, spt_mma_allo_1D_lim, spt_mma_free_1D
#define _TYPE_ type(spintable)
#  define _FUNC_NAME_ spt_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ spt_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'spt_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module citrans
