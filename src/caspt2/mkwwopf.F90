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

subroutine MKWWOPF(IVEC,JVEC,NOP2,OP2)
! Presently symmetry blocking is disregarded, but index pair
! permutation symmetry is used.
! NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)

! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation cases BVAT(+) and BVAT(-), i.e. cases 8 and 9, to
! construct the zero-, one-, and two-body
! expansions of the product (Op in IVEC conjugated)(Op in JVEC)
! as operating on the CASSCF space.
! Formulae used:
! For the F+ case (i.e. case 8)
! W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy + 2 Etyux)
! For the F- case (i.e. case 9)
! W1(tu,ab)(conj)*W2(xy,cd) = (dac*dbd)*(2 Etxuy - 2 Etyux)

use SUPERINDEX, only: MTGEU, MTGTU
use EQSOLV, only: MODVEC
use caspt2_module, only: NASHT, NASUP, NINDEP, NISUP, NSYM, NTGEUES, NTGTUES
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NOP2
real(kind=wp), intent(inout) :: OP2(NOP2)
integer(kind=iwp) :: ICASE, IIEND, IISTA, ISCT, ISYM, ITABS, ITU, ITUABS, ITX, ITY, IUABS, IUX, IUY, IW1, IW2, IWPROD, IXABS, IXY, &
                     IXYABS, IYABS, JTXUY, JTYUX, MDVEC, NAS, NCOL, NIS, NWPROD
real(kind=wp) :: W_PROD
real(kind=wp), allocatable :: WPROD(:)
real(kind=wp), allocatable, target :: W1(:), W2_H(:)
real(kind=wp), pointer :: W2(:)

! FIRST THE F+ i.e. BVAT+ i.e. CASE 8 -----------------------------
ICASE = 8
! Loop over symmetry ISYM
do ISYM=1,NSYM
  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  if (NINDEP(ISYM,ICASE) == 0) cycle
  ! Allocate space for one section of excitation amplitudes:
  ! Pick up a symmetry block of W1 and W2
  MDVEC = MODVEC(ISYM,ICASE)
  call mma_allocate(W1,NAS*MDVEC,Label='W1')
  if (JVEC == IVEC) then
    W2 => W1
  else
    call mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
    W2 => W2_H
  end if
  NWPROD = NAS**2
  ! Allocate space for the contraction:
  call mma_allocate(WPROD,NWPROD,Label='WPROD')
  WPROD(:) = Zero
  ! Sectioning loop added:
  ISCT = 0
  do IISTA=1,NIS,MDVEC
    ISCT = ISCT+1
    IIEND = min(IISTA-1+MDVEC,NIS)
    NCOL = 1+IIEND-IISTA
    call RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
    call RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
    ! Multiply WProd = (W1 sect )*(W2 sect transpose)
    call DGEMM_('N','T',NAS,NAS,NCOL,One,W1,NAS,W2,NAS,One,WPROD,NAS)
  end do
  ! Deallocate W1 and W2
  call mma_deallocate(W1)
  if (JVEC /= IVEC) call mma_deallocate(W2_H)
  nullify(W2)

  ! Loop over (TU)
  do ITU=1,NAS
    IW1 = ITU
    ITUABS = ITU+NTGEUES(ISYM)
    ITABS = MTGEU(1,ITUABS)
    IUABS = MTGEU(2,ITUABS)
    ! Loop over (XY)
    do IXY=1,NAS
      IW2 = IXY
      IXYABS = IXY+NTGEUES(ISYM)
      IXABS = MTGEU(1,IXYABS)
      IYABS = MTGEU(2,IXYABS)
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      ITY = ITABS+NASHT*(IYABS-1)
      IUX = IUABS+NASHT*(IXABS-1)
      IWPROD = IW1+NAS*(IW2-1)
      W_PROD = WPROD(IWPROD)
      ! Remember: C For the F+ case (i.e. case 8)
      ! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy + 2 Etyux)
      ! Contrib to 2-particle operator, from 2 Etxuy:
      if (ITX >= IUY) then
        JTXUY = (ITX*(ITX-1))/2+IUY
      else
        JTXUY = (IUY*(IUY-1))/2+ITX
      end if
      OP2(JTXUY) = OP2(JTXUY)+Two*W_PROD
      ! Contrib to 2-particle operator, from 2 Etyux:
      if (ITY >= IUX) then
        JTYUX = (ITY*(ITY-1))/2+IUX
      else
        JTYUX = (IUX*(IUX-1))/2+ITY
      end if
      OP2(JTYUX) = OP2(JTYUX)+Two*W_PROD
    end do
  end do
  ! Deallocate matrix product:
  call mma_deallocate(WPROD)
end do
! THEN THE F- i.e. BVAT- i.e. CASE 9 -----------------------------
ICASE = 9
! Loop over symmetry ISYM
do ISYM=1,NSYM
  NAS = NASUP(ISYM,ICASE)
  NIS = NISUP(ISYM,ICASE)
  if (NINDEP(ISYM,ICASE) == 0) cycle
  ! Allocate space for one section of excitation amplitudes:
  ! Pick up a symmetry block of W1 and W2
  MDVEC = MODVEC(ISYM,ICASE)
  call mma_allocate(W1,NAS*MDVEC,Label='W1')
  if (JVEC == IVEC) then
    W2 => W1
  else
    call mma_allocate(W2_H,NAS*MDVEC,Label='W2_H')
    W2 => W2_H
  end if
  NWPROD = NAS**2
  ! Allocate space for the contraction:
  call mma_allocate(WPROD,NWPROD,Label='WPROD')
  WPROD(:) = Zero
  ! Sectioning loop added:
  ISCT = 0
  do IISTA=1,NIS,MDVEC
    ISCT = ISCT+1
    IIEND = min(IISTA-1+MDVEC,NIS)
    NCOL = 1+IIEND-IISTA
    call RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
    call RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
    ! Multiply WProd = (W1 sect )*(W2 sect transpose)
    call DGEMM_('N','T',NAS,NAS,NCOL,One,W1,NAS,W2,NAS,One,WPROD,NAS)
  end do
  ! Deallocate W1 and W2
  call mma_deallocate(W1)
  if (JVEC /= IVEC) call mma_deallocate(W2_H)
  nullify(W2)

  ! Loop over (TU)
  do ITU=1,NAS
    IW1 = ITU
    ITUABS = ITU+NTGTUES(ISYM)
    ITABS = MTGTU(1,ITUABS)
    IUABS = MTGTU(2,ITUABS)
    ! Loop over (XY)
    do IXY=1,NAS
      IW2 = IXY
      IXYABS = IXY+NTGTUES(ISYM)
      IXABS = MTGTU(1,IXYABS)
      IYABS = MTGTU(2,IXYABS)
      ITX = ITABS+NASHT*(IXABS-1)
      IUY = IUABS+NASHT*(IYABS-1)
      ITY = ITABS+NASHT*(IYABS-1)
      IUX = IUABS+NASHT*(IXABS-1)
      IWPROD = IW1+NAS*(IW2-1)
      W_PROD = WPROD(IWPROD)
      ! Remember: C For the F- case (i.e. case 9)
      ! W1(tu,ij)(conj)*W2(xy,kl) = (dik*djl)*(2 Etxuy - 2 Etyux)
      ! Contrib to 2-particle operator, from 2 Etxuy:
      if (ITX >= IUY) then
        JTXUY = (ITX*(ITX-1))/2+IUY
      else
        JTXUY = (IUY*(IUY-1))/2+ITX
      end if
      OP2(JTXUY) = OP2(JTXUY)+Two*W_PROD
      ! Contrib to 2-particle operator, from -2 Etyux:
      if (ITY >= IUX) then
        JTYUX = (ITY*(ITY-1))/2+IUX
      else
        JTYUX = (IUX*(IUX-1))/2+ITY
      end if
      OP2(JTYUX) = OP2(JTYUX)-Two*W_PROD
    end do
  end do
  ! Deallocate matrix product:
  call mma_deallocate(WPROD)
end do

end subroutine MKWWOPF
