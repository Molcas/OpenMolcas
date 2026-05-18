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

subroutine MKWWOPE(IVEC,JVEC,OP0,OP1)
! Presently symmetry blocking is disregarded.

! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation case VJAI, i.e. cases 6 and 7, to express the
! product (Op in IVEC conjugated)(Op in JVEC) as a one-body
! operator on the CASSCF space.
! Formula used:
!  (W1(t,aij) conj)*(W2(x,bkl)) = dik*djl*dab*(2*dtx - Etx)
! the same for both cases 6 and 7.

use definitions, only: iwp, wp
use constants, only: Zero, One, Two
use EQSOLV, only: MODVEC
use stdalloc, only: mma_allocate, mma_deallocate
use caspt2_module, only: NASHT, NSYM, NASUP, NISUP, NINDEP, NAES

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC
real(kind=wp), intent(inout) :: OP0, OP1(NASHT,NASHT)
real(kind=wp), allocatable :: W1(:), W2(:), WPROD(:)
integer(kind=iwp) ICASE, ISYM, NAS, NIS, MDVEC, NWPROD, ISCT, IISTA, IIEND, NCOL, IT, IW1, ITABS, IX, IW2, IWPROD, IXABS
real(kind=wp) W_PROD

! Loop over cases
do ICASE=6,7
  ! Loop over symmetry ISYM
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NINDEP(ISYM,ICASE) == 0) cycle
    ! Allocate space for one section of excitation amplitudes:
    MDVEC = MODVEC(ISYM,ICASE)
    call mma_allocate(W1,NAS*MDVEC,Label='W1')
    call mma_allocate(W2,NAS*MDVEC,Label='W2')
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
      ! Multiply WProd = (W1)*(W2 transpose)
      call DGEMM_('N','T',NAS,NAS,NCOL,One,W1,NAS,W2,NAS,One,WPROD,NAS)
    end do
    ! Deallocate space for this block of excitation amplitudes:
    call mma_deallocate(W1)
    call mma_deallocate(W2)

    ! Loop over (T)
    do IT=1,NAS
      IW1 = IT
      ITABS = IT+NAES(ISYM)
      ! Loop over (X)
      do IX=1,NAS
        IW2 = IX
        IXABS = IX+NAES(ISYM)
        IWPROD = IW1+NAS*(IW2-1)
        W_PROD = WPROD(IWPROD)
        OP1(ITABS,IXABS) = OP1(ITABS,IXABS)-W_PROD
        if (ITABS == IXABS) OP0 = OP0+Two*W_PROD
      end do
    end do
    ! Deallocate matrix product
    call mma_deallocate(WPROD)
  end do
  ! End of loop over cases.
end do

end subroutine MKWWOPE
