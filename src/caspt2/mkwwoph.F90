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

subroutine MKWWOPH(IVEC,JVEC,OP0)
! Given the coefficients for two excitation operators, available in
! vectors numbered IVEC and JVEC on file, use the blocks for
! excitation case BJAI, i.e. cases 12 and 13, to express the
! product (Op in IVEC conjugated)(Op in JVEC) as a zero-body
! operator, i.e. a scalar factor, in the CASSCF space.
! Formula used:
!  (W1(ij,ab) conj)*(W2(kl,cd)) = dik*djl*dac*dbd
! the same for both cases 10 and 11.

use EQSOLV, only: MODVEC
use caspt2_module, only: NASUP, NINDEP, NISUP, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC
real(kind=wp), intent(inout) :: OP0
integer(kind=iwp) :: ICASE, IIEND, IISTA, ISCT, ISYM, MDVEC, NAS, NCOL, NIS, NSCT
real(kind=wp), allocatable :: W1(:), W2(:)
real(kind=wp), external :: DDot_

! Loop over cases
do ICASE=12,13
  ! Loop over symmetry ISYM
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NINDEP(ISYM,ICASE) == 0) cycle
    ! Allocate space for one section of excitation amplitudes:
    MDVEC = MODVEC(ISYM,ICASE)
    call mma_allocate(W1,NAS*MDVEC,Label='W1')
    call mma_allocate(W2,NAS*MDVEC,Label='W2')
    ! Sectioning loop added:
    ISCT = 0
    do IISTA=1,NIS,MDVEC
      ISCT = ISCT+1
      IIEND = min(IISTA-1+MDVEC,NIS)
      NCOL = 1+IIEND-IISTA
      NSCT = NAS*NCOL
      call RDSCTC(ISCT,ISYM,ICASE,IVEC,W1,NAS*MDVEC)
      call RDSCTC(ISCT,ISYM,ICASE,JVEC,W2,NAS*MDVEC)
    ! Pick up a symmetry block of W1 and W2
      OP0 = OP0+DDOT_(NSCT,W1,1,W2,1)
    end do
    call mma_deallocate(W1)
    call mma_deallocate(W2)
  end do
  ! End of loop over cases.
end do

end subroutine MKWWOPH
