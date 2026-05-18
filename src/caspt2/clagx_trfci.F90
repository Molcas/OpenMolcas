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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CLagX_TrfCI(NCONF,CI)

use caspt2_global, only: TAT, TORB
use caspt2_module, only: NISH, NRAS1, NRAS2, NRAS3, NSSH, NSYM, STSYM
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCONF
real(kind=wp), intent(inout) :: CI(NCONF)
integer(kind=iwp) :: I, IJ, IOFF1, IOFF2, ISYM, J, JI, NI, NISH_SAVE(8), NR1, NR2, NR3, NS

TAT(:) = Zero

IOFF1 = 0
IOFF2 = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NR1 = NRAS1(ISYM)
  NR2 = NRAS2(ISYM)
  NR3 = NRAS3(ISYM)
  NS = NSSH(ISYM)
  ! Skip inactive transformation matrix:
  IOFF1 = IOFF1+NI**2
  ! Copy RAS1 transformation matrix transposed to TAT:
  do I=1,NR1
    do J=1,NR1
      IJ = I+NR1*(J-1)
      JI = J+NR1*(I-1)
      TAT(IOFF2+JI) = TORB(IOFF1+IJ)
    end do
  end do
  IOFF1 = IOFF1+NR1**2
  IOFF2 = IOFF2+NR1**2
  ! Copy RAS2 transformation matrix transposed to TAT:
  do I=1,NR2
    do J=1,NR2
      IJ = I+NR2*(J-1)
      JI = J+NR2*(I-1)
      TAT(IOFF2+JI) = TORB(IOFF1+IJ)
    end do
  end do
  IOFF1 = IOFF1+NR2**2
  IOFF2 = IOFF2+NR2**2
  ! Copy RAS2 transformation matrix transposed to TAT:
  do I=1,NR3
    do J=1,NR3
      IJ = I+NR3*(J-1)
      JI = J+NR3*(I-1)
      TAT(IOFF2+JI) = TORB(IOFF1+IJ)
    end do
  end do
  IOFF1 = IOFF1+NR3**2
  IOFF2 = IOFF2+NR3**2
  ! Skip virtual transformation matrix:
  IOFF1 = IOFF1+NS**2
end do
! Transform SGM to use original MO:

! In difference to TORB TAT only includes the active orbitals. We trick
! MkTraCI to work with this by temporarily setting the number of inactive
! orbitals to zero.
NISH_SAVE(:) = NISH(:)
NISH(:) = 0
call mkTraCI(size(TAT),TAT,STSYM,nConf,CI)
NISH(:) = NISH_SAVE(:)

end subroutine CLagX_TrfCI
