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

subroutine NRCNF2(NORB,ISM,NCNF2)
! Returns the array NCNF2, which contains the number of
! (sub-)configurations with NCLS closed-shell and NOPN open-shell
! orbitals and having symmetry label LSYM, stored as
!          NCNF2(LSYM,IPOS)
! with IPOS=(NOCC*(NOCC+1))/2+NOPN+1, NOCC=NCLS+NOPN,
! provided that 0<=NCLS, 0<=NOPN, and NOCC<=NORB.
! Prerequisite: The orbital symmetry labels stored in ISM.
! Method: Induction

use Symmetry_Info, only: MUL, nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NORB, ISM(NORB)
integer(kind=iwp), intent(out) :: NCNF2(nIrrep,(NORB+1)*(NORB+2)/2)
integer(kind=iwp) :: IPOS1, IPOS2, IPOS3, ISYM, JSYM, L, NCLS, NEW, NOCC, NOPN

NCNF2(:,:) = 0
NCNF2(1,1) = 1
do L=1,NORB
  do NOCC=L,1,-1
    do NOPN=0,NOCC
      NCLS = NOCC-NOPN
      IPOS1 = (NOCC*(NOCC+1))/2+NOPN+1
      IPOS2 = IPOS1-NOCC
      IPOS3 = IPOS2-1
      do ISYM=1,nIrrep
        NEW = NCNF2(ISYM,IPOS1)
        if (NCLS > 0) NEW = NEW+NCNF2(ISYM,IPOS2)
        JSYM = MUL(ISM(L),ISYM)
        if (NOPN > 0) NEW = NEW+NCNF2(JSYM,IPOS3)
        NCNF2(ISYM,IPOS1) = NEW
      end do
    end do
  end do
end do

end subroutine NRCNF2
