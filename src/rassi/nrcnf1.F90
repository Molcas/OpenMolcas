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

subroutine NRCNF1(MAXEL,NORB,NGAS,NGASLIM,NGASORB,NCNF1,MXTMP,NCNF2)
! Returns the array NCNF1, which contains the number of
! configurations with the following criteria:
!    Orbital indices range from 1..NORB
!    GAS restrictions described by NGASLIM, NGASORB
!    Total number of electrons is at most MAXEL
!    NCLS closed-shell and NOPN open-shell orbitals
!    Symmetry label LSYM
! for all possible values of NCLS,NOPN and LSYM, stored as
!          NCNF1(LSYM,IPOS)
! with IPOS=(NOCC*(NOCC+1))/2+NOPN+1, NOCC=NCLS+NOPN,
! provided that 0<=NCLS, 0<=NOPN, and NOCC<=MIN(NORB,MAXEL).
! MXTMP=Max nr of orbitals in one GAS partition.
! Prerequisite: The orbital symmetry labels stored in ISM.
!               The GAS restriction arrays
! Method: Induction over GAS partitions.

use Symmetry_Info, only: MUL, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: MaxEl, NORB, NGAS, NGASLIM(2,NGAS), NGASORB(nIrrep,NGAS), MXTMP
integer(kind=iwp), intent(out) :: NCNF1(nIrrep,((MAXEL+1)*(MAXEL+2))/2), NCNF2(nIrrep,((MXTMP+1)*(MXTMP+2))/2)
integer(kind=iwp) :: IGAS, II, IPOS, IPOSNW, IPOSOLD, ISYM, ISYMNW, ISYMOLD, MAXOCC, MXOCCOLD, NCLS, NCLSNW, NCLSOLD, NELMN, &
                     NELMX, NEW, NG, NO, NOCC, NOCCMX, NOCCNW, NOCCOLD, NOPN, NOPNNW, NOPNOLD, NX, NY
integer(kind=iwp), allocatable :: ISM(:)

MAXOCC = min(MAXEL,NORB)
! Initialize:
NCNF1(:,:) = 0
NCNF1(1,1) = 1
call mma_allocate(ISM,NORB,Label='ISM')
! Max nr of occupied orbitals so far:
NOCCMX = 0
do IGAS=1,NGAS
  MXOCCOLD = NOCCMX
  ! Nr of orbitals in this partition
  II = 0
  do ISYM=1,nIrrep
    NG = NGASORB(ISYM,IGAS)
    ISM(II+1:II+NG) = ISYM
    II = II+NG
  end do
  NO = sum(NGASORB(1:nIrrep,IGAS))
  NELMN = max(0,NGASLIM(1,IGAS))
  NELMX = min(2*NO,NGASLIM(2,IGAS))
  call NRCNF2(NO,ISM,NCNF2)
  do NOCCNW=min(MAXOCC,MXOCCOLD+NELMX),0,-1
    do NOPNNW=0,NOCCNW
      NCLSNW = NOCCNW-NOPNNW
      IPOSNW = (NOCCNW*(NOCCNW+1))/2+NOPNNW+1
      do ISYMNW=1,nIrrep
        NEW = 0
        do NOCC=NELMN/2,min(NELMX,NO)
          do NOPN=max(0,2*NOCC-NELMX),min(2*NOCC-NELMN,NOCC,NELMX)
            NCLS = NOCC-NOPN
            IPOS = (NOCC*(NOCC+1))/2+NOPN+1
            do ISYM=1,nIrrep
              NY = NCNF2(ISYM,IPOS)
              if (NY == 0) cycle
              NCLSOLD = NCLSNW-NCLS
              if (NCLSOLD < 0) cycle
              NOPNOLD = NOPNNW-NOPN
              if (NOPNOLD < 0) cycle
              NOCCOLD = NCLSOLD+NOPNOLD
              if (NOCCOLD > MXOCCOLD) cycle
              ISYMOLD = MUL(ISYM,ISYMNW)
              IPOSOLD = (NOCCOLD*(NOCCOLD+1))/2+NOPNOLD+1
              NX = NCNF1(ISYMOLD,IPOSOLD)
              if (NX == 0) cycle
              NEW = NEW+NCNF1(ISYMOLD,IPOSOLD)*NCNF2(ISYM,IPOS)
              NOCCMX = max(NOCCNW,NOCCMX)
            end do
          end do
        end do
        NCNF1(ISYMNW,IPOSNW) = NEW
      end do
    end do
  end do

end do
call mma_deallocate(ISM)

end subroutine NRCNF1
