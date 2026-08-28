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

!#define _DEBUGPRINT_
subroutine ZOOS(ISM,IBLTP,MAXSYM,IOCOC,NSSOA,NSSOB,NOCTPA,NOCTPB,IDC,IOOS,NOOS,NCOMB,IXPND)
! Generate offsets for CI vector for RAS CI expansion of given symmetry
! Combination type is defined by IDC
! Total number of combinations NCOMB is also obtained
!
! Symmetry is defined through ISM
!
! ICBLTP gives typo of symmetry block :
! = 0 : symmetry block is not included
! = 1 : symmetry block is included, all OO types
! = 2 : symmetry block is included, lower OO types
!
! If IXPND /= 0, the diagonal blocks are always choosen expanded
!
! ========
!  Output
! ========
!
! IOOS(IOCA,IOCB,ISYM) : Start of block with alpha strings of
!                        symmetry ISYM and type IOCA, and
!                        betastrings of type IOCB
! NOOS(IOCA,IOCB,ISYM) : Number of combinations
! The ordering used for the CI vector is
!
!    SYMMETRY  FOR ALPHA STRINGS..(GIVES SYMMETRY OF BETA STRING)
!         OCCUPATION TYPE  FOR ALPHA STRING
!            OCCUPATION TYPE FOR    BETA STRING
!                BETA STRING (COLUMN INDEX)
!                ALPHA STRINGS (ROW INDEX)
!    END OF LOOPS

use Symmetry_Info, only: Mul
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: ISM, IBLTP(*), MAXSYM, NOCTPA, NOCTPB, IOCOC(NOCTPA,NOCTPB), NSSOA(MAXSYM,NOCTPA), &
                                 NSSOB(MAXSYM,NOCTPB), IDC, IXPND
integer(kind=iwp), intent(out) :: IOOS(NOCTPA,NOCTPB,MAXSYM), NOOS(NOCTPA,NOCTPB,MAXSYM), NCOMB
integer(kind=iwp) :: IAOCC, IASYM, IBOCC, IBSYM, IREST1, MXBOCC

IOOS(:,:,:) = 0
NOOS(:,:,:) = 0
!ICBLTP(1:MAXSYM) = 0
NCOMB = 0
do IASYM=1,MAXSYM

  IBSYM = Mul(IASYM,ISM)
  if (IBSYM == 0) cycle
  ! Allowed combination symmetry block?
  if ((IDC /= 1) .and. (IBLTP(IASYM) == 0)) cycle
  ! Allowed occupation combinations
  do IAOCC=1,NOCTPA
    if (IBLTP(IASYM) == 1) then
      MXBOCC = NOCTPB
      IREST1 = 0
    else
      MXBOCC = IAOCC
      IREST1 = 1
    end if
    do IBOCC=1,MXBOCC
      ! Is this block allowed
      if (IOCOC(IAOCC,IBOCC) == 1) then
        IOOS(IAOCC,IBOCC,IASYM) = NCOMB+1
        if ((IXPND == 0) .and. (IREST1 == 1) .and. (IAOCC == IBOCC)) then
          NCOMB = NCOMB+(NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
          NOOS(IAOCC,IBOCC,IASYM) = (NSSOA(IASYM,IAOCC)+1)*NSSOB(IBSYM,IBOCC)/2
        else
          NCOMB = NCOMB+NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
          NOOS(IAOCC,IBOCC,IASYM) = NSSOA(IASYM,IAOCC)*NSSOB(IBSYM,IBOCC)
        end if
      end if
      !write(u6,*) ' NOOS(IA,IB,ISM) ',NOOS(IAOCC,IBOCC,IASYM)
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' =============='
write(u6,*) ' ZOOS reporting'
write(u6,*) ' =============='
write(u6,*)
write(u6,*) ' NSSOA, NSSOB (input)'
call IWRTMA(NSSOA,MAXSYM,NOCTPA,MAXSYM,NOCTPA)
call IWRTMA(NSSOB,MAXSYM,NOCTPB,MAXSYM,NOCTPB)
write(u6,*)
write(u6,*) ' Number of combinations obtained ',NCOMB
write(u6,*) ' Offsets for combination OOS blocks'
do IASYM=1,MAXSYM
  write(u6,'(A,I2)') '  Symmetry ',IASYM
  call IWRTMA(IOOS(:,:,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end do
write(u6,*) ' Number of  combinations per OOS blocks'
do IASYM=1,MAXSYM
  write(u6,'(A,I2)') '  Symmetry ',IASYM
  call IWRTMA(NOOS(:,:,IASYM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
end do
#endif

return

end subroutine ZOOS
