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

subroutine MKSGNUM(STSYM,SGS,CIS,EXS)
! PURPOSE: FOR ALL UPPER AND LOWER WALKS
!          COMPUTE THE DIRECT ARC WEIGHT SUM AND THE
!          REVERSE ARC WEIGHT SUM, RESPECTIVELY.
!          STORE THE DATA IN THE TABLES USGN AND LSGN

use Symmetry_Info, only: Mul
use gugx, only: CIStruct, EXStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: STSYM
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IC, ICODE, ICONF, IDAWSUM, ILOFF, ILW, IPOS, IRAWSUM, ISTEP, ISYM, IUOFF, IUW, JPOS, JSYM, LEV, LV, MIDV, &
                     NLW, NUW
integer(kind=iwp), allocatable :: ISTEPVEC(:)

call mma_allocate(EXS%USGN,SGS%MxUp,CIS%nMidV,Label='EXS%USGN')
call mma_allocate(EXS%LSGN,SGS%MxDwn,CIS%nMidV,Label='EXS%LSGN')
call mma_allocate(ISTEPVEC,SGS%nLev,Label='ISTEPVEC')

! INITIALIZE NUMBERING TABLES

EXS%USGN(:,:) = 0
EXS%LSGN(:,:) = 0

! MAIN LOOP RUNS OVER MIDVERTICES AND SYMMETRIES

ICONF = 0
do MIDV=1,CIS%nMidV
  do ISYM=1,SGS%nSym
    IUOFF = 1+CIS%IOW(1,ISYM,MIDV)
    NUW = CIS%NOW(1,ISYM,MIDV)
    JSYM = Mul(ISYM,STSYM)
    ILOFF = 1+CIS%IOW(2,JSYM,MIDV)
    NLW = CIS%NOW(2,JSYM,MIDV)
    if ((NUW == 0) .or. (NLW == 0)) cycle

    ! LOOP OVER ALL UPPER WALKS

    do IUW=1,NUW
      IPOS = IUOFF+CIS%nIpWlk*(IUW-1)
      ! UNPACK THE UPPER WALK STEP VECTOR
      ICODE = CIS%iCase(IPOS)
      JPOS = 0
      do LEV=SGS%MidLev+1,SGS%nLev
        JPOS = JPOS+1
        if (JPOS == 16) then
          JPOS = 1
          IPOS = IPOS+1
          ICODE = CIS%iCase(IPOS)
        end if
        ISTEP = mod(ICODE,4)
        ISTEPVEC(LEV) = ISTEP
        ICODE = ICODE/4
      end do
      ! GET REVERSE ARC WEIGHT FOR UPPER WALK
      IRAWSUM = 1
      LV = 1
      do LEV=SGS%nLev,SGS%MidLev+1,-1
        IC = ISTEPVEC(LEV)
        LV = SGS%Down(LV,IC)
        IRAWSUM = IRAWSUM+SGS%RAW(LV,IC)
      end do
      EXS%USGN(IRAWSUM,MIDV) = IUW
    end do

    ! LOOP OVER ALL LOWER WALKS

    do ILW=1,NLW
      IPOS = ILOFF+CIS%nIpWlk*(ILW-1)
      ! UNPACK WALK STEP VECTOR
      ICODE = CIS%iCase(IPOS)
      JPOS = 0
      do LEV=1,SGS%MidLev
        JPOS = JPOS+1
        if (JPOS == 16) then
          JPOS = 1
          IPOS = IPOS+1
          ICODE = CIS%iCase(IPOS)
        end if
        ISTEP = mod(ICODE,4)
        ISTEPVEC(LEV) = ISTEP
        ICODE = ICODE/4
      end do
      ! GET DIRECT ARC WEIGHT FOR THE LOWER WALK
      IDAWSUM = 1
      LV = SGS%nVert
      do LEV=1,SGS%MidLev
        IC = ISTEPVEC(LEV)
        LV = SGS%Up(LV,IC)
        IDAWSUM = IDAWSUM+SGS%DAW(LV,IC)
      end do
      EXS%LSGN(IDAWSUM,MIDV) = ICONF
      ICONF = ICONF+NUW
    end do

  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' LSGN IN SUBROUTINE MKSGNUM'
do MIDV=1,CIS%nMidV
  write(u6,'(1X,''MIDV='',I3,/,(20I6))') MIDV,EXS%LSGN(:,MIDV)
end do
write(u6,*)
write(u6,*) ' USGN IN SUBROUTINE MKSGNUM'
do MIDV=1,CIS%nMidV
  write(u6,'(1X,''MIDV='',I3,/,(20I6))') MIDV,EXS%USGN(:,MIDV)
end do
write(u6,*)
#endif

call mma_deallocate(ISTEPVEC)

end subroutine MKSGNUM
