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

subroutine MKCOT(SGS,CIS)
! PURPOSE: SET UP COUNTER AND OFFSET TABLES FOR WALKS AND CSFS
! NOTE:    TO GET GET VARIOUS COUNTER AND OFFSET TABLES
!          THE DOWN-CHAIN TABLE IS SCANNED TO PRODUCE ALL POSSIBLE
!          WALKS. POSSIBLY, THERE ARE MORE EFFICIENT WAYS, BUT
!          SINCE ONLY UPPER AND LOWER WALKS ARE REQUIRED
!          THEIR NUMBER IS VERY LIMITTED, EVEN FOR LARGE CASES.

use Symmetry_Info, only: Mul
use gugx, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp) :: IHALF, ILND, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, LEV, LEV1, LEV2, MV, NUW
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IS
#endif
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

CIS%nIpWlk = 1+(SGS%MidLev-1)/15
CIS%nIpWlk = max(CIS%nIpWlk,1+(SGS%nLev-SGS%MidLev-1)/15)
call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF')
call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF')
call mma_allocate(SGS%Scr,[1,3],[0,SGS%nLev],Label='SGS%Scr')

! CLEAR ARRAYS IOW AND NOW

CIS%NOW(:,:,:) = 0
CIS%IOW(:,:,:) = 0

! CLEAR ARRAYS IOCSF AND NOCSF

CIS%IOCSF(:,:,:) = 0
CIS%NOCSF(:,:,:) = 0

! START MAIN LOOP OVER UPPER AND LOWER WALKS, RESPECTIVELY.

do IHALF=1,2
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = SGS%nLev
    LEV2 = SGS%MidLev
  else
    IVTSTA = SGS%MVSta
    IVTEND = SGS%MVEnd
    LEV1 = SGS%MidLev
    LEV2 = 0
  end if

  ! LOOP OVER VERTICES STARTING AT TOP OF SUBGRAPH

  do IVTOP=IVTSTA,IVTEND
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH
    LEV = LEV1
    SGS%Scr(IVERT,LEV) = IVTOP
    SGS%Scr(ISYM,LEV) = 1
    SGS%Scr(ISTEP,LEV) = -1
    do while (LEV <= LEV1)
      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX
      IVT = SGS%Scr(IVERT,LEV)
      do ISTP=SGS%Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) exit
      end do
      ! NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
      if (ISTP > 3) then
        SGS%Scr(ISTEP,LEV) = -1
        LEV = LEV+1
        cycle
      end if
      ! SUCH AN ARC WAS FOUND. WALK DOWN:
      SGS%Scr(ISTEP,LEV) = ISTP
      ISML = 1
      if ((ISTP == 1) .or. (ISTP == 2)) ISML = SGS%ISm(LEV)
      LEV = LEV-1
      SGS%Scr(ISYM,LEV) = Mul(ISML,SGS%Scr(ISYM,LEV+1))
      SGS%Scr(IVERT,LEV) = IVB
      SGS%Scr(ISTEP,LEV) = -1
      if (LEV > LEV2) cycle
      ! WE HAVE REACHED THE BOTTOM LEVEL. THE WALK IS COMPLETE.
      ! FIND MIDVERTEX NUMBER ORDERING NUMBER AND SYMMETRY OF THIS WALK
      MV = SGS%Scr(IVERT,SGS%MidLev)+1-SGS%MVSta
      IWSYM = SGS%Scr(ISYM,LEV2)
      ILND = 1+CIS%NOW(IHALF,IWSYM,MV)
      ! SAVE THE MAX WALK NUMBER FOR GIVEN SYMMETRY AND MIDVERTEX
      CIS%NOW(IHALF,IWSYM,MV) = ILND
      ! BACK UP ONE LEVEL AND TRY AGAIN:
      LEV = LEV+1
    end do
  end do
end do

call CSFCOUNT(CIS,SGS%nSym,NUW)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)
write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,'(A,I2,A,8I6)') '  MV=',MV,'    UPPER WALKS:',(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,'(A,8I6)') '           LOWER WALKS:',(CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do ISTP=1,SGS%nSym
    write(u6,'(A,I2,A,8I6)') ' ISTP=',ISTP,'  CONFIGURATIONS:',(CIS%NOCSF(IS,MV,ISTP),IS=1,SGS%nSym)
  end do
end do
#endif

end subroutine MKCOT
