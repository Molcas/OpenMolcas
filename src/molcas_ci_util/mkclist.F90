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

subroutine MKCLIST(SGS,CIS)
! PURPOSE: CONSTRUCT THE COMPRESSED CASE-LIST, I.E.,
!          STORE THE STEP VECTOR FOR ALL POSSIBLE WALKS
!          IN THE ARRAY ICASE. GROUPS OF 15 CASES ARE PACKED
!          INTO ONE INTEGER WORD.

use Symmetry_Info, only: Mul
use gugx, only: CIStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
type(SGStruct), intent(inout) :: SGS
type(CIStruct), intent(inout) :: CIS
integer(kind=iwp) :: IC, IHALF, ILND, IPOS, ISML, ISTP, IVB, IVT, IVTEND, IVTOP, IVTSTA, IWSYM, L, LEV, LEV1, LEV2, LL, MV
logical(kind=iwp) :: Found
integer(kind=iwp), parameter :: IVERT = 1, ISYM = 2, ISTEP = 3

call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')
call mma_allocate(SGS%Scr,[1,3],[0,SGS%nLev],Label='SGS%Scr',safe='*')

! CLEAR ARRAY NOW. IT WILL BE RESTORED FINALLY

CIS%NOW(:,:,:) = 0

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
    ! SET CURRENT LEVEL=TOP LEVEL OF SUBGRAPH:
    LEV = LEV1
    SGS%Scr(IVERT,LEV) = IVTOP
    SGS%Scr(ISYM,LEV) = 1
    SGS%Scr(ISTEP,LEV) = -1
    do while (LEV <= LEV1)
      ! FIND FIRST POSSIBLE UNTRIED ARC DOWN FROM CURRENT VERTEX:
      IVT = SGS%Scr(IVERT,LEV)
      Found = .false.
      do ISTP=SGS%Scr(ISTEP,LEV)+1,3
        IVB = SGS%Down(IVT,ISTP)
        if (IVB /= 0) then
          Found = .true.
          exit
        end if
      end do
      if (Found) then
        ! ALT A -- SUCH AN ARC WAS FOUND. WALK DOWN:
        SGS%Scr(ISTEP,LEV) = ISTP
        ISML = 1
        if ((ISTP == 1) .or. (ISTP == 2)) ISML = SGS%ISm(LEV)
        LEV = LEV-1
        SGS%Scr(ISYM,LEV) = Mul(ISML,SGS%Scr(ISYM,LEV+1))
        SGS%Scr(IVERT,LEV) = IVB
        SGS%Scr(ISTEP,LEV) = -1
        if (LEV > LEV2) cycle
        ! WE HAVE REACHED THE LOWER LEVEL. THE WALK IS COMPLETE.
        ! MIDVERTEX NUMBER:
        MV = SGS%Scr(IVERT,SGS%MidLev)+1-SGS%MVSta
        ! SYMMETRY LABEL OF THIS WALK:
        IWSYM = SGS%Scr(ISYM,LEV2)
        ! ITS ORDERING NUMBER WITHIN THE SAME BATCH OF (IHALF,IWSYM,MV):
        ILND = 1+CIS%NOW(IHALF,IWSYM,MV)
        CIS%NOW(IHALF,IWSYM,MV) = ILND
        ! CONSEQUENTLY, THE POSITION IMMEDIATELY BEFORE THIS COMPRESSED WALK:
        IPOS = CIS%IOW(IHALF,IWSYM,MV)+(ILND-1)*CIS%nIpWlk
        ! PACK THE STEPS IN GROUPS OF 15 LEVELS PER INTEGER:
        do LL=LEV2+1,LEV1,15
          IC = 0
          do L=min(LL+14,LEV1),LL,-1
            IC = 4*IC+SGS%Scr(ISTEP,L)
          end do
          IPOS = IPOS+1
          CIS%ICase(IPOS) = IC
        end do
        ! FINISHED WITH THIS WALK. BACK UP ONE LEVEL AND TRY AGAIN:
        LEV = LEV+1
      else
        ! ALT B -- NO SUCH ARC WAS POSSIBLE. GO UP ONE STEP AND TRY AGAIN.
        SGS%Scr(ISTEP,LEV) = -1
        LEV = LEV+1
      end if
    end do
  end do
end do

call mma_deallocate(SGS%Scr)

end subroutine MKCLIST
