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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************

!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
subroutine MKSEG(SGS,CIS,EXS)
! PURPOSE: CONSTRUCT THE TABLES ISGM AND VSGM.
! ISGM(IVLT,ISGT) REFERS TO A SEGMENT OF THE SEGMENT TYPE
!    ISGT=1,..,26, WHOSE TOP LEFT VERTEX IS IVLT. ISGM GIVES
!    ZERO IF THE SEGMENT IS IMPOSSIBLE IN THE GRAPH DEFINED BY
!    THE PALDUS TABLE DRT. ELSE IT IS THE BOTTOM LEFT VERTEX
!    NUMBER OF THE SEGMENT. THE SEGMENT VALUE IS THEN VSGM.

use gugx, only: CIStruct, EXStruct, SGStruct
use segtab, only: IC1, IC2, ISVC, ITVPT
use stdalloc, only: mma_allocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IA, IAL, IB, IBL, ISGT, ITT, IV, IV1, IV2, IVL, IVLB, IVLT, IVRB, IVRT, LEV, MV, MVLL, MVRR
real(kind=wp) :: V
integer(kind=iwp), parameter :: IATAB = 3, IBTAB = 4
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ICL, ICR, ID
character(len=20) :: TEXT
character(len=*), parameter :: FRML(7) = ['        1.0         ','       -1.0         ','        1/(B+1)     ', &
                                          '       -1/(B+1)     ','  SQRT(   B /(B+1)) ','  SQRT((B+2)/(B+1)) ', &
                                          '  SQRT(B(B+2))/(B+1)']
#endif

call mma_allocate(CIS%IVR,SGS%nVert,2,Label='CIS%IVR')
call mma_allocate(CIS%ISGM,SGS%nVert,26,Label='CIS%ISGM')
call mma_allocate(CIS%VSGM,SGS%nVert,26,Label='CIS%VSGM')
call mma_allocate(EXS%MVL,CIS%nMidV,2,Label='EXS%MVL')
call mma_allocate(EXS%MVR,CIS%nMidV,2,Label='EXS%MVR')

CIS%IVR(:,:) = 0
do LEV=1,SGS%nLev
  IV1 = SGS%LTV(LEV)
  IV2 = SGS%LTV(LEV-1)-1
  do IVL=IV1,IV2
    IAL = SGS%DRT(IVL,IATAB)
    IBL = SGS%DRT(IVL,IBTAB)
    do IV=IVL+1,IV2
      IA = SGS%DRT(IV,IATAB)
      if (IA == IAL) then
        IB = SGS%DRT(IV,IBTAB)
        if (IB == (IBL-1)) CIS%IVR(IVL,1) = IV
      else if (IA == (IAL-1)) then
        IB = SGS%DRT(IV,IBTAB)
        if (IB == (IBL+1)) CIS%IVR(IVL,2) = IV
      end if
    end do
  end do
end do
! CONSTRUCT THE MVL AND MVR TABLES:
do IVL=SGS%MVSta,SGS%MVEnd
  MVLL = IVL-SGS%MVSta+1
  MVRR = 0
  if (CIS%IVR(IVL,1) /= 0) MVRR = CIS%IVR(IVL,1)-SGS%MVSta+1
  EXS%MVR(MVLL,1) = MVRR
  MVRR = 0
  if (CIS%IVR(IVL,2) /= 0) MVRR = CIS%IVR(IVL,2)-SGS%MVSta+1
  EXS%MVR(MVLL,2) = MVRR
  EXS%MVL(MVLL,1) = 0
  EXS%MVL(MVLL,2) = 0
end do
do MV=1,CIS%nMidV
  if (EXS%MVR(MV,1) /= 0) EXS%MVL(EXS%MVR(MV,1),1) = MV
  if (EXS%MVR(MV,2) /= 0) EXS%MVL(EXS%MVR(MV,2),2) = MV
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' MIDVERT PAIR TABLES MVL,MVR IN MKSEG:'
write(u6,*) ' MVL TABLE:'
write(u6,1234) (MV,EXS%MVL(MV,1),EXS%MVL(MV,2),MV=1,CIS%nMidV)
write(u6,*) ' MVR TABLE:'
write(u6,1234) (MV,EXS%MVR(MV,1),EXS%MVR(MV,2),MV=1,CIS%nMidV)
write(u6,*)
write(u6,*) ' VERTEX PAIR TABLE IVR IN MKSEG:'
write(u6,1234) (IVL,CIS%IVR(IVL,1),CIS%IVR(IVL,2),IVL=1,SGS%nVert)
#endif

! INITIALIZE SEGMENT TABLES, AND MARK VERTICES AS UNUSABLE:
do IVLT=1,SGS%nVert
  do ISGT=1,26
    CIS%ISGM(IVLT,ISGT) = 0
    CIS%VSGM(IVLT,ISGT) = Zero
  end do
end do
do IVLT=1,SGS%nVert
  do ISGT=1,26
    ITT = ITVPT(ISGT)
    IVRT = IVLT
    if ((ITT == 1) .or. (ITT == 2)) IVRT = CIS%IVR(IVLT,ITT)
    if (IVRT == 0) cycle
    IVLB = SGS%Down(IVLT,IC1(ISGT))
    if (IVLB == 0) cycle
    IVRB = SGS%Down(IVRT,IC2(ISGT))
    if (IVRB == 0) cycle
    ! SEGMENT IS NOW ACCEPTED AS POSSIBLE.

    CIS%ISGM(IVLT,ISGT) = IVLB
    IB = SGS%DRT(IVLT,IBTAB)
    select case (ISVC(ISGT))
      case (1)
        V = One
      case (2)
        V = -One
      case (3)
        V = One/real(1+IB,kind=wp)
      case (4)
        V = -One/real(1+IB,kind=wp)
      case (5)
        V = sqrt(real(IB,kind=wp)/real(1+IB,kind=wp))
      case (6)
        V = sqrt(real(2+IB,kind=wp)/real(1+IB,kind=wp))
      case (7)
        V = sqrt(real(IB*(2+IB),kind=wp))/real(1+IB,kind=wp)
      case default
        V = Zero ! Dummy assignment
        call Abend()
    end select
    CIS%VSGM(IVLT,ISGT) = V
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' SEGMENT TABLE IN MKSEG.'
write(u6,*) ' VLT SGT ICL ICR VLB       SEGMENT TYPE         FORMULA'
do IV=1,SGS%nVert
  do ISGT=1,26
    ID = CIS%ISGM(IV,ISGT)
    if (ID == 0) cycle
    ICL = IC1(ISGT)
    ICR = IC2(ISGT)
    if (ISGT <= 4) TEXT = '  WALK SECTION.'
    if ((ISGT >= 5) .and. (ISGT <= 8)) TEXT = ' TOP SEGMENT.'
    if ((ISGT >= 9) .and. (ISGT <= 18)) TEXT = ' MID-SEGMENT.'
    if ((ISGT >= 19) .and. (ISGT <= 22)) TEXT = ' BOTTOM SEGMENT.'
    if (ISGT > 22) TEXT = ' DOWN-WALK SECTION.'
    write(u6,2345) IV,ISGT,ICL,ICR,ID,TEXT,FRML(ISVC(ISGT))
  end do
end do
1234 format(3(3(1X,I4),4X))
2345 format(1X,5I4,5X,A20,5X,A20)
# endif

end subroutine MKSEG
