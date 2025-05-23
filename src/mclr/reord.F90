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
! Copyright (C) 1990, Markus P. Fuelscher                              *
!               1990, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine REORD(SGS,CIS,EXS,NCONF,IMODE,ICONF,ISPIN,kSym,CIOLD)
! AUTHOR:  M.P. FUELSCHER AND J. OLSEN
!          UNIV. OF LUND, SWEDEN 1990
!
! PURPOSE: CONSTRUCT THE REINDEXING ARRAY WHICH REORDERS
!          THE CSFS GENERATED BY THE DETERMINANT CODE INTO
!          THE SPLIT GRAPH GUGA ORDER.
!          FOR THAT CONSTRUCT FOR EACH CSF THE CORRESPONDING
!          STEP VECTOR AND PASS IT TO THE FUNCTIONS IPHASE
!          AND ISGNUM WHICH COMPUTES THE THE PHASE FACTOR
!          INVOLVED WHEN GOING FROM THE SYMMETRIC TO THE
!          UNITARY GROUP AND THE SPLIT ORDERING NUMBER.
!          FINALLY, THE VECTORS OF CI COEFFICIENTS ARE
!          REORDERED AND, TWO MODE ARE POSSIBLE:
!          IMODE=0 : FROM SYMMETRIC GROUP TO SPLIT GRAPH UGA ORDER
!          IMODE=1 : FROM SPLIT GRAPH UGA TO SYMMETRIC GROUP ORDER

use MCLR_data, only: minop, NCNATS, NCPCNT, NTYP
use gugx, only: CIStruct, EXStruct, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(in) :: EXS
integer(kind=iwp), intent(in) :: nCONF, IMODE, ICONF(*), ISPIN(*), kSym
real(kind=wp), intent(inout) :: CIOLD(NCONF)
integer(kind=iwp) :: IC, ICL, ICNBS, ICNBS0, ICSBAS, ICSFJP, IICSF, IOPEN, IP, IPBAS, ISG, ITYP, IWALK(50)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif
real(kind=wp) :: PHASE
real(kind=wp), allocatable :: CINEW(:)
integer(kind=iwp), external :: IPHASE, ISGNUM

call mma_allocate(CINEW,NCONF,Label='CINEW')

! LOOP OVER CONFIGURATIONS TYPES

ICSFJP = 0
ICNBS0 = 0 ! dummy initialize
IPBAS = 0 ! dummy initialize
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1 ! MINOP IS NOT INITIALIZED TEOEAW
  ICL = (SGS%nActEl-IOPEN)/2
  ! BASE ADDRESS FOR CONFIGURATION OF THIS TYPE
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNATS(ITYP-1,kSym)*(SGS%nActEl+IOPEN-1)/2
  end if
  ! BASE ADDRESS FOR PROTOTYPE SPIN COUPLINGS
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NCPCNT(ITYP-1)*(IOPEN-1)
  end if

  ! LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
  ! SPIN COUPLINGS

  do IC=1,NCNATS(ITYP,kSym)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    do IICSF=1,NCPCNT(ITYP)
      ICSFJP = ICSFJP+1
      ICSBAS = IPBAS+(IICSF-1)*IOPEN
      ! COMPUTE STEP VECTOR
      call STEPVEC(ICONF(ICNBS),ICONF(ICNBS+ICL),ICL,IOPEN,ISPIN(ICSBAS),SGS%nLev,IWALK)
      ! GET SPLIT GRAPH ORDERING NUMBER
      ! FUNCTION ISGNUM
      ISG = ISGNUM(SGS%nLev,SGS%nVert,SGS%MidLev,SGS%MVSta,CIS%nMidV,SGS%MxUp,SGS%MxDwn,SGS%Down,SGS%Up,SGS%DAW,SGS%RAW,EXS%USGN, &
                   EXS%LSGN,IWALK)
      ! GET PHASE PHASE FACTOR
      IP = IPHASE(SGS%nLev,SGS%nVert,SGS%DRT,SGS%Up,IWALK)
      ! NOW REORDER THIS ELEMENT OF THE CI-VECTOR
      PHASE = real(IP,kind=wp)
      select case (iMode)
        case (0)
          CINEW(ISG) = CIOLD(ICSFJP)*PHASE
        case (1)
          CINEW(ICSFJP) = CIOLD(ISG)*PHASE
      end select
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' OLD CI-VECTORS IN SUBROUTINE REORD'
write(u6,'(10F12.8)') (CIOLD(I),I=1,NCONF)
write(u6,*) ' NEW CI-VECTORS IN SUBROUTINE REORD'
write(u6,'(10F12.8)') (CINEW(I),I=1,NCONF)
write(u6,*)
#endif

CIOLD(:) = CINEW(:)
call mma_deallocate(CINEW)

end subroutine REORD
