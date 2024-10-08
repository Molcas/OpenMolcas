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
subroutine NRCOUP(SGS,CIS,EXS)

use Symmetry_Info, only: Mul
use gugx, only: CIStruct, EXStruct, SGStruct
use segtab, only: IBVPT, IC1, ITVPT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: IBSYM, ICL, INDEO, INDEOB, INDEOT, IP, IPQ, IQ, ISGT, ISYDS1, ISYM, ISYUS1, ITSYM, IVLB, IVLT, LEV, LFTSYM, &
                     MV, MV1, MV2, MV3, MV4, MV5, MXDWN, MXUP, N, NDWNS1, NSGMX, NSGTMP, NT1TMP, NT2TMP, NT3TMP, NT4TMP, NT5TMP, &
                     NUPS1, NUW
integer(kind=iwp), allocatable :: NRL(:,:,:)
integer(kind=iwp), parameter :: LTAB = 1
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IS, IST, NCP
#endif

call mma_allocate(CIS%NOW,2,SGS%nSym,CIS%nMidV,Label='CIS%NOW',safe='*')
call mma_allocate(CIS%IOW,2,SGS%nSym,CIS%nMidV,Label='CIS%IOW',safe='*')
call mma_allocate(CIS%NCSF,SGS%nSym,Label='CIS%NCSF',safe='*')
call mma_allocate(CIS%NOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%NOCSF',safe='*')
call mma_allocate(CIS%IOCSF,SGS%nSym,CIS%nMidV,SGS%nSym,Label='CIS%IOCSF',safe='*')

EXS%MxEO = (SGS%nLev*(SGS%nLev+5))/2
call mma_allocate(EXS%NOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%NOCP')
call mma_allocate(EXS%IOCP,EXS%MxEO,SGS%nSym,CIS%nMidV,Label='EXS%IOCP')
call mma_allocate(NRL,[1,SGS%nSym],[1,SGS%nVert],[0,EXS%MxEO],Label='NRL')

NRL(:,1:SGS%MVEnd,:) = 0
NRL(1,1,0) = 1

do IVLT=1,SGS%MVSta-1
  LEV = SGS%DRT(IVLT,LTAB)
  do ISGT=1,26
    IVLB = CIS%ISGM(IVLT,ISGT)
    if (IVLB == 0) cycle
    ICL = IC1(ISGT)
    ISYM = 1
    if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
    do ITSYM=1,SGS%nSym
      IBSYM = Mul(ITSYM,ISYM)

      select case (ISGT)
        case (:4)
          ! THIS IS AN UPPER WALK.
          NRL(IBSYM,IVLB,0) = NRL(IBSYM,IVLB,0)+NRL(ITSYM,IVLT,0)
        case (5:8)
          ! THIS IS AN TOP SEGMENT.
          INDEO = LEV+(IBVPT(ISGT)-1)*SGS%nLev
          NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,0)
        case (9:18)
          ! THIS IS A MID-SEGMENT.
          do IP=LEV+1,SGS%nLev
            INDEOT = IP+(ITVPT(ISGT)-1)*SGS%nLev
            INDEOB = IP+(IBVPT(ISGT)-1)*SGS%nLev
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
          end do
        case (19:22)
          ! THIS IS A BOTTOM SEGMENT.
          do IP=LEV+1,SGS%nLev
            INDEOT = IP+(ITVPT(ISGT)-1)*SGS%nLev
            IPQ = (IP*(IP-1))/2+LEV
            INDEOB = IPQ+2*SGS%nLev
            NRL(IBSYM,IVLB,INDEOB) = NRL(IBSYM,IVLB,INDEOB)+NRL(ITSYM,IVLT,INDEOT)
          end do
        case default
          ! THIS IS A LOWER WALK.
          do INDEO=2*SGS%nLev+1,EXS%MxEO
            NRL(IBSYM,IVLB,INDEO) = NRL(IBSYM,IVLB,INDEO)+NRL(ITSYM,IVLT,INDEO)
          end do
      end select
    end do
  end do
end do

MXUP = 0
do MV=1,CIS%nMidV
  IVLT = MV+SGS%MVSta-1
  do LFTSYM=1,SGS%nSym
    CIS%NOW(1,LFTSYM,MV) = NRL(LFTSYM,IVLT,0)
    MXUP = max(MXUP,CIS%NOW(1,LFTSYM,MV))
    do INDEO=1,EXS%MxEO
      EXS%NOCP(INDEO,LFTSYM,MV) = NRL(LFTSYM,IVLT,INDEO)
    end do
  end do
end do

NRL(:,SGS%MVSta:,:) = 0
NRL(1,SGS%nVert,0) = 1
do IVLT=SGS%nVert-1,SGS%MVSta,-1
  LEV = SGS%DRT(IVLT,LTAB)
  do ISGT=1,26
    IVLB = CIS%ISGM(IVLT,ISGT)
    if (IVLB == 0) cycle
    ICL = IC1(ISGT)
    ISYM = 1
    if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
    do ITSYM=1,SGS%nSym
      IBSYM = Mul(ITSYM,ISYM)
      select case (ISGT)
        case (23:)
          ! THIS IS A LOWER WALK.
          NRL(ITSYM,IVLT,0) = NRL(ITSYM,IVLT,0)+NRL(IBSYM,IVLB,0)
        case (19:22)
          ! THIS IS AN BOTTOM SEGMENT.
          INDEO = LEV+(ITVPT(ISGT)-1)*SGS%nLev
          NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,0)
        case (9:18)
          ! THIS IS A MID-SEGMENT.
          do IQ=1,LEV-1
            INDEOT = IQ+(ITVPT(ISGT)-1)*SGS%nLev
            INDEOB = IQ+(IBVPT(ISGT)-1)*SGS%nLev
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
          end do
        case (5:8)
          ! THIS IS AN TOP SEGMENT.
          do IQ=1,LEV-1
            INDEOB = IQ+(IBVPT(ISGT)-1)*SGS%nLev
            IPQ = (LEV*(LEV-1))/2+IQ
            INDEOT = IPQ+2*SGS%nLev
            NRL(ITSYM,IVLT,INDEOT) = NRL(ITSYM,IVLT,INDEOT)+NRL(IBSYM,IVLB,INDEOB)
          end do
        case default
          ! THIS IS AN UPPER WALK.
          do IPQ=1,(LEV*(LEV-1))/2
            INDEO = 2*SGS%nLev+IPQ
            NRL(ITSYM,IVLT,INDEO) = NRL(ITSYM,IVLT,INDEO)+NRL(IBSYM,IVLB,INDEO)
          end do
      end select
    end do
  end do
end do

MXDWN = 0
do MV=1,CIS%nMidV
  IVLT = MV+SGS%MVSta-1
  do LFTSYM=1,SGS%nSym
    CIS%NOW(2,LFTSYM,MV) = NRL(LFTSYM,IVLT,0)
    MXDWN = max(MXDWN,CIS%NOW(2,LFTSYM,MV))
    do INDEO=1,EXS%MxEO
      N = NRL(LFTSYM,IVLT,INDEO)
      if (N /= 0) EXS%NOCP(INDEO,LFTSYM,MV) = N
    end do
  end do
end do

EXS%nICOup = 0
do INDEO=1,EXS%MxEO
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      EXS%IOCP(INDEO,LFTSYM,MV) = EXS%nICOup
      EXS%nICOup = EXS%nICOup+EXS%NOCP(INDEO,LFTSYM,MV)
    end do
  end do
end do

call CSFCOUNT(CIS,SGS%nSym,NUW)

!AR INSERT FOR US IN SIGMA ROUTINE

NSGMX = 1
NDWNS1 = 0  ! Dummy initialization
NSGTMP = max(MXUP,MXDWN)
do MV3=1,CIS%nMidV
  MV1 = EXS%MVL(MV3,2)
  MV2 = EXS%MVL(MV3,1)
  MV4 = EXS%MVR(MV3,1)
  MV5 = EXS%MVR(MV3,2)
  do ISYUS1=1,SGS%nSym
    NUPS1 = CIS%NOW(1,ISYUS1,MV3)
    do ISYDS1=1,SGS%nSym
      NDWNS1 = CIS%NOW(2,ISYDS1,MV3)
      NSGMX = max(NSGMX,CIS%NOCSF(ISYUS1,MV3,ISYDS1))

      if (MV1 /= 0) then
        NT4TMP = NUPS1*CIS%NOW(2,ISYDS1,MV1)
        NSGTMP = max(NSGTMP,NT4TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV1)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV2 /= 0) then
        NT3TMP = NUPS1*CIS%NOW(2,ISYDS1,MV2)
        NSGTMP = max(NSGTMP,NT3TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV2)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV4 /= 0) then
        NT1TMP = NUPS1*CIS%NOW(2,ISYDS1,MV4)
        NSGTMP = max(NSGTMP,NT1TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV4)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

      if (MV5 /= 0) then
        NT2TMP = NUPS1*CIS%NOW(2,ISYDS1,MV5)
        NSGTMP = max(NSGTMP,NT2TMP)
        NT5TMP = NDWNS1*CIS%NOW(1,ISYUS1,MV5)
        NSGTMP = max(NSGTMP,NT5TMP)
      end if

    end do
  end do
end do

call mma_allocate(EXS%SGTMP,NSGTMP,Label='EXS%SGTMP')

#ifdef _DEBUGPRINT_
write(u6,600) MXUP,MXDWN,NSGMX,NSGMX,NSGTMP

!AR END OF INSERT
write(u6,*)
write(u6,*) ' TOTAL NR OF WALKS: UPPER ',NUW
write(u6,*) '                    LOWER ',CIS%nWalk-NUW
write(u6,*) '                     SUM  ',CIS%nWalk
write(u6,*) ' TOTAL NR OF COUPL COEFFS ',EXS%nICOup
INDEO = 2*SGS%nLev+1
write(u6,*) '         OF TYPE 1&2 ONLY:',EXS%IOCP(INDEO,1,1)
write(u6,*)
write(u6,*) ' NR OF CONFIGURATIONS/SYMM:'
write(u6,'(8(1X,I8))') (CIS%NCSF(IS),IS=1,SGS%nSym)
write(u6,*)

write(u6,*)
write(u6,*) ' NR OF WALKS AND CONFIGURATIONS IN NRCOUP'
write(u6,*) ' BY MIDVERTEX AND SYMMETRY.'
do MV=1,CIS%nMidV
  write(u6,*)
  write(u6,1234) MV,(CIS%NOW(1,IS,MV),IS=1,SGS%nSym)
  write(u6,1235) (CIS%NOW(2,IS,MV),IS=1,SGS%nSym)
  do IST=1,SGS%nSym
    write(u6,1236) IST,(CIS%NOCSF(IS,MV,IST),IS=1,SGS%nSym)
  end do
end do
write(u6,*)
write(u6,*) ' NR OF COUPLING COEFFICIENTS:'
write(u6,*) ' 1. OPEN LOOPS TYPE 1:'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      NCP = EXS%NOCP(IP,IS,MV)
      if (NCP /= 0) write(u6,2345) IP,MV,IS,NCP
    end do
  end do
end do
write(u6,*)
write(u6,*) ' 2. OPEN LOOPS TYPE 2:'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      NCP = EXS%NOCP(SGS%nLev+IP,IS,MV)
      if (NCP /= 0) write(u6,2345) IP,MV,IS,NCP
    end do
  end do
end do
write(u6,*)
write(u6,*) ' 3. CLOSED LOOPS:'
do IP=2,SGS%nLev
  do IQ=1,IP-1
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do IS=1,SGS%nSym
        NCP = EXS%NOCP(INDEO,IS,MV)
        if (NCP /= 0) write(u6,2346) IP,IQ,MV,IS,NCP
      end do
    end do
  end do
end do
600 format(/,' MAXIMUM NUMBER OF WALKS', &
           /,' UPPER ',I6,' LOWER ',I6, &
           /,' LENGTH OF LARGEST WORK ARRAYS IN SUBROUTINE SIGMA', &
           /,' TEMPORARY SGM1 ',I7, &
           /,' TEMPORARY SGM2 ',I7, &
           /,' NSGTMP         ',I7)
1234 format('  MV=',I2,'    UPPER WALKS:',8I6)
1235 format('           LOWER WALKS:',8I6)
1236 format(' IST=',I2,'  CONFIGURATIONS:',8I6)
2345 format(' P=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
2346 format(' P=',I2,'  Q=',I2,'  MV=',I2,' SYMM ',I1,' NOCP=',I4)
#endif

call mma_deallocate(NRL)

end subroutine NRCOUP
