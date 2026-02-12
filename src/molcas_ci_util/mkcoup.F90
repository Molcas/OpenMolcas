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

subroutine MKCOUP(SGS,CIS,EXS)
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
! Purpose: Compute and return the table ICOUP(1..3,ICOP).
! The number of coupling coeffs is obtained from NOCP, the offset to
! the ICOP numbering is given by IOCP. The numbers ICOUP(1..3,ICOP) are
! then the ket and bra half-walks, enumerated by the Lund scheme,
! and the index into the VTAB table (the pool of possible values of
! coupling coefficients).

! Any loop is regarded as a segment path from top to midlevel, or
! from midlevel to bottom.
! The segment path is described by the table ISGPTH. It is
! essentially a list of which one of segments nr 1..26 that are
! used at each level. The segments are of three types:
! Type 0: Upwalk segment or top loop segment.
! Type 1: Mid segment.
! Type 2: Bottom segment.
! Type 3: Downwalk segment.
! ISGPTH(IVLFT,LEV)=Left upper vertex.
! ISGPTH(ITYPE,LEV)=Type of segment, (0..3).
! ISGPTH(IAWSL,LEV)=Left arc weight sum (from top, or from bottom).
! ISGPTH(IAWSR,LEV)=Similar, right.
! ISGPTH(ILS  ,LEV)=Left symmetry label (accum from top or bottom).
! ISGPTH(ICS  ,LEV)=Left coupling case number (0..3).
! ISGPTH(ISEG ,LEV)=Segment type (1..26).
! These indices are used to denote the columns of table ISGPTH.

use Symmetry_Info, only: Mul
use gugx, only: CIStruct, EXStruct, SGStruct
use segtab, only: IBVPT, IC1, IC2, ITVPT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(inout) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp) :: i, i1, i2, IAWS, IC, ICL, ICOP, ICR, IHALF, iLnd, IndEO, iP, iPos, iQ, iS, iSg, iSgt, iSym, iT, iTyp, iTypMx, &
                     iTypT, iVlb, iVlt, iVrt, iVrTop, iVTab, iVTEnd, iVTop, iVTSta, L, Lev, Lev1, Lev2, LftSym, LL, MV, nCheck, &
                     nVTab_Final
real(kind=wp) :: C
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I3, ICOP1, ICOP2, ICP1, ICP2, N, NRC, NRCPQ
real(kind=wp) :: CP
#endif
integer(kind=iwp), allocatable :: ILNDW(:), ISGPTH(:,:)
real(kind=wp), allocatable :: val(:), VTab(:)
integer(kind=iwp), parameter :: IVLFT = 1, ITYPE = 2, IAWSL = 3, IAWSR = 4, ILS = 5, ICS = 6, ISEG = 7, nVTab = 5000

call mma_allocate(ILNDW,CIS%nWalk,Label='ILNDW')
call mma_allocate(ISGPTH,[1,7],[0,SGS%nLev],Label='ISGPTH')
call mma_allocate(val,[0,SGS%nLev],Label='val')
call mma_allocate(VTab,nVTab,Label='VTab')

call mma_allocate(EXS%ICoup,3,EXS%nICoup,Label='EXS%ICoup')
call mma_allocate(CIS%ICase,CIS%nWalk*CIS%nIpWlk,Label='CIS%ICase',safe='*')

!CIS%nIpWlk = 1+(SGS%MidLev-1)/15
!CIS%nIpWlk = max(CIS%nIpWlk,1+(nLev-SGS%MidLev-1)/15)
! NOW IS ZEROED AND WILL BE USED AS AN ARRAY OF COUNTERS, BUT WILL BE RESTORED FINALLY.
do IHALF=1,2
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      CIS%NOW(IHALF,IS,MV) = 0
    end do
  end do
end do
! SIMILAR FOR THE COUPLING COEFFICIENT TABLE:
do INDEO=1,EXS%MxEO
  do MV=1,CIS%nMidV
    do IS=1,SGS%nSym
      EXS%NOCP(INDEO,IS,MV) = 0
    end do
  end do
end do

! COUPLING COEFFICIENT VALUE TABLE:
NVTAB_FINAL = 2
VTab(1) = One
VTab(2) = -One

NCHECK = 0

do IHALF=1,2
  if (IHALF == 1) then
    IVTSTA = 1
    IVTEND = 1
    LEV1 = SGS%nLev
    LEV2 = SGS%MidLev
    ITYPMX = 0
  else
    IVTSTA = SGS%MVSta
    IVTEND = SGS%MVEnd
    LEV1 = SGS%MidLev
    LEV2 = 0
    ITYPMX = 2
  end if
  do IVTOP=IVTSTA,IVTEND
    do ITYP=0,ITYPMX
      IVRTOP = IVTOP
      if (ITYP > 0) IVRTOP = CIS%IVR(IVTOP,ITYP)
      if (IVRTOP == 0) cycle
      LEV = LEV1
      ISGPTH(IVLFT,LEV) = IVTOP
      ISGPTH(ITYPE,LEV) = ITYP
      ISGPTH(IAWSL,LEV) = 0
      ISGPTH(IAWSR,LEV) = 0
      ISGPTH(ILS,LEV) = 1
      ISGPTH(ISEG,LEV) = 0
      val(LEV) = One
      do while (LEV <= LEV1)
        ITYPT = ISGPTH(ITYPE,LEV)
        IVLT = ISGPTH(IVLFT,LEV)
        do ISGT=ISGPTH(ISEG,LEV)+1,26
          IVLB = CIS%ISGM(IVLT,ISGT)
          if (IVLB == 0) cycle
          if (ITYPT == ITVPT(ISGT)) exit
        end do
        if (ISGT > 26) then
          ISGPTH(ISEG,LEV) = 0
          LEV = LEV+1
        else
          ISGPTH(ISEG,LEV) = ISGT
          ICL = IC1(ISGT)
          ISYM = 1
          if ((ICL == 1) .or. (ICL == 2)) ISYM = SGS%ISm(LEV)
          IVRT = IVLT
          if ((ITYPT == 1) .or. (ITYPT == 2)) IVRT = CIS%IVR(IVLT,ITYPT)
          ICR = IC2(ISGT)
          ISGPTH(ICS,LEV) = ICL
          LEV = LEV-1
          ISGPTH(IAWSL,LEV) = ISGPTH(IAWSL,LEV+1)+SGS%MAW(IVLT,ICL)
          ISGPTH(IAWSR,LEV) = ISGPTH(IAWSR,LEV+1)+SGS%MAW(IVRT,ICR)
          val(LEV) = val(LEV+1)*CIS%VSGM(IVLT,ISGT)
          ISGPTH(ILS,LEV) = Mul(ISYM,ISGPTH(ILS,LEV+1))
          ISGPTH(IVLFT,LEV) = IVLB
          ISGPTH(ITYPE,LEV) = IBVPT(ISGT)
          ISGPTH(ISEG,LEV) = 0
          if (LEV > LEV2) cycle

          MV = ISGPTH(IVLFT,SGS%MidLev)+1-SGS%MVSta
          LFTSYM = ISGPTH(ILS,LEV2)
          IT = ISGPTH(ITYPE,SGS%MidLev)
          if (IT == 0) IT = 3
          if (ISGPTH(ITYPE,LEV2) == 0) IT = 0

          if (IT == 0) then
            ILND = 1+CIS%NOW(IHALF,LFTSYM,MV)
            IAWS = ISGPTH(IAWSL,LEV2)
            ILNDW(IAWS) = ILND
            CIS%NOW(IHALF,LFTSYM,MV) = ILND
            IPOS = CIS%IOW(IHALF,LFTSYM,MV)+(ILND-1)*CIS%nIpWlk
            do LL=LEV2+1,LEV1,15
              IC = 0
              do L=min(LL+14,LEV1),LL,-1
                IC = 4*IC+ISGPTH(ICS,L)
              end do
              IPOS = IPOS+1
              CIS%ICase(IPOS) = IC
            end do
          else
            IP = 0
            IQ = 0
            do L=LEV2+1,LEV1
              ISG = ISGPTH(ISEG,L)
              if ((ISG >= 5) .and. (ISG <= 8)) IP = L
              if ((ISG >= 19) .and. (ISG <= 22)) IQ = L
            end do
            if (IP == 0) IP = IQ
            INDEO = SGS%nLev*(IT-1)+IP
            if (IT == 3) INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
            ICOP = 1+EXS%NOCP(INDEO,LFTSYM,MV)
            EXS%NOCP(INDEO,LFTSYM,MV) = ICOP
            ICOP = EXS%IOCP(INDEO,LFTSYM,MV)+ICOP
            NCHECK = NCHECK+1
            if (ICOP > EXS%nICoup) then
              write(u6,*) ' ERROR: NICOUP=',EXS%nICoup
              write(u6,*) ' NR OF COUPS PRODUCED:',NCHECK
              write(u6,*) '           TYPE NR IT:',IT
              write(u6,*) '            IP,IQ    :',IP,IQ
              write(u6,*) '            INDEO    :',INDEO
              write(u6,*) '        MIDVERTEX MV :',MV
              write(u6,*) ' LEFT SYMMETRY LFTSYM:',LFTSYM
              write(u6,*) ' COUP OFFSET IOCP    :',EXS%IOCP(INDEO,LFTSYM,MV)
              write(u6,*) ' COUP SERIAL NR ICOP :',ICOP
              write(u6,*) ' D:O, WITHOUT OFFSET :',ICOP-EXS%IOCP(INDEO,LFTSYM,MV)
              write(u6,*) ' CURRENT NOCP NUMBER :',EXS%NOCP(INDEO,LFTSYM,MV)
              call ABEND()
            end if

            C = val(LEV2)
            do I=1,NVTAB_FINAL
              IVTAB = I
              if (abs(C-VTab(I)) < 1.0e-10_wp) exit
            end do
            if (I > NVTAB_FINAL) then
              NVTAB_FINAL = NVTAB_FINAL+1
              if (NVTAB_FINAL > nVTab) then
                write(u6,*) 'MKCOUP: NVTAB_FINAL=',NVTAB_FINAL
                write(u6,*) 'NVTAB_FINAL should not be allowed to grow'
                write(u6,*) 'beyond nVTab which was set provisionally'
                write(u6,*) 'in subroutine GINIT in file ginit.f.'
                write(u6,*) 'Now nVTab=',nVTab
                write(u6,*) 'This may indicate a problem with your input.'
                write(u6,*) 'If you do want to do this big calculation, try'
                write(u6,*) 'increasing nVTab in GINIT and recompile.'
                call ABEND()
              end if
              VTab(NVTAB_FINAL) = C
              IVTAB = NVTAB_FINAL
            end if
            EXS%ICoup(1,ICOP) = ISGPTH(IAWSL,LEV2)
            EXS%ICoup(2,ICOP) = ISGPTH(IAWSR,LEV2)
            EXS%ICoup(3,ICOP) = IVTAB
            if (ICOP > EXS%nICoup) then
              write(u6,*) 'MKCOUP: ICOP>NICOUP!'
              call ABEND()
            end if
          end if

          LEV = LEV+1
        end if
      end do

    end do
  end do
end do
! RENUMBER THE COUPLING COEFFICIENT INDICES BY LUND SCHEME:
do ICOP=1,EXS%nICoup
  I1 = EXS%ICoup(1,ICOP)
  I2 = EXS%ICoup(2,ICOP)
  EXS%ICoup(1,ICOP) = ILNDW(I1)
  EXS%ICoup(2,ICOP) = ILNDW(I2)
end do

call mma_deallocate(val)
call mma_deallocate(ISGPTH)
call mma_deallocate(ILNDW)

#ifdef _DEBUGPRINT_
ICOP1 = 0
ICOP2 = 0
write(u6,*) ' NR OF DIFFERENT VALUES OF COUP:',NVTAB_FINAL
do ICOP=1,EXS%nICoup
  I3 = EXS%ICoup(3,ICOP)
  if (I3 == 1) ICOP1 = ICOP1+1
  if (I3 == 2) ICOP2 = ICOP2+1
end do
write(u6,*)
write(u6,*) ' NR OF COUPS WITH VALUE  1.0:',ICOP1
write(u6,*) ' NR OF COUPS WITH VALUE -1.0:',ICOP2
write(u6,*)
write(u6,*) ' COUPLING COEFFICIENTS:'
write(u6,*) '    IP    IQ    MV LFTSYM NOCP'
write(u6,*) ' 1. OPEN LOOPS TYPE 1.'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(IP,LFTSYM,MV)
      write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
    end do
  end do
end do
write(u6,*) ' 2. OPEN LOOPS TYPE 2.'
do IP=1,SGS%nLev
  INDEO = SGS%nLev+IP
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(INDEO,LFTSYM,MV)
      write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,N
    end do
  end do
end do
write(u6,*) ' 3. CLOSED LOOPS.'
do IP=1,SGS%nLev
  do IQ=1,IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do LFTSYM=1,SGS%nSym
        N = EXS%NOCP(INDEO,LFTSYM,MV)
        write(u6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,N
      end do
    end do
  end do
end do
write(u6,*)
write(u6,*) ' COUPLING COEFFICIENTS:'
write(u6,*) '    IP    IQ    MV LFTSYM ICOP ICOUP1&2   COUP'
write(u6,*) ' 1. OPEN LOOPS TYPE 1.'
do IP=1,SGS%nLev
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(IP,LFTSYM,MV)
      ICOP = EXS%IOCP(IP,LFTSYM,MV)
      do I=1,N
        ICOP = ICOP+1
        ICP1 = EXS%ICoup(1,ICOP)
        ICP2 = EXS%ICoup(2,ICOP)
        CP = VTab(EXS%ICoup(3,ICOP))
        write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
      end do
    end do
  end do
end do
write(u6,*) ' 2. OPEN LOOPS TYPE 2.'
do IP=1,SGS%nLev
  INDEO = SGS%nLev+IP
  do MV=1,CIS%nMidV
    do LFTSYM=1,SGS%nSym
      N = EXS%NOCP(INDEO,LFTSYM,MV)
      ICOP = EXS%IOCP(INDEO,LFTSYM,MV)
      do I=1,N
        ICOP = ICOP+1
        ICP1 = EXS%ICoup(1,ICOP)
        ICP2 = EXS%ICoup(2,ICOP)
        CP = VTab(EXS%ICoup(3,ICOP))
        write(u6,'(6X,6(1X,I5),F10.7)') IP,MV,LFTSYM,ICOP,ICP1,ICP2,CP
      end do
    end do
  end do
end do
write(u6,*) ' 3. CLOSED LOOPS.'
do IP=1,SGS%nLev
  do IQ=1,IP
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do MV=1,CIS%nMidV
      do LFTSYM=1,SGS%nSym
        N = EXS%NOCP(INDEO,LFTSYM,MV)
        ICOP = EXS%IOCP(INDEO,LFTSYM,MV)
        do I=1,N
          ICOP = ICOP+1
          ICP1 = EXS%ICoup(1,ICOP)
          ICP2 = EXS%ICoup(2,ICOP)
          CP = VTab(EXS%ICoup(3,ICOP))
          write(u6,'(7(1X,I5),F10.7)') IP,IQ,MV,LFTSYM,ICOP,ICP1,ICP2,CP
        end do
      end do
    end do
  end do
end do
write(u6,*)
write(u6,*) ' CONVENTIONAL NR OF COUPLING COEFFS, BY PAIR:'
NRC = 0
do IP=2,SGS%MidLev
  do IQ=1,IP-1
    NRCPQ = 0
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(INDEO,LFTSYM,MV)*CIS%NOW(1,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
do IP=SGS%MidLev+1,SGS%nLev
  do IQ=1,SGS%MidLev
    NRCPQ = 0
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(IP,LFTSYM,MV)*EXS%NOCP(IQ,ISYM,MV)
        NRCPQ = NRCPQ+EXS%NOCP(SGS%nLev+IP,LFTSYM,MV)*EXS%NOCP(SGS%nLev+IQ,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
do IP=SGS%MidLev+2,SGS%nLev
  do IQ=SGS%MidLev+1,IP-1
    NRCPQ = 0
    INDEO = 2*SGS%nLev+(IP*(IP-1))/2+IQ
    do LFTSYM=1,SGS%nSym
      ISYM = LFTSYM
      do MV=1,CIS%nMidV
        NRCPQ = NRCPQ+EXS%NOCP(INDEO,LFTSYM,MV)*CIS%NOW(2,ISYM,MV)
      end do
    end do
    write(u6,'(1X,2I5,5X,I9)') IP,IQ,NRCPQ
    NRC = NRC+NRCPQ
  end do
end do
write(u6,*)
write(u6,*) ' TOTAL CONVENTIONAL COUPLING COEFFS:',NRC
#endif

call mma_allocate(EXS%VTab,nVTab_Final,Label='EXS%VTab')
EXS%VTab(1:nVTab_final) = VTab(1:nVTab_final)
call mma_deallocate(VTab)

end subroutine MKCOUP
