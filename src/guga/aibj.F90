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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine AIBJ(L0,L1,L2,L3,ITAI)

use guga_global, only: IADD10, BS3, BS4, COUP, COUP1, IA, IB, ICASE, ICH, ICOUP, ICOUP1, IJ, IOUT, IRC, IVF0, IWAY, IX, J1, J2, &
                       JM, JM1, JNDX, JRC, LN, Lu_10, MXVERT, NBUF, NMAT
use guga_util_global, only: COP, ICOP1, nCOP
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: I, IABIJ, IC1, IC11, IC2, IC22, ICP1, ICP2, IDIF, IFAB, IFAI, II, IID, IJJ, IJM, IJS, IN_, IN2, IND1, IND2, &
                     IND3, ISTOP, IT1, IT2, ITAIL, ITT1, ITT2, ITURN, ITYP, J, JJ, JJ1, JJD, JND1, JND2, JOJ, JTURN, KM, KM1, &
                     LTYP, NI, NJ, NUMM(7)
real(kind=wp) :: COPL, COPLA, COPLA0
logical(kind=iwp) :: first1, first2, skip1, skip2
integer(kind=iwp), external :: ICUNP

IC1 = 0      ! dummy initialize
IC2 = 0      ! dummy initialize
COPLA = Zero ! dummy initialize
do I=1,7
  NUMM(I) = 0
end do
IOUT = 0
NMAT = 0
do NI=1,LN
  do NJ=1,NI
    I = ICH(NI)
    J = ICH(NJ)
    if (I <= J) then
      I = ICH(NJ)
      J = ICH(NI)
    end if
    LTYP = 0
    IOUT = IOUT+1
    ICOP1(IOUT) = 0
    if (IOUT >= NBUF) then
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end if
    IOUT = IOUT+1
    !PAM96 ICOP1(IOUT) = ior(I,ishft(J,10))
    !ICOP1(IOUT) = I+2**10*J
    ICOP1(IOUT) = ior(I,ishft(J,10))
    if (IOUT >= NBUF) then
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end if
    IJS = IJ(I+1)+1
    IJM = IJ(I)
    ! FIRST ORDER INTERACTION
    ! TRIPLET-VALENCE INTERACTIONS
    JTURN = 1
    ITURN = 0
    ITT1 = 2
    ITT2 = 0
    do
      IT1 = ITT1*MXVERT
      !ulf IT2 = ITT2*300
      IT2 = ITT2*MXVERT
      if (ITT2 == 0) then
        II = 0
        IID = 0
      else
        II = IRC(ITT2)
        IID = JRC(ITT2)
      end if
      JJ = IRC(ITT1)
      JJD = JRC(ITT1)
      ITYP = JTURN
      if (ITURN /= 0) ITYP = JTURN-2
      do IJJ=IJS,IJM
        ITAIL = IX(IT2+IJJ)
        if (IT1 /= IT2) call TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
        IWAY(I) = 1
        loop_1: do
          KM = I
          J2(KM+1) = IJJ
          J1(KM+1) = IJJ
          first1 = .true.
          if (I == J) then
            skip1 = .true.
          else
            skip1 = .false.
            call LOOP1(KM,ISTOP,IT1,IT2)
            if (ISTOP == 1) exit loop_1
          end if
          loop_2: do
            if (first1 .and. (.not. skip1)) then
              KM = KM-1
              IWAY(KM) = 1
              if (KM == J) skip1 = .true.
            end if
            first1 = .false.
            if (.not. skip1) then
              call LOOP5(KM,ISTOP,IT1,IT2)
              if (ISTOP == 0) then
                first1 = .true.
                cycle loop_2
              end if
              KM = KM+1
              if (KM == I) cycle loop_1
              cycle loop_2
            end if
            skip1 = .false.
            IWAY(J) = 1
            do
              KM = J
              JM(KM) = IVF0+1
              JM1(KM) = IVF0+1
              IABIJ = 0
              if (I /= J) then
                if (ITURN == 0) call LOOP10(KM,ISTOP,IT1,IT2)
                IFAI = 1
                if (ITURN == 1) call LOOP13(KM,ISTOP,IFAI,IT1,IT2)
                if (ISTOP == 1) then
                  if (I == J) exit
                  KM = KM+1
                  if (KM == I) cycle loop_1
                  cycle loop_2
                else if (J1(KM) == J2(KM)) then
                  IABIJ = 1
                  IC1 = ICOUP(KM)
                  IC2 = ICOUP1(KM)
                  KM1 = KM+1
                  IDIF = IA(J1(KM1))-IA(J2(KM1))
                  if (IWAY(J) == 2) COPLA = COUP(J+1)
                  if ((IDIF == 0) .and. (IWAY(J) == 5)) COPLA = COUP(J+1)*BS4(IB(J2(J+1))+1)
                  if ((IDIF == 1) .and. (IWAY(J) == 4)) COPLA = COUP(J+1)*BS3(IB(J2(J+1))+1)
                end if
              else
                if (ITURN == 0) call LOOP7(KM,ISTOP,IT1,IT2)
                IFAI = 1
                if (ITURN == 1) call LOOP26(KM,ISTOP,IFAI,IT1,IT2)
                LTYP = 1
                if (IWAY(I) == 5) LTYP = 0
                if (ISTOP == 1) then
                  if (I == J) exit
                  KM = KM+1
                  if (KM == I) cycle loop_1
                  cycle loop_2
                end if
              end if
              first2 = .true.
              do
                if (first2) then
                  KM = KM-1
                  if (KM /= 0) IWAY(KM) = 1
                end if
                first2 = .false.
                if (KM /= 0) then
                  JM(KM) = IVF0+1
                  JM1(KM) = IVF0+1
                  if (ITURN == 0) call LOOP17(KM,ISTOP,IT1,IT2)
                  IFAI = 0
                  if (J1(KM+1) == J2(KM+1)) then
                    if (I == J) then
                      if (ICOUP(KM+1) == ICOUP1(KM+1)) IFAI = 1
                    else if (IABIJ /= 0) then
                      IC11 = ICOUP(KM+1)-IC1
                      IC22 = ICOUP1(KM+1)-IC2
                      if (IC11 == IC22) IFAI = 1
                    end if
                  end if
                  if (ITURN == 1) call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
                  if (ISTOP == 0) then
                    first2 = .true.
                    cycle
                  end if
                else
                  COPL = COUP(1)
                  if ((JTURN >= 6) .or. (JTURN == 3)) COPL = COUP1(1)
                  IFAB = 0
                  skip2 = .false.
                  if (JTURN == 5) then
                    if ((I == J) .and. (LTYP == 1)) skip2 = .true.
                  end if
                  if ((.not. skip2) .and. (ITT1 == ITT2)) then
                    if (I == J) then
                      if ((LTYP == 1) .and. (ICOUP(1) > ICOUP1(1))) skip2 = .true.
                    else if (IABIJ /= 0) then
                      IC11 = ICOUP(1)-IC1
                      IC22 = ICOUP1(1)-IC2
                      if (IC11 == IC22) IFAB = 1
                    end if
                  end if
                  if (.not. skip2) then
                    COPLA0 = COPLA
                    if (IFAB == 0) COPLA0 = Zero
                    do IN_=1,ITAIL
                      ICP1 = ICOUP(1)+IN_
                      JND1 = JNDX(II+ICP1)
                      if (JND1 == 0) cycle
                      ICP1 = JND1-IID
                      if (ITT1 == ITT2) then
                        IN2 = IN_
                      else
                        IN2 = ITAI(IN_)
                        if (IN2 == 0) cycle
                      end if
                      ICP2 = ICOUP1(1)+IN2
                      JND2 = JNDX(JJ+ICP2)
                      if (JND2 == 0) cycle
                      ICP2 = JND2-JJD
                      if ((ITT1 == ITT2) .and. (ICP1 == ICP2)) then
                        JJ1 = (JJD+ICP1-1)*LN+I
                        JOJ = ICUNP(ICASE,JJ1)
                        if (JOJ > 1) JOJ = JOJ-1
                        COPLA0 = JOJ-2
                        IFAB = 1
                      end if
                      IOUT = IOUT+1
                      NUMM(JTURN) = NUMM(JTURN)+1
                      COP(IOUT) = COPL
                      !PAM96 IND1 = ior(IFAB,ishft(ITURN,1))
                      !PAM96 IND2 = ior(IND1,ishft(ITYP,2))
                      !PAM96 IND3 = ior(IND2,ishft(ICP1,5))
                      !PAM96 ICOP1(IOUT) = ior(IND3,ishft(ICP2,18))
                      !IND1 = IFAB+2**1*ITURN
                      !IND2 = IND1+2**2*ITYP
                      !IND3 = IND2+2**5*ICP1
                      !ICOP1(IOUT) = IND3+2**18*ICP2
                      IND1 = ior(IFAB,ishft(ITURN,1))
                      IND2 = ior(IND1,ishft(ITYP,2))
                      IND3 = ior(IND2,ishft(ICP1,5))
                      ICOP1(IOUT) = ior(IND3,ishft(ICP2,18))
                      if (IOUT >= NBUF) then
                        ICOP1(nCOP+1) = NBUF
                        call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
                        call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
                        NMAT = NMAT+NBUF
                        IOUT = 0
                      end if
                      if (IFAB == 0) cycle
                      IOUT = IOUT+1
                      NUMM(JTURN) = NUMM(JTURN)+1
                      COP(IOUT) = COPLA0
                      ICOP1(IOUT) = 1
                      if (IOUT < NBUF) cycle
                      ICOP1(nCOP+1) = NBUF
                      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
                      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
                      NMAT = NMAT+NBUF
                      IOUT = 0
                    end do
                  end if
                end if
                KM = KM+1
                if (KM == J) exit
              end do
            end do
            if (.true.) exit loop_2
          end do loop_2
          if (.true.) exit loop_1
        end do loop_1
      end do
      select case (JTURN)
        case default !(1)
          ! SINGLET-VALENCE INTERACTIONS
          JTURN = 2
          ITT1 = 3
        case (2)
          ! TRIPLET-TRIPLET INTERACTIONS
          JTURN = 3
          ITURN = 1
          ITT1 = 2
          ITT2 = 2
        case (3)
          ! SINGLET-SINGLET INTERACTIONS
          JTURN = 4
          ITT1 = 3
          ITT2 = 3
        case (4)
          ! TRIPLET-SINGLET INTERACTIONS
          JTURN = 5
          ITT1 = 2
          ITT2 = 3
        case (5)
          ! SINGLET-TRIPLET INTERACTIONS
          JTURN = 6
          ITT1 = 3
          ITT2 = 2
        case (6)
          ! DOUBLET-DOUBLET INTERACTIONS
          JTURN = 7
          ITT1 = 1
          ITT2 = 1
        case (7)
          exit
      end select
    end do
  end do
end do
ICOP1(nCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(nCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
write(u6,600) NMAT
write(u6,610) (NUMM(I),I=1,7)

return

600 format(/,6X,'COEFFICIENTS FOR AIBJ',I9)
610 format(6X,'DIFFERENT TYPES',7I9)

end subroutine AIBJ
