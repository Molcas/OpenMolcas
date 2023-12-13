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

subroutine AI(JTYP,ITAI,L0,L1,L2,L3)

use guga_global, only: COUP, IADD10, ICOUP, ICOUP1, ICH, IJ, IOUT, IRC, IWAY, IX, J1, J2, JNDX, JRC, LN, Lu_10, MXVERT, NBUF, NMAT
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JTYP, L0(*), L1(*), L2(*), L3(*)
integer(kind=iwp), intent(_OUT_) :: ITAI(*)
integer(kind=iwp) :: I, ICP1, ICP2, II, IID, IJJ, IJM, IJS, IN_, IN2, IND, ISTOP, IT1, IT2, ITAIL, ITT1, ITT2, ITURN, ITYP, JJ, &
                     JJD, JMAX, JND1, JND2, JOUT, KM, NI
real(kind=wp) :: CHKSUM
logical(kind=iwp) :: first

IOUT = 0
NMAT = 0
JMAX = 0
do NI=1,LN
  IOUT = IOUT+1
  ICOP1(IOUT) = 0
  JOUT = 0
  if (IOUT >= NBUF) then
    ICOP1(NCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
  end if
  I = ICH(NI)
  IOUT = IOUT+1
  ICOP1(IOUT) = I
  if (IOUT >= NBUF) then
    ICOP1(NCOP+1) = NBUF
    call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
    call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
    NMAT = NMAT+NBUF
    IOUT = 0
  end if
  IJS = IJ(I+1)+1
  IJM = IJ(I)
  if (JTYP == 1) then
    ITURN = 2
    ITT1 = 2
    ITT2 = 1
  else
    ! DOUBLET-VALENCE INTERACTIONS
    ITURN = 1
    ITT1 = 1
    ITT2 = 0
  end if
  do
    IT1 = ITT1*MXVERT
    IT2 = ITT2*MXVERT
    if (ITT1 == 0) then
      JJ = 0
      JJD = 0
    else
      JJ = IRC(ITT1)
      JJD = JRC(ITT1)
    end if
    if (ITT2 == 0) then
      II = 0
      IID = 0
    else
      II = IRC(ITT2)
      IID = JRC(ITT2)
    end if
    ITYP = ITURN
    do IJJ=IJS,IJM
      ITAIL = IX(IT2+IJJ)
      call TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
      IWAY(I) = 1
      do
        KM = I
        J2(KM+1) = IJJ
        J1(KM+1) = IJJ
        call LOOP1(KM,ISTOP,IT1,IT2)
        if (ISTOP == 1) exit
        first = .true.
        do
          if (first) then
            KM = KM-1
            if (KM /= 0) IWAY(KM) = 1
            first = .false.
          end if
          if (KM /= 0) then
            call LOOP5(KM,ISTOP,IT1,IT2)
            if (ISTOP == 0) then
              first = .true.
              cycle
            end if
            KM = KM+1
            if (KM == I) exit
          else
            do IN_=1,ITAIL
              ICP1 = ICOUP(1)+IN_
              JND1 = JNDX(II+ICP1)
              if (JND1 == 0) cycle
              ICP1 = JND1-IID
              IN2 = ITAI(IN_)
              if (IN2 == 0) cycle
              ICP2 = ICOUP1(1)+IN2
              JND2 = JNDX(JJ+ICP2)
              if (JND2 == 0) cycle
              ICP2 = JND2-JJD
              IOUT = IOUT+1
              JOUT = JOUT+1
              if (JOUT > JMAX) JMAX = JOUT
              COP(IOUT) = COUP(1)
              !PAM96 IND = ior(ITYP,ishft(ICP2,6))
              !IND=ITYP+2**6*ICP2
              IND = ior(ITYP,ishft(ICP2,6))
              !PAM96 ICOP1(IOUT) = ior(IND,ishft(ICP1,19))
              !ICOP1(IOUT)=IND+2**19*ICP1
              ICOP1(IOUT) = ior(IND,ishft(ICP1,19))
              if (IOUT < NBUF) cycle
              ICOP1(NCOP+1) = NBUF
              call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
              call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
              NMAT = NMAT+NBUF
              IOUT = 0
            end do
            KM = KM+1
            if (KM == I) exit
          end if
        end do
      end do
    end do
    select case (ITURN)
      case default !(1)
        ! TRIPLET-DOUBLET INTERACTIONS
        ITURN = 2
        ITT1 = 2
        ITT2 = 1
      case (2)
        ! SINGLET-DOUBLET INTERACTIONS
        ITURN = 3
        ITT1 = 3
        ITT2 = 1
      case (3)
        exit
    end select
  end do
end do
ICOP1(NCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(NCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
CHKSUM = Zero
do I=1,NCOP
  CHKSUM = CHKSUM+COP(I)
end do
call ADD_INFO('GUGA_CHKSUM',[CHKSUM],1,8)
if (JTYP == 0) write(u6,600) NMAT
if (JTYP == 0) then
  return
end if
write(u6,601) NMAT
IAD10(1) = JMAX
write(u6,602) JMAX

return

600 format(/,6X,'COEFFICIENTS FOR AI',I11)
601 format(/,6X,'COEFFICIENTS FOR ABCI',I9)
602 format(6X,'MAXIMUM NUMBER OF ELEMENTS',I6)

end subroutine AI
