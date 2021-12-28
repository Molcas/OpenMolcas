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

subroutine AIAI(BUFOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB,NBINS)

use guga_global, only: COUP, COUP1, IADD11, ICH, ICOUP, IJ, ILIM, IRC, IV0, IVF0, IWAY, IX, J1, J2, JM, JM1, JNDX, LN, LNP, Lu_11, &
                       MXVERT
use Definitions, only: wp, iwp, RtoI

implicit none
real(kind=wp), intent(inout) :: BUFOUT(*)
integer(kind=iwp), intent(inout) :: INDOUT(*), ICAD(*), IBUFL(*)
integer(kind=iwp), intent(in) :: KBUF, NTPB, NBINS
integer(kind=iwp) :: I, IAD110, ICP, ICPP, ICQ, IDIV, IFAI, IIJ, IJJ, IJM, IJS, IN_, IPOS, ISTOP, ISU, ISUM, IT1, IT2, ITAIL, ITT, &
                     IVL, JND1, KBUF0, KBUF1, KBUF2, KM, NBN, NI
real(kind=wp) :: CP
logical(kind=iwp) :: do_loop

KBUF0 = RtoI*KBUF
KBUF1 = KBUF0+KBUF+1
KBUF2 = KBUF1+1
IDIV = RtoI
do NI=1,LN
  I = ICH(NI)
  IIJ = (I*(I+1))/2
  IJS = IJ(I+1)+1
  IJM = IJ(I)
  do ITT=2,ILIM
    IT1 = (ITT-1)*MXVERT
    IT2 = IT1
    do IJJ=IJS,IJM
      ITAIL = IX(IT1+IJJ)
      IWAY(I) = 1
      loop_1: do
        KM = I
        J2(KM+1) = IJJ
        J1(KM+1) = IJJ
        JM(KM) = IVF0+1
        JM1(KM) = IVF0+1
        IFAI = 0
        call LOOP26(KM,ISTOP,IFAI,IT1,IT2)
        if (ISTOP == 1) exit
        if (J1(KM) /= J2(KM)) cycle
        do_loop = .true.
        do while (do_loop)
          KM = KM-1
          if (KM == 0) then
            IVL = J2(1)
            ISUM = IV0-IVL
            ISU = IRC(ISUM)
            if (ISUM == 3) CP = COUP(1)
            if (ISUM /= 3) CP = COUP1(1)
            do IN_=1,ITAIL
              JND1 = JNDX(ISU+ICOUP(1)+IN_)
              if (JND1 == 0) cycle
              IPOS = (JND1-1)*LNP+IIJ
              NBN = (IPOS-1)/NTPB+1
              IBUFL(NBN) = IBUFL(NBN)+1
              ICQ = ICAD(NBN)
              ICP = ICQ/IDIV+IBUFL(NBN)
              BUFOUT(ICP) = CP
              ICPP = ICQ+KBUF0+IBUFL(NBN)
              INDOUT(ICPP) = IPOS
              if (IBUFL(NBN) < KBUF) cycle
              INDOUT(ICQ+KBUF1) = KBUF
              IAD110 = IADD11
              call iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
              INDOUT(ICQ+KBUF2) = IAD110
              IBUFL(NBN) = 0
            end do
            if (I == 1) cycle loop_1
            KM = 1
          else
            IWAY(KM) = 1
          end if
          do
            JM(KM) = IVF0+1
            JM1(KM) = IVF0+1
            IFAI = 0
            call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
            if (ISTOP == 1) then
              KM = KM+1
              if (KM == I) cycle loop_1
            else if (J1(KM) == J2(KM)) then
              exit
            end if
          end do
        end do
        do_loop = .false.
      end do loop_1
    end do
  end do
end do
! EMPTY LAST BUFFERS
do I=1,NBINS
  ICQ = ICAD(I)
  INDOUT(ICQ+KBUF1) = IBUFL(I)
  IAD110 = IADD11
  call iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
  ICAD(I) = IAD110
end do

return

end subroutine AIAI
