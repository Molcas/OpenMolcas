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
!***********************************************************************

subroutine AIAI(BUFOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB,NBINS)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
dimension BUFOUT(*), INDOUT(*), ICAD(*), IBUFL(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"

KBUF0 = RTOI*KBUF
KBUF1 = KBUF0+KBUF+1
KBUF2 = KBUF1+1
IDIV = RTOI
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
32    KM = I
      J2(KM+1) = IJJ
      J1(KM+1) = IJJ
      JM(KM) = IVF0+1
      JM1(KM) = IVF0+1
      IFAI = 0
      call LOOP26(KM,ISTOP,IFAI,IT1,IT2)
      if (ISTOP == 1) GO TO 30
      if (J1(KM) /= J2(KM)) GO TO 32
41    KM = KM-1
      if (KM == 0) GO TO 51
      IWAY(KM) = 1
42    JM(KM) = IVF0+1
      JM1(KM) = IVF0+1
      IFAI = 0
      call LOOP23(KM,ISTOP,IFAI,IT1,IT2)
      if (ISTOP == 1) GO TO 43
      if (J1(KM) /= J2(KM)) GO TO 42
      GO TO 41
43    KM = KM+1
      if (KM == I) GO TO 32
      GO TO 42
51    IVL = J2(1)
      ISUM = IV0-IVL
      ISU = IRC(ISUM)
      if (ISUM == 3) CP = COUP(1)
      if (ISUM /= 3) CP = COUP1(1)
      do IN=1,ITAIL
        JND1 = JNDX(ISU+ICOUP(1)+IN)
        if (JND1 == 0) GO TO 62
        IPOS = (JND1-1)*LNP+IIJ
        NBN = (IPOS-1)/NTPB+1
        IBUFL(NBN) = IBUFL(NBN)+1
        ICQ = ICAD(NBN)
        ICP = ICQ/IDIV+IBUFL(NBN)
        BUFOUT(ICP) = CP
        ICPP = ICQ+KBUF0+IBUFL(NBN)
        INDOUT(ICPP) = IPOS
        if (IBUFL(NBN) < KBUF) GO TO 62
        INDOUT(ICQ+KBUF1) = KBUF
        IAD110 = IADD11
        call iDAFILE(Lu_11,1,INDOUT(ICQ+1),KBUF2,IADD11)
        INDOUT(ICQ+KBUF2) = IAD110
        IBUFL(NBN) = 0
62      continue
      end do
      if (I == 1) GO TO 32
      KM = 1
      GO TO 42
30    continue
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
