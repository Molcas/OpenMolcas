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

subroutine AI(JTYP,ITAI,L0,L1,L2,L3)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
dimension ITAI(*), L0(*), L1(*), L2(*), L3(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"

IOUT = 0
NMAT = 0
JMAX = 0
do NI=1,LN
  IOUT = IOUT+1
  ICOP1(IOUT) = 0
  JOUT = 0
  if (IOUT < NBUF) GO TO 460
  ICOP1(NCOP+1) = NBUF
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
460 I = ICH(NI)
  IOUT = IOUT+1
  ICOP1(IOUT) = I
  if (IOUT < NBUF) GO TO 11
  ICOP1(NCOP+1) = NBUF
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
11 IJS = IJ(I+1)+1
  IJM = IJ(I)
  if (JTYP == 1) GO TO 101
  ! DOUBLET-VALENCE INTERACTIONS
  ITURN = 1
  ITT1 = 1
  ITT2 = 0
150 IT1 = ITT1*MXVERT
  IT2 = ITT2*MXVERT
  JJ = 0
  if (ITT1 /= 0) JJ = IRC(ITT1)
  JJD = 0
  if (ITT1 /= 0) JJD = JRC(ITT1)
  II = 0
  if (ITT2 /= 0) II = IRC(ITT2)
  IID = 0
  if (ITT2 /= 0) IID = JRC(ITT2)
  ITYP = ITURN
  do IJJ=IJS,IJM
    ITAIL = IX(IT2+IJJ)
    call TAIL(I,IJJ,ITAI,ITAIL,L0,L1,L2,L3,IT1,IT2)
    IWAY(I) = 1
32  KM = I
    J2(KM+1) = IJJ
    J1(KM+1) = IJJ
    call LOOP1(KM,ISTOP,IT1,IT2)
    if (ISTOP == 1) GO TO 30
41  KM = KM-1
    if (KM == 0) GO TO 51
    IWAY(KM) = 1
42  call LOOP5(KM,ISTOP,IT1,IT2)
    if (ISTOP == 0) GO TO 41
53  KM = KM+1
    if (KM == I) GO TO 32
    GO TO 42
51  do IN=1,ITAIL
      ICP1 = ICOUP(1)+IN
      JND1 = JNDX(II+ICP1)
      if (JND1 == 0) GO TO 80
      ICP1 = JND1-IID
      IN2 = ITAI(IN)
      if (IN2 == 0) GO TO 80
      ICP2 = ICOUP1(1)+IN2
      JND2 = JNDX(JJ+ICP2)
      if (JND2 == 0) GO TO 80
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
      if (IOUT < NBUF) GO TO 80
      ICOP1(NCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
80    continue
    end do
    GO TO 53
30  continue
  end do
  GO TO(101,102,10),ITURN
  ! TRIPLET-DOUBLET INTERACTIONS
101 ITURN = 2
  ITT1 = 2
  ITT2 = 1
  GO TO 150
  ! SINGLET-DOUBLET INTERACTIONS
102 ITURN = 3
  ITT1 = 3
  ITT2 = 1
  GO TO 150
10 continue
end do
ICOP1(NCOP+1) = IOUT
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
NMAT = NMAT+IOUT
ICOP1(NCOP+1) = -1
call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
CHKSUM = 0.0d0
do I=1,NCOP
  CHKSUM = CHKSUM+COP(I)
end do
call ADD_INFO('GUGA_CHKSUM',[CHKSUM],1,8)
if (JTYP == 0) write(IW,600) NMAT
600 format(/,6X,'COEFFICIENTS FOR AI',I11)
if (JTYP == 0) then
  return
end if
write(IW,601) NMAT
601 format(/,6X,'COEFFICIENTS FOR ABCI',I9)
IAD10(1) = JMAX
write(IW,602) JMAX
602 format(6X,'MAXIMUM NUMBER OF ELEMENTS',I6)

return

end subroutine AI
