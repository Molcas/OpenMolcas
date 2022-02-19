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
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine PRWF_CPF(ICASE,JSY,INDEX,C)

implicit real*8(A-H,O-Z)
dimension C(*), index(*), JSY(*)
dimension ICASE(*)
#include "SysDef.fh"
#include "cpfmcpf.fh"
dimension IOC(57), IORB(57), ISP(57), ILSYM(57)
external DNRM2_
! Statement functions
!PAM97      EXTERNAL UNPACK
!PAM97      INTEGER UNPACK
!PAM97      JO(L)=UNPACK(QOCC((L+29)/30), 2*L-(2*L-1)/60*60, 2)
!PAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
JO(L) = ICUNP(ICASE,L)
JSYM(L) = JSUNP_CPF(JSY,L)

NA = 0 ! dummy initialized
NB = 0 ! dummy initialized
ILIM = 4
if (IFIRST /= 0) ILIM = 2
NCONF = JSC(ILIM)
if (ISDCI == 1) call DSCAL_(NCONF,1.0d0/DNRM2_(NCONF,C,1),C,1)
JCONF = JSC(1)
THRC = CTRSH
write(6,5) THRC
call XFLUSH(6)
if (ISDCI == 0) write(6,6)

do J=1,LN
  IORB(J+2) = J
  ILSYM(J+2) = NSM(J)
end do

do I=1,NCONF
  JJ = I
  IJ = I
  if (I == IREF0) then
    write(6,105) I,C(I),'REFERENCE'
    call XFLUSH(6)
    GO TO 26
  end if
  CI = C(I)
  if (abs(CI) < THRC) GO TO 10
  if (I <= JCONF) then
    write(6,105) I,CI,'VALANCE'
    call XFLUSH(6)
    GO TO 26
  end if
  if (I <= JSC(2)) then
    JMIN = IRC(1)+1
    write(6,105) I,CI,'DOUBLET'
    call XFLUSH(6)
  else if (I <= JSC(3)) then
    JMIN = IRC(2)+1
  else
    JMIN = IRC(3)+1
  end if
  IX1 = IRC(ILIM)
  do J=JMIN,IX1
    JJ = J-1
    if (index(J) >= IJ) GO TO 25
  end do
25 continue
26 continue
  NSJ = MUL(JSYM(JJ),LSYM)
  JVIR = I-index(JJ)
  if (I > JCONF) JVIR = IJ-index(JJ)
  II1 = (JJ-1)*LN
  do II=1,LN
    II1 = II1+1
    ISP(II+2) = JO(II1)
    JOJ = ISP(II+2)
    if (JOJ > 1) JOJ = JOJ-1
    IOC(II+2) = JOJ
  end do
  if (JJ <= IRC(1)) then
    IORB(1) = 0
    IOC(1) = 0
    ISP(1) = 0
    ILSYM(1) = 0
    IORB(2) = 0
    IOC(2) = 0
    ISP(2) = 0
    ILSYM(2) = 0
    GO TO 100
  end if
  if (JJ <= IRC(2)) then
    IORB(2) = JVIR+NSYS(NSJ)+LN
    IOC(2) = 1
    ISP(2) = 1
    ILSYM(2) = NSJ
    IORB(1) = 0
    IOC(1) = 0
    ISP(1) = 0
    ILSYM(1) = 0
    GO TO 100
  end if
  IN = 0
  do II=1,NVIRT
    NA = II
    NSI = MUL(NSJ,NSM(LN+II))
    J1 = NSYS(NSI)+1
    J2 = NSYS(NSI+1)
    if (J2 > II) J2 = II
    if (J2 < J1) GO TO 46
    do J=J1,J2
      NB = J
      IN = IN+1
      if (IN == JVIR) GO TO 48
    end do
46  continue
  end do
48 continue
  IORB(1) = LN+NB
  IOC(1) = 1
  ISP(1) = 1
  ILSYM(1) = NSM(IORB(1))
  if (NA == NB) then
    IORB(2) = IORB(1)
    IOC(2) = 2
    ISP(2) = 3
    ILSYM(2) = NSM(IORB(2))
    IORB(1) = 0
    IOC(1) = 0
    ISP(1) = 0
    ILSYM(1) = 0
    CI = CI/SQ2
    if (abs(CI) < THRC) GO TO 10
  else
    IORB(2) = LN+NA
    IOC(2) = 1
    ISP(2) = 2
    ILSYM(2) = NSM(IORB(2))
  end if
  if (JJ <= IRC(3)) write(6,105) I,CI,'TRIPLET'
  if (JJ > IRC(3)) write(6,105) I,CI,'SINGLET'
100 continue
  write(6,*)
  call XFLUSH(6)
  if (LN+2 <= 36) then
    write(6,120) 'ORBITALS     ',(IORB(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,120) 'OCCUPATION   ',(IOC(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,120) 'SPIN-COUPLING',(ISP(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,120) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
    call XFLUSH(6)
  else
    write(6,121) 'ORBITALS     ',(IORB(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,121) 'OCCUPATION   ',(IOC(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,121) 'SPIN-COUPLING',(ISP(J),J=1,LN+2)
    call XFLUSH(6)
    write(6,121) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
    call XFLUSH(6)
  end if
10 continue
end do

return

5 format(//6X,'PRINTOUT OF CI-COEFFICIENTS LARGER THAN',F10.2,/)
6 format(/6X,'WAVE FUNCTION NOT NORMALIZED',/)
105 format(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,3X,A)
120 format(6X,A,36I3)
121 format(6X,A,55I2)

end subroutine PRWF_CPF
