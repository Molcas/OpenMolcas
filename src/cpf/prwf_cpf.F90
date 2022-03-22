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

subroutine PRWF_CPF(ICASE,JSY,INDX,C)

use cpf_global, only: CTRSH, ILIM, IRC, IREF0, ISDCI, JSC, LN, LSYM, NCONF, NSM, NSYS, NVIRT, SQ2
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), INDX(*), JSY(*)
real(kind=wp), intent(inout) :: C(*)
integer(kind=iwp) :: I, II, II1, IIN, IJ, ILSYM(57), IOC(57), IORB(57), ISP(57), IX1, J, J1, J2, JCONF, JJ, JMIN, JOJ, JVIR, NA, &
                     NB, NSI, NSJ
real(kind=wp) :: CI, CNRM, THRC
integer(kind=iwp), external :: ICUNP, JSUNP
real(kind=wp), external :: DNRM2_

NA = 0 ! dummy initialized
NB = 0 ! dummy initialized
NCONF = JSC(ILIM)
if (ISDCI == 1) then
  CNRM = DNRM2_(NCONF,C,1)
  C(1:NCONF) = C(1:NCONF)/CNRM
end if
JCONF = JSC(1)
THRC = CTRSH
write(u6,5) THRC
if (ISDCI == 0) write(u6,6)

do J=1,LN
  IORB(J+2) = J
  ILSYM(J+2) = NSM(J)
end do

do I=1,NCONF
  JJ = I
  IJ = I
  if (I == IREF0) then
    write(u6,105) I,C(I),'REFERENCE'
  else
    CI = C(I)
    if (abs(CI) < THRC) cycle
    if (I <= JCONF) then
      write(u6,105) I,CI,'VALANCE'
    else
      if (I <= JSC(2)) then
        JMIN = IRC(1)+1
        write(u6,105) I,CI,'DOUBLET'
      else if (I <= JSC(3)) then
        JMIN = IRC(2)+1
      else
        JMIN = IRC(3)+1
      end if
      IX1 = IRC(ILIM)
      do J=JMIN,IX1
        JJ = J-1
        if (INDX(J) >= IJ) exit
      end do
    end if
  end if
  NSJ = MUL(JSUNP(JSY,JJ),LSYM)
  JVIR = I-INDX(JJ)
  if (I > JCONF) JVIR = IJ-INDX(JJ)
  II1 = (JJ-1)*LN
  do II=1,LN
    II1 = II1+1
    ISP(II+2) = ICUNP(ICASE,II1)
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
  else if (JJ <= IRC(2)) then
    IORB(2) = JVIR+NSYS(NSJ)+LN
    IOC(2) = 1
    ISP(2) = 1
    ILSYM(2) = NSJ
    IORB(1) = 0
    IOC(1) = 0
    ISP(1) = 0
    ILSYM(1) = 0
  else
    IIN = 0
    outer: do II=1,NVIRT
      NA = II
      NSI = MUL(NSJ,NSM(LN+II))
      J1 = NSYS(NSI)+1
      J2 = NSYS(NSI+1)
      if (J2 > II) J2 = II
      do J=J1,J2
        NB = J
        IIN = IIN+1
        if (IIN == JVIR) exit outer
      end do
    end do outer
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
      if (abs(CI) < THRC) cycle
    else
      IORB(2) = LN+NA
      IOC(2) = 1
      ISP(2) = 2
      ILSYM(2) = NSM(IORB(2))
    end if
    if (JJ <= IRC(3)) write(u6,105) I,CI,'TRIPLET'
    if (JJ > IRC(3)) write(u6,105) I,CI,'SINGLET'
  end if
  write(u6,*)
  if (LN+2 <= 36) then
    write(u6,120) 'ORBITALS     ',(IORB(J),J=1,LN+2)
    write(u6,120) 'OCCUPATION   ',(IOC(J),J=1,LN+2)
    write(u6,120) 'SPIN-COUPLING',(ISP(J),J=1,LN+2)
    write(u6,120) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
  else
    write(u6,121) 'ORBITALS     ',(IORB(J),J=1,LN+2)
    write(u6,121) 'OCCUPATION   ',(IOC(J),J=1,LN+2)
    write(u6,121) 'SPIN-COUPLING',(ISP(J),J=1,LN+2)
    write(u6,121) 'SYMMETRY     ',(ILSYM(J),J=1,LN+2)
  end if
end do

return

5 format(//6X,'PRINTOUT OF CI-COEFFICIENTS LARGER THAN',F10.2,/)
6 format(/6X,'WAVE FUNCTION NOT NORMALIZED',/)
105 format(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,3X,A)
120 format(6X,A,36I3)
121 format(6X,A,55I2)

end subroutine PRWF_CPF
