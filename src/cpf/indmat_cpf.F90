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

subroutine INDMAT_CPF(JSY,INDX,ISAB,ISMAX,JREFX)

use cpf_global, only: IFIRST, ILIM, IPRINT, IRC, IREF0, ISC, JJS, JSC, LN, LSYM, NDIAG, NNS, NSM, NSYM, NSYS, NVIR, NVIRT
use Symmetry_Info, only: Mul
use Definitions, only: iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), JREFX(*)
integer(kind=iwp), intent(_OUT_) :: INDX(*), ISAB(*), ISMAX
integer(kind=iwp) :: I, ICOUS(8), II, IIN, IN0, IN2, IND, IR, IR1, IR2, IX1, IX2, IX3, IX4, JCONF, JJM, JSCI, NA, NB, NSAB, NSS
integer(kind=iwp), external :: JSUNP

! DETERMINE REFERENCE STATE
JCONF = ISC(1)
do IR=1,JCONF
  if (JREFX(IR) == 1) IREF0 = IR
end do
if (IPRINT > 5) write(u6,999) IREF0,(JREFX(IR),IR=1,JCONF)

NSYS(1) = 0
do I=2,NSYM
  NSYS(I) = NSYS(I-1)+NVIR(I-1)
end do
NSYS(NSYM+1) = NVIRT
do I=1,NSYM
  ICOUS(I) = 0
  NNS(I) = 0
end do
ISMAX = 0
IN0 = -NVIRT
do NA=1,NVIRT
  IN0 = IN0+NVIRT
  IIN = IN0
  IN2 = -NVIRT+NA
  do NB=1,NA
    IIN = IIN+1
    IN2 = IN2+NVIRT
    NSAB = MUL(NSM(LN+NA),NSM(LN+NB))
    ICOUS(NSAB) = ICOUS(NSAB)+1
    ISAB(IIN) = ICOUS(NSAB)
    if (ISMAX < ISAB(IIN)) ISMAX = ISAB(IIN)
    ISAB(IN2) = ISAB(IIN)
    if (ISAB(IIN) > NNS(NSAB)) NNS(NSAB) = ISAB(IIN)
  end do
  NDIAG(NA) = ISAB(IIN)
end do
IND = 0
IR = IRC(1)
do II=1,IR
  IND = IND+1
  INDX(II) = IND
end do
JSC(1) = IND
IR1 = IR+1
IR2 = IRC(2)
do II=IR1,IR2
  INDX(II) = IND
  NSS = MUL(JSUNP(JSY,II),LSYM)
  IND = IND+NVIR(NSS)
end do
JSC(2) = IND
if (IFIRST == 0) then
  IR1 = IR2+1
  IR2 = IRC(4)
  JSC(3) = IND
  do II=IR1,IR2
    INDX(II) = IND
    NSS = MUL(JSUNP(JSY,II),LSYM)
    IND = IND+ICOUS(NSS)
    if (II == IRC(3)) JSC(3) = IND
  end do
  JSC(4) = IND
end if
IX1 = JSC(1)
IX2 = JSC(2)-JSC(1)
write(u6,213)
if (IFIRST == 0) then
  JJM = (JJS(LSYM+1)-JJS(LSYM))*NVIRT
  IX3 = JSC(3)-JSC(2)-JJM
  IX4 = JSC(4)-JSC(3)
  write(u6,215) IX1,IX2,IX3,IX4
else
  write(u6,216) IX1,IX2
  JJM = 0
end if
JSCI = JSC(ILIM)-JJM
write(u6,50) ISC(ILIM),JSCI

return

50 format(//6X,'FORMAL NUMBER OF CONFIGURATIONS',I8,/8X,'REAL NUMBER OF CONFIGURATIONS',I8)
213 format(//,6X,'FULL-SPACE CONFIGURATIONS (REAL)')
215 format(/,6X,'NUMBER OF VALENCE STATES',I16,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7, &
           /,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,/,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
216 format(/,6X,'NUMBER OF VALENCE STATES',I14,/,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
999 format(2X,I3,2X,'JREFX',10I5)

end subroutine INDMAT_CPF
