!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine PRWF_MRCI(ICSPCK,INTSYM,INDX,C,JREFX)

implicit real*8(A-H,O-Z)
dimension C(*), INDX(*), ICSPCK(*), INTSYM(*), JREFX(*)
character*12 CSFTYP
character*14 FORM0, FORM1, FORM2, FORM
character*14 FORM00, FORM01, FORM02
#include "SysDef.fh"
#include "mrci.fh"
dimension IOC(32), IORBI(32), ISP(32), ILSYM(32)
!PAM97 external UNPACK
!PAM97 integer UNPACK
data FORM00/'(1X,A,6X,53I2)'/
data FORM01/'(1X,A,3X,54I2)'/
data FORM02/'(1X,A,55I2)   '/
data FORM0/'(1X,A,6X,34I3)'/
data FORM1/'(1X,A,3X,35I3)'/
data FORM2/'(1X,A,36I3)   '/
! STATEMENT FUNCTIONS FOR RETRIEVING GUGA CASE NUMBERS AND INTERNAL
! SYMMETRY LABEL:
!PAM97 JO(L) = UNPACK(CSPCK((L+29)/30),2*L-(2*L-1)/60*60,2)
JO(L) = ICUNP(ICSPCK,L)
!PAM96 JSYM(L) = UNPACK(INTSYM((L+9)/10),3*mod(L-1,10)+1,3)+1
JSYM(L) = JSUNP(INTSYM,L)

NA = 0
NB = 0
ILIM = 4
LN2 = LN+2
if (IFIRST /= 0) ILIM = 2
SCALE = 1.0/sqrt(DDOT_(NCONF,C,1,C,1))
call DSCAL_(NCONF,SCALE,C,1)
do J=1,LN
  IORBI(J+2) = IORB(J)
  ILSYM(J+2) = NSM(J)
end do
JCONF = JSC(1)
write(6,'(A,F5.3)') '      CI-COEFFICIENTS LARGER THAN ',CTRSH
call XFLUSH(6)
do IS=1,NSYM
  if (NFMO(IS) > 0) then
    write(6,*) ' NOTE: THE FOLLOWING ORBITALS WERE FROZEN'
    call XFLUSH(6)
    write(6,*) ' ALREADY AT THE INTEGRAL TRANSFORMATION STEP'
    call XFLUSH(6)
    write(6,*) ' AND DO NOT EXPLICITLY APPEAR:'
    call XFLUSH(6)
    write(6,'(6X,A,8I4)') '  SYMMETRY:',(I,I=1,NSYM)
    call XFLUSH(6)
    write(6,'(6X,A,8I4)') 'PRE-FROZEN:',(NFMO(I),I=1,NSYM)
    call XFLUSH(6)
    goto 6
  end if
end do
6 continue
write(6,*) ' ORDER OF SPIN-COUPLING: (PRE-FROZEN, NOT SHOWN)'
call XFLUSH(6)
write(6,*) '                         (FROZEN, NOT SHOWN)'
call XFLUSH(6)
write(6,*) '                          VIRTUAL'
call XFLUSH(6)
write(6,*) '                          ADDED VALENCE'
call XFLUSH(6)
write(6,*) '                          INACTIVE'
call XFLUSH(6)
write(6,*) '                          ACTIVE'
call XFLUSH(6)
write(6,*)
call XFLUSH(6)
write(6,*) ' ORBITALS ARE NUMBERED WITHIN EACH SEPARATE SYMMETRY.'
call XFLUSH(6)
do I=1,NCONF
  CI = C(I)
  ACI = abs(C(I))
  if (I <= JCONF) then
    JMIN = 1
    NREXT = 0
    if (JREFX(I) /= 0) then
      CSFTYP = '   REFERENCE'
      CLIM = 0.0d00
    else
      CSFTYP = '     VALENCE'
      CLIM = CTRSH
    end if
  else if (I <= JSC(2)) then
    NREXT = 1
    JMIN = 1+IRC(1)
    CSFTYP = '     DOUBLET'
    CLIM = CTRSH
  else if (I <= JSC(3)) then
    NREXT = 2
    JMIN = 1+IRC(2)
    CSFTYP = '     TRIPLET'
    CLIM = CTRSH
  else
    NREXT = 2
    JMIN = 1+IRC(3)
    CSFTYP = '     SINGLET'
    CLIM = CTRSH
  end if
  if (ACI < CLIM) goto 10
  JJ = I
  IJ = I
  if (NREXT > 0) then
    JMAX = IRC(ILIM)
    do J=JMIN,JMAX
      JJ = J
      if (INDX(J) < IJ) GO TO 20
      JJ = JJ-1
      goto 25
20    continue
    end do
25  continue
  end if
  NSJ = MUL(JSYM(JJ),LSYM)
  JVIR = IJ-INDX(JJ)
  II1 = (JJ-1)*LN
  do II=1,LN
    II1 = II1+1
    ISP(II+2) = JO(II1)
    IOC(II+2) = (1+ISP(II+2))/2
  end do
  if (NREXT == 0) then
    FORM = FORM0
    if (LN2 > 36) FORM = FORM00
    LN1 = 3
  else if (NREXT == 1) then
    IO = JVIR+NVIRP(NSJ)+LN
    IORBI(2) = IORB(IO)
    IOC(2) = 1
    ISP(2) = 1
    ILSYM(2) = NSJ
    FORM = FORM1
    if (LN2 > 36) FORM = FORM01
    LN1 = 2
  else
    IN = 0
    do II=1,NVIRT
      NA = II
      NSI = MUL(NSJ,NSM(LN+II))
      J1 = NVIRP(NSI)+1
      J2 = NVIRP(NSI)+NVIR(NSI)
      if (J2 > II) J2 = II
      if (J2 < J1) GO TO 46
      do J=J1,J2
        NB = J
        IN = IN+1
        if (IN == JVIR) GO TO 48
      end do
46    continue
    end do
48  continue
    IOC(1) = 1
    ISP(1) = 1
    ILSYM(1) = NSM(LN+NB)
    IO = LN+NB
    IORBI(1) = IORB(IO)
    if (NA == NB) then
      IORBI(2) = IORBI(1)
      IOC(2) = 2
      ISP(2) = 3
      ILSYM(2) = NSM(IO)
      CI = CI*sqrt(0.5d00)
      if (abs(CI) < CTRSH) GO TO 10
      FORM = FORM1
      if (LN2 > 36) FORM = FORM01
      LN1 = 2
    else
      IOC(2) = 1
      ISP(2) = 2
      if (CSFTYP == '     TRIPLET') ISP(2) = 1
      IO = LN+NA
      IORBI(2) = IORB(IO)
      ILSYM(2) = NSM(IO)
      FORM = FORM2
      if (LN2 > 36) FORM = FORM02
      LN1 = 1
    end if
  end if
  write(6,*)
  call XFLUSH(6)
  write(6,105) I,C(I),CSFTYP
  call XFLUSH(6)
105 format(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,A)
  write(6,FORM) 'SYMMETRY     ',(ILSYM(J),J=LN1,LN2)
  call XFLUSH(6)
  write(6,FORM) 'ORBITALS     ',(IORBI(J),J=LN1,LN2)
  call XFLUSH(6)
  write(6,FORM) 'OCCUPATION   ',(IOC(J),J=LN1,LN2)
  call XFLUSH(6)
  write(6,FORM) 'SPIN-COUPLING',(ISP(J),J=LN1,LN2)
  call XFLUSH(6)
10 continue
end do

return

end subroutine PRWF_MRCI
