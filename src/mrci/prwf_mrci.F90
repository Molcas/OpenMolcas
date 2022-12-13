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

use mrci_global, only: CTRSH, IFIRST, IORB, IRC, JSC, LN, LSYM, NCONF, NFMO, NSM, NSYM, NVIR, NVIRP, NVIRT
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ICSPCK(*), INTSYM(*), INDX(*), JREFX(*)
real(kind=wp), intent(inout) :: C(*)
integer(kind=iwp) :: I, II, II1, IIN, IJ, ILIM, ILSYM(32), IO, IOC(32), IORBI(32), IS, ISP(32), J, J1, J2, JCONF, JJ, JMAX, JMIN, &
                     JVIR, LN1, LN2, NA, NB, NREXT, NSI, NSJ
real(kind=wp) :: ACI, CI, CLIM, SCL
character(len=14) :: FRMT
character(len=12) :: CSFTYP
character(len=*), parameter :: FORM00 = '(1X,A,6X,53I2)', FORM01 = '(1X,A,3X,54I2)', FORM02 = '(1X,A,55I2)', &
                               FORM0 = '(1X,A,6X,34I3)', FORM1 = '(1X,A,3X,35I3)', FORM2 = '(1X,A,36I3)'
integer(kind=iwp), external :: ICUNP, JSUNP
real(kind=wp), external :: DDOT_

NA = 0
NB = 0
ILIM = 4
LN2 = LN+2
if (IFIRST /= 0) ILIM = 2
SCL = One/sqrt(DDOT_(NCONF,C,1,C,1))
C(1:NCONF) = SCL*C(1:NCONF)
do J=1,LN
  IORBI(J+2) = IORB(J)
  ILSYM(J+2) = NSM(J)
end do
JCONF = JSC(1)
write(u6,'(A,F5.3)') '      CI-COEFFICIENTS LARGER THAN ',CTRSH
do IS=1,NSYM
  if (NFMO(IS) > 0) then
    write(u6,*) ' NOTE: THE FOLLOWING ORBITALS WERE FROZEN'
    write(u6,*) ' ALREADY AT THE INTEGRAL TRANSFORMATION STEP'
    write(u6,*) ' AND DO NOT EXPLICITLY APPEAR:'
    write(u6,'(6X,A,8I4)') '  SYMMETRY:',(I,I=1,NSYM)
    write(u6,'(6X,A,8I4)') 'PRE-FROZEN:',(NFMO(I),I=1,NSYM)
    exit
  end if
end do
write(u6,*) ' ORDER OF SPIN-COUPLING: (PRE-FROZEN, NOT SHOWN)'
write(u6,*) '                         (FROZEN, NOT SHOWN)'
write(u6,*) '                          VIRTUAL'
write(u6,*) '                          ADDED VALENCE'
write(u6,*) '                          INACTIVE'
write(u6,*) '                          ACTIVE'
write(u6,*)
write(u6,*) ' ORBITALS ARE NUMBERED WITHIN EACH SEPARATE SYMMETRY.'
do I=1,NCONF
  CI = C(I)
  ACI = abs(C(I))
  if (I <= JCONF) then
    JMIN = 1
    NREXT = 0
    if (JREFX(I) /= 0) then
      CSFTYP = '   REFERENCE'
      CLIM = Zero
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
  if (ACI < CLIM) exit
  JJ = I
  IJ = I
  if (NREXT > 0) then
    JMAX = IRC(ILIM)
    do J=JMIN,JMAX
      JJ = J
      if (INDX(J) >= IJ) then
        JJ = JJ-1
        exit
      end if
    end do
  end if
  NSJ = MUL(JSUNP(INTSYM,JJ),LSYM)
  JVIR = IJ-INDX(JJ)
  II1 = (JJ-1)*LN
  do II=1,LN
    II1 = II1+1
    ISP(II+2) = ICUNP(ICSPCK,II1)
    IOC(II+2) = (1+ISP(II+2))/2
  end do
  if (NREXT == 0) then
    FRMT = FORM0
    if (LN2 > 36) FRMT = FORM00
    LN1 = 3
  else if (NREXT == 1) then
    IO = JVIR+NVIRP(NSJ)+LN
    IORBI(2) = IORB(IO)
    IOC(2) = 1
    ISP(2) = 1
    ILSYM(2) = NSJ
    FRMT = FORM1
    if (LN2 > 36) FRMT = FORM01
    LN1 = 2
  else
    IIN = 0
    outer: do II=1,NVIRT
      NA = II
      NSI = MUL(NSJ,NSM(LN+II))
      J1 = NVIRP(NSI)+1
      J2 = NVIRP(NSI)+NVIR(NSI)
      if (J2 > II) J2 = II
      if (J2 >= J1) then
        do J=J1,J2
          NB = J
          IIN = IIN+1
          if (IIN == JVIR) exit outer
        end do
      end if
    end do outer
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
      CI = CI*sqrt(Half)
      if (abs(CI) < CTRSH) cycle
      FRMT = FORM1
      if (LN2 > 36) FRMT = FORM01
      LN1 = 2
    else
      IOC(2) = 1
      ISP(2) = 2
      if (CSFTYP == '     TRIPLET') ISP(2) = 1
      IO = LN+NA
      IORBI(2) = IORB(IO)
      ILSYM(2) = NSM(IO)
      FRMT = FORM2
      if (LN2 > 36) FRMT = FORM02
      LN1 = 1
    end if
  end if
  write(u6,*)
  write(u6,105) I,C(I),CSFTYP
  write(u6,FRMT) 'SYMMETRY     ',(ILSYM(J),J=LN1,LN2)
  write(u6,FRMT) 'ORBITALS     ',(IORBI(J),J=LN1,LN2)
  write(u6,FRMT) 'OCCUPATION   ',(IOC(J),J=LN1,LN2)
  write(u6,FRMT) 'SPIN-COUPLING',(ISP(J),J=LN1,LN2)
end do

return

105 format(/6X,'CONFIGURATION',I7,3X,'COEFFICIENT',F10.6,A)

end subroutine PRWF_MRCI
