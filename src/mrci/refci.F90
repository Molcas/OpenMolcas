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

subroutine REFCI(HREF,AREF,EREF,ICSPCK,CISEL,PLEN)

use mrci_global, only: CSEL, IREFX, IROOT, LN, NCOMP, NCSPCK, NREF, NRROOT, NSEL, SSEL
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: HREF((NREF*(NREF+1))/2)
real(kind=wp), intent(out) :: AREF(NREF,NREF), EREF(NREF), PLEN(NREF)
integer(kind=iwp), intent(in) :: ICSPCK(NCSPCK)
real(kind=wp), intent(in) :: CISEL(NREF,NSEL)
integer(kind=iwp) :: I, IC, IOFF, IREF, J, JJ, K, K1, K2, NC, NPRT
character(len=48) :: FORM1, FORM2, FORM3, FORM4
character(len=2) :: STR
integer(kind=iwp), external :: ICUNP

write(u6,*)
write(u6,*) repeat('-',60)
write(u6,*) '   REFERENCE CI CALCULATION.'
write(u6,*) repeat('-',60)
if (NSEL == 0) then
  write(u6,*) ' ROOT SELECTION BY ENERGY ORDERING.'
  if (NRROOT == 1) then
    write(u6,'(A,I8)') '  ONE SINGLE ROOT, NUMBER.....: ',IROOT(1)
  else
    write(u6,*) ' THE FOLLOWING ROOTS WILL BE SELECTED:'
    write(u6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
  end if
else
  write(u6,*) ' ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
  write(u6,*) ' THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
  write(u6,*) ' THEIR PROJECTIONS ONTO A SELECTION SPACE.'
  if (NRROOT == 1) then
    write(u6,*) ' SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
  else
    write(u6,'(A,I2,A)') ' SELECT THE ',NRROOT,' EIGENVECTORS WITH LARGEST PROJECTION.'
  end if
  write(u6,*) ' THE SELECTION SPACE IS SPANNED BY THE FOLLOWING VECTORS (NONZERO COMPONENTS ONLY):'
  JJ = 0
  do I=1,NSEL
    write(u6,'(A,I2)') ' VECTOR NR. ',I
    NC = NCOMP(I)
    write(u6,'(5X,I2,5X,A20,F12.8)') (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
    JJ = JJ+NC
  end do
end if
write(u6,*)
call JACSCF(HREF,AREF,EREF,NREF,-1,1.0e-11_wp)
call ORDER(AREF,EREF,NREF)
call CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)
if (NSEL > 0) then
  write(u6,*) ' THE FOLLOWING ROOTS WERE SELECTED:'
  write(u6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
end if
write(STR,'(I2)') LN
FORM1 = '(2X,'//STR//'X,A,I7,2(8X,I7))'
FORM2 = '(2X,'//STR//'X,A,3F15.8)'
FORM3 = '('' CSF NR'',I5,'' CASE '','//STR//'I1,3(F13.6,2X))'
FORM4 = '(''       '',I5,''      '','//STR//'I1,3(F13.6,2X))'
write(u6,*)
write(u6,*) '        LOWEST REFERENCE CI ROOTS:'
NPRT = min(NREF,IROOT(NRROOT)+2)
do K1=1,NPRT,3
  K2 = min(NPRT,K1+2)
  write(u6,FORM1) '            ROOT',(K,K=K1,K2)
  if (NSEL > 0) write(u6,FORM2) 'SELECTION WEIGHT',(PLEN(K),K=K1,K2)
  write(u6,FORM2) '          ENERGY',(EREF(K),K=K1,K2)
  do IREF=1,NREF
    IC = IREFX(IREF)
    IOFF = LN*(IC-1)
    if (IREF == 1) then
      write(u6,FORM3) IC,(ICUNP(ICSPCK,IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
    else
      write(u6,FORM4) IC,(ICUNP(ICSPCK,IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
    end if
  end do
  write(u6,*)
end do
write(u6,*)

return

end subroutine REFCI
