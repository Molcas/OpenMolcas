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

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension HREF((NREF*(NREF+1))/2), AREF(NREF,NREF), EREF(NREF)
dimension PLEN(NREF), CISEL(NREF,NSEL), ICSPCK(NCSPCK)
character*2 STR
character*48 FORM1, FORM2, FORM3, FORM4
!PAM97 external UNPACK
!PAM97 integer UNPACK
!PAM97 JCASE(L) = UNPACK(CSPCK((L+29)/30),2*L-(2*L-1)/60*60,2)
JCASE(L) = ICUNP(ICSPCK,L)

write(6,*)
call XFLUSH(6)
write(6,*) ('-',I=1,60)
call XFLUSH(6)
write(6,*) '   REFERENCE CI CALCULATION.'
call XFLUSH(6)
write(6,*) ('-',I=1,60)
call XFLUSH(6)
if (NSEL == 0) then
  write(6,*) ' ROOT SELECTION BY ENERGY ORDERING.'
  call XFLUSH(6)
  if (NRROOT == 1) then
    write(6,'(A,I8)') '  ONE SINGLE ROOT, NUMBER.....: ',IROOT(1)
    call XFLUSH(6)
  else
    write(6,*) ' THE FOLLOWING ROOTS WILL BE SELECTED:'
    call XFLUSH(6)
    write(6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
    call XFLUSH(6)
  end if
else
  write(6,*) ' ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
  call XFLUSH(6)
  write(6,*) ' THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
  call XFLUSH(6)
  write(6,*) ' THEIR PROJECTIONS ONTO A SELECTION SPACE.'
  call XFLUSH(6)
  if (NRROOT == 1) then
    write(6,*) ' SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
    call XFLUSH(6)
  else
    write(6,'(A,I2,A)') ' SELECT THE ',NRROOT,' EIGENVECTORS WITH LARGEST PROJECTION.'
    call XFLUSH(6)
  end if
  write(6,*) ' THE SELECTION SPACE IS SPANNED BY THE FOLLOWING VECTORS (NONZERO COMPONENTS ONLY):'
  call XFLUSH(6)
  JJ = 0
  do I=1,NSEL
    write(6,'(A,I2)') ' VECTOR NR. ',I
    call XFLUSH(6)
    NC = NCOMP(I)
    write(6,'(5X,I2,5X,A20,F12.8)') (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
    call XFLUSH(6)
    JJ = JJ+NC
  end do
end if
write(6,*)
call XFLUSH(6)
call JACSCF(HREF,AREF,EREF,NREF,-1,1.0D-11)
call ORDER(AREF,EREF,NREF)
call CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)
if (NSEL > 0) then
  write(6,*) ' THE FOLLOWING ROOTS WERE SELECTED:'
  call XFLUSH(6)
  write(6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
  call XFLUSH(6)
end if
write(STR,'(I2)') LN
call XFLUSH(6)
FORM1 = '(2X,'//STR//'X,A,I7,2(8X,I7))'
FORM2 = '(2X,'//STR//'X,A,3F15.8)'
FORM3 = '('' CSF NR'',I5,'' CASE '','//STR//'I1,3(F13.6,2X))'
FORM4 = '(''       '',I5,''      '','//STR//'I1,3(F13.6,2X))'
write(6,*)
call XFLUSH(6)
write(6,*) '        LOWEST REFERENCE CI ROOTS:'
call XFLUSH(6)
NPRT = min(NREF,IROOT(NRROOT)+2)
do K1=1,NPRT,3
  K2 = min(NPRT,K1+2)
  write(6,FORM1) '            ROOT',(K,K=K1,K2)
  call XFLUSH(6)
  if (NSEL > 0) write(6,FORM2) 'SELECTION WEIGHT',(PLEN(K),K=K1,K2)
  write(6,FORM2) '          ENERGY',(EREF(K),K=K1,K2)
  call XFLUSH(6)
  do IREF=1,NREF
    IC = IREFX(IREF)
    IOFF = LN*(IC-1)
    if (IREF == 1) then
      write(6,FORM3) IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      call XFLUSH(6)
    else
      write(6,FORM4) IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      call XFLUSH(6)
    end if
  end do
  write(6,*)
  call XFLUSH(6)
end do
write(6,*)
call XFLUSH(6)

return

end subroutine REFCI
