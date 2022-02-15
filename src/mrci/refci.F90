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
real(kind=wp) :: HREF((NREF*(NREF+1))/2), AREF(NREF,NREF), EREF(NREF), CISEL(NREF,NSEL), PLEN(NREF)
integer(kind=iwp) :: ICSPCK(NCSPCK)
integer(kind=iwp) :: I, IC, IOFF, IREF, J, JJ, K, K1, K2, NC, NPRT
character(len=48) :: FORM1, FORM2, FORM3, FORM4
character(len=2) :: STR
integer(kind=iwp), external :: ICUNP
!Statement function
integer(kind=iwp) :: JCASE, L
!PAM97 integer(kind=iwp), external :: UNPACK
!PAM97 JCASE(L) = UNPACK(CSPCK((L+29)/30),2*L-(2*L-1)/60*60,2)
JCASE(L) = ICUNP(ICSPCK,L)

write(u6,*)
call XFLUSH(u6)
write(u6,*) ('-',I=1,60)
call XFLUSH(u6)
write(u6,*) '   REFERENCE CI CALCULATION.'
call XFLUSH(u6)
write(u6,*) ('-',I=1,60)
call XFLUSH(u6)
if (NSEL == 0) then
  write(u6,*) ' ROOT SELECTION BY ENERGY ORDERING.'
  call XFLUSH(u6)
  if (NRROOT == 1) then
    write(u6,'(A,I8)') '  ONE SINGLE ROOT, NUMBER.....: ',IROOT(1)
    call XFLUSH(u6)
  else
    write(u6,*) ' THE FOLLOWING ROOTS WILL BE SELECTED:'
    call XFLUSH(u6)
    write(u6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
    call XFLUSH(u6)
  end if
else
  write(u6,*) ' ROOT SELECTION BY PROJECTION: THE EIGENVECTORS OF'
  call XFLUSH(u6)
  write(u6,*) ' THE REFERENCE CI ARE ORDERED BY DECREASING SIZE OF'
  call XFLUSH(u6)
  write(u6,*) ' THEIR PROJECTIONS ONTO A SELECTION SPACE.'
  call XFLUSH(u6)
  if (NRROOT == 1) then
    write(u6,*) ' SELECT THE EIGENVECTOR WITH LARGEST PROJECTION.'
    call XFLUSH(u6)
  else
    write(u6,'(A,I2,A)') ' SELECT THE ',NRROOT,' EIGENVECTORS WITH LARGEST PROJECTION.'
    call XFLUSH(u6)
  end if
  write(u6,*) ' THE SELECTION SPACE IS SPANNED BY THE FOLLOWING VECTORS (NONZERO COMPONENTS ONLY):'
  call XFLUSH(u6)
  JJ = 0
  do I=1,NSEL
    write(u6,'(A,I2)') ' VECTOR NR. ',I
    call XFLUSH(u6)
    NC = NCOMP(I)
    write(u6,'(5X,I2,5X,A20,F12.8)') (J,SSEL(JJ+J),CSEL(JJ+J),J=1,NC)
    call XFLUSH(u6)
    JJ = JJ+NC
  end do
end if
write(u6,*)
call XFLUSH(u6)
call JACSCF(HREF,AREF,EREF,NREF,-1,1.0e-11_wp)
call ORDER(AREF,EREF,NREF)
call CI_SELECT_MRCI(NREF,AREF,PLEN,NSEL,CISEL,NRROOT,IROOT)
if (NSEL > 0) then
  write(u6,*) ' THE FOLLOWING ROOTS WERE SELECTED:'
  call XFLUSH(u6)
  write(u6,'(12(A,I2))') ' ROOTS NR ',IROOT(1),(',',IROOT(I),I=2,NRROOT-1),', AND ',IROOT(NRROOT)
  call XFLUSH(u6)
end if
write(STR,'(I2)') LN
call XFLUSH(u6)
FORM1 = '(2X,'//STR//'X,A,I7,2(8X,I7))'
FORM2 = '(2X,'//STR//'X,A,3F15.8)'
FORM3 = '('' CSF NR'',I5,'' CASE '','//STR//'I1,3(F13.6,2X))'
FORM4 = '(''       '',I5,''      '','//STR//'I1,3(F13.6,2X))'
write(u6,*)
call XFLUSH(u6)
write(u6,*) '        LOWEST REFERENCE CI ROOTS:'
call XFLUSH(u6)
NPRT = min(NREF,IROOT(NRROOT)+2)
do K1=1,NPRT,3
  K2 = min(NPRT,K1+2)
  write(u6,FORM1) '            ROOT',(K,K=K1,K2)
  call XFLUSH(u6)
  if (NSEL > 0) write(u6,FORM2) 'SELECTION WEIGHT',(PLEN(K),K=K1,K2)
  write(u6,FORM2) '          ENERGY',(EREF(K),K=K1,K2)
  call XFLUSH(u6)
  do IREF=1,NREF
    IC = IREFX(IREF)
    IOFF = LN*(IC-1)
    if (IREF == 1) then
      write(u6,FORM3) IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      call XFLUSH(u6)
    else
      write(u6,FORM4) IC,(JCASE(IOFF+J),J=1,LN),(AREF(IREF,K),K=K1,K2)
      call XFLUSH(u6)
    end if
  end do
  write(u6,*)
  call XFLUSH(u6)
end do
write(u6,*)
call XFLUSH(u6)

return

end subroutine REFCI
