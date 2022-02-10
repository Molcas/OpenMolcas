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

subroutine INDMAT(ICSPCK,INTSYM,INDX,ISAB,JREFX,CISEL)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INDX(*), ICSPCK(*), INTSYM(*)
dimension ISAB(NVIRT,NVIRT), JREFX(*)
dimension CISEL(NREF,NSEL)
character*20 STR20
!PAM97 external UNPACK
!PAM97 integer UNPACK
!Statement functions
!PAM97 JCASE(L)=UNPACK(CSPCK((L+29)/30), 2*L-(2*L-1)/60*60, 2)
JCASE(L) = ICUNP(ICSPCK,L)
!PAM96 JSYM(L)=UNPACK(INTSYM((L+9)/10),3*mod(L-1,10)+1,3)+1
JSYM(L) = JSUNP(INTSYM,L)

ILIM = 4
if (IFIRST /= 0) ILIM = 2
do I=1,NSYM
  NVPAIR(I) = 0
end do
ISMAX = 0
do NA=1,NVIRT
  do NB=1,NA
    NSAB = MUL(NSM(LN+NA),NSM(LN+NB))
    NVPAIR(NSAB) = NVPAIR(NSAB)+1
    ISAB(NA,NB) = NVPAIR(NSAB)
    ISMAX = max(ISMAX,ISAB(NA,NB))
    ISAB(NB,NA) = ISAB(NA,NB)
  end do
  NDIAG(NA) = ISAB(NA,NA)
end do
! INDX - STARTING POINT IN CI VECTOR OF EACH BLOCK WITH A
! COMMON INTERNAL WALK.
! VALENCE CONFIGURATIONS:
IR = IRC(1)
do II=1,IR
  INDX(II) = II
end do
JSC(1) = IR
! SINGLY EXCITED CONFIGURATIONS:
IR1 = IR+1
IR2 = IRC(2)
IND = IR
do II=IR1,IR2
  INDX(II) = IND
  NSS = MUL(JSYM(II),LSYM)
  IND = IND+NVIR(NSS)
end do
JSC(2) = IND
NCDOUB = IND-JSC(1)
if (IFIRST == 0) then
  ! DOUBLY EXCITED CONFIGURATIONS:
  IR1 = IR2+1
  IR2 = IRC(4)
  JSC(3) = JSC(2)
  do II=IR1,IR2
    INDX(II) = IND
    NSS = MUL(JSYM(II),LSYM)
    IND = IND+NVPAIR(NSS)
    if (II == IRC(3)) JSC(3) = IND
  end do
  JSC(4) = IND
  JJM = (JJS(LSYM+1)-JJS(LSYM))*NVIRT
  NCTRIP = JSC(3)-JSC(2)-JJM
  NCSING = JSC(4)-JSC(3)
else
  JJM = 0
  NCTRIP = 0
  NCSING = 0
end if
NCONF = JSC(ILIM)
! LIST THE REFERENCE CONFIGURATIONS, AND AT THE SAME TIME,
! IDENTIFY CSFS GIVEN IN SELECTION VECTOR INPUT:
write(6,*)
call XFLUSH(6)
write(6,*) '      LIST OF REFERENCE CONFIGURATIONS.'
call XFLUSH(6)
write(6,*) '     CONF NR:    GUGA CASE NUMBERS OF ACTIVE ORBITALS:'
call XFLUSH(6)
do I=1,IRC(1)
  IREF = JREFX(I)
  if (IREF == 0) goto 135
  IREFX(IREF) = I
  IOFF = LN*(I-1)
  write(6,'(5X,I6,7X,30I1)') I,(JCASE(IOFF+J),J=1,LN)
  call XFLUSH(6)
  JJ = 0
  do ISEL=1,NSEL
    CISEL(IREF,ISEL) = 0.0d00
    NC = NCOMP(ISEL)
    do IC=1,NC
      STR20 = SSEL(JJ+IC)
      do ILEV=1,LN
        JCAS1 = JCASE(IOFF+ILEV)
        read(STR20(ILEV:ILEV),'(I1)') JCAS2
        if (JCAS1 /= JCAS2) goto 133
      end do
      CISEL(IREF,ISEL) = CSEL(JJ+IC)
      JJ = JJ+NC
      goto 134
133   continue
    end do
    JJ = JJ+NC
134 continue
  end do
135 continue
end do
! ORTHONORMALIZE THE SELECTION VECTORS:
do ISEL=1,NSEL
  do JSEL=1,ISEL-1
    X = DDOT_(NREF,CISEL(1,JSEL),1,CISEL(1,ISEL),1)
    call DAXPY_(NREF,-X,CISEL(1,JSEL),1,CISEL(1,ISEL),1)
  end do
  X = 1.0d00/DDOT_(NREF,CISEL(1,ISEL),1,CISEL(1,ISEL),1)
  call DSCAL_(NREF,X,CISEL(1,ISEL),1)
end do
write(6,*)
call XFLUSH(6)
write(6,*) '      REAL CONFIGURATIONS:'
call XFLUSH(6)
if (IFIRST == 0) then
  write(6,215) NREF,NCVAL-NREF,NCDOUB,NCTRIP,NCSING
  call XFLUSH(6)
215 format(/,6X,'               REFERENCE ',I8,/,6X,'           OTHER VALENCE ',I8,/,6X,' DOUBLET COUPLED SINGLES ',I8, &
           /,6X,' TRIPLET COUPLED DOUBLES ',I8,/,6X,' SINGLET COUPLED DOUBLES ',I8)
else
  write(6,216) NREF,NCVAL-NREF,NCDOUB
  call XFLUSH(6)
216 format(/,6X,'               REFERENCE ',I8,/,6X,'           OTHER VALENCE ',I8,/,6X,' DOUBLET COUPLED SINGLES ',I8)
end if
JSCI = JSC(ILIM)-JJM
write(6,'(6X,A,I8)') '                  TOTAL ',JSCI
call XFLUSH(6)

return

end subroutine INDMAT
