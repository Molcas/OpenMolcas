************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SGMOP(IMODE,IORBTAB,ISSTAB,IFSBTAB1,IFSBTAB2,
     &                   COEFF,SGM,PSI)
      IMPLICIT NONE
      REAL*8 COEFF(*),PSI(*),SGM(*)
      REAL*8 CFFPHS,SCALE
      INTEGER IORBTAB(*),NASPRT,NASORB
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER IFSB1,IBLKPOS1,ISST1,KSTARR1,NSBS1
      INTEGER IFSB2,IBLKPOS2,ISST2,KSTARR2,NSBS2
      INTEGER NSSTP,NHSH2,KHSH2,NFSB1,ISST,NSBS
      INTEGER IPH,NPOP1,IMODE,KSSTOP,KSBSOP
      INTEGER NDI,NDJ,KPOS,ISPART,KSSTTB
      INTEGER I,J,IPOS1,IPOS2,ISUM,LSBSET
      INTEGER ISSTARR(50),NDIARR(50),NDJARR(50),IPHARR(50)
      INTEGER ISBS1,ISBS2,ISORB,KOINFO,KSBSCR,KSBSAN
      INTEGER KSBS1,KSBS2,KSORB,KSSTCR,KSSTAN,MORSBITS
      INTEGER NDETS1,NDETS2,IERR
#include "WrkSpc.fh"

C Purpose: Add to the wave function SGM the result of applying
C an operator to PSI. The operator is a sum of creators (IMODE=1)
C or annihilators (IMODE=-1) multiplied by coefficients. Only the
C sector of SGM described by IFSBTAB1 will be updated.
C The orbital table:
      NASPRT= IORBTAB(9)
      NASORB= IORBTAB(4)
      KOINFO=19
C The substring table:
      MORSBITS=ISSTAB(6)
      NSSTP   =ISSTAB(7)
      KSSTTB=15
      KSSTAN=ISSTAB( 9)
      KSSTCR=ISSTAB(10)
      KSBSAN=ISSTAB(13)
      KSBSCR=ISSTAB(14)
      IF(IMODE.EQ.1) THEN
        KSSTOP=KSSTAN
        KSBSOP=KSBSAN
      ELSE
        KSSTOP=KSSTCR
        KSBSOP=KSBSCR
      END IF
C The FS blocks of the SGM wave function:
      NFSB1=IFSBTAB1(3)
      NDETS1=IFSBTAB1(5)
      KSTARR1=8
C The FS blocks of the PSI wave function:
      NDETS2=IFSBTAB2(5)
      NHSH2=IFSBTAB2(6)
      KHSH2=IFSBTAB2(7)
      KSTARR2=8
C Make an array with nr of earlier substrings for each
C substring type:
      CALL GETMEM('NSBSET','Allo','Inte',LSBSET,NSSTP)
      ISUM=0
      DO ISST=1,NSSTP
        IWORK(LSBSET-1+ISST)=ISUM
        NSBS=ISSTAB(KSSTTB+5*(ISST-1))
        ISUM=ISUM+NSBS
      END DO

C Loop over FS blocks of the SGM wave function
      DO IFSB1=1,NFSB1
        KPOS=KSTARR1+(NASPRT+2)*(IFSB1-1)
        DO ISPART=1,NASPRT
          ISSTARR(ISPART)=IFSBTAB1(KPOS-1+ISPART)
        END DO
        IBLKPOS1 =IFSBTAB1(KPOS+NASPRT+1)
CTEST      write(*,'(1x,a,8I8)')'SGM FSB1,IBLKPOS1=',IBLKPOS1
CTEST      write(*,'(1x,a,8I8)')'SGM FSB1=',(ISSTARR(I),I=1,NASPRT)
C Initial values for lower and higher dimensions.
C Also, extra phase factor due to spin orbitals in higher substrings.
        NDI=1
        DO ISPART=1,NASPRT
          NDIARR(ISPART)=NDI
          ISST1=ISSTARR(ISPART)
          NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))
          NDI=NDI*NSBS1
        END DO
        NDJ=1
        IPH=1
        DO ISPART=NASPRT,1,-1
          NDJARR(ISPART)=NDJ
          IPHARR(ISPART)=IPH
          ISST1=ISSTARR(ISPART)
          NSBS1=ISSTAB(KSSTTB+0+5*(ISST1-1))
          NPOP1=ISSTAB(KSSTTB+1+5*(ISST1-1))
          NDJ=NDJ*NSBS1
          IF(NPOP1.NE.2*(NPOP1/2)) IPH=-IPH
        END DO
C Loop over active orbitals:
        DO ISORB=1,NASORB
          IF(COEFF(ISORB).EQ.0.0D0) GOTO 200
          ISPART=IORBTAB(KOINFO+6+8*(ISORB-1))
          CFFPHS=DBLE(IPHARR(ISPART))*COEFF(ISORB)
          KSORB =IORBTAB(KOINFO+7+8*(ISORB-1))
CTEST      write(*,'(1x,a,8I8)')'ISORB,ISPART,KSORB:',
CTEST     &                      ISORB,ISPART,KSORB
          ISST1=ISSTARR(ISPART)
          NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))

C Modify the bra substring type by annih or creating ISORB
          ISST2=ISSTAB(KSSTOP-1+KSORB+MORSBITS*(ISST1-1))
          IF(ISST2.EQ.0) GOTO 200

CTEST      write(*,'(1x,a,8I8)')'ISST1,ISST2:',ISST1,ISST2
C Determine dimensions for multiple daxpy:
C Dimension for earlier subpartitions is NDI
C Dimension for later   subpartitions is NDJ
C Dimensions for present subpartition are NSBS1,NSBS2
          NDI=NDIARR(ISPART)
          NDJ=NDJARR(ISPART)

          NSBS2=ISSTAB(KSSTTB+5*(ISST2-1))
          ISSTARR(ISPART)=ISST2
C Get the corresponding FS block number
CTEST      write(*,'(1x,a,8I8)')'PSI FSB2=',(ISSTARR(I),I=1,NASPRT)
          CALL HSHGET(ISSTARR,NASPRT,NASPRT+2,IFSBTAB2(KSTARR2),
     &                NHSH2,IFSBTAB2(KHSH2),IFSB2)
CTEST      write(*,'(1x,a,8I8)')'IFSB1,IFSB2:',IFSB1,IFSB2
          ISSTARR(ISPART)=ISST1
          IF(IFSB2.EQ.0) GOTO 200
          KPOS=KSTARR2+(NASPRT+2)*(IFSB2-1)
          IBLKPOS2 =IFSBTAB2(KPOS+NASPRT+1)
C Now loop over ket substrings in this subpartition
          DO KSBS1=1,NSBS1
            ISBS1=KSBS1+IWORK(LSBSET-1+ISST1)
            ISBS2=ISSTAB(KSBSOP-1+KSORB+MORSBITS*(ISBS1-1))
            IF(ISBS2.EQ.0) GOTO 100
            IF(ISBS2.GT.0) THEN
              SCALE= CFFPHS
              ISBS2= ISBS2
            ELSE
              SCALE=-CFFPHS
              ISBS2=-ISBS2
            END IF
            KSBS2=ISBS2-IWORK(LSBSET-1+ISST2)

C CALL some multiple daxpy...
            DO I=0,NDI-1
             DO J=0,NDJ-1
              IPOS1=IBLKPOS1+I+NDI*(KSBS1-1+NSBS1*J)
              IPOS2=IBLKPOS2+I+NDI*(KSBS2-1+NSBS2*J)
              SGM(IPOS1)=SGM(IPOS1)+SCALE*PSI(IPOS2)
              IERR=0
              IF(IPOS1.LT.1 .OR. IPOS1.GT.NDETS1) IERR=1
              IF(IPOS2.LT.1 .OR. IPOS2.GT.NDETS2) IERR=1
              IF(IERR.NE.0) THEN
                WRITE(6,*)' SGMOP addressing error.'
                WRITE(6,'(1x,a,8I8)')' SGM dimension NDETS1=',NDETS1
                WRITE(6,'(1x,a,8I8)')' Position IPOS1=',IPOS1
                WRITE(6,'(1x,a,8I8)')' PSI dimension NDETS2=',NDETS2
                WRITE(6,'(1x,a,8I8)')' Position IPOS2=',IPOS2
                CALL ABEND()
              END IF
C Temporary test print:
CTEST              if(PSI(IPOS2).ne.0.0d0) then
CTEST                WRITE(*,'(1x,f16.8,2i8)')SCALE,IPOS1,IPOS2
CTEST                write(*,'(1x,a,8I8)')'IFSB1,IFSB2:',IFSB1,IFSB2
CTEST                write(*,'(1x,a,8I8)')'IBLKPOS1,NDI,NSBS1:',
CTEST     &                                IBLKPOS1,NDI,NSBS1
CTEST                write(*,'(1x,a,8I8)')'IBLKPOS2,NDI,NSBS2:',
CTEST     &                                IBLKPOS2,NDI,NSBS2
CTEST                write(*,'(1x,a,8I8)')'I,KSBS1,J:',I,KSBS1,J
CTEST                write(*,'(1x,a,8I8)')'I,KSBS2,J:',I,KSBS2,J
CTEST              end if
C End of test prints
             END DO
            END DO

 100        CONTINUE
          END DO

 200      CONTINUE
C End of loop over orbitals
        END DO
C End of loop over FS blocks
      END DO
      CALL GETMEM('NSBSET','Free','Inte',LSBSET,NSSTP)
      RETURN
      END
