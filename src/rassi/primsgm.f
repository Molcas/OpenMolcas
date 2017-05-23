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
      SUBROUTINE PRIMSGM(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB1,
     &                   IFSBTAB2,COEFF,SGM,PSI)
      IMPLICIT NONE
      REAL*8 PSI(*),SGM(*)
      REAL*8 COEFF,CFFPHS,SCALE
      INTEGER IORBTAB(*),NASPRT
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB1(*),IFSBTAB2(*)
      INTEGER IFSB1,IBLKPOS1,ISST1,KSTARR1,NSBS1
      INTEGER IFSB2,IBLKPOS2,ISST2,KSTARR2,NSBS2
      INTEGER NSSTP,NHSH2,KHSH2,NFSB1,ISST,NSBS
      INTEGER NPOP1,IMODE,KSSTOP,KSBSOP
      INTEGER NDI,NDJ,KPOS,ISPART,KSSTTB
      INTEGER I,J,IPOS1,IPOS2,ISUM,LSBSET,ISP
      INTEGER ISSTARR(50)
      INTEGER ISBS1,ISBS2,ISORB,KOINFO,KSBSCR,KSBSAN
      INTEGER KSBS1,KSBS2,KSORB,KSSTCR,KSSTAN,MORSBITS
C     INTEGER NBLKDET1,NBLKDET2,NDETS1,NDETS2,IERR
      INTEGER NBLKDET1,NBLKDET2,NDETS1,NDETS2
C     INTEGER LDUM,NDUM
#include "WrkSpc.fh"
C Purpose: Add to the wave function SGM the result of applying
C an operator to PSI. The operator is a single creator (IMODE=1)
C or annihilator (IMODE=-1) multiplied by COEFF. Only the
C sector of SGM described by IFSBTAB1 will be updated.
      IF(COEFF.EQ.0.0D0) RETURN
C The orbital table:
      NASPRT= IORBTAB(9)
      KOINFO=19
      ISPART=IORBTAB(KOINFO+6+8*(ISORB-1))
      KSORB =IORBTAB(KOINFO+7+8*(ISORB-1))
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
        DO ISP=1,NASPRT
          ISSTARR(ISP)=IFSBTAB1(KPOS-1+ISP)
        END DO
        NBLKDET1 =IFSBTAB1(KPOS+NASPRT  )
        IBLKPOS1 =IFSBTAB1(KPOS+NASPRT+1)
C Initial values for lower and higher dimensions.
C Also, extra phase factor due to spin orbitals in higher substrings.
        NDI=1
        DO ISP=1,ISPART-1
          ISST1=ISSTARR(ISP)
          NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))
          NDI=NDI*NSBS1
        END DO
        NDJ=1
        CFFPHS=COEFF
        DO ISP=ISPART+1,NASPRT
          ISST1=ISSTARR(ISP)
          NSBS1=ISSTAB(KSSTTB+0+5*(ISST1-1))
          NPOP1=ISSTAB(KSSTTB+1+5*(ISST1-1))
          NDJ=NDJ*NSBS1
          IF(NPOP1.NE.2*(NPOP1/2)) CFFPHS=-CFFPHS
        END DO
        ISST1=ISSTARR(ISPART)
        NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))

C Modify the bra substring type by annih or creating ISORB
        ISST2=ISSTAB(KSSTOP-1+KSORB+MORSBITS*(ISST1-1))
        IF(ISST2.EQ.0) GOTO 200

C Determine dimensions for multiple daxpy:
C Dimension for earlier subpartitions is NDI
C Dimension for later   subpartitions is NDJ
C Dimensions for present subpartition are NSBS1,NSBS2

        NSBS2=ISSTAB(KSSTTB+5*(ISST2-1))
        ISSTARR(ISPART)=ISST2
C Get the corresponding FS block number
        CALL HSHGET(ISSTARR,NASPRT,NASPRT+2,IFSBTAB2(KSTARR2),
     &                NHSH2,IFSBTAB2(KHSH2),IFSB2)
        ISSTARR(ISPART)=ISST1
        IF(IFSB2.EQ.0) GOTO 200
        KPOS=KSTARR2+(NASPRT+2)*(IFSB2-1)
        NBLKDET2 =IFSBTAB2(KPOS+NASPRT  )
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
          IF (NDI.EQ.1) THEN
           IF (NDJ.EQ.1) THEN
            IPOS1=IBLKPOS1+(KSBS1-1)
            IPOS2=IBLKPOS2+(KSBS2-1)
            SGM(IPOS1)=SGM(IPOS1)+SCALE*PSI(IPOS2)
           ELSE
            DO J=0,NDJ-1
             IPOS1=IBLKPOS1+(KSBS1-1+NSBS1*J)
             IPOS2=IBLKPOS2+(KSBS2-1+NSBS2*J)
             SGM(IPOS1)=SGM(IPOS1)+SCALE*PSI(IPOS2)
            END DO
           END IF
          ELSE
           IF (NDJ.EQ.1) THEN
            DO I=0,NDI-1
             IPOS1=IBLKPOS1+I+NDI*(KSBS1-1)
             IPOS2=IBLKPOS2+I+NDI*(KSBS2-1)
             SGM(IPOS1)=SGM(IPOS1)+SCALE*PSI(IPOS2)
            END DO
           ELSE
            DO J=0,NDJ-1
             DO I=0,NDI-1
              IPOS1=IBLKPOS1+I+NDI*(KSBS1-1+NSBS1*J)
              IPOS2=IBLKPOS2+I+NDI*(KSBS2-1+NSBS2*J)
              SGM(IPOS1)=SGM(IPOS1)+SCALE*PSI(IPOS2)
             END DO
            END DO
           END IF
          END IF

 100      CONTINUE
        END DO

C End of loop over FS blocks
 200    CONTINUE
      END DO
      CALL GETMEM('NSBSET','Free','Inte',LSBSET,NSSTP)
      RETURN
      END
