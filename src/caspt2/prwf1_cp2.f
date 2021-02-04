************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PRWF1_CP2(NOCSF,IOCSF,NOW,IOW,ISYCI,CI,THR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION CI(*)
      CHARACTER(LEN=256) LINE
      CHARACTER(LEN=1) CODE(0:3)

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
      DIMENSION ICS(MXLEV)
      DATA CODE /'0','u','d','2'/

C -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
C -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
C -- THE   U P P E R   PART OF THE WALK.

C SVC: set up a CSF string length as LENCSF
      LINE=' '
      LENCSF=0
      ISY=0
      DO LEV=1,NLEV
        IF(ISY.NE.ISM(LEV)) THEN
          ISY=ISM(LEV)
          LENCSF=LENCSF+1
        END IF
        LENCSF=LENCSF+1
      END DO
      LENCSF=MIN(LENCSF,256)
      LENCSF=MAX(LENCSF,10)

 100  FORMAT(2X,A10,2X,A16,2X,A,2(2X,A13))
 200  FORMAT(2X,I10,2X,'(',I2,':',I1,':',I4,'/',I4,')',
     &       2X,A,2(2X,F13.6))

C Size of occup/spin coupling part of line:
      WRITE(6,*)' Occupation of active orbitals, and spin coupling'
      WRITE(6,*)' of open shells. (u,d: Spin up or down).'
      WRITE(6,*)' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
      LINE(1:10)='Occupation'
      WRITE(6,100)
     & 'Conf','SGUGA info      ',LINE(1:LENCSF),
     & 'Coefficient','Weight'

C     SVC2010:
C     allocate scratch memory for determinant expansion
      IF (PRSD) THEN
        CALL GETMEM ('LEX','ALLO','INTEGER',LLEX,NLEV)
      END IF

      LINE=' '

C -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
      DO 40 MV=1,NMIDV
        DO 41 ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,ISYCI)
          IF(NCI.EQ.0) GOTO 41
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=MUL(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,ISYCI)
          IUW0=LICASE-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=LICASE-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO 30 IDWN=1,NDWN
            DO 31 IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF).LT.THR) GOTO  31
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=IWORK(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO 10 LEV=1,MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=IWORK(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
  10            CONTINUE
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=IWORK(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=IWORK(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
C -- PRINT IT!
              K=0
              ISY=0
              DO LEV=1,NLEV
                IF(ISY.NE.ISM(LEV)) THEN
                  ISY=ISM(LEV)
                  K=K+1
                  LINE(K:K)=' '
                END IF
                K=K+1
                LINE(K:K)=CODE(ICS(LEV))
              END DO
              WRITE(6,200)
     &               ICONF,MV,ISYUP,IUP,IDWN,
     &               LINE(1:LENCSF),COEF,COEF**2
C     SVC2010 experimental: add determinant expansion
              IF (PRSD) THEN
c     Specify projected spin in half integer units
C     Default: use maximum spin projection
               IMS = ISPIN-1
               WRITE(6,*)
               CALL EXPCSF (ICS, NLEV, IMS, IWORK(LLEX))
               WRITE(6,*)
              ENDIF
  31        CONTINUE
  30      CONTINUE
  41    CONTINUE
  40  CONTINUE
C     SVC2010: free scratch for determinant expansion
      IF (PRSD) THEN
        CALL GETMEM ('LEX','FREE','INTEGER',LLEX,NLEV)
      END IF
      WRITE(6,*)
      RETURN
      END
