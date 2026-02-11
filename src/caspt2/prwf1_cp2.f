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
      SUBROUTINE PRWF1_CP2(NOCSF,IOCSF,NOW,IOW,ISYCI,CI,THR,nMidV)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp, u6
      use gugx, only: SGS, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, ISPIN, PRSD
      use gugx, only: MxLev
      IMPLICIT None
      Integer(kind=iwp), Intent(In):: ISYCI, nMidV
      integer(kind=iwp), intent(in):: NOCSF(NSYM,NMIDV,NSYM),
     &                                IOCSF(NSYM,NMIDV,NSYM)
      integer(kind=iwp), intent(in):: NOW(2,NSYM,NMIDV),
     &                                IOW(2,NSYM,NMIDV)
      real(kind=wp), intent(in):: CI(*), THR

      CHARACTER(LEN=256) LINE
      CHARACTER(LEN=1) :: CODE(0:3)=['0','u','d','2']
      integer(kind=iwp) ICS(MXLEV)
      integer(kind=iwp) :: nLev, nIpWlk
      integer(kind=iwp), ALLOCATABLE:: LEX(:)
      real(kind=wp) COEF
      integer(kind=iwp) IC1, ICDPOS, ICDWN, ICONF, ICUP, ICUPOS, IDW0,
     &                  IDWN, IDWNSV, IMS, ISY, ISYDWN, ISYUP, IUP,
     &                  IUW0, K, LENCSF, LEV, MV, NCI, NDWN, NNN, NUP

      nLev  = SGS%nLev
      nIpWlk= CIS%nIpWlk

C -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
C -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
C -- THE   U P P E R   PART OF THE WALK.

C SVC: set up a CSF string length as LENCSF
      LINE=' '
      LENCSF=0
      ISY=0
      DO LEV=1,NLEV
        IF(ISY/=SGS%ISM(LEV)) THEN
          ISY=SGS%ISM(LEV)
          LENCSF=LENCSF+1
        END IF
        LENCSF=LENCSF+1
      END DO
      LENCSF=MIN(LENCSF,256)
      LENCSF=MAX(LENCSF,10)


C Size of occup/spin coupling part of line:
      WRITE(u6,*)' Occupation of active orbitals, and spin coupling'
      WRITE(u6,*)' of open shells. (u,d: Spin up or down).'
      WRITE(u6,*)' SGUGA info is (Midvert:IsyUp:UpperWalk/LowerWalk)'
      LINE(1:10)='Occupation'
      WRITE(u6,'(2X,A10,2X,A16,2X,A,2(2X,A13))')
     & 'Conf','SGUGA info      ',LINE(1:LENCSF),
     & 'Coefficient','Weight'

C     SVC2010:
C     allocate scratch memory for determinant expansion
      IF (PRSD) THEN
        CALL mma_allocate(LEX,NLEV,LABEL='LEX')
      END IF

      LINE=' '

C -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
      DO MV=1,NMIDV
        DO ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,ISYCI)
          IF(NCI==0) Cycle
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=Mul(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,ISYCI)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO IDWN=1,NDWN
            DO IUP=1,NUP
              ICONF=ICONF+1
              COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
              IF(ABS(COEF)<THR) Cycle
              IF(IDWNSV/=IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=CIS%ICASE(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO LEV=1,SGS%MIDLEV
                  NNN=NNN+1
                  IF(NNN==16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=CIS%ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
                End Do
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=CIS%ICASE(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=SGS%MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN==16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=CIS%ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
C -- PRINT IT!
              K=0
              ISY=0
              DO LEV=1,NLEV
                IF(ISY/=SGS%ISM(LEV)) THEN
                  ISY=SGS%ISM(LEV)
                  K=K+1
                  LINE(K:K)=' '
                END IF
                K=K+1
                LINE(K:K)=CODE(ICS(LEV))
              END DO
              WRITE(u6,"(2X,I10,2X,'(',I2,':',I1,':',I4,'/',I4,')',
     &       2X,A,2(2X,F13.6))")
     &               ICONF,MV,ISYUP,IUP,IDWN,
     &               LINE(1:LENCSF),COEF,COEF**2
C     SVC2010 experimental: add determinant expansion
              IF (PRSD) THEN
c     Specify projected spin in half integer units
C     Default: use maximum spin projection
               IMS = ISPIN-1
               WRITE(u6,*)
               CALL EXPCSF (ICS, NLEV, IMS, LEX, coef, 0)
               WRITE(u6,*)
              ENDIF
            End Do
          End Do
        End Do
      End Do

C     SVC2010: free scratch for determinant expansion
      IF (PRSD) CALL mma_deallocate(LEX)
      WRITE(u6,*)

      END SUBROUTINE PRWF1_CP2
