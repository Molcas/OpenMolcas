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
      SUBROUTINE HAM3(OP0,OP1,NOP2,OP2,NOP3,OP3,ISYCI,CI,SGM)
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: SGS, CIS, EXS
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION OP1(NASHT,NASHT),OP2(NOP2),OP3(NOP3)
      DIMENSION CI(*),SGM(*)
C Local arrays:
      DIMENSION IATOG(MXLEV)
      Real*8, Allocatable:: SGM1(:), SGM2(:)
      Integer nLev
      nLev = SGS%nLev

C Purpose: Compute and add a contribution to SGM which is
C obtained from a sum of zero- one- two- and three-electron
C operators acting on wave function CI.

C Note that the coefficients in OP1 and OP2 must have been
C modified by adding elements from OP2 and OP3, as done in
C subroutine MODOP.

C Presently symmetry blocking is disregarded for OP2, OP3, but
C index pair C permutation symmetry is used.
C NOP2=(NASHT**2+1 over 2)  (Binomial coefficient)
C NOP3=(NASHT**2+2 over 3)  (Binomial coefficient)

      IF(NCONF.EQ.0) RETURN
      IF(ABS(OP0).GT.1.0D-15) THEN
        CALL DAXPY_(NCONF,OP0,CI,1,SGM,1)
      END IF
      IF(NACTEL.EQ.0) RETURN

C Unless this is a special-case wave function, reserve space
C for intermediate results of elementary excitations.
      IF(ISCF.EQ.0) THEN
        CALL mma_allocate(SGM1,MXCI,Label='SGM1')
        IF(NACTEL.GE.2) CALL mma_allocate(SGM2,MXCI,Label='SGM2')
      END IF
C Special cases:
      OCCNO=0.0d0
      IF(ISCF.EQ.1) OCCNO=2.0D0
      IF(ISCF.EQ.2) OCCNO=1.0D0

C Create reorder table giving the GUGA level, i.e. CI-coupling
C ordinal number of each active orbital.
      ITABS=0
      DO ISYM=1,NSYM
        DO I=1,NLEV
          IF(SGS%ISM(I).EQ.ISYM) THEN
            ITABS=ITABS+1
            IATOG(ITABS)=I
          END IF
        END DO
      END DO

      DO IZ=1,NASHT
       DO IY=1,NASHT
        IYZ=IY+(IZ-1)*NASHT
        ISYZ=MUL(IASYM(IY),IASYM(IZ))
        ISYM1=MUL(ISYZ,ISYCI)
        NSGM1=CIS%NCSF(ISYM1)
        IF(NSGM1.EQ.0) CYCLE
        IF(ISCF.EQ.0) THEN
C The general case:
C Compute SGM1:=E(IY,IZ) PSI
          CALL DCOPY_(NSGM1,[0.0D0],0,SGM1,1)
          LEVY=IATOG(IY)
          LEVZ=IATOG(IZ)
          CALL SIGMA1(SGS,CIS,EXS,
     &                LEVY,LEVZ,1.0D00,ISYCI,CI,SGM1)
C Add non-zero 1-el contribution to SGM:
          IF(ISYZ.EQ.1) THEN
            X=OP1(IY,IZ)
            IF(ABS(X).GT.1.0D-15) THEN
              CALL DAXPY_(NCONF,X,SGM1,1,SGM,1)
CTEST      WRITE(*,*)' op1:',X
CTEST      WRITE(*,*)' iyz, sgm(1):',iyz,sgm(1)
            END IF
          END IF
        ELSE
C Closed-shell or hi-spin case:
          IF(IY.NE.IZ) CYCLE
          X=OCCNO*OP1(IY,IZ)
          SGM(1)=SGM(1)+X*CI(1)
        END IF
        IF(NACTEL.EQ.1) CYCLE
        DO IX=IZ,NASHT
         IVMIN=1
         IF(IX.EQ.IZ) IVMIN=IY
         DO IV=IVMIN,NASHT
          IVX=IV+(IX-1)*NASHT
          ISVX=MUL(IASYM(IV),IASYM(IX))
          ISVXYZ=MUL(ISVX,ISYZ)
          IVXYZ=(IVX*(IVX-1))/2+IYZ
          ISYM2=MUL(ISVX,ISYM1)
          NSGM2=CIS%NCSF(ISYM2)
          IF(NSGM2.EQ.0) CYCLE
          IF(ISCF.EQ.0) THEN
C The general case:
C Compute SGM2:=E(IV,IX) SGM1
            CALL DCOPY_(NSGM2,[0.0D0],0,SGM2,1)
            LEVV=IATOG(IV)
            LEVX=IATOG(IX)
            CALL SIGMA1(SGS,CIS,EXS,
     &                  LEVV,LEVX,1.0D00,ISYM1,SGM1,SGM2)
C Add non-zero 2-el contribution to SGM:
            IF(ISVXYZ.EQ.1) THEN
              X=OP2(IVXYZ)
              IF(ABS(X).GT.1.0D-15) THEN
                CALL DAXPY_(NCONF,X,SGM2,1,SGM,1)
CTEST      WRITE(*,*)' op2:',X
CTEST      WRITE(*,*)' ivxyz, sgm(1):',ivxyz,sgm(1)
              END IF
            END IF
          ELSE
C Closed-shell or hi-spin case:
            IF(IY.NE.IZ) CYCLE
            IF(IV.NE.IX) CYCLE
            X=(OCCNO**2)*OP2(IVXYZ)
            SGM(1)=SGM(1)+X*CI(1)
          END IF
          IF(NACTEL.EQ.2) CYCLE
          DO IU=IX,NASHT
           ITMIN=1
           IF(IU.EQ.IX) ITMIN=IV
           DO IT=ITMIN,NASHT
            ITU=IT+(IU-1)*NASHT
            ISTU=MUL(IASYM(IT),IASYM(IU))
            IF(ISTU.NE.ISVXYZ) CYCLE
            ITUVXYZ=((ITU+1)*ITU*(ITU-1))/6+IVXYZ
            X=OP3(ITUVXYZ)
            IF(ABS(X).LT.1.0D-15) CYCLE
C Add non-zero 3-el contribution to SGM:
            IF(ISCF.EQ.0) THEN
              LEVT=IATOG(IT)
              LEVU=IATOG(IU)
              CALL SIGMA1(SGS,CIS,EXS,
     &                    LEVT,LEVU,X,ISYM2,SGM2,SGM)
CTEST      WRITE(*,*)' op3:',X
CTEST      WRITE(*,*)' ituvxyz, sgm(1):',ituvxyz,sgm(1)
            ELSE
C Closed-shell or hi-spin case:
              IF(IT.NE.IU) CYCLE
              IF(IV.NE.IX) CYCLE
              IF(IY.NE.IZ) CYCLE
              X=(OCCNO**3)*OP3(ITUVXYZ)
              SGM(1)=SGM(1)+X*CI(1)
            END IF
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO

C Deallocate temporary arrays, if any:
      IF(ISCF.EQ.0) THEN
        CALL mma_deallocate(SGM1)
        IF(NACTEL.GE.2) CALL mma_deallocate(SGM2)
      END IF

      END SUBROUTINE HAM3
