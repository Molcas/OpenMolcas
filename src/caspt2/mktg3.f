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
      SUBROUTINE MKTG3(LSYM1,LSYM2,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
      IMPLICIT REAL*8 (a-h,o-z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"

#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
      DIMENSION TG1(NASHT,NASHT),TG2(NASHT,NASHT,NASHT,NASHT)
      DIMENSION TG3(NTG3)
      DIMENSION CI1(MXCI),CI2(MXCI)
C Procedure for computing 1-body, 2-body, and 3-body transition
C density elements with active indices only.

C In: Wave functions CI1, with symmetry LSYM1, and CI2, with
C  symmetry LSYM2.
C
C Out: Transition density matrices, denoted here TG1, TG2 and TG3.
C Storage: TG1 and TG2 are simple two- and four-index arrays, and
C includes also such zeroes that are implied by symmetry.
C But TG3 is quite large, and while it is stored with zeroes, it
C is made more compact by the following addressing:

C <Psi1|E_tuvxyz|Psi2> is stored in TG3(ITG3) where
C    ITG3= ((i+1)*i*(i-1))/6 + (j*(j-1))/2 + k
C     i  = max(tu,vx,yz)
C     j  = mid(tu,vx,yz)
C     k  = min(tu,vx,yz)
C tu stands for the pair index tu= t + NASHT*(u-1), etc., and t is
C the usual active orbital number, when they are enumerated across
C all the symmetries (The ''absolute'' active index).

      CALL QENTER('MKTG3')

C Put in zeroes. Recognize special cases:
      OVL=1.0D0
      IF(NASHT.EQ.0) GOTO 999
      IF(LSYM1.NE.LSYM2) OVL=0.0D0
      CALL DCOPY_(NASHT**2,0.0D0,0,TG1,1)
      CALL DCOPY_(NASHT**4,0.0D0,0,TG2,1)
      CALL DCOPY_(NTG3,0.0D0,0,TG3,1)
      IF(NACTEL.EQ.0) GOTO 999

      IF(ISCF.EQ.0) GOTO 100

C -Special code for the closed-shell or hi-spin cases:
C ISCF=1 for closed-shell, =2 for hispin
      OCC=2.0D0
      IF(ISCF.EQ.2) OCC=1.0D0
      DO IT=1,NASHT
        TG1(IT,IT)=OCC
      END DO
      IF(NACTEL.EQ.1) GOTO 999
      DO IT=1,NASHT
       DO IU=1,NASHT
        TG2(IT,IT,IU,IU)=TG1(IT,IT)*TG1(IU,IU)
        IF(IU.EQ.IT) THEN
         TG2(IT,IT,IU,IU)=TG2(IT,IT,IU,IU)-TG1(IT,IU)
         ELSE
          TG2(IT,IU,IU,IT)=-TG1(IT,IT)
         END IF
        END DO
       END DO
      IF(NACTEL.EQ.2) GOTO 999
       DO IT1=1,NLEV
        DO IU1=1,NLEV
         IND1=IT1+NASHT*(IU1-1)
         DO IT2=1,NLEV
          DO IU2=1,IU1
           IND2=IT2+NASHT*(IU2-1)
           IF(IND2.GT.IND1) GOTO 199
           DO IT3=1,NLEV
            DO IU3=1,IU2
             IND3=IT3+NASHT*(IU3-1)
             IF(IND3.GT.IND2) GOTO 198
             VAL=TG1(IT1,IU1)*TG1(IT2,IU2)*TG1(IT3,IU3)

C Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Add here the necessary Kronecker deltas times 2-body matrix
C elements and lower, so we get a true normal-ordered density matrix
C element.

C <PSI1|E(T1,U1,T2,U2,T3,U3)|PSI2>
C = <PSI1|E(T1,U1)E(T2,U2)E(T3,U3)|PSI2>
C -D(T3,U2)*(TG2(T1,U1,T2,U3)+D(T2,U1)*TG1(T1,U3))
C -D(T2,U1)*TG2(T1,U2,T3,U3)
C -D(T3,U1)*TG2(T2,U2,T1,U3)

      IF(IT3.EQ.IU2) THEN
        VAL=VAL-TG2(IT1,IU1,IT2,IU3)
        IF(IT2.EQ.IU1) THEN
          VAL=VAL-TG1(IT1,IU3)
        END IF
      END IF
      IF(IT2.EQ.IU1) THEN
        VAL=VAL-TG2(IT1,IU2,IT3,IU3)
      END IF
      IF(IT3.EQ.IU1) THEN
        VAL=VAL-TG2(IT2,IU2,IT1,IU3)
      END IF

C VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
      ITG3=((IND1+1)*IND1*(IND1-1))/6+(IND2*(IND2-1))/2+IND3
      TG3(ITG3)=VAL


 198        CONTINUE
           END DO
          END DO
 199      CONTINUE
         END DO
        END DO
       END DO
      END DO
      GOTO 999

 100  CONTINUE
C Here, for regular CAS or RAS cases.

C Special pair index allows true RAS cases to be handled:
      CALL GETMEM('IPTOLEV','ALLO','INTE',LP2LEV,2*NASHT**2)
      LP2LEV1=LP2LEV
      LP2LEV2=LP2LEV+NASHT**2
      IP=0
C First, IL < JL pairs.
      DO IL=1,NLEV-1
       DO JL=IL+1,NLEV
        IP=IP+1
        IWORK(LP2LEV1-1+IP)=IL
        IWORK(LP2LEV2-1+IP)=JL
       END DO
      END DO
C Then, IL = JL pairs.
      DO IL=1,NLEV
        IP=IP+1
        IWORK(LP2LEV1-1+IP)=IL
        IWORK(LP2LEV2-1+IP)=IL
      END DO
C Last, IL > JL pairs.
      DO IL=2,NLEV
       DO JL=1,IL-1
        IP=IP+1
        IWORK(LP2LEV1-1+IP)=IL
        IWORK(LP2LEV2-1+IP)=JL
       END DO
      END DO
C If now any matrix element E(t1u1)E(t2u2)..E(tnun) is arranged
C such that the pair indices are non-decreasing, then the matrix
C element can be correctly computed by performing explicit
C excitations within the RAS space.
C But we also need the 'usual' pair index in order to use the
C packed addressing.

      NCI1=NCSF(LSYM1)
C Overlap:
      IF(LSYM1.EQ.LSYM2) OVL=DDOT_(NCI1,CI1,1,CI2,1)
C Allocate as many vectors as possible:
C Wishful thinking:
      NVECS=2*NASHT**2+1
C But what is really available?
      CALL GETMEM('DUMMY','MAX ','REAL',L,NTG3WRK)
      NTG3WRK=MIN(MXCI*NVECS,NTG3WRK)
      NVECS=NTG3WRK/MXCI
      NTG3WRK=NVECS*MXCI
C Find optimal subdivision of available vectors:
      NYZBUF=NINT(DBLE(NVECS-1)/DBLE(NASHT))
      NYZBUF=MAX(1,NYZBUF)
      NTUBUF=MIN(NASHT**2,NVECS-1-NYZBUF)
      NYZBUF=NVECS-1-NTUBUF
C Insufficient memory?
      IF(NTUBUF.LE.0) THEN
        WRITE(6,*)' Too little memory left for MKTG3.'
        WRITE(6,*)' Need at least 3 vectors of length MXCI=',MXCI
        CALL ABEND()
      END IF
      IF(NTUBUF.LE.(NASHT**2)/5) THEN
        WRITE(6,*)' WARNING: MKTG3 will be inefficient owing to'
        WRITE(6,*)' small memory.'
      END IF
      CALL GETMEM('TG3WRK','ALLO','REAL',LTG3WRK,NTG3WRK)
C And divide it up:
      LSGM1=LTG3WRK
      LTAU=LSGM1+NTUBUF*MXCI
      LSGM2=LTAU+MXCI

C Sectioning loops over pair indices IP3 (ket side):
      DO IP3STA=1,NASHT**2,NYZBUF
       IP3END=MIN(NASHT**2,IP3STA-1+NYZBUF)
C Compute a section of sigma vectors E(YZ)*PSI2 to memory:
       LTO=LSGM2
       DO IP3=IP3STA,IP3END
C Translate to levels in the SGUGA coupling order:
        IL=IWORK(LP2LEV1-1+IP3)
        JL=IWORK(LP2LEV2-1+IP3)
        IY=L2ACT(IL)
        IZ=L2ACT(JL)
        IYS=IASYM(IY)
        IZS=IASYM(IZ)
        ISSG2=MUL(MUL(IYS,IZS),LSYM2)
        CALL DCOPY_(MXCI,0.0D0,0,WORK(LTO),1)
C LTO is first element of Sigma2 = E(YZ) Psi2
        CALL SIGMA1_CP2(IL,JL,1.0D00,LSYM2,CI2,WORK(LTO),
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        IF(ISSG2.EQ.LSYM1) THEN
          TG1(IY,IZ)=DDOT_(NCI1,CI1,1,WORK(LTO),1)
        END IF
        LTO=LTO+MXCI
       END DO
C Sectioning loops over pair indices IP1 (bra side):
       DO IP1STA=IP3STA,NASHT**2,NTUBUF
        IP1END=MIN(NASHT**2,IP1STA-1+NTUBUF)
C Compute a section of sigma vectors E(UT)*PSI1 to memory:
        LTO=LSGM1
        DO IP1=IP1STA,IP1END
C Translate to levels:
         JL=IWORK(LP2LEV1-1+IP1)
         IL=IWORK(LP2LEV2-1+IP1)
         IT=L2ACT(IL)
         IU=L2ACT(JL)
         ITS=IASYM(IT)
         IUS=IASYM(IU)
         ISSG1=MUL(MUL(ITS,IUS),LSYM1)
         CALL DCOPY_(MXCI,0.0D0,0,WORK(LTO),1)
         CALL SIGMA1_CP2(IL,JL,1.0D00,LSYM1,CI1,WORK(LTO),
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
         LTO=LTO+MXCI
        END DO
C Now compute as many elements as possible:
        LFROM=LSGM2
        DO IP3=IP3STA,IP3END
         IY=L2ACT(IWORK(LP2LEV1-1+IP3))
         IZ=L2ACT(IWORK(LP2LEV2-1+IP3))
C LFROM will be start element of Sigma2=E(YZ) Psi2
         IYZ=IY+NASHT*(IZ-1)
         IYS=IASYM(IY)
         IZS=IASYM(IZ)
         ISSG2=MUL(MUL(IYS,IZS),LSYM2)
         DO IP2=IP3,IP1END
          IL=IWORK(LP2LEV1-1+IP2)
          JL=IWORK(LP2LEV2-1+IP2)
          IV=L2ACT(IL)
          IX=L2ACT(JL)
          IVX=IV+NASHT*(IX-1)
          IVS=IASYM(IV)
          IXS=IASYM(IX)
          ISTAU=MUL(MUL(IVS,IXS),ISSG2)
          NTAU=NCSF(ISTAU)
          CALL DCOPY_(MXCI,0.0D0,0,WORK(LTAU),1)
C LTAU  will be start element of Tau=E(VX) Sigma2=E(VX) E(YZ) Psi2
          CALL SIGMA1_CP2(IL,JL,1.0D00,ISSG2,WORK(LFROM),WORK(LTAU),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          IF(ISTAU.EQ.LSYM1) THEN
           TG2(IV,IX,IY,IZ)=DDOT_(NTAU,WORK(LTAU),1,CI1,1)
          END IF
          DO IP1=MAX(IP2,IP1STA),IP1END
           IT=L2ACT(IWORK(LP2LEV1-1+IP1))
           IU=L2ACT(IWORK(LP2LEV2-1+IP1))
           ITS=IASYM(IT)
           IUS=IASYM(IU)
           ISSG1=MUL(MUL(ITS,IUS),LSYM1)
           IF(ISSG1.EQ.ISTAU) THEN
            L=LSGM1+MXCI*(IP1-IP1STA)
            VAL=DDOT_(NTAU,WORK(LTAU),1,WORK(L),1)
            ITU=IT+NASHT*(IU-1)
C Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Code to put it in correct place:
            IF(ITU.LT.IVX) THEN
              IF(ITU.GE.IYZ) THEN
                JTU=IVX
                JVX=ITU
                JYZ=IYZ
              ELSE IF(IVX.LT.IYZ) THEN
                  JTU=IYZ
                  JVX=IVX
                  JYZ=ITU
              ELSE
                  JTU=IVX
                  JVX=IYZ
                  JYZ=ITU
              END IF
            ELSE
              IF(ITU.LT.IYZ) THEN
                JTU=IYZ
                JVX=ITU
                JYZ=IVX
              ELSE IF (IVX.GE.IYZ) THEN
                JTU=ITU
                JVX=IVX
                JYZ=IYZ
              ELSE
                JTU=ITU
                JVX=IYZ
                JYZ=IVX
              END IF
            END IF
            JTUVXYZ=((JTU+1)*JTU*(JTU-1))/6+(JVX*(JVX-1))/2+JYZ
            TG3(JTUVXYZ)=VAL

C End of symmetry requirement IF-clause:
           END IF
C End of IP1 loop.
          END DO
C End of IP2 loop.
         END DO
         LFROM=LFROM+MXCI
C End of IP3 loop.
        END DO
C End of IP1STA sectioning loop
       END DO
C End of IP3STA sectioning loop
      END DO
      CALL GETMEM('TG3WRK','FREE','REAL',LTG3WRK,NTG3WRK)
C Now the computed elements of TG2 contain <PSI1|E(IT1,IU1)E(IT2,IU2)|PSI2>
C and TG3 contains <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
C Add here the necessary Kronecker deltas times 2-body matrix
C elements and lower, so we get a true normal-ordered density matrix
C element.

C First, the 2-particle density matrix:
C <PSI1|E(T,U,V,X)|PSI2>  = <PSI1|E(TU)E(VX)|PSI2> - D(V,U)*TG2(T,U,V,X)
      DO IP1=1,NASHT**2
       IT=L2ACT(IWORK(LP2LEV1-1+IP1))
       IU=L2ACT(IWORK(LP2LEV2-1+IP1))
       DO IP2=1,IP1
        IV=L2ACT(IWORK(LP2LEV1-1+IP2))
        IX=L2ACT(IWORK(LP2LEV2-1+IP2))
        IF(IV.EQ.IU) TG2(IT,IU,IV,IX)=TG2(IT,IU,IV,IX)-TG1(IT,IX)
        TG2(IV,IX,IT,IU)=TG2(IT,IU,IV,IX)
       END DO
      END DO
C and then the 3-particle density matrix:
C <PSI1|E(T,U,V,X,Y,Z)|PSI2>  = <PSI1|E(TU)E(VX)E(YZ)|PSI2>
C -D(Y,X)*(TG2(T,U,V,Z)+D(V,U)*TG1(T,Z))
C -D(V,U)*TG2(T,X,Y,Z) C -D(Y,U)*TG2(V,X,T,Z)
      DO IP1=1,NASHT**2
       IT=L2ACT(IWORK(LP2LEV1-1+IP1))
       IU=L2ACT(IWORK(LP2LEV2-1+IP1))
       ITU=IT+NASHT*(IU-1)
       ITS=IASYM(IT)
       IUS=IASYM(IU)
       IS1=MUL(MUL(ITS,IUS),LSYM1)
       DO IP2=1,IP1
        IV=L2ACT(IWORK(LP2LEV1-1+IP2))
        IX=L2ACT(IWORK(LP2LEV2-1+IP2))
        IVX=IV+NASHT*(IX-1)
        IVS=IASYM(IV)
        IXS=IASYM(IX)
        IS2=MUL(MUL(IVS,IXS),IS1)
        DO IP3=1,IP2
         IY=L2ACT(IWORK(LP2LEV1-1+IP3))
         IZ=L2ACT(IWORK(LP2LEV2-1+IP3))
         IYS=IASYM(IY)
         IZS=IASYM(IZ)
         IS3=MUL(MUL(IYS,IZS),IS2)
         IF(IS3.EQ.LSYM2) THEN
          IYZ=IY+NASHT*(IZ-1)
          IF(ITU.LT.IVX) THEN
            IF(ITU.GE.IYZ) THEN
              JTU=IVX
              JVX=ITU
              JYZ=IYZ
            ELSE IF(IVX.LT.IYZ) THEN
                JTU=IYZ
                JVX=IVX
                JYZ=ITU
            ELSE
                JTU=IVX
                JVX=IYZ
                JYZ=ITU
            END IF
          ELSE
            IF(ITU.LT.IYZ) THEN
              JTU=IYZ
              JVX=ITU
              JYZ=IVX
            ELSE IF (IVX.GE.IYZ) THEN
              JTU=ITU
              JVX=IVX
              JYZ=IYZ
            ELSE
              JTU=ITU
              JVX=IYZ
              JYZ=IVX
            END IF
          END IF
          JTUVXYZ=((JTU+1)*JTU*(JTU-1))/6+(JVX*(JVX-1))/2+JYZ
          VAL=TG3(JTUVXYZ)
          IF(IY.EQ.IX) THEN
           VAL=VAL-TG2(IT,IU,IV,IZ)
           IF(IV.EQ.IU) THEN
            VAL=VAL-TG1(IT,IZ)
           END IF
          END IF
          IF(IV.EQ.IU) THEN
           VAL=VAL-TG2(IT,IX,IY,IZ)
          END IF
          IF(IY.EQ.IU) THEN
           VAL=VAL-TG2(IV,IX,IT,IZ)
          END IF
          TG3(JTUVXYZ)=VAL
         END IF
        END DO
       END DO
      END DO
      CALL GETMEM('IPTOLEV','FREE','INTE',LP2LEV,2*NASHT**2)

 999  CONTINUE
      CALL QEXIT('MKTG3')
      RETURN
      END
