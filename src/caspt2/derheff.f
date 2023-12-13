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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine DerHEff(CLag,VECROT)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "pt2_guga.fh"

      INTEGER IST,JST

      INTEGER I
      INTEGER NTG1,NTG2,NTG3
      INTEGER IDCI,LCI1,LCI2
      REAL*8 OVL,DUMMY(1)
C
      Dimension CLag(nConf,nState),VECROT(*)
C     return

C We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
C Note: Need proper allocation even if unused, sinced allocated
C arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      CALL GETMEM('DTG1','ALLO','REAL',LDTG1,NTG1)
      CALL GETMEM('DTG2','ALLO','REAL',LDTG2,NTG2)
      CALL GETMEM('DTG3','ALLO','REAL',LDTG3,NTG3)
      CALL DCOPY_(NTG1,[0.0D+00],0,WORK(LDTG1),1)
      CALL DCOPY_(NTG2,[0.0D+00],0,WORK(LDTG2),1)
      CALL DCOPY_(NTG3,[0.0D+00],0,WORK(LDTG3),1)
C
      !! OVL will contain the derivative contribution?
      !! It should be ignored
      OVL = 0.0D+00
      CALL DerHeffX(IVECW,IVECC,OVL,WORK(LDTG1),WORK(LDTG2),WORK(LDTG3))
C
      CALL GETMEM('MCCI1','ALLO','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','ALLO','REAL',LCI2,MXCI)
      CALL GETMEM('MCCI3','ALLO','REAL',LCI3,MXCI)
C
      IF(ISCF.EQ.0) THEN
        IDCI=IDTCEX
        JST = jState
        CALL DCOPY_(NCONF,[0.0D+00],0,WORK(LCI1),1)
        DO I=1,NSTATE
          IST = I
          IF (IST.EQ.JST) THEN
            CALL DDAFILE(LUCIEX,2,WORK(LCI2),NCONF,IDCI)
          Else If (ABS(VECROT(IST)).le.1.0d-12) Then
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,2,WORK(LCI3),NCONF,IDCI)
            CALL DAXPY_(NCONF,VECROT(IST),WORK(LCI3),1,WORK(LCI1),1)
          END IF
        END DO
C
        CALL DCOPY_(NCONF,[0.0D+00],0,WORK(LCI3),1)
        CALL DERTG3(.TRUE.,STSYM,STSYM,WORK(LCI1),WORK(LCI2),OVL,
     &              WORK(LDTG1),WORK(LDTG2),NTG3,WORK(LDTG3),
     &              WORK(LCI3),CLag(1,JST))
C
        DO I=1,NSTATE
          IST = I
          IF (IST.EQ.JST) THEN
            CYCLE
          Else If (ABS(VECROT(IST)).le.1.0d-12) Then
            CYCLE
          ELSE
            CALL DAXPY_(NCONF,VECROT(IST),WORK(LCI3),1,CLag(1,IST),1)
          END IF
        END DO
      END IF
C
      CALL GETMEM('MCCI1','FREE','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','FREE','REAL',LCI2,MXCI)
      CALL GETMEM('MCCI3','FREE','REAL',LCI3,MXCI)
C
      CALL GETMEM('DTG1','FREE','REAL',LDTG1,NTG1)
      CALL GETMEM('DTG2','FREE','REAL',LDTG2,NTG2)
      CALL GETMEM('DTG3','FREE','REAL',LDTG3,NTG3)
C
C
      RETURN
C
      End Subroutine DerHEff
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DerHeffX(IVEC,JVEC,OVL,DTG1,DTG2,DTG3)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
C Compute the coupling Hamiltonian element defined as
C     HEL = < ROOT1 | H * OMEGA | ROOT2 >
C assuming that IVEC contains a contravariant representation of
C H|ROOT1>, JVEC contains a contravariant representation of
C OMEGA|ROOT2>, and OVL, TG1, TG2, TG3 contain the overlap (normally
C expected to be 0 or 1) and active transition density matrices of ROOT1
C and ROOT2. See also subroutine TSVEC for explanations.

C SVC (March 2014): modification of original code to handle distributed
C RHS arrays. There is now a main HCOUP subroutine that loops over cases
C and irreps and gets access to the process-specific block of the RHS.
C The coupling for that block is computed by the subroutine HCOUP_BLK.

#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      Dimension DTG1(NASHT,NASHT)
      Dimension DTG2(NASHT,NASHT,NASHT,NASHT)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)
      Dimension DTG3(*)

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


C Sketch of procedure:
C  HEL=0.0D0
C  Loop over every (case/symmetry)-block.
C           If (No such vector block) Skip to end of loop
C           Allocate two places for this block, VEC1 and VEC2
C           Read VEC1 as IVEC component from file.
C           Read VEC2 as JVEC component from file.
C           Loop nest, computing
C              HEL := HEL + VEC1*GOM*VEC2
C           End of loop nest
C           Deallocate VEC1 and VEC2
C  End of loop.

      ! HEL=0.0D0
C     HECOMP=0.0D0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          ! HEBLK=0.0D0

          IF(NAS*NIS.EQ.0) GOTO 1
          IF(NIN.EQ.0) GOTO 1

          CALL RHS_ALLO (NAS,NIS,lg_V1)
          CALL RHS_ALLO (NAS,NIS,lg_V2)
          CALL RHS_READ (NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
          CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
          CALL RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
          CALL RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

          IF ((iLo1.NE.iLo2) .OR.
     &        (iHi1.NE.iHi2) .OR.
     &        (jLo1.NE.jLo2) .OR.
     &        (jHi1.NE.jHi2)) THEN
            WRITE(6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
            CALL ABEND()
          END IF

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                      DBL_MB(MV1),DBL_MB(MV2),OVL,
     &                      DTG1,DTG2,DTG3)
          ELSE
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                      WORK(MV1),WORK(MV2),OVL,
     &                      DTG1,DTG2,DTG3)
          END IF
#else
          CALL DerHEffX_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                    WORK(MV1),WORK(MV2),OVL,
     &                    DTG1,DTG2,DTG3)
#endif
C
          CALL RHS_RELEASE (lg_V1,IASTA1,IAEND1,IISTA1,IIEND1)
          CALL RHS_RELEASE (lg_V2,IASTA2,IAEND2,IISTA2,IIEND2)
          CALL RHS_FREE (NAS,NIS,lg_V1)
          CALL RHS_FREE (NAS,NIS,lg_V2)

 1        CONTINUE
C         HECOMP(ICASE,ISYM)=HEBLK
C         HEL=HEL+HEBLK
        END DO
      END DO

C Sum-reduce the per-process contributions
C     CALL GADGOP_SCAL(HEL,'+')
C     NHECOMP=14*9
C     CALL GADGOP(HECOMP,NHECOMP,'+')

C     IF(IPRGLB.GE.DEBUG) THEN
C       DO ICASE=1,13
C         SUMSYM=0.0D0
C         DO ISYM=1,NSYM
C           SUMSYM=SUMSYM+HECOMP(ICASE,ISYM)
C         END DO
C         HECOMP(ICASE,NSYM+1)=SUMSYM
C       END DO

C       DO ISYM=1,NSYM+1
C         SUMCASE=0.0D0
C         DO ICASE=1,13
C           SUMCASE=SUMCASE+HECOMP(ICASE,ISYM)
C         END DO
C         HECOMP(14,ISYM)=SUMCASE
C       END DO

C       WRITE(6,'(20a4)')('----',i=1,20)
C       WRITE(6,*)'HCOUP: The contributions to the Hamiltonian coupling'
C       WRITE(6,*)' elements, by case and by symmetry label.'
C       DO IC=1,13
C         WRITE(6,'(1X,A8,9F12.8)')
C    &      CASES(IC),(HECOMP(IC,IS),IS=1,NSYM+1)
C       END DO
C       CALL XFLUSH(6)
C       WRITE(6,'(1X,A8,9F12.8)')
C    &    'Summed: ', (HECOMP(14,IS),IS=1,NSYM+1)
C       WRITE(6,*)
C     END IF


      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DerHEffX_BLK(ICASE,ISYM,NAS,IISTA,IIEND,V1,V2,OVL,
     &                        DTG1,DTG2,DTG3)
      USE SUPERINDEX
C Compute a contribution to the coupling Hamiltonian element (HEL)
C defined as HEL = < ROOT1 | H * OMEGA | ROOT2 >. The contribution
C arises from the block V_(A,I), with A=1,NAS and I=IISTA,IIEND,
C with A the active superindex and I the inactive superindex. Since
C the inactive superindex is partitioned over processes, each process
C only computes part of the HEL value, which is then sum reduced in the
C calling subroutine.
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "eqsolv.fh"

      DIMENSION V1(*), V2(*)

      Dimension DTG1(NASHT,NASHT)
      Dimension DTG2(NASHT,NASHT,NASHT,NASHT)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)
      Dimension DTG3(*)


      ! HEBLK=0.0D0

      IF (IISTA.LE.0) RETURN

      NISBLK=IIEND-IISTA+1
      SELECT CASE (ICASE)
************************************************************************
      CASE (1)
        DO IAS=1,NAS
          IASABS=NTUVES(ISYM)+IAS
          ITABS=MTUV(1,IASABS)
          IUABS=MTUV(2,IASABS)
          IVABS=MTUV(3,IASABS)
          DO JAS=1,NAS
            JASABS=NTUVES(ISYM)+JAS
            IXABS=MTUV(1,JASABS)
            IYABS=MTUV(2,JASABS)
            IZABS=MTUV(3,JASABS)
C Compute and use SA(ITABS IUABS IVABS, IXABS IYABS IZABS)
C Formulae used:
C  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz -
C         - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
C Gvuxtyz is stored using full permutation symmetry of three pairs
C (vu),(xt), and (yz):
            IND1=IVABS+NASHT*(IUABS-1)
            IND2=IXABS+NASHT*(ITABS-1)
            IND3=IYABS+NASHT*(IZABS-1)
            IF(IND2.GT.IND3) THEN
              IF(IND1.GT.IND2) THEN
                JND1=IND1
                JND2=IND2
                JND3=IND3
              ELSE IF(IND1.GT.IND3) THEN
                JND1=IND2
                JND2=IND1
                JND3=IND3
              ELSE
                JND1=IND2
                JND2=IND3
                JND3=IND1
              END IF
            ELSE
              IF(IND1.GT.IND3) THEN
                JND1=IND1
                JND2=IND3
                JND3=IND2
              ELSE IF(IND1.GT.IND2) THEN
                JND1=IND3
                JND2=IND1
                JND3=IND2
              ELSE
                JND1=IND3
                JND2=IND2
                JND3=IND1
              END IF
            END IF
            ITG3=((JND1+1)*JND1*(JND1-1))/6+(JND2*(JND2-1))/2+JND3
C  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz -
C         - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
C Compute TMP=Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
C           TMP=TG3(ITG3)
C           IF(IYABS.EQ.IUABS) THEN
C             TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
C           END IF
C           IF(IYABS.EQ.ITABS) THEN
C             TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
C             IF(IXABS.EQ.IUABS) THEN
C               TMP=TMP+TG1(IVABS,IZABS)
C             END IF
C           END IF
C           IF(IXABS.EQ.IUABS) THEN
C             TMP=TMP+TG2(IVABS,ITABS,IYABS,IZABS)
C           END IF
C SA is the negative of this, and then some correction:
C           SA=-TMP
C           IF(IXABS.EQ.ITABS) THEN
C             SA=SA+2.0D0*TG2(IVABS,IUABS,IYABS,IZABS)
C             IF(IYABS.EQ.IUABS) THEN
C               SA=SA+2.0D0*TG1(IVABS,IZABS)
C             END IF
C           END IF
C SA has been computed.

C           HEBLK=HEBLK+SA*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            IF(IXABS.EQ.ITABS) THEN
              DTG2(IVABS,IUABS,IYABS,IZABS)
     *          = DTG2(IVABS,IUABS,IYABS,IZABS) + 2.0D+00*VAL
              IF(IYABS.EQ.IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + 2.0D+00*VAL
              END IF
            END IF
            VAL = -VAL
            DTG3(ITG3) = DTG3(ITG3) + VAL
            IF(IYABS.EQ.IUABS) THEN
              DTG2(IVABS,IZABS,IXABS,ITABS)
     *          = DTG2(IVABS,IZABS,IXABS,ITABS) + VAL
            END IF
            IF(IYABS.EQ.ITABS) THEN
              DTG2(IVABS,IUABS,IXABS,IZABS)
     *          = DTG2(IVABS,IUABS,IXABS,IZABS) + VAL
              IF(IXABS.EQ.IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + VAL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              DTG2(IVABS,ITABS,IYABS,IZABS)
     *          = DTG2(IVABS,ITABS,IYABS,IZABS) + VAL
            END IF
          END DO
        END DO
************************************************************************
      CASE(4)
        DO IAS=1,NAS
          IASABS=NTUVES(ISYM)+IAS
          IXABS=MTUV(1,IASABS)
          IUABS=MTUV(2,IASABS)
          IVABS=MTUV(3,IASABS)
          DO JAS=1,NAS
            JASABS=NTUVES(ISYM)+JAS
            ITABS=MTUV(1,JASABS)
            IYABS=MTUV(2,JASABS)
            IZABS=MTUV(3,JASABS)
C Compute and use SC(IXABS IUABS IVABS, ITABS IYABS IZABS)
C In SBMAT, the formula is written as SC(tuv,xyz)
C    = Gvutxyz +dyu Gvztx + dyx Gvutz + dtu Gvxyz + dtu dyx Gvz
C Rewritten, in order to reuse same quantities as in SA:
C  SC(xuv,tyz)
C    = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
C Gvuxtyz is stored using full permutation symmetry of three pairs
C (vu),(xt), and (yz):
            IND1=IVABS+NASHT*(IUABS-1)
            IND2=IXABS+NASHT*(ITABS-1)
            IND3=IYABS+NASHT*(IZABS-1)
            IF(IND2.GT.IND3) THEN
              IF(IND1.GT.IND2) THEN
                JND1=IND1
                JND2=IND2
                JND3=IND3
              ELSE IF(IND1.GT.IND3) THEN
                JND1=IND2
                JND2=IND1
                JND3=IND3
              ELSE
                JND1=IND2
                JND2=IND3
                JND3=IND1
              END IF
            ELSE
              IF(IND1.GT.IND3) THEN
                JND1=IND1
                JND2=IND3
                JND3=IND2
              ELSE IF(IND1.GT.IND2) THEN
                JND1=IND3
                JND2=IND1
                JND3=IND2
              ELSE
                JND1=IND3
                JND2=IND2
                JND3=IND1
              END IF
            END IF
            ITG3=((JND1+1)*JND1*(JND1-1))/6+(JND2*(JND2-1))/2+JND3
C  SC(xuv,tyz) (rewritten, swapping x and t)
C    = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
C           TMP=TG3(ITG3)
C           IF(IYABS.EQ.IUABS) THEN
C             TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
C           END IF
C           IF(IYABS.EQ.ITABS) THEN
C             TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
C             IF(IXABS.EQ.IUABS) THEN
C               TMP=TMP+TG1(IVABS,IZABS)
C             END IF
C           END IF
C           IF(IXABS.EQ.IUABS) THEN
C             TMP=TMP+TG2(IVABS,ITABS,IYABS,IZABS)
C           END IF
C           SC= TMP

C           HEBLK=HEBLK+SC*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG3(ITG3) = DTG3(ITG3) + VAL
            IF(IYABS.EQ.IUABS) THEN
              DTG2(IVABS,IZABS,IXABS,ITABS)
     *          = DTG2(IVABS,IZABS,IXABS,ITABS) + VAL
            END IF
            IF(IYABS.EQ.ITABS) THEN
              DTG2(IVABS,IUABS,IXABS,IZABS)
     *          = DTG2(IVABS,IUABS,IXABS,IZABS) + VAL
              IF(IXABS.EQ.IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + VAL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              DTG2(IVABS,ITABS,IYABS,IZABS)
     *          = DTG2(IVABS,ITABS,IYABS,IZABS) + VAL
            END IF
          END DO
        END DO
************************************************************************
      CASE(2)
        DO IAS=1,NAS
          IASABS=NTGEUES(ISYM)+IAS
          ITABS=MTGEU(1,IASABS)
          IUABS=MTGEU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGEUES(ISYM)+JAS
            IXABS=MTGEU(1,JASABS)
            IYABS=MTGEU(2,JASABS)
C Formulae used:
C    SB(tu,xy)=
C    = 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu
C      -4dxu dyt + 2dxu Gyt
C    SB(tu,yx)=
C    = 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu
C      -4dyu dxt + 2dyu Gxt
C    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
C    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)
C           SBtuxy=2.0d0*TG2(IXABS,ITABS,IYABS,IUABS)
C           SBtuyx=2.0d0*TG2(IYABS,ITABS,IXABS,IUABS)
C           IF(IXABS.EQ.ITABS) THEN
C             SBtuxy=SBtuxy-4.0d0*TG1(IYABS,IUABS)
C             SBtuyx=SBtuyx+2.0d0*TG1(IYABS,IUABS)
C             IF(IYABS.EQ.IUABS) THEN
C               SBtuxy=SBtuxy+8.0d0*OVL
C               SBtuyx=SBtuyx-4.0d0*OVL
C             END IF
C           END IF
C           IF(IYABS.EQ.IUABS) THEN
C             SBtuxy=SBtuxy-4.0d0*TG1(IXABS,ITABS)
C             SBtuyx=SBtuyx+2.0d0*TG1(IXABS,ITABS)
C           END IF
C           IF(IYABS.EQ.ITABS) THEN
C             SBtuxy=SBtuxy+2.0d0*TG1(IXABS,IUABS)
C             SBtuyx=SBtuyx-4.0d0*TG1(IXABS,IUABS)
C             IF(IXABS.EQ.IUABS) THEN
C               SBtuxy=SBtuxy-4.0d0*OVL
C               SBtuyx=SBtuyx+8.0d0*OVL
C             END IF
C           END IF
C           IF(IXABS.EQ.IUABS) THEN
C             SBtuxy=SBtuxy+2.0d0*TG1(IYABS,ITABS)
C             SBtuyx=SBtuyx-4.0d0*TG1(IYABS,ITABS)
C           END IF

C           SBP=SBtuxy + SBtuyx

C           HEBLK=HEBLK+SBP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(IXABS,ITABS,IYABS,IUABS)
     *        = DTG2(IXABS,ITABS,IYABS,IUABS) + 2.0D+00*VAL
            DTG2(IYABS,ITABS,IXABS,IUABS)
     *        = DTG2(IYABS,ITABS,IXABS,IUABS) + 2.0D+00*VAL
            IF(IXABS.EQ.ITABS) THEN
              DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS)
     *          - 4.0D+00*VAL + 2.0D+00*VAL
              IF(IYABS.EQ.IUABS) THEN
                OVL = OVL + 8.0D+00*VAL - 4.0D+00*VAL
              END IF
            END IF
            IF(IYABS.EQ.IUABS) THEN
              DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)
     *          - 4.0D+00*VAL + 2.0D+00*VAL
            END IF
            IF(IYABS.EQ.ITABS) THEN
              DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS)
     *          + 2.0D+00*VAL - 4.0D+00*VAL
              IF(IXABS.EQ.IUABS) THEN
                OVL = OVL - 4.0D+00*VAL + 8.0D+00*VAL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS)
     *          + 2.0D+00*VAL - 4.0D+00*VAL
            END IF
          END DO
        END DO
************************************************************************
      CASE(3)
        DO IAS=1,NAS
          IASABS=NTGTUES(ISYM)+IAS
          ITABS=MTGTU(1,IASABS)
          IUABS=MTGTU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGTUES(ISYM)+JAS
            IXABS=MTGTU(1,JASABS)
            IYABS=MTGTU(2,JASABS)
C Formulae used:
C    SB(tu,xy)=
C    = 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu
C      -4dxu dyt + 2dxu Gyt
C    SB(tu,yx)=
C    = 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu
C      -4dyu dxt + 2dyu Gxt
C    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
C    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)
C           SBtuxy=2.0d0*TG2(IXABS,ITABS,IYABS,IUABS)
C           SBtuyx=2.0d0*TG2(IYABS,ITABS,IXABS,IUABS)
C           IF(IXABS.EQ.ITABS) THEN
C             SBtuxy=SBtuxy-4.0d0*TG1(IYABS,IUABS)
C             SBtuyx=SBtuyx+2.0d0*TG1(IYABS,IUABS)
C             IF(IYABS.EQ.IUABS) THEN
C               SBtuxy=SBtuxy+8.0d0*OVL
C               SBtuyx=SBtuyx-4.0d0*OVL
C             END IF
C           END IF
C           IF(IYABS.EQ.IUABS) THEN
C             SBtuxy=SBtuxy-4.0d0*TG1(IXABS,ITABS)
C             SBtuyx=SBtuyx+2.0d0*TG1(IXABS,ITABS)
C           END IF
C           IF(IYABS.EQ.ITABS) THEN
C             SBtuxy=SBtuxy+2.0d0*TG1(IXABS,IUABS)
C             SBtuyx=SBtuyx-4.0d0*TG1(IXABS,IUABS)
C             IF(IXABS.EQ.IUABS) THEN
C               SBtuxy=SBtuxy-4.0d0*OVL
C               SBtuyx=SBtuyx+8.0d0*OVL
C             END IF
C           END IF
C           IF(IXABS.EQ.IUABS) THEN
C             SBtuxy=SBtuxy+2.0d0*TG1(IYABS,ITABS)
C             SBtuyx=SBtuyx-4.0d0*TG1(IYABS,ITABS)
C           END IF

C           SBM=SBtuxy - SBtuyx

C           HEBLK=HEBLK+SBM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(IXABS,ITABS,IYABS,IUABS)
     *        = DTG2(IXABS,ITABS,IYABS,IUABS) + 2.0D+00*VAL
            DTG2(IYABS,ITABS,IXABS,IUABS)
     *        = DTG2(IYABS,ITABS,IXABS,IUABS) - 2.0D+00*VAL
            IF(IXABS.EQ.ITABS) THEN
              DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS)
     *          - 4.0D+00*VAL - 2.0D+00*VAL
              IF(IYABS.EQ.IUABS) THEN
                OVL = OVL + 8.0D+00*VAL + 4.0D+00*VAL
              END IF
            END IF
            IF(IYABS.EQ.IUABS) THEN
              DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)
     *          - 4.0D+00*VAL - 2.0D+00*VAL
            END IF
            IF(IYABS.EQ.ITABS) THEN
              DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS)
     *          + 2.0D+00*VAL + 4.0D+00*VAL
              IF(IXABS.EQ.IUABS) THEN
                OVL = OVL - 4.0D+00*VAL - 8.0D+00*VAL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS)
     *          + 2.0D+00*VAL + 4.0D+00*VAL
            END IF
          END DO
        END DO
************************************************************************
      CASE(5)
        NAS1=NAS/2
        DO IAS1=1,NAS1
          IAS2=IAS1+NAS1
          IASABS=NTUES(ISYM)+IAS1
          ITABS=MTU(1,IASABS)
          IUABS=MTU(2,IASABS)
          DO JAS1=1,NAS1
            JAS2=JAS1+NAS1
            JASABS=NTUES(ISYM)+JAS1
            IXABS=MTU(1,JASABS)
            IYABS=MTU(2,JASABS)
C Formulae used:
C    SD11(tu1,xy1)=2*(Gutxy + dtx Guy)
C    SD12(tu2,xy1)= -(Gutxy + dtx Guy)
C    SD21(tu2,xy1)= -(Gutxy + dtx Guy)
C    SD22(tu2,xy2)= -Gxtuy +2*dtx Guy
C           GUTXY= TG2(IUABS,ITABS,IXABS,IYABS)
C           SD11=2.0D0*GUTXY
C           SD12= -GUTXY
C           SD21= -GUTXY
C           SD22= -TG2(IXABS,ITABS,IUABS,IYABS)
C           IF(ITABS.EQ.IXABS) THEN
C             GUY=TG1(IUABS,IYABS)
C             SD11=SD11+2.0D0*GUY
C             SD12=SD12 -GUY
C             SD21=SD21 -GUY
C             SD22=SD22+2.0D0*GUY
C           END IF

C           HEBLK=HEBLK+SD11*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS1),NAS)
C           HEBLK=HEBLK+SD12*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS1),NAS)
C           HEBLK=HEBLK+SD21*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS2),NAS)
C           HEBLK=HEBLK+SD22*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS2),NAS)
            VAL11 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS1),NAS)
            VAL12 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS2),NAS)
            VAL21 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS1),NAS)
            VAL22 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS2),NAS)
            DTG2(IUABS,ITABS,IXABS,IYABS)
     *        = DTG2(IUABS,ITABS,IXABS,IYABS)
     *        + 2.0D+00*VAL11 - VAL12 - VAL21
            DTG2(IXABS,ITABS,IUABS,IYABS)
     *        = DTG2(IXABS,ITABS,IUABS,IYABS) - VAL22
            IF(ITABS.EQ.IXABS) THEN
              DTG1(IUABS,IYABS) = DTG1(IUABS,IYABS)
     *          + 2.0D+00*VAL11 - VAL12 - VAL21 + 2.0D+00*VAL22
            END IF
          END DO
        END DO
************************************************************************
      CASE(6)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SE(t,x)=2*dxt - Dxt
C           SE=-TG1(IXABS,ITABS)
C           IF(IXABS.EQ.ITABS) SE=SE+2.0d0*OVL
C           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - VAL
            IF(IXABS.EQ.ITABS) OVL=OVL+2.0D+00*VAL
          END DO
        END DO
************************************************************************
      CASE(7)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SE(t,x)=2*dxt - Dxt
C           SE=-TG1(IXABS,ITABS)
C           IF(IXABS.EQ.ITABS) SE=SE+2.0d0*OVL
C           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - VAL
            IF(IXABS.EQ.ITABS) OVL=OVL+2.0D+00*VAL
          END DO
        END DO
************************************************************************
      CASE(8)
C ========================================================
C Compute and use SFP(ITABS IUABS , IXABS IYABS)
C and (later, similar) SFM(ITABS IUABS , IXABS IYABS)
        DO IAS=1,NAS
          IASABS=NTGEUES(ISYM)+IAS
          ITABS=MTGEU(1,IASABS)
          IUABS=MTGEU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGEUES(ISYM)+JAS
            IXABS=MTGEU(1,JASABS)
            IYABS=MTGEU(2,JASABS)
C Formulae used:
C    SF(tu,xy)= 2 Gtxuy
C    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
C    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
C           SFtuxy=2.0d0*TG2(ITABS,IXABS,IUABS,IYABS)
C           SFtuyx=2.0d0*TG2(ITABS,IYABS,IUABS,IXABS)

C           SFP=SFtuxy + SFtuyx
C           HEBLK=HEBLK+SFP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(ITABS,IXABS,IUABS,IYABS)
     *        = DTG2(ITABS,IXABS,IUABS,IYABS) + 2.0D+00*VAL
            DTG2(ITABS,IYABS,IUABS,IXABS)
     *        = DTG2(ITABS,IYABS,IUABS,IXABS) + 2.0D+00*VAL
          END DO
        END DO
************************************************************************
      CASE(9)
C ========================================================
C Compute and use SFM(ITABS IUABS, IXABS ,IYABS)
        DO IAS=1,NAS
          IASABS=NTGTUES(ISYM)+IAS
          ITABS=MTGTU(1,IASABS)
          IUABS=MTGTU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGTUES(ISYM)+JAS
            IXABS=MTGTU(1,JASABS)
            IYABS=MTGTU(2,JASABS)
C Formulae used:
C    SF(tu,xy)= 4 Ptxuy
C    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
C    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
C           SFtuxy=2.0d0*TG2(ITABS,IXABS,IUABS,IYABS)
C           SFtuyx=2.0d0*TG2(ITABS,IYABS,IUABS,IXABS)

C           SFM=SFtuxy - SFtuyx
C           HEBLK=HEBLK+SFM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(ITABS,IXABS,IUABS,IYABS)
     *        = DTG2(ITABS,IXABS,IUABS,IYABS) + 2.0D+00*VAL
            DTG2(ITABS,IYABS,IUABS,IXABS)
     *        = DTG2(ITABS,IYABS,IUABS,IXABS) - 2.0D+00*VAL
          END DO
        END DO
************************************************************************
C CASES GP, GM
C Compute and use SG(ITABS , IXABS) (Same for cases GP and GM)
************************************************************************
      CASE(10)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SG(t,x)= Gtx
C           SG= TG1(ITABS,IXABS)

C           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)
     *        + DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(11)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SG(t,x)= Gtx
C           SG= TG1(ITABS,IXABS)

C           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)
     *        + DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(12)
        OVL = OVL + DDOT_(NAS*NISBLK,V2,1,V1,1)
        IF(ABS(OVL).GE.1.0D-12) THEN
C         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      CASE(13)
        OVL = OVL + DDOT_(NAS*NISBLK,V2,1,V1,1)
        IF(ABS(OVL).GE.1.0D-12) THEN
C         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      END SELECT
      Return
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DERTG3(DOG3,LSYM1,LSYM2,CI1,CI2,OVL,DTG1,DTG2,NTG3,
     *                  DTG3,CLAG1,CLAG2)
      IMPLICIT REAL*8 (a-h,o-z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
      DIMENSION DTG1(NASHT,NASHT),DTG2(NASHT,NASHT,NASHT,NASHT)
      DIMENSION DTG3(NTG3)
      DIMENSION CI1(MXCI),CI2(MXCI)
      DIMENSION CLAG1(MXCI),CLAG2(MXCI)
      LOGICAL   DOG3
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


C Put in zeroes. Recognize special cases:
C     OVL=1.0D0
      IF(NASHT.EQ.0) GOTO 999
C     IF(LSYM1.NE.LSYM2) OVL=0.0D0
C     CALL DCOPY_(NASHT**2,[0.0D0],0,TG1,1)
C     CALL DCOPY_(NASHT**4,[0.0D0],0,TG2,1)
C     CALL DCOPY_(NTG3,[0.0D0],0,TG3,1)
      IF(NACTEL.EQ.0) GOTO 999

      IF(ISCF.EQ.0) GOTO 100
      write (6,*) "Here is the special case"
      write (6,*) "not yet"
      call abend

C -Special code for the closed-shell or hi-spin cases:
C ISCF=1 for closed-shell, =2 for hispin
      ! OCC=2.0D0
      ! IF(ISCF.EQ.2) OCC=1.0D0
      DO IT=1,NASHT
      ! TG1(IT,IT)=OCC
      END DO
      IF(NACTEL.EQ.1) GOTO 999
      DO IT=1,NASHT
       DO IU=1,NASHT
      ! TG2(IT,IT,IU,IU)=TG1(IT,IT)*TG1(IU,IU)
      ! IF(IU.EQ.IT) THEN
      !  TG2(IT,IT,IU,IU)=TG2(IT,IT,IU,IU)-TG1(IT,IU)
      !  ELSE
      !   TG2(IT,IU,IU,IT)=-TG1(IT,IT)
      !  END IF
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
C            VAL=TG1(IT1,IU1)*TG1(IT2,IU2)*TG1(IT3,IU3)

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
C       VAL=VAL-TG2(IT1,IU1,IT2,IU3)
        IF(IT2.EQ.IU1) THEN
C         VAL=VAL-TG1(IT1,IU3)
        END IF
      END IF
      IF(IT2.EQ.IU1) THEN
C       VAL=VAL-TG2(IT1,IU2,IT3,IU3)
      END IF
      IF(IT3.EQ.IU1) THEN
C       VAL=VAL-TG2(IT2,IU2,IT1,IU3)
      END IF

C VAL is now =<PSI1|E(IT1,IU1,IT2,IU2,IT3,IU3)|PSI2>
      ! ITG3=((IND1+1)*IND1*(IND1-1))/6+(IND2*(IND2-1))/2+IND3
C     TG3(ITG3)=VAL


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

C First, the 3-particle density matrix:
C <PSI1|E(T,U,V,X,Y,Z)|PSI2>  = <PSI1|E(TU)E(VX)E(YZ)|PSI2>
C -D(Y,X)*(TG2(T,U,V,Z)+D(V,U)*TG1(T,Z))
C -D(V,U)*TG2(T,X,Y,Z) C -D(Y,U)*TG2(V,X,T,Z)
      IF (DOG3) THEN
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
          VAL=DTG3(JTUVXYZ)
          IF(IY.EQ.IX) THEN
           DTG2(IT,IU,IV,IZ)=DTG2(IT,IU,IV,IZ)-VAL
           IF(IV.EQ.IU) THEN
            DTG1(IT,IZ)=DTG1(IT,IZ)-VAL
           END IF
          END IF
          IF(IV.EQ.IU) THEN
           DTG2(IT,IX,IY,IZ)=DTG2(IT,IX,IY,IZ)-VAL
          END IF
          IF(IY.EQ.IU) THEN
           DTG2(IV,IX,IT,IZ)=DTG2(IV,IX,IT,IZ)-VAL
          END IF
         END IF
        END DO
       END DO
      END DO
      END IF
C
C Then, the 2-particle density matrix:
C <PSI1|E(T,U,V,X)|PSI2>  = <PSI1|E(TU)E(VX)|PSI2> - D(V,U)*TG2(T,U,V,X)
      DO IP1=1,NASHT**2
       IT=L2ACT(IWORK(LP2LEV1-1+IP1))
       IU=L2ACT(IWORK(LP2LEV2-1+IP1))
       DO IP2=1,IP1
        IV=L2ACT(IWORK(LP2LEV1-1+IP2))
        IX=L2ACT(IWORK(LP2LEV2-1+IP2))
        IF (IP1.ne.IP2) Then
          DTG2(IT,IU,IV,IX)=DTG2(IT,IU,IV,IX)+DTG2(IV,IX,IT,IU)
          DTG2(IV,IX,IT,IU) = 0.0D+00
        End If
        IF(IV.EQ.IU) DTG1(IT,IX)=DTG1(IT,IX)-DTG2(IT,IU,IV,IX)
       END DO
      END DO
C
C If now any matrix element E(t1u1)E(t2u2)..E(tnun) is arranged
C such that the pair indices are non-decreasing, then the matrix
C element can be correctly computed by performing explicit
C excitations within the RAS space.
C But we also need the 'usual' pair index in order to use the
C packed addressing.

      NCI1=NCSF(LSYM1)
C Overlap:
C     IF(LSYM1.EQ.LSYM2) OVL=DDOT_(NCI1,CI1,1,CI2,1)
      IF(LSYM1.EQ.LSYM2) THEN
        Call DaXpY_(NCI1,OVL,CI1,1,CLAG2,1)
        Call DaXpY_(NCI1,OVL,CI2,1,CLAG1,1)
      END IF
C     write (*,*) "overlap = ",DDOT_(NCI1,CI1,1,CI2,1)
C Allocate as many vectors as possible:
C Wishful thinking:
      NVECS=2*NASHT**2+1
C But what is really available?
      CALL GETMEM('DUMMY','MAX ','REAL',L,NTG3WRK)
      NTG3WRK=NTG3WRK/2
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
        WRITE(6,*)' Need at least 6 vectors of length MXCI=',MXCI
        CALL ABEND()
      END IF
      IF(NTUBUF.LE.(NASHT**2)/5) THEN
        WRITE(6,*)' WARNING: MKTG3 will be inefficient owing to'
        WRITE(6,*)' small memory.'
      END IF
      CALL GETMEM('TG3WRK','ALLO','REAL',LTG3WRK,NTG3WRK)
      CALL GETMEM('BUF1','ALLO','REAL',LBUF1,MXCI)
C
      CALL GETMEM('DTU','ALLO','REAL',LDTU,MXCI*NTUBUF)
      CALL GETMEM('DYZ','ALLO','REAL',LDYZ,MXCI*NYZBUF)
C
C And divide it up:
      !! LSGM1: NTUBUF vectors
      !! LTAU : 1 vector
      !! LSGM2: NYZBUF vectors
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
        CALL DCOPY_(MXCI,[0.0D0],0,WORK(LTO),1)
C LTO is first element of Sigma2 = E(YZ) Psi2
        CALL SIGMA1_CP2(IL,JL,1.0D00,LSYM2,CI2,WORK(LTO),
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        IF(ISSG2.EQ.LSYM1.AND.DTG1(IY,IZ).NE.0.0D+00) THEN
          !! It is possible to calculate the contribution using
          !! DGEMV, but DAXPY seems to be faster than DGEMV
          Call DaXpY_(NCI1,DTG1(IY,IZ),WORK(LTO),1,CLAG1,1)
        END IF
        LTO=LTO+MXCI
       END DO
C
       CALL DCopy_(MXCI*NYZBUF,[0.0D+00],0,WORK(LDYZ),1)
C Sectioning loops over pair indices IP1 (bra side):
       DO IP1STA=IP3STA,NASHT**2,NTUBUF
        IP1END=MIN(NASHT**2,IP1STA-1+NTUBUF)
C Compute a section of sigma vectors E(UT)*PSI1 to memory:
        LTO=LSGM1
        !! <Psi1|E(TU)
        DO IP1=IP1STA,IP1END
C Translate to levels:
         JL=IWORK(LP2LEV1-1+IP1)
         IL=IWORK(LP2LEV2-1+IP1)
         IT=L2ACT(IL)
         IU=L2ACT(JL)
         ITS=IASYM(IT)
         IUS=IASYM(IU)
         ISSG1=MUL(MUL(ITS,IUS),LSYM1)
         CALL DCOPY_(MXCI,[0.0D0],0,WORK(LTO),1)
         CALL SIGMA1_CP2(IL,JL,1.0D00,LSYM1,CI1,WORK(LTO),
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
         IF (ISSG1.EQ.LSYM1.AND.DTG1(IU,IT).NE.0.0D+00
     &       .AND.IP3STA.EQ.1) THEN
          Call DaXpY_(NCI1,DTG1(IU,IT),WORK(LTO),1,CLAG2,1)
         END IF
         LTO=LTO+MXCI
        END DO
C
        CALL DCopy_(MXCI*NTUBUF,[0.0D+00],0,WORK(LDTU),1)
C Now compute as many elements as possible:
        LFROM=LSGM2
        LFROMD=LDYZ
        DO IP3=IP3STA,IP3END
         IY=L2ACT(IWORK(LP2LEV1-1+IP3))
         IZ=L2ACT(IWORK(LP2LEV2-1+IP3))
C LFROM will be start element of Sigma2=E(YZ) Psi2
         IYZ=IY+NASHT*(IZ-1)
         IYS=IASYM(IY)
         IZS=IASYM(IZ)
         ISSG2=MUL(MUL(IYS,IZS),LSYM2)
         IM=IWORK(LP2LEV1-1+IP3)
         JM=IWORK(LP2LEV2-1+IP3)
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
          CALL DCOPY_(MXCI,[0.0D0],0,WORK(LTAU),1)
C LTAU  will be start element of Tau=E(VX) Sigma2=E(VX) E(YZ) Psi2
          !! LTAU = EvxEyz|Psi2>
          CALL SIGMA1_CP2(IL,JL,1.0D00,ISSG2,WORK(LFROM),WORK(LTAU),
     &     IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &     IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &     WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          IF(ISTAU.EQ.LSYM1.AND.DTG2(IV,IX,IY,IZ).NE.0.0D+00) THEN
C          DTG2(IV,IX,IY,IZ)=DDOT_(NTAU,WORK(LTAU),1,CI1,1)
           !! For left derivative: <I|Evx Eyz|Psi2>
           Call DaXpY_(NTAU,DTG2(IV,IX,IY,IZ),WORK(LTAU),1,CLAG1,1)
           !! For right derivative: <Psi1|Evx Eyz|I>
           IF (IP2.GE.IP1STA.AND.IP2.LE.IP1END) THEN
              ibuf = lsgm1+mxci*(ip2-ip1sta)
              Call DaXpY_(MXCI,DTG2(IV,IX,IY,IZ),WORK(IBUF),1,
     *                    WORK(LFROMD),1)
           ELSE
         CALL SIGMA1_CP2(JL,IL,DTG2(IV,IX,IY,IZ),ISSG2,CI1,WORK(LFROMD),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
           END IF
           DTG2(IV,IX,IY,IZ) = 0.0D+00
          END IF
          IF (DOG3) THEN
          CALL DCopy_(MXCI,[0.0D+00],0,WORK(LBUF1),1)
          DO IP1=MAX(IP2,IP1STA),IP1END
           IT=L2ACT(IWORK(LP2LEV1-1+IP1))
           IU=L2ACT(IWORK(LP2LEV2-1+IP1))
           ITS=IASYM(IT)
           IUS=IASYM(IU)
           ISSG1=MUL(MUL(ITS,IUS),LSYM1)
           IF(ISSG1.EQ.ISTAU) THEN
C           L=LSGM1+MXCI*(IP1-IP1STA)
C           VAL=DDOT_(NTAU,WORK(LTAU),1,WORK(L),1)
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
            IF (DTG3(JTUVXYZ).NE.0.0D+00) THEN
              !! For left derivative: <I|Evx Eyz|Psi2> * Dtuvxyz
              !! I don't understand, but this is much faster than
              !! processing all possible vectors at once with DGEMM
              !! (and DGER) after finishing the IP1 loop below
              Call DaXpY_(MXCI,DTG3(JTUVXYZ),
     *                    WORK(LTAU),1,WORK(LDTU+MXCI*(IP1-IP1STA)),1)
              !! For right derivative: <Psi1|Etu|I> * Dtuvxyz
              !! This is also (slightly) faster than DGEMV, apparently
              Call DaXpY_(MXCI,DTG3(JTUVXYZ),
     *                    WORK(LSGM1+MXCI*(IP1-IP1STA)),1,
     *                    WORK(LBUF1),1)
              !! Prepare for the right derivative
C             WORK(LBUF2+IP1-MAX(IP2,IP1STA)) = DTG3(JTUVXYZ)
             END IF
C End of symmetry requirement IF-clause:
           END IF
C End of IP1 loop.
          END DO
C         !! For left derivative: <I|Evx Eyz|Psi2> * Dtuvxyz
C         CALL DGEMM_('N','T',MXCI,IP1END-MAX(IP2,IP1STA)+1,1,
C    &                1.0D+00,WORK(LTAU),MXCI,
C    &                        WORK(LBUF2),IP1END-MAX(IP2,IP1STA)+1,
C    &                1.0D+00,WORK(LDTU+MXCI*(MAX(IP2,IP1STA)-1)),MXCI)
C         !! For right derivative: <Psi1|Etu|I> * Dtuvxyz
C         CALL DGEMV_('N',MXCI,IP1END-MAX(IP2,IP1STA)+1,
C    &                1.0D+00,WORK(LSGM1+MXCI*(MAX(IP2,IP1STA)-1)),MXCI,
C    &                        WORK(LBUF2),1,
C    &                0.0D+00,WORK(LBUF1),1)
          !! Second operator for the right derivative:
          !! <Psi1|Etu Evx|I> * Dtuvxyz
          CALL SIGMA1_CP2(JL,IL,1.0D+00,ISTAU,WORK(LBUF1),WORK(LFROMD),
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
          END IF !! End of DOG3 clause
C End of IP2 loop.
         END DO
         LFROM=LFROM+MXCI
         LFROMD=LFROMD+MXCI
C End of IP3 loop.
        END DO
C
        LTO=LDTU
        !! <I|Etu Evx Eyz|Psi2> * Dtuvxyz
        DO IP1=IP1STA,IP1END
C Translate to levels:
         IL=IWORK(LP2LEV1-1+IP1)
         JL=IWORK(LP2LEV2-1+IP1)
         IT=L2ACT(IL)
         IU=L2ACT(JL)
         ITS=IASYM(IT)
         IUS=IASYM(IU)
         ISSG1=MUL(MUL(ITS,IUS),LSYM1)
         CALL SIGMA1_CP2(IL,JL,1.0D00,LSYM1,WORK(LTO),CLAG1,
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
         LTO=LTO+MXCI
        END DO
C End of IP1STA sectioning loop
       END DO
C
       LTO=LDYZ
       DO IP3=IP3STA,IP3END
        IY=L2ACT(IWORK(LP2LEV1-1+IP3))
        IZ=L2ACT(IWORK(LP2LEV2-1+IP3))
C LFROM will be start element of Sigma2=E(YZ) Psi2
        IYZ=IY+NASHT*(IZ-1)
        IYS=IASYM(IY)
        IZS=IASYM(IZ)
        ISSG2=MUL(MUL(IYS,IZS),LSYM2)
        IM=IWORK(LP2LEV1-1+IP3)
        JM=IWORK(LP2LEV2-1+IP3)
C LTO is first element of Sigma2 = E(YZ) Psi2
        CALL SIGMA1_CP2(JM,IM,1.0D00,LSYM2,WORK(LTO),CLAG2,
     &    IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &    IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &    WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
        LTO=LTO+MXCI
       END DO
C End of IP3STA sectioning loop
      END DO
C
      CALL GETMEM('TG3WRK','FREE','REAL',LTG3WRK,NTG3WRK)
      CALL GETMEM('BUF1','FREE','REAL',LBUF1,MXCI)
      CALL GETMEM('DTU','FREE','REAL',LDTU,MXCI*NTUBUF)
      CALL GETMEM('DYZ','FREE','REAL',LDYZ,MXCI*NYZBUF)
      CALL GETMEM('IPTOLEV','FREE','INTE',LP2LEV,2*NASHT**2)

 999  CONTINUE
      RETURN
      END
