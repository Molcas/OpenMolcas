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
* Copyright (C) 2014, Steven Vancoillie                                *
************************************************************************
      SUBROUTINE HCOUP(IVEC,JVEC,OVL,TG1,TG2,TG3,HEL)
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
#include "output.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      Dimension TG1(NASHT,NASHT)
      Dimension TG2(NASHT,NASHT,NASHT,NASHT)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)
      Dimension TG3(*)

      DIMENSION HECOMP(14,9)

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

      HEL=0.0D0
      HECOMP=0.0D0
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          HEBLK=0.0D0

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
            CALL HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                     DBL_MB(MV1),DBL_MB(MV2),OVL,HEBLK,
     &                     TG1,TG2,TG3)
          ELSE
            CALL HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                     WORK(MV1),WORK(MV2),OVL,HEBLK,
     &                     TG1,TG2,TG3)
          END IF
#else
          CALL HCOUP_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                   WORK(MV1),WORK(MV2),OVL,HEBLK,
     &                   TG1,TG2,TG3)
#endif
          CALL RHS_RELEASE (lg_V1,IASTA1,IAEND1,IISTA1,IIEND1)
          CALL RHS_RELEASE (lg_V2,IASTA2,IAEND2,IISTA2,IIEND2)
          CALL RHS_FREE (NAS,NIS,lg_V1)
          CALL RHS_FREE (NAS,NIS,lg_V2)

 1        CONTINUE
          HECOMP(ICASE,ISYM)=HEBLK
          HEL=HEL+HEBLK
        END DO
      END DO

C Sum-reduce the per-process contributions
      CALL GADGOP_SCAL(HEL,'+')
      NHECOMP=14*9
      CALL GADGOP(HECOMP,NHECOMP,'+')

      IF(IPRGLB.GE.DEBUG) THEN
        DO ICASE=1,13
          SUMSYM=0.0D0
          DO ISYM=1,NSYM
            SUMSYM=SUMSYM+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(ICASE,NSYM+1)=SUMSYM
        END DO

        DO ISYM=1,NSYM+1
          SUMCASE=0.0D0
          DO ICASE=1,13
            SUMCASE=SUMCASE+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(14,ISYM)=SUMCASE
        END DO

        WRITE(6,'(20a4)')('----',i=1,20)
        WRITE(6,*)'HCOUP: The contributions to the Hamiltonian coupling'
        WRITE(6,*)' elements, by case and by symmetry label.'
        DO IC=1,13
          WRITE(6,'(1X,A8,9F12.8)')
     &      CASES(IC),(HECOMP(IC,IS),IS=1,NSYM+1)
        END DO
        CALL XFLUSH(6)
        WRITE(6,'(1X,A8,9F12.8)')
     &    'Summed: ', (HECOMP(14,IS),IS=1,NSYM+1)
        WRITE(6,*)
      END IF


      END

      SUBROUTINE HCOUP_BLK(ICASE,ISYM,NAS,IISTA,IIEND,V1,V2,OVL,HEBLK,
     &                     TG1,TG2,TG3)
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
#include "output.fh"
#include "SysDef.fh"
#include "eqsolv.fh"

      DIMENSION V1(*), V2(*)

      Dimension TG1(NASHT,NASHT)
      Dimension TG2(NASHT,NASHT,NASHT,NASHT)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)
      Dimension TG3(*)


      HEBLK=0.0D0

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
            TMP=TG3(ITG3)
            IF(IYABS.EQ.IUABS) THEN
              TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
            END IF
            IF(IYABS.EQ.ITABS) THEN
              TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
              IF(IXABS.EQ.IUABS) THEN
                TMP=TMP+TG1(IVABS,IZABS)
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              TMP=TMP+TG2(IVABS,ITABS,IYABS,IZABS)
            END IF
C SA is the negative of this, and then some correction:
            SA=-TMP
            IF(IXABS.EQ.ITABS) THEN
              SA=SA+2.0D0*TG2(IVABS,IUABS,IYABS,IZABS)
              IF(IYABS.EQ.IUABS) THEN
                SA=SA+2.0D0*TG1(IVABS,IZABS)
              END IF
            END IF
C SA has been computed.

            HEBLK=HEBLK+SA*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            TMP=TG3(ITG3)
            IF(IYABS.EQ.IUABS) THEN
              TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
            END IF
            IF(IYABS.EQ.ITABS) THEN
              TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
              IF(IXABS.EQ.IUABS) THEN
                TMP=TMP+TG1(IVABS,IZABS)
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              TMP=TMP+TG2(IVABS,ITABS,IYABS,IZABS)
            END IF
            SC= TMP

            HEBLK=HEBLK+SC*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            SBtuxy=2.0d0*TG2(IXABS,ITABS,IYABS,IUABS)
            SBtuyx=2.0d0*TG2(IYABS,ITABS,IXABS,IUABS)
            IF(IXABS.EQ.ITABS) THEN
              SBtuxy=SBtuxy-4.0d0*TG1(IYABS,IUABS)
              SBtuyx=SBtuyx+2.0d0*TG1(IYABS,IUABS)
              IF(IYABS.EQ.IUABS) THEN
                SBtuxy=SBtuxy+8.0d0*OVL
                SBtuyx=SBtuyx-4.0d0*OVL
              END IF
            END IF
            IF(IYABS.EQ.IUABS) THEN
              SBtuxy=SBtuxy-4.0d0*TG1(IXABS,ITABS)
              SBtuyx=SBtuyx+2.0d0*TG1(IXABS,ITABS)
            END IF
            IF(IYABS.EQ.ITABS) THEN
              SBtuxy=SBtuxy+2.0d0*TG1(IXABS,IUABS)
              SBtuyx=SBtuyx-4.0d0*TG1(IXABS,IUABS)
              IF(IXABS.EQ.IUABS) THEN
                SBtuxy=SBtuxy-4.0d0*OVL
                SBtuyx=SBtuyx+8.0d0*OVL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              SBtuxy=SBtuxy+2.0d0*TG1(IYABS,ITABS)
              SBtuyx=SBtuyx-4.0d0*TG1(IYABS,ITABS)
            END IF

            SBP=SBtuxy + SBtuyx

            HEBLK=HEBLK+SBP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            SBtuxy=2.0d0*TG2(IXABS,ITABS,IYABS,IUABS)
            SBtuyx=2.0d0*TG2(IYABS,ITABS,IXABS,IUABS)
            IF(IXABS.EQ.ITABS) THEN
              SBtuxy=SBtuxy-4.0d0*TG1(IYABS,IUABS)
              SBtuyx=SBtuyx+2.0d0*TG1(IYABS,IUABS)
              IF(IYABS.EQ.IUABS) THEN
                SBtuxy=SBtuxy+8.0d0*OVL
                SBtuyx=SBtuyx-4.0d0*OVL
              END IF
            END IF
            IF(IYABS.EQ.IUABS) THEN
              SBtuxy=SBtuxy-4.0d0*TG1(IXABS,ITABS)
              SBtuyx=SBtuyx+2.0d0*TG1(IXABS,ITABS)
            END IF
            IF(IYABS.EQ.ITABS) THEN
              SBtuxy=SBtuxy+2.0d0*TG1(IXABS,IUABS)
              SBtuyx=SBtuyx-4.0d0*TG1(IXABS,IUABS)
              IF(IXABS.EQ.IUABS) THEN
                SBtuxy=SBtuxy-4.0d0*OVL
                SBtuyx=SBtuyx+8.0d0*OVL
              END IF
            END IF
            IF(IXABS.EQ.IUABS) THEN
              SBtuxy=SBtuxy+2.0d0*TG1(IYABS,ITABS)
              SBtuyx=SBtuyx-4.0d0*TG1(IYABS,ITABS)
            END IF

            SBM=SBtuxy - SBtuyx

            HEBLK=HEBLK+SBM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            GUTXY= TG2(IUABS,ITABS,IXABS,IYABS)
            SD11=2.0D0*GUTXY
            SD12= -GUTXY
            SD21= -GUTXY
            SD22= -TG2(IXABS,ITABS,IUABS,IYABS)
            IF(ITABS.EQ.IXABS) THEN
              GUY=TG1(IUABS,IYABS)
              SD11=SD11+2.0D0*GUY
              SD12=SD12 -GUY
              SD21=SD21 -GUY
              SD22=SD22+2.0D0*GUY
            END IF

            HEBLK=HEBLK+SD11*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS1),NAS)
            HEBLK=HEBLK+SD12*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS1),NAS)
            HEBLK=HEBLK+SD21*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS2),NAS)
            HEBLK=HEBLK+SD22*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS2),NAS)
          END DO
        END DO
************************************************************************
      CASE(6)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SE(t,x)=2*dxt - Dxt
            SE=-TG1(IXABS,ITABS)
            IF(IXABS.EQ.ITABS) SE=SE+2.0d0*OVL
            HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(7)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SE(t,x)=2*dxt - Dxt
            SE=-TG1(IXABS,ITABS)
            IF(IXABS.EQ.ITABS) SE=SE+2.0d0*OVL
            HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            SFtuxy=2.0d0*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=2.0d0*TG2(ITABS,IYABS,IUABS,IXABS)

            SFP=SFtuxy + SFtuyx
            HEBLK=HEBLK+SFP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            SFtuxy=2.0d0*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=2.0d0*TG2(ITABS,IYABS,IUABS,IXABS)

            SFM=SFtuxy - SFtuyx
            HEBLK=HEBLK+SFM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
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
            SG= TG1(ITABS,IXABS)

            HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(11)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
C Formula used: SG(t,x)= Gtx
            SG= TG1(ITABS,IXABS)

            HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(12)
        IF(ABS(OVL).GE.1.0D-12) THEN
          HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      CASE(13)
        IF(ABS(OVL).GE.1.0D-12) THEN
          HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      END SELECT
      Return
      END
