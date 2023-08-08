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
!       SUBROUTINE MSGradPrep(UEFF)
! C
!       USE CHOVEC_IO
!       use caspt2_gradient, only: iRoot1, iRoot2
! C
!       IMPLICIT REAL*8 (A-H,O-Z)
! C
! #include "rasdim.fh"
! #include "caspt2.fh"
! #include "eqsolv.fh"
! #include "WrkSpc.fh"
! #include "sigma.fh"

! #include "caspt2_grad.fh"
! #include "csfbas.fh"
! C
! #include "pt2_guga.fh"
! C
! #include "chocaspt2.fh"
! C
!       DIMENSION UEFF(nState,nState)
!       DIMENSION VECROT(nState,nState)
! C
!       iRoot1 = iRlxRoot
!       iRoot2 = iRlxRoot
!       DO ilStat = 1, nState
!         DO jlStat = 1, nState
!           TMP = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
!      *        + UEFF(ilStat,iRoot2)*UEFF(jlStat,iRoot1)
!           VECROT(ilStat,jlStat) = TMP*0.5D+00
!         END DO
!       END DO
! C
!       END SUBROUTINE MSGradPrep
C
C-----------------------------------------------------------------------
C
      Subroutine MS_Res(MODE,IST,JST,Scal)
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
      REAL*8 DVALUE

      INTEGER I
      INTEGER NTG1,NTG2,NTG3,LTG1,LTG2,LTG3
      INTEGER IDCI,LCI1,LCI2
      REAL*8 OVL,DUMMY(1)

C We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
C Note: Need proper allocation even if unused, sinced allocated
C arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      CALL GETMEM('TG1','ALLO','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','ALLO','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','ALLO','REAL',LTG3,NTG3)
      WORK(LTG1)=0.0D0
      WORK(LTG2)=0.0D0
      WORK(LTG3)=0.0D0

      CALL GETMEM('MCCI1','ALLO','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','ALLO','REAL',LCI2,MXCI)
      IF(ISCF.EQ.0) THEN
C Read root vectors nr. IST and JST from LUCI.
        IDCI=IDTCEX
        DO I=1,NSTATE
          IF(I.EQ.IST) THEN
            CALL DDAFILE(LUCIEX,2,WORK(LCI1),NCONF,IDCI)
            IF(I.EQ.JST) THEN
              CALL DCOPY_(NCONF,WORK(LCI1),1,WORK(LCI2),1)
            END IF
          ELSE IF(I.EQ.JST) THEN
            CALL DDAFILE(LUCIEX,2,WORK(LCI2),NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          END IF
        END DO
      END IF

      CALL MKTG3(STSYM,STSYM,WORK(LCI1),WORK(LCI2),OVL,
     &           WORK(LTG1),WORK(LTG2),NTG3,WORK(LTG3))
      CALL GETMEM('MCCI1','FREE','REAL',LCI1,MXCI)
      CALL GETMEM('MCCI2','FREE','REAL',LCI2,MXCI)

      !! Do similar to RHS_STRANS. Multiply the solution vector (T) with
      !! the overlap-like term constructe with transition density
      !! matrices. The output is IVECC (MODE=1) or IVECC2 (MODE=2)
      IF (MODE.EQ.1) THEN
        CALL MS_STRANS(IVECW,IVECC,OVL,WORK(LTG1),WORK(LTG2),
     &                 WORK(LTG3),DVALUE,SCAL)
      ELSE IF (MODE.EQ.2) THEN
        CALL MS_STRANS(IVECC,IVECC2,OVL,WORK(LTG1),WORK(LTG2),
     &                 WORK(LTG3),DVALUE,SCAL)
      ELSE IF (MODE.EQ.3) THEN
        CALL MS_STRANS(IVECW,IVECC,OVL,WORK(LTG1),WORK(LTG2),
     &                 WORK(LTG3),DVALUE,SCAL)
      END IF

      CALL GETMEM('TG1','FREE','REAL',LTG1,NTG1)
      CALL GETMEM('TG2','FREE','REAL',LTG2,NTG2)
      CALL GETMEM('TG3','FREE','REAL',LTG3,NTG3)

      RETURN
      End Subroutine MS_Res
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MS_STRANS(IVEC,JVEC,OVL,TG1,TG2,TG3,HEL,SCAL)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_output, only:iPrGlb,debug
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
C     if (icase.ne.12.and.icase.ne.13) cycle ! H
C     if (icase.ne.10.and.icase.ne.11) cycle ! G
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
            CALL MS_STRANS_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                      DBL_MB(MV1),DBL_MB(MV2),OVL,
     &                      TG1,TG2,TG3,SCAL)
          ELSE
            CALL MS_STRANS_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                      WORK(MV1),WORK(MV2),OVL,
     &                      TG1,TG2,TG3,SCAL)
          END IF
#else
          CALL MS_STRANS_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                    WORK(MV1),WORK(MV2),OVL,
     &                    TG1,TG2,TG3,SCAL)
#endif
          !! Save T*S
          CALL RHS_SAVE (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
C
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
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MS_STRANS_BLK(ICASE,ISYM,NAS,IISTA,IIEND,V1,V2,OVL,
     &                         TG1,TG2,TG3,SCAL)
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

      Dimension TG1(NASHT,NASHT)
      Dimension TG2(NASHT,NASHT,NASHT,NASHT)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)
      Dimension TG3(*)


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

C           HEBLK=HEBLK+SA*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SA*SCAL,V1(IAS),NAS,V2(JAS),NAS)
C           Do i = 1, NAS
C             V2(JAS+i-1) = SA*V1(iAS-1)
C           End Do
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

C           HEBLK=HEBLK+SC*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SC*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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

C           HEBLK=HEBLK+SBP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SBP*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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

C           HEBLK=HEBLK+SBM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SBM*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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

C           HEBLK=HEBLK+SD11*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS1),NAS)
C           HEBLK=HEBLK+SD12*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS1),NAS)
C           HEBLK=HEBLK+SD21*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS2),NAS)
C           HEBLK=HEBLK+SD22*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS2),NAS)
            Call DaXpY_(NISBLK,SD11*SCAL,V1(IAS1),NAS,V2(JAS1),NAS)
            Call DaXpY_(NISBLK,SD12*SCAL,V1(IAS1),NAS,V2(JAS2),NAS)
            Call DaXpY_(NISBLK,SD21*SCAL,V1(IAS2),NAS,V2(JAS1),NAS)
            Call DaXpY_(NISBLK,SD22*SCAL,V1(IAS2),NAS,V2(JAS2),NAS)
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
C           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SE*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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
C           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SE*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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
C           HEBLK=HEBLK+SFP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SFP*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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
C           HEBLK=HEBLK+SFM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SFM*SCAL,V1(IAS),NAS,V2(JAS),NAS)
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

C           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
C           Do i = 1, NISBLK
C             V2(JAS+NAS*(i-1)) = SG*V1(IAS+NAS*(i-1))
C           End Do
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

C           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
C           Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
            Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
C           Do i = 1, NISBLK
C             V2(JAS+NAS*(i-1)) = SG*V1(IAS+NAS*(i-1))
C           End Do
          END DO
        END DO
************************************************************************
      CASE(12)
        IF(ABS(OVL).GE.1.0D-12) THEN
C         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      CASE(13)
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
      Subroutine LoadCI_XMS(Bas,Mode,CI,Istate,U0)
      implicit real(8) (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      character(len=1) Bas
      real(8) CI(Nconf),U0(nState,nState)
      integer ID, Istate
C
C     MODE=0 is equivalent to LoadCI (XMS basis)
C     MODE=1 constructs the CI vector in CASSCF basis (back-transformed)
C     CSF in natural (Bas=N) or quasi-canonical (Bas=C) orbital basis
C
      If (Bas.eq.'N'.or.Bas.eq.'n') Then
        ID=IDCIEX !! natural
      Else If (Bas.eq.'C'.or.Bas.eq.'c') Then
        ID=IDTCEX !! quasi-canonical
      ELse
        write (6,*) "the first argument in LoadCI_XMS should be either",
     *              "N (natural) or C (quasi-canonical)"
        call abend
      End If
C
      If (Mode.eq.0) Then
        do I=1,Istate-1
          call ddafile(LUCIEX,0,CI,Nconf,ID)
        end do
        call ddafile(LUCIEX,2,CI,Nconf,ID)
      Else If (Mode.eq.1) Then
        Call DCopy_(Nconf,[0.0D+00],0,CI,1)
        Call GetMem('WRK','ALLO','REAL',ipWRK,Nconf)
        do I=1,Nstate
          call ddafile(LUCIEX,2,Work(ipWRK),Nconf,ID)
          Call DaXpY_(Nconf,U0(Istate,I),Work(ipWRK),1,CI,1)
        end do
        Call GetMem('WRK','FREE','REAL',ipWRK,Nconf)
      End If

      return
      End Subroutine LoadCI_XMS
C
C-----------------------------------------------------------------------
C
      Subroutine XMS_Grad(CLag,H0,U0,UEFF,OMGDER)
C
      use caspt2_gradient, only: do_nac, do_csf, iRoot1, iRoot2
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Dimension CLag(nConf,nState),H0(nState,nState),U0(nState,nState),
     *          UEFF(nState,nState)
      Dimension SLag(nState*nState)
      Dimension OMGDER(nState,nState)
C
C     The XMS rotation applies to any variants: XMS-CASPT2, XDW-CASPT2,
C     and RMS-CASPT2.
C
      If (IFXMS .or. IFRMS) Then
C
        Call GetMem('CI1   ','ALLO','REAL',ipCI1 ,nConf)
        Call GetMem('CI2   ','ALLO','REAL',ipCI2 ,nConf)
        Call GetMem('SGM1  ','ALLO','REAL',ipSGM1,nConf)
        Call GetMem('SGM2  ','ALLO','REAL',ipSGM2,nConf)
        Call GetMem('TG1   ','ALLO','REAL',ipTG1 ,nAshT**2)
        Call GetMem('TG2   ','ALLO','REAL',ipTG2 ,nAshT**4)
C
        Call GetMem('DG1   ','ALLO','REAL',ipDG1 ,nAshT**2)
        Call GetMem('DG2   ','ALLO','REAL',ipDG2 ,nAshT**4)
C
!       ----- Construct pseudo-density matrix -----

        !! First, we need to consider the derivative of the rotation
        !! of the CASSCF energy (Heff[1]).

        !! Forward transformation of UEFF (CASSCF basis to XMS basis)
        !! SLag is used as a working array
        Call DGEMM_('T','N',nState,nState,nState,
     &              1.0D+00,U0,nState,UEFF,nState,
     &              0.0D+00,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C       write (6,*) "ueff in xms basis"
C       call sqprt(ueff,nstate)
C
        !! Then actual calculation
C       write (*,*) "ipfimo"
C       call sqprt(work(ipfimo),nbast)
        Call DCopy_(nCLag ,[0.0D+00],0,Work(ipCLag),1)
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,Work(ipCI1),iStat,U0)
          Do jStat = 1, nState
            Call LoadCI_XMS('C',0,Work(ipCI2),jStat,U0)
C
            Call Dens2T_RPT2(Work(ipCI1),Work(ipCI2),
     *                       Work(ipSGM1),Work(ipSGM2),
     *                       Work(ipTG1),Work(ipTG2))
            Call DScal_(nAshT**2,0.5D+00,Work(ipTG1),1)
            Call DScal_(nAshT**4,0.5D+00,Work(ipTG2),1)
C
            Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipDG1),1)
            Call DCopy_(nAshT**4,[0.0D+00],0,Work(ipDG2),1)
C         Call Cnst_SA_CLag(iStat.eq.jStat,Work(ipTG1),Work(ipTG2),
C    *                      Work(ipDG1),Work(ipDG2),EEE)
C     write (*,*) istat,jstat,eee

C         call sqprt(work(ipdg1),nasht)
C         call sqprt(work(ipdg2),nasht**2)
            Call CnstInt(1,Work(ipDG1),Work(ipDG2))
C         call sqprt(work(ipdg1),nasht)
C         call sqprt(work(ipdg2),nasht**2)
C         call abend
C
C         Scal = UEFF(iStat,iRlxRoot)*UEFF(jStat,iRlxRoot)
            If (do_nac) Then
              Scal =(UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
     *             + UEFF(jStat,iRoot1)*UEFF(iStat,iRoot2))*0.5d+00
            Else
              Scal = UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
            End If
C       write (*,*) " scal in xms"
C       write (*,*) istat,jstat,scal
            if (IFDW .and. zeta >= 0.0d0) then
              scal = scal + OMGDER(iStat,jStat)
            end if
            Call DScal_(nAshT**2,Scal,Work(ipDG1),1)
            Call DScal_(nAshT**4,Scal,Work(ipDG2),1)
C         call sqprt(work(ipdg1),nasht)
C         call sqprt(work(ipdg2),nasht**2)
C
            Call GetMem('DG3   ','ALLO','REAL',ipDG3 ,nAshT**6)
C         Call GetMem('DF1   ','ALLO','REAL',ipDF1 ,nAshT**2)
C         Call GetMem('DF2   ','ALLO','REAL',ipDF2 ,nAshT**4)
C         Call GetMem('DF3   ','ALLO','REAL',ipDF3 ,nAshT**6)
            Call DCopy_(nAshT**6,[0.0D+00],0,Work(ipDG3),1)
C         Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipDF1),1)
C         Call DCopy_(nAshT**4,[0.0D+00],0,Work(ipDF2),1)
C         Call DCopy_(nAshT**6,[0.0D+00],0,Work(ipDF3),1)
C         IFF=1
C         JSTATE = iStat
C         If (iStat.eq.jStat) Then
C         Call CLagSym(nAshT,Work(ipDG1),Work(ipDG2),
C    *                 Work(ipDF1),Work(ipDF2),0)
C         Call CnstCLag(IFF,Work(ipCLag+nConf*(iStat-1)),
C    *                  Work(ipDG1),Work(ipDG2),Work(ipDG3),
C    *                  Work(ipDF1),Work(ipDF2),Work(ipDF3),
C    *                  Work(ipTG1),
C    *                  Work(ipDF1),Work(ipDF2),Work(ipDF3),
C    *     U0)
C         Else
            STSYM=1
            NTG1=NASHT**2
            NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
            OVL=0.0D+00
            CALL DERTG3(.False.,STSYM,STSYM,WORK(ipCI1),WORK(ipCI2),OVL,
     &                  WORK(ipDG1),WORK(ipDG2),NTG3,WORK(ipDG3),
     &                  Work(ipCLag+nConf*(iStat-1)),
     &                  Work(ipCLag+nConf*(jStat-1)))
C         End If
            Call GetMem('DG3   ','FREE','REAL',ipDG3 ,nAshT**6)
C         Call GetMem('DF1   ','FREE','REAL',ipDF1 ,nAshT**2)
C         Call GetMem('DF2   ','FREE','REAL',ipDF2 ,nAshT**4)
C         Call GetMem('DF3   ','FREE','REAL',ipDF3 ,nAshT**6)
          End Do
        End Do
C
        !! Back transformation of UEFF (XMS basis to CASSCF basis)
        !! SLag is used as a working array
        Call DGEMM_('N','N',nState,nState,nState,
     &              1.0D+00,U0,nState,UEFF,nState,
     &              0.0D+00,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C     write (6,*) "ueff in casscf basis"
C     call sqprt(ueff,nstate)
C
        Call GetMem('DG1   ','FREE','REAL',ipDG1 ,nAshT**2)
        Call GetMem('DG2   ','FREE','REAL',ipDG2 ,nAshT**4)
C
      !! Add to the full CI Lagrangian
      !! Work(ipCLag) is in quasi-canonical basis, so transformation
      !! to natural basis is required.
        IF(ORBIN.EQ.'TRANSFOR') Then
          Do iState = 1, nState
            Call CLagX_TrfCI(Work(ipCLag+nConf*(iState-1)))
          End Do
        End If
        Call DaXpY_(nCLag,1.0D+00,Work(ipCLag),1,CLag,1)
C
        !! Compute the Lagrange multiplier for XMS
        !! The diagonal element is always zero.
        !! The code has an additional scaling with 0.5,
        !! because some contributions are doubled.
C     write (*,*) "istat,jstat,scal"
        Call DCopy_(nState*nState,[0.0d+00],0,SLag,1)
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,Work(ipCI1),iStat,U0)
          EEI = H0(iStat,iStat)
          Do jStat = 1, iStat
            If (iStat.eq.jStat) Then
              SLag(iStat+nState*(jStat-1)) = 0.0D+00
            Else
              Call LoadCI_XMS('N',0,Work(ipCI2),jStat,U0)
              EEJ = H0(jStat,jStat)
C
              Scal = DDOT_(nConf,Work(ipCI1),1,CLag(1,jStat),1)
     *             - DDOT_(nConf,Work(ipCI2),1,CLag(1,iStat),1)
C      write (*,*) "original scal = ", scal
              If (do_csf) Then
                !! JCTC 2017, 13, 2561: eq.(66)
                !! iStat and jStat: XMS
                !! kStat and lStat: CASSCF
                fact = 0.0d+00
                Do kStat = 1, nState
                  Do lStat = 1, nState
                    fact = fact
     *                   + (UEFF(kStat,iRoot1)*UEFF(lStat,iRoot2)
     *                   - UEFF(kStat,iRoot2)*UEFF(lStat,iRoot1))*0.5d0
     *                   * U0(kStat,iStat)*U0(lStat,jStat)
                  End Do
                End Do
                Scal = Scal+fact*(ENERGY(iRoot2)-ENERGY(iRoot1))*2.0d0
C      write (*,*) "scal after= ", scal
C      write (*,*) fact,energy(iroot2)-energy(iroot1)
              End If
              Scal = 0.25D+00*Scal/(EEJ-EEI)
C           write (*,*) istat,jstat,scal
              SLag(iStat+nState*(jStat-1)) = Scal
              SLag(jStat+nState*(iStat-1)) = Scal
            End If
          End Do
        End Do
C
      !! Back transformation of UEFF (XMS basis to CASSCF basis)
      !! SLag is used as a working array
C     Call DGEMM_('N','N',nState,nState,nState,
C    &            1.0D+00,U0,nState,UEFF,nState,
C    &            0.0D+00,Work(ipTG2),nState)
C     Call DCopy_(nState**2,Work(ipTG2),1,UEFF,1)
C
      !! Either subtract the CI derivative of Heff(1)
      !! or add off-diagonal couplings in rhs_sa (Z-vector)
      !Call DaXpY_(nCLag,-1.0D+00,Work(ipCLag),1,CLag,1)
C
        !! Finally, construct the pseudo-density matrix
        !! d = \sum_{ST} w_{ST} * d_{ST}
        Call GetMem('G1    ','ALLO','REAL',ipG1  ,nAshT**2)
        Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipG1),1)
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,Work(ipCI1),iStat,U0)
          Do jStat = 1, nState
            If (ABS(SLag(iStat+nState*(jStat-1))).le.1.0d-10) Cycle
            Call LoadCI_XMS('C',0,Work(ipCI2),jStat,U0)
C
C         Call Dens2T_RPT2(Work(ipCI1),Work(ipCI2),
C    *                     Work(ipSGM1),Work(ipSGM2),
C    *                     Work(ipTG1),Work(ipTG2))
            Call Dens1T_RPT2(Work(ipCI1),Work(ipCI2),
     *                       Work(ipSGM1),Work(ipTG1))
C         call sqprt(work(iptg1),nasht)
            Scal = SLag(iStat+nState*(jStat-1))*2.0d+00
C         write (*,*) "istat,jstat=",istat,jstat
C         write (*,*) "scal = ", scal
            Call DaXpY_(nAshT**2,Scal,Work(ipTG1),1,Work(ipG1),1)
          End Do
        End Do
C     write (*,*) "ipG1"
C     call sqprt(work(ipG1),nasht)
C
        Call GetMem('CI1   ','FREE','REAL',ipCI1 ,nConf)
        Call GetMem('CI2   ','FREE','REAL',ipCI2 ,nConf)
        Call GetMem('SGM1  ','FREE','REAL',ipSGM1,nConf)
        Call GetMem('SGM2  ','FREE','REAL',ipSGM2,nConf)
        Call GetMem('TG1   ','FREE','REAL',ipTG1 ,nAshT**2)
        Call GetMem('TG2   ','FREE','REAL',ipTG2 ,nAshT**4)
C
C       ----- Calculate orbital derivatives -----
C
        Call GetMem('RDMEIG','ALLO','REAL',ipRDMEIG,nAshT**2)
        Call GetMem('DPT   ','ALLO','REAL',ipDPT  ,nBasSq)
C     Call GetMem('FIFA  ','ALLO','REAL',ipFIFA ,nBasSq)
C
        Call DCopy_(nBasSq,[0.0D+00],0,Work(ipDPT),1)

        Call GetMem('TRFMAT','ALLO','REAL',ipTrf   ,nBasSq)
        Call GetMem('RDMSA ','ALLO','REAL',ipRDMSA ,nAshT**2)
        Call GetMem('WRK1  ','ALLO','REAL',ipWRK1  ,nBasSq)
        Call GetMem('WRK2  ','ALLO','REAL',ipWRK2  ,nBasSq)
C
        !! Construct always state-averaged density; XMS basis is always
        !! generated with the state-averaged density.
        Call DCopy_(nDRef,[0.0D+00],0,Work(ipWRK1),1)
        Call GetMem('LCI','ALLO','REAL',LCI,nConf)
        Wgt  = 1.0D+00/nState
        Do iState = 1, nState
C       Call DaXpY_(nDRef,Wgt,Work(LDMix+nDRef*(iState-1)),1,
C    *              Work(ipWRK1),1)
          Call LoadCI(WORK(LCI),iState)
          call POLY1(WORK(LCI))
          call GETDREF(WORK(ipWRK2))
          Call DaXpY_(nDRef,Wgt,Work(ipWRK2),1,Work(ipWRK1),1)
        End Do
        Call GetMem('LCI','FREE','REAL',LCI,nConf)
        Call SQUARE(Work(ipWRK1),Work(ipRDMSA),1,nAshT,nAshT)
C
        nOrbI = nBas(1) - nDel(1) !! nOrb(1)
        nBasI = nBas(1)
C       Call SQUARE(Work(LFIFA),Work(ipFIFA),1,nOrbI,nOrbI)
        !! ipFIFASA is in natural orbital basis
        Call DCopy_(nBsqT,[0.0D+0],0,Work(ipTrf),1)
        Call CnstTrf(Work(LTOrb),Work(ipTrf))
C
        !! FIFA: natural -> quasi-canonical
C       write (*,*) "fifasa"
C       call sqprt(work(ipfifasa),norbi)
        If (IFDW .or. IFRMS) Then
          Call DGemm_('T','N',nOrbI,nOrbI,nOrbI,
     *                1.0D+00,Work(ipTrf),nBasI,Work(ipFIFASA),nOrbI,
     *                0.0D+00,Work(ipWRK1),nOrbI)
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,
     *                1.0D+00,Work(ipWRK1),nOrbI,Work(ipTrf),nBasI,
     *                0.0D+00,Work(ipFIFASA),nOrbI)
        End If
        Call DCopy_(nBasSq,Work(ipFIFASA),1,Work(ipFIFA),1)
C     Call SQUARE(Work(LFIFA),Work(ipFIFA),1,nOrbI,nOrbI)
C     write (*,*) "fifa in quasi-canonical correct"
C     call sqprt(work(ipfifa),norbi)
C
        !! Orbital derivatives of FIFA
        !! Both explicit and implicit orbital derivatives are computed
        !! Also, compute G(D) and put the active contribution to RDMEIG
        Call DCopy_(nOLag ,[0.0D+00],0,Work(ipOLag),1)
        Call EigDer2(Work(ipRDMEIG),Work(ipTrf),Work(ipFIFA),
     *               Work(ipRDMSA),Work(ipG1),
     *               Work(ipWRK1),Work(ipWRK2))
C
        !! Add to PT2 density
        !! No inactive contributions. Correct as long as CASSCF CI
        !! vector are orthogonal.
        Call AddDEPSA(Work(ipDPT),Work(ipG1))
        Call DPT2_TrfStore(1.0D+00,Work(ipDPT),Work(ipDPT2),
     *                   Work(ipTrf),Work(ipWRK1))
C
        !! Finalize OLag (anti-symetrize) and construct WLag
        Call OLagFinal(Work(ipOLag),Work(ipTrf))
C
        Call GetMem('RDMSA ','FREE','REAL',ipRDMSA ,nAshT**2)
        Call GetMem('WRK1  ','FREE','REAL',ipWRK1 ,nBasSq)
        Call GetMem('WRK2  ','FREE','REAL',ipWRK2 ,nBasSq)
C
        Call GetMem('DPT   ','FREE','REAL',ipDPT  ,nBasSq)
C
C       ----- Calculate CI derivatives -----
C
        !! use quasi-canonical CSF rather than natural CSF
C     ISAV = IDCIEX
C     IDCIEX = IDTCEX
        Call DCopy_(nCLag,[0.0D+00],0,Work(ipCLag),1)
C
        !! 1) Explicit CI derivative
        !! a: Extract FIFA in the AS for explicit CI derivative
        !!    Note that this FIFA uses state-averaged density matrix
C       call sqprt(work(ipfifa),nbast)
        Do iAsh = 1, nAshT
          Do jAsh = 1, nAshT
            Work(ipG1+iAsh-1+nAshT*(jAsh-1)) = Work(ipFIFA + nFro(1)
     &          + nIsh(1) + iAsh-1 + nBas(1)*(nFro(1)+nIsh(1)+jAsh-1))
          End Do
        End Do
C       call sqprt(work(ipg1),nasht)
        !! Transform quasi-canonical to natural
        nCor = nFro(1)+nIsh(1)
C     write (*,*) nfro(1),nish(1)
C     call sqprt(work(iptrf),nbast)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipTrf+nBasT*nCor+nCor),nBasT,
     *                      Work(ipG1),nAshT,
     *              0.0D+00,Work(ipFIFA),nAshT)
C     call sqprt(work(ipfifa),nasht)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipFIFA),nAshT,
     *                      Work(ipTrf+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipG1),nAshT)
C     call sqprt(work(ipg1),nasht)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipTrf+nBasT*nCor+nCor),nBasT,
     *                      Work(ipRDMEIG),nAshT,
     *              0.0D+00,Work(ipFIFA),nAshT)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,
     *              1.0D+00,Work(ipFIFA),nAshT,
     *                      Work(ipTrf+nBasT*nCor+nCor),nBasT,
     *              0.0D+00,Work(ipRDMEIG),nAshT)
        Call GetMem('TRFMAT','FREE','REAL',ipTrf   ,nBasSq)
C     Call GetMem('FIFA  ','FREE','REAL',ipFIFA ,nBasSq)
C
        !! b: FIFA in the inactive space
        TRC=0.0D0
C     DO ISYM=1,NSYM
        DO I=1,NISH(1) ! ISYM)
C         II=IOFF(ISYM)+(I*(I+1))/2
          II=(I*(I+1))/2
          TRC=TRC+WORK(LFIFA+II-1)
        END DO
C     END DO
* Contribution from inactive orbitals:
        EINACT=2.0D0*TRC
        !! This EINACT may be wrong. Perhaps, WORK(LFIFA) has to be
        !! back-transformed to natural orbital basis. However, this does
        !! not contribute to the final gradient as long as all the
        !! (internal) CI vectors are orthogonal.
C
        !! c: Finally, compute explicit CI derivative
        !! y_{I,T} = w_{ST} <I|f|S>
        !! Here, ipG1 is FIFA = ftu
C     write (*,*) einact
C     call sqprt(work(ipg1),nasht)
C     do i = 1, nstate*(nstate+1)/2
C     write (*,'(i3,f20.10)') i,slag(i)
C     end do
        Call CLagEigT(Work(ipCLag),Work(ipG1),SLag,EINACT)
C     do i = 1, nconf*nstate
C     write (*,'(i3,f20.10)') i,work(ipclag+i-1)
C     end do
C
        !! 2) Implicit CI derivative
        Call CLagEig(.False.,Work(ipCLag),Work(ipRDMEIG))
C     do i = 1, nconf*nstate
C     write (*,'(i3,f20.10)') i,work(ipclag+i-1)
C     end do
C
        Call GetMem('RDMEIG','FREE','REAL',ipRDMEIG,nAshT**2)
        Call GetMem('G1    ','FREE','REAL',ipG1    ,nAshT**2)
C
      End If
C
      If (do_csf) Then
        !! Eq (68)
        !! I think this contributes even in XMS-CASPT2.
        !! However, maybe I implemented other terms in a different way
        !! from BAGEL.
        Call GetMem('CI1   ','ALLO','REAL',ipCI1 ,nConf)
        !! Do in the natural orbital in the XMS basis
        !! because ipCLag is like that
        !! CASSCF -> XMS
        Call DGEMM_('T','N',nState,nState,nState,
     *              1.0D+00,U0,nState,UEFF,nState,
     *              0.0D+00,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C
        EDIFF = ENERGY(iRoot2)-ENERGY(iRoot1)
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,Work(ipCI1),iStat,U0)
          Do jStat = 1, nState
            Scal = (UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
     *           -  UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1))*0.5d+00
            Scal = Scal * EDIFF
            Call DaXpY_(nConf,Scal,Work(ipCI1),1,
     *                  Work(ipCLag+nConf*(jStat-1)),1)
          End Do
        End Do
C
        !! XMS -> CASSCF
        Call DGEMM_('N','N',nState,nState,nState,
     *              1.0D+00,U0,nState,UEFF,nState,
     *              0.0D+00,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C
        Call GetMem('CI1   ','FREE','REAL',ipCI1 ,nConf)
      End If
C
      !! Finally, add to the total CI derivative array
      Call DaXpY_(nConf*nState,1.0D+00,Work(ipCLag),1,CLag,1)
C     IDCIEX = ISAV
C
      End Subroutine XMS_Grad
C
C-----------------------------------------------------------------------
C
      !! From poly3
      SUBROUTINE CLagEigT(CLag,RDMEIG,SLag,EINACT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"
C
      DIMENSION CLag(nConf,nState),RDMEIG(*),SLag(*)
C
      !! RDMEIG
      Call GetMem('LCI1','ALLO','REAL',LCI1,nConf)
      Call GetMem('LCI2','ALLO','REAL',LCI2,nConf)
      Do iStat = 1, nState
        If (ISCF.EQ.0) Then
          Call LoadCI(Work(LCI1),iStat)
        Else
          Work(LCI1) = 1.0D+00
        End If
        !! Skip iStat = jStat because SLag is zero
        Do jStat = 1, nState !! iStat-1
          If (ISCF.EQ.0) Then
            Call LoadCI(Work(LCI2),jStat)
          Else
            Work(LCI2) = 1.0D+00
          End If
          !! One of doubling is due to the scaling factor
          !! One of doubling is due to the symmetry of iStat and jStat
C         Scal = SLag(iStat+nState*(jStat-1))*4.0d+00
          Scal = SLag(iStat+nState*(jStat-1))*2.0d+00
          If (ABS(Scal).le.1.0D-09) Cycle
C
          Call Poly1_CLagT(Work(LCI1),Work(LCI2),
     *                     CLag(1,iStat),CLag(1,jStat),RDMEIG,Scal)
          !! Inactive terms
          Call DaXpY_(nConf,Scal*EINACT,Work(LCI1),1,CLag(1,jStat),1)
          Call DaXpY_(nConf,Scal*EINACT,Work(LCI2),1,CLag(1,iStat),1)
        End Do
      End Do
C
      Call GetMem('LCI1','FREE','REAL',LCI1,nConf)
      Call GetMem('LCI2','FREE','REAL',LCI2,nConf)
C
      Return
C
      End Subroutine CLagEigT
C
C-----------------------------------------------------------------------
C
      SUBROUTINE POLY1_CLagT(CI1,CI2,CLag1,CLag2,RDMEIG,Scal)
      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES THE 1-EL DENSITY
* MATRIX FOR A CASSCF WAVE FUNCTION.
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"

      REAL*8, INTENT(IN) :: CI1(NCONF),CI2(NCONF)

      INTEGER LSGM1

      INTEGER I
      REAL*8 :: CLag1(*), CLag2(*), RDMEIG(*),Scal


      IF(NLEV.GT.0) THEN
        CALL GETMEM('LSGM1','ALLO','REAL',LSGM1 ,MXCI)
        CALL DENS1T_RPT2_CLag(CI1,CI2,WORK(LSGM1),
     *                        CLag1,CLag2,RDMEIG,Scal)
      END IF

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in common block in pt2_guga.fh
* CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV.GT.0) THEN
        CALL GETMEM('LSGM1','FREE','REAL',LSGM1 ,MXCI)
      END IF


      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DENS1T_RPT2_CLag (CI1,CI2,SGM1,CLag1,CLag2,RDMEIG,SCAL)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      LOGICAL RSV_TSK

      REAL*8 CI1(MXCI),CI2(MXCI),SGM1(MXCI)
      !! Symmetry?
      REAL*8 CLag1(nConf),CLag2(nConf),RDMEIG(NLEV,NLEV),Scal

C     REAL*8 GTU

      INTEGER ID
      INTEGER IST,ISU,ISTU
      INTEGER IT,IU,LT,LU

      INTEGER ITASK,LTASK,LTASK2T,LTASK2U,NTASKS

      INTEGER ISSG,NSGM

* Purpose: Compute the 1-electron density matrix array G1.


* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.fh

* SVC20100311: set up a task table with LT,LU
* SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks=(nLev**2+nLev)/2

      nTasks = nLev**2
      CALL GETMEM ('Tasks','ALLO','INTE',lTask,2*nTasks)
      lTask2T=lTask
      lTask2U=lTask+nTasks

      iTask=0
C     DO LT=1,nLev
C       DO LU=1,LT
C         iTask=iTask+1
C         iWork(lTask2T+iTask-1)=LT
C         iWork(lTask2U+iTask-1)=LU
C       ENDDO
C     ENDDO
      ! First, IL < JL pairs.
      Do LT = 1, nLev-1
        Do LU = LT+1, nLev
          iTask = iTask + 1
          iWork(lTask2T+iTask-1) = LT
          iWork(lTask2U+iTask-1) = LU
        End Do
      End Do
      ! Then, IL = JL pairs.
      Do LT = 1, nLev
        iTask = iTask + 1
        iWork(lTask2T+iTask-1) = LT
        iWork(lTask2U+iTask-1) = LT
      End Do
      ! Last, IL > JL pairs.
      Do LT = 2, nLev
        Do LU = 1, LT-1
          iTask = iTask + 1
          iWork(lTask2T+iTask-1) = LT
          iWork(lTask2U+iTask-1) = LU
        End Do
      End Do
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

* SVC20100311: BEGIN SEPARATE TASK EXECUTION
 500  If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
      LT=iWork(lTask2T+iTask-1)
        IST=ISM(LT)
        IT=L2ACT(LT)
        LU=iWork(lTask2U+iTask-1)
          ! LTU=iTask
          ISU=ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,STSYM)
          NSGM=NCSF(ISSG)
          IF(NSGM.EQ.0) GOTO 500
* GETSGM2 computes E_UT acting on CI and saves it on SGM1
          CALL GETSGM2(LU,LT,STSYM,CI1,SGM1)
          IF(ISTU.EQ.1) THEN
            ! Symmetry not yet
C            write(6,*) "it,iu = ", it,iu
            Call DaXpY_(NSGM,RDMEIG(IT,IU)*SCAL,SGM1,1,CLag2,1)
C           if (IT.ne.IU)
C    *        Call DaXpY_(NSGM,2.0d+00*RDMEIG(IT,IU),SGM1,1,CLag,1)
          END IF
          CALL GETSGM2(LU,LT,STSYM,CI2,SGM1)
          IF(ISTU.EQ.1) THEN
            Call DaXpY_(NSGM,RDMEIG(IT,IU)*SCAL,SGM1,1,CLag1,1)
          END IF

* SVC: The master node now continues to only handle task scheduling,
*      needed to achieve better load balancing. So it exits from the task
*      list. It has to do it here since each process gets at least one
*      task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GOTO 500
 501  CONTINUE

      CALL Free_Tsk(ID)

      CALL GETMEM ('Tasks','FREE','INTE',lTask,2*nTasks)

C 99  CONTINUE


      RETURN
      END
C
C-----------------------------------------------------------------------
C
      Subroutine CnstInt(Mode,INT1,INT2)
C
      Use CHOVEC_IO
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "chocaspt2.fh"
#include "caspt2_grad.fh"
C
      Real*8 INT1(nAshT,nAshT),INT2(nAshT,nAshT,nAshT,nAshT)
C
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
C
      Call DCopy_(nAshT**2,[0.0D+00],0,INT1,1)
      Call DCopy_(nAshT**4,[0.0D+00],0,INT2,1)
C
      iSym=1
      nFroI = nFro(iSym)
      nIshI = nIsh(iSym)
      nCorI = nFroI+nIshI
      nBasI = nBas(iSym)
      ! nOrbI = nOrb(iSym)

C     CALL GETMEM('FIMO  ','ALLO','REAL',ipFIMO,nBasSq)
      CALL GETMEM('WRK1  ','ALLO','REAL',ipWRK1,nBasSq)
      CALL GETMEM('WRK2  ','ALLO','REAL',ipWRK2,nBasSq)
C     Call SQUARE(Work(LFIMO),Work(ipFIMO),1,nOrbI,nOrbI)
C
C     --- One-Electron Integral
C
      !! Read H_{\mu \nu}
C     IRC=-1
C     IOPT=6
C     ICOMP=1
C     ISYLBL=1
C     CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,WRK2,ISYLBL)
C     !! triangular -> square transformation
C     Call Square(WRK2,WRK1,1,nBasT,nBasT)
C     !! AO -> MO transformation
C     Call DGemm_('T','N',nBasT,nBasT,nBasT,
C    *            1.0D+00,Work(LCMOPT2),nBasT,WRK1,nBasT,
C    *            0.0D+00,WRK2,nBasT)
C     Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *            1.0D+00,WRK2,nBasT,Work(LCMOPT2),nBasT,
C    *            0.0D+00,WRK1,nBasT)
      !! Inactive energy
C     Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C       RIn_Ene = RIn_Ene + 2.0d+00*WRK1(iCorI,iCorI)
C     End Do
      !! Put in INT1
C     Do iAshI = 1, nAsh(iSym)
C       Do jAshI = 1, nAsh(iSym)
C         Val = WRK1(nCorI+iAshI,nCorI+jAshI)
C         INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
C       End Do
C     End Do
      Do iAshI = 1, nAsh(iSym)
        Do jAshI = 1, nAsh(iSym)
C         Val = FIMO(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
          Val = Work(ipFIMO+nCorI+iAshI-1+nBasI*(nCorI+jAshI-1))
          INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
        End Do
      End Do
C
C     --- Two-Electron Integral
C
      iSymA = 1
      iSymI = 1
      iSymB = 1
      iSymJ = 1
C     If (.not.IfChol) Then
C       Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C         iOrb = iCorI
C         jOrb = iCorI
C         Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
C         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
C           RIn_Ene = RIn_Ene + 2.0d+00*WRK1(jCorI,jCorI)
C         End Do
C         Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
C         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
C           RIn_Ene = RIn_Ene - WRK1(jCorI,jCorI)
C         End Do
C       End Do
C     End If
C
      If (IfChol) Then
        Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
        Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
        Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
        DO JSYM=1,NSYM
          IB1=NBTCHES(JSYM)+1
          IB2=NBTCHES(JSYM)+NBTCH(JSYM)
C
          MXBGRP=IB2-IB1+1
          IF (MXBGRP.LE.0) CYCLE
          CALL GETMEM('BGRP','ALLO','INTE',LBGRP,2*MXBGRP)
          IBGRP=1
          DO IB=IB1,IB2
           IWORK(LBGRP  +2*(IBGRP-1))=IB
           IWORK(LBGRP+1+2*(IBGRP-1))=IB
           IBGRP=IBGRP+1
          END DO
          NBGRP=MXBGRP
C
          CALL MEMORY_ESTIMATE(JSYM,IWORK(LBGRP),NBGRP,
     &                         NCHOBUF,MXPIQK,NADDBUF)
          CALL GETMEM('KETBUF','ALLO','REAL',LKET,NCHOBUF)
C         write(6,*) "nchobuf= ", nchobuf
C         write(6,*) "nbgrp= ", nbgrp
C         write(6,*) "nbtch= ", nbtch(jsym)
          Do IBGRP=1,NBGRP
C
            IBSTA=IWORK(LBGRP  +2*(IBGRP-1))
            IBEND=IWORK(LBGRP+1+2*(IBGRP-1))
C           write(6,*) ibsta,ibend
C
            NV=0
            DO IB=IBSTA,IBEND
              NV=NV+NVLOC_CHOBATCH(IB)
            END DO
C
            !! int2(tuvx) = (tu|vx)/2
            !! This can be computed without frozen orbitals
            Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                                Work(LKET),nKet,
     &                                IBSTA,IBEND)
C
C           CALL GETMEM('WRKCHO','ALLO','REAL',LWRKCHO,NV)
C           Call DCopy_(NV,[0.0D+00],0,Work(LWRKCHO),1)
            If (IBGRP.EQ.1) SCAL = 0.0D+00
            If (IBGRP.NE.1) SCAL = 1.0D+00
            Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,
     *                  0.5D+00,Work(LKET),NASH(JSYM)**2,
     *                          Work(LKET),NASH(JSYM)**2,
     *                  SCAL   ,INT2,NASH(JSYM)**2)
C           CALL GETMEM('WRKCHO','FREE','REAL',LWRKCHO,NV)
          End Do
          CALL GETMEM('KETBUF','FREE','REAL',LKET,NCHOBUF)
          CALL GETMEM('BGRP','FREE','INTE',LBGRP,2*MXBGRP)
        End Do
      Else
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,
     *                Work(ipWRK1),Work(ipWRK2))
            !! Put in INT1
C           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C             INT1(iAshI,jAshI) = INT1(iAshI,jAshI)
C    *          + 2.0d+00*WRK1(iCorI,iCorI)
C           End Do
            !! Put in INT2
            Do kAshI = 1, nAsh(iSym)
              Do lAshI = 1, nAsh(iSym)
                INT2(iAshI,jAshI,kAshI,lAshI)
     *        = INT2(iAshI,jAshI,kAshI,lAshI)
C    *        + WRK1(nCorI+kAshI,nCorI+lAshI)*0.5d+00
     *        + Work(ipWRK1+nCorI+kAshI-1+nBasT*(nCorI+lAshI-1))*0.5d+00
              End Do
            End Do
C
C           Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT1
C           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C             INT1(iAshI,jAshI) = INT1(iAshI,jAshI) - WRK1(iCorI,iCorI)
C           End Do
          End Do
        End Do
      End If
C     CALL GETMEM('FIMO  ','FREE','REAL',ipFIMO,nBasSq)
      CALL GETMEM('WRK1  ','FREE','REAL',ipWRK1,nBasSq)
      CALL GETMEM('WRK2  ','FREE','REAL',ipWRK2,nBasSq)
C     write(6,*) "int2"
C     call sqprt(int2,25)
C     call sqprt(int1,5)
C     call sqprt(fimo,12)
      If (Mode.eq.0) Then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          iTU = iT + nAshT*(iU-1)
          Do iV = 1, nAshT
            Do iX = 1, nAshT
              iVX = iV + nAshT*(iX-1)
              If (iVX.gt.iTU) Then
               INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX) + INT2(iV,iX,iT,iU)
               INT2(iV,iX,iT,iU) = 0.0D+00
              End If
            End Do
          End Do
        End Do
      End Do
      End If
C
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          Do iX = 1, nAshT
C           INT1(IT,IU) = INT1(IT,IU) - INT2(IT,IX,IX,IU)
          End Do
        End Do
      End Do
C
      Return
C
      End Subroutine CnstInt
C
C-----------------------------------------------------------------------
C
      Subroutine CnstAntiC(DPT2Canti,UEFF,U0)
C
      use caspt2_gradient, only: iRoot1, iRoot2
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Dimension DPT2Canti(*),UEFF(nState,nState),U0(*)
C
      Call GetMem('CI1   ','ALLO','REAL',ipCI1 ,nConf)
      Call GetMem('CI2   ','ALLO','REAL',ipCI2 ,nConf)
      Call GetMem('SGM1  ','ALLO','REAL',ipSGM1,nConf)
      Call GetMem('SGM2  ','ALLO','REAL',ipSGM2,nConf)
      Call GetMem('TG1   ','ALLO','REAL',ipTG1 ,nAshT**2)
      Call GetMem('G1    ','ALLO','REAL',ipG1  ,nAshT**2)
C
      Call DCopy_(nAshT**2,[0.0D+00],0,Work(ipG1),1)
      Do iStat = 1, nState
        !! UEFF is in the CASSCF basis, so the CI coefficient
        !! has to be transformed(-back) accordingly
        Call LoadCI_XMS('N',1,Work(ipCI1),iStat,U0)
        Do jStat = 1, nState
          If (iStat.eq.jStat) Cycle
          Call LoadCI_XMS('N',1,Work(ipCI2),jStat,U0)
C
          Call Dens1T_RPT2(Work(ipCI1),Work(ipCI2),
     *                     Work(ipSGM1),Work(ipTG1))
          Scal = UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1)
     *         - UEFF(jStat,iRoot2)*UEFF(iStat,iRoot1)
          Scal = Scal*0.5d+00
C         write (*,*) istat,jstat,scal
C         call sqprt(work(iptg1),5)
          Call DaXpY_(nAshT**2,Scal,Work(ipTG1),1,Work(ipG1),1)
        End Do
      End Do
C     call sqprt(work(ipg1),5)
C
      Call GetMem('CI1   ','FREE','REAL',ipCI1 ,nConf)
      Call GetMem('CI2   ','FREE','REAL',ipCI2 ,nConf)
      Call GetMem('SGM1  ','FREE','REAL',ipSGM1,nConf)
      Call GetMem('SGM2  ','FREE','REAL',ipSGM2,nConf)
      Call GetMem('TG1   ','FREE','REAL',ipTG1 ,nAshT**2)
C
      !! The DPT2Canti computed so far has been doubled
      Call DScal_(nBasT**2,0.5d+00,DPT2Canti,1)
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2.gt.0) Then
          !! Add active orbital density
          !! Probably incorrect if symmetry
          Do iOrb0 = 1, nAsh(iSym)
            ! iOrb1 = nIsh(iSym)+iOrb0
            iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
            Do jOrb0 = 1, nAsh(iSym)
              ! jOrb1 = nIsh(iSym)+jOrb0
              jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
              DPT2Canti(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          = DPT2Canti(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))
     *          + Work(ipG1+iOrb0-1+nAsh(iSym)*(jOrb0-1))
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
C     write (*,*) "dpt2anti after"
C     call sqprt(dpt2canti,nbast)
C
      Call GetMem('G1    ','FREE','REAL',ipG1  ,nAshT**2)
C
      !! Add orbital response
      Call GetMem('WRK1  ','ALLO','REAL',ipWRK1,nBasT**2)
      Call GetMem('WRK2  ','ALLO','REAL',ipWRK2,nBasT**2)
C
      Call DCopy_(nBasT**2,DPT2Canti,1,Work(ipWRK1),1)
      !! Scale with the CASPT2 energy difference
      Scal = ENERGY(iRoot2)-ENERGY(iRoot1)
      Call DScal_(nBasT**2,Scal,Work(ipWRK1),1)
      !! anti-symmetrize the orbital response
      Call DGeSub(Work(ipWRK1),nBas(1),'N',
     &            Work(ipWRK1),nBas(1),'T',
     &            Work(ipWRK2),nBas(1),
     &            nBas(1),nBas(1))
C      write (*,*) "ipwrk1"
C      call sqprt(work(ipwrk1),nbast)
C
      !! Probably, the way CSF term is computed in MOLCAS is different
      !! from in BAGEL; see Eqs.(51)--(53). In the active space, i.e.
      !! for SA-CASSCF, the orbital response in the active space cancel
      !! (see JCP 2004, 120, 7322), but not in off-diagonal blocks.
      !! The way MOLCAS computes adds more than the one BAGEL does, so
      !! the orbital response has to be subtracted?
      Call DaXpY_(nBasT**2,-1.0d+00,Work(ipWRK1),1,Work(ipOLagFull),1)
C
      Call GetMem('WRK1  ','FREE','REAL',ipWRK1,nBasT**2)
      Call GetMem('WRK2  ','FREE','REAL',ipWRK2,nBasT**2)
C
      !! S[x]^A:D = S[x]:D^A, according to alaska_util/csfgrad.f
      !! so antisymmetrize D
      Do i = 1, nBasT
        Do j = 1, i-1
          Scal = DPT2Canti(i+nBasT*(j-1))
     *         - DPT2Canti(j+nBasT*(i-1))
          DPT2Canti(i+nBasT*(j-1)) =  Scal*0.5D+00
          DPT2Canti(j+nBasT*(i-1)) = -Scal*0.5D+00
        End Do
      End Do
C     write (*,*) "dpt2anti sym"
C     call sqprt(dpt2canti,nbast)
C     write (*,*) "dpt2c"
C     call sqprt(work(ipdpt2c),nbast)
C
      Return
C
      End Subroutine CnstAntiC
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DENS1T_RPT2 (CI1,CI2,SGM1,G1)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_output, only:iPrGlb,debug
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      LOGICAL RSV_TSK

      REAL*8 CI1(MXCI),CI2(MXCI),SGM1(MXCI)
      REAL*8 G1(NLEV,NLEV)

      REAL*8 GTU

      INTEGER ID
      INTEGER IST,ISU,ISTU
      INTEGER IT,IU,LT,LU

      INTEGER ITASK,LTASK,LTASK2T,LTASK2U,NTASKS

      INTEGER ISSG,NSGM

      REAL*8, EXTERNAL :: DDOT_,DNRM2_

c Purpose: Compute the 1- and 2-electron density matrix
c arrays G1 and G2.


      CALL DCOPY_(NG1,[0.0D0],0,G1,1)

C For the special cases, there is no actual CI-routines involved:
c Special code for hi-spin case:
      IF(ISCF.EQ.2) THEN
        DO IT=1,NASHT
          G1(IT,IT)=1.0D00
        END DO
        GOTO 99
      END IF
c Special code for closed-shell:
      IF(ISCF.EQ.1 .AND. NACTEL.GT.0) THEN
        DO IT=1,NASHT
          G1(IT,IT)=2.0D00
        END DO
        GOTO 99
      END IF

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.fh

C-SVC20100311: set up a task table with LT,LU
      nTasks=(nLev**2+nLev)/2
      nTasks= nLev**2

      CALL GETMEM ('Tasks','ALLO','INTE',lTask,2*nTasks)
      lTask2T=lTask
      lTask2U=lTask+nTasks

      iTask=0
      DO LT=1,nLev
        DO LU=1,nLev!LT
          iTask=iTask+1
          iWork(lTask2T+iTask-1)=LT
          iWork(lTask2U+iTask-1)=LU
        ENDDO
      ENDDO
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

C-SVC20100311: BEGIN SEPARATE TASK EXECUTION
 500  If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
C     LTU=0
C     DO 140 LT=1,NLEV
      LT=iWork(lTask2T+iTask-1)
        IST=ISM(LT)
        IT=L2ACT(LT)
C       DO 130 LU=1,LT
        LU=iWork(lTask2U+iTask-1)
C         LTU=LTU+1
          ! LTU=iTask
          ISU=ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,STSYM)
          NSGM=NCSF(ISSG)
          IF(NSGM.EQ.0) GOTO 500
          CALL GETSGM2(LU,LT,STSYM,CI1,SGM1)
          IF(ISTU.EQ.1) THEN
            GTU=DDOT_(NSGM,CI2,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
          END IF

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
#if defined (_MOLCAS_MPP_) && !defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GOTO 500
 501  CONTINUE
      CALL Free_Tsk(ID)

      CALL GETMEM ('Tasks','FREE','INTE',lTask,2*nTasks)

      CALL GAdSUM (G1,NG1)

  99  CONTINUE

      IF(iPrGlb.GE.DEBUG) THEN
        WRITE(6,'("DEBUG> ",A)')
     &   "DENS2_RPT2: norms of the density matrices:"
        WRITE(6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
      ENDIF


      RETURN
      END
C
C-----------------------------------------------------------------------
C
      Subroutine DWDER(OMGDER,HEFF,SLag)
C
      use definitions, only:wp
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      Dimension OMGDER(nState,nState),HEFF(nState,nState),
     *          SLag(nState,nState)
C
      Do ilStat = 1, nState
        Ebeta = HEFF(ilStat,ilStat)
        Do jlStat = 1, nState
          Ealpha = HEFF(jlStat,jlStat)
C
          Factor = 0.0d+00
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWType.EQ.1) Then
              Factor = Factor + exp(-zeta*(Ealpha - Egamma)**2)
            Else If (DWType.EQ.2) Then
              Factor = Factor
     *          + exp(-zeta*(Ealpha/HEFF(jlStat,klStat))**2)
            Else If (DWType.EQ.3) Then
              Dag = abs(Ealpha - Egamma) + 1.0e-9_wp
              Hag = abs(HEFF(jlStat,klStat))
              Factor = Factor
     *          + exp(-zeta*Dag/(sqrt(Hag)+tiny(Hag)))
            End If
          End Do
C
          DEROMG = OMGDER(ilStat,jlStat)
C
          !! derivative of alpha-beta
          If (DWType.EQ.1) Then
            DERAB = EXP(-ZETA*(Ealpha-Ebeta)**2)/Factor
            Scal = -2.0D+00*ZETA*DERAB*(Ealpha-Ebeta)*DEROMG
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
          Else If (DWType.EQ.2) Then
            DERAB = EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
            DERAB = DERAB*2.0D+00*DEROMG*ZETA
     *            *(Ealpha/HEFF(jlStat,ilStat))**2
            Scal  = -DERAB/Ealpha
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            Scal  = +DERAB/HEFF(jlStat,ilStat)
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          Else If (DWType.EQ.3) Then
            Dag = abs(Ealpha - Ebeta) + 1.0e-9_wp
            Hag = abs(HEFF(jlStat,ilStat))
            DERAB = EXP(-ZETA*Dag/(sqrt(Hag)+Tiny(Hag)))/Factor
            DERAB = -ZETA*DERAB*Dag/(sqrt(Hag)+tiny(Hag))
     *              *DEROMG
            Scal = DERAB/Dag
            If (Ealpha-Ebeta.le.0.0d+00) Scal = -Scal
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
            Scal  = -DERAB/(sqrt(Hag)+tiny(Hag))/sqrt(Hag)*0.5d+00
            If (HEFF(jlStat,ilStat).le.0.0d+00) Scal = -Scal
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          End If
C
          !! derivative of alpha-gamma
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWtype.EQ.1) Then
              DERAC = EXP(-ZETA*(Ealpha-Ebeta)**2)/(Factor*Factor)
     *               *EXP(-ZETA*(Ealpha-Egamma)**2)
              Scal = 2.0D+00*ZETA*DERAC*(Ealpha-Egamma)*DEROMG
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
            Else If (DWType.EQ.2) Then
              DERAC = EXP(-ZETA*(Ealpha/HEFF(jlStat,klStat))**2)/Factor
     *              * EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
              DERAC =-DERAC*2.0D+00*DEROMG*ZETA
     *              *(Ealpha/HEFF(jlStat,klStat))**2
              Scal  = -DERAC/Ealpha
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              Scal  = +DERAC/HEFF(jlStat,klStat)
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            Else IF (DWType.EQ.3) Then
              Dbg = abs(Ealpha - Egamma) + 1.0e-9_wp
              Hbg = abs(HEFF(jlStat,klStat))
              DERAC =EXP(-ZETA*Dag/(sqrt(Hag)+tiny(Hag)))/Factor
     *              *EXP(-ZETA*Dbg/(sqrt(Hbg)+tiny(Hbg)))/Factor
              DERAC = ZETA*DERAC*Dbg/(sqrt(Hbg)+tiny(Hbg))
     *                *DEROMG
              Scal = DERAC/Dbg
              If (Ealpha-Egamma.le.0.0d+00) Scal = -Scal
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
              Scal  = -DERAC/(sqrt(Hbg)+tiny(Hbg))/sqrt(Hbg)*0.5d+00
              If (HEFF(jlStat,klStat).le.0.0d+00) Scal = -Scal
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            End If
          End Do
        End Do
      End Do
C
      Return
C
      End Subroutine DWDER
