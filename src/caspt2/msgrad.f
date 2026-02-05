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
C
      Subroutine MS_Res(MODE,IST,JST,Scal)
      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV, only: IVECC, IVECC2, IVECW
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: STSYM, NCONF, NASHT, ISCF, NSTATE
      use pt2_guga, only: MXCI
      use Constants, only: Zero
C
C     Compute the derivative of E^PT2 with respct to the T amplitude
C
      implicit none

      integer(kind=iwp), intent(in) :: MODE, IST, JST
      real(kind=wp), intent(in) :: Scal

      integer(kind=iwp) :: I, NTG1, NTG2, NTG3, IDCI
      real(kind=wp) :: DVALUE, OVL, DUMMY(1)

      real(kind=wp),allocatable :: TG1(:),TG2(:),TG3(:),CI1(:),CI2(:)

C We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
C Note: Need proper allocation even if unused, sinced allocated
C arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      call mma_allocate(TG1,NTG1,Label='TG1')
      call mma_allocate(TG2,NTG2,Label='TG2')
      call mma_allocate(TG3,NTG3,Label='TG3')
      TG1(:) = Zero
      TG2(:) = Zero
      TG3(:) = Zero

      call mma_allocate(CI1,MXCI,Label='MCCI1')
      call mma_allocate(CI2,MXCI,Label='MCCI2')
      IF(ISCF == 0) THEN
C Read root vectors nr. IST and JST from LUCI.
        IDCI=IDTCEX
        DO I=1,NSTATE
          IF(I == IST) THEN
            CALL DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
            IF(I == JST) THEN
              CI2(1:NCONF) = CI1(1:NCONF)
            END IF
          ELSE IF(I == JST) THEN
            CALL DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          END IF
        END DO
      END IF

      CALL MKTG3(STSYM,STSYM,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)

!! Do similar to RHS_STRANS. Multiply the solution vector (T) with
!! the overlap-like term constructe with transition density
!! matrices. The output is IVECC (MODE=1) or IVECC2 (MODE=2)
      IF (MODE == 1) THEN
        CALL MS_STRANS(IVECW,IVECC,OVL,TG1,TG2,TG3,DVALUE,SCAL)
      ELSE IF (MODE == 2) THEN
        CALL MS_STRANS(IVECC,IVECC2,OVL,TG1,TG2,TG3,DVALUE,SCAL)
      ELSE IF (MODE == 3) THEN
        CALL MS_STRANS(IVECW,IVECC,OVL,TG1,TG2,TG3,DVALUE,SCAL)
      END IF

      call mma_deallocate(TG1)
      call mma_deallocate(TG2)
      call mma_deallocate(TG3)

      End Subroutine MS_Res
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MS_STRANS(IVEC,JVEC,OVL,TG1,TG2,TG3,HEL,SCAL)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM, NASHT, NASUP, NISUP, NINDEP, CASES
      use Constants, only: Zero
      use definitions, only: wp, iwp, u6

      implicit none
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

      integer(kind=iwp), intent(in) :: IVEC, JVEC
      real(kind=wp), intent(inout) :: OVL, TG1(NASHT,NASHT),
     &  TG2(NASHT,NASHT,NASHT,NASHT), TG3(*)
      real(kind=wp), intent(out) :: HEL
      real(kind=wp), intent(in) :: SCAL
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

      real(kind=wp) :: HECOMP(14,9), HEBLK, SUMSYM, SUMCASE
      integer(kind=iwp) :: ICASE, ISYM, NAS, NIN, NIS, lg_V1, iLo1,
     &  iHi1, jLo1, jHi1, MV1, lg_V2, iLo2, iHi2, jLo2, jHi2, MV2,
     &  NHECOMP, i, IC, IS

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

      HEL=Zero
      HECOMP=Zero
      DO ICASE=1,13
C     if (icase /= 12.and.icase /= 13) cycle ! H
C     if (icase /= 10.and.icase /= 11) cycle ! G
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          HEBLK=Zero

          if (NAS*NIS /= 0 .and. NIN /= 0) then
            CALL RHS_ALLO (NAS,NIS,lg_V1)
            CALL RHS_ALLO (NAS,NIS,lg_V2)
            CALL RHS_READ (NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
            CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
            CALL RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
            CALL RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

            IF ((iLo1 /= iLo2) .OR. (iHi1 /= iHi2) .OR.
     &          (jLo1 /= jLo2) .OR. (jHi1 /= jHi2)) THEN
              WRITE(u6,'(1X,A)')'HCOUP: Error: block mismatch, abort...'
              CALL ABEND()
            END IF

#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              CALL MS_STRANS_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                        DBL_MB(MV1),DBL_MB(MV2),OVL,
     &                        TG1,TG2,TG3,SCAL)
            ELSE
#endif
              CALL MS_STRANS_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                        GA_Arrays(MV1)%A,
     &                        GA_Arrays(MV2)%A,OVL,
     &                        TG1,TG2,TG3,SCAL)
#ifdef _MOLCAS_MPP_
            END IF
#endif
            !! Save T*S
            CALL RHS_SAVE (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
C
            CALL RHS_RELEASE (lg_V1,iLo1,iHi1,jLo1,jHi1)
            CALL RHS_RELEASE (lg_V2,iLo2,iHi2,jLo2,jHi2)
            CALL RHS_FREE (lg_V1)
            CALL RHS_FREE (lg_V2)
          end if
          HECOMP(ICASE,ISYM)=HEBLK
          HEL=HEL+HEBLK
        END DO
      END DO

C Sum-reduce the per-process contributions
      CALL GADGOP_SCAL(HEL,'+')
      NHECOMP=14*9
      CALL GADGOP(HECOMP,NHECOMP,'+')

      IF(IPRGLB >= DEBUG) THEN
        DO ICASE=1,13
          SUMSYM=Zero
          DO ISYM=1,NSYM
            SUMSYM=SUMSYM+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(ICASE,NSYM+1)=SUMSYM
        END DO

        DO ISYM=1,NSYM+1
          SUMCASE=Zero
          DO ICASE=1,13
            SUMCASE=SUMCASE+HECOMP(ICASE,ISYM)
          END DO
          HECOMP(14,ISYM)=SUMCASE
        END DO

        WRITE(u6,'(20a4)')('----',i=1,20)
       WRITE(u6,*)'HCOUP: The contributions to the Hamiltonian coupling'
        WRITE(u6,*)' elements, by case and by symmetry label.'
        DO IC=1,13
          WRITE(u6,'(1X,A8,9F12.8)')
     &      CASES(IC),(HECOMP(IC,IS),IS=1,NSYM+1)
        END DO
        CALL XFLUSH(u6)
        WRITE(u6,'(1X,A8,9F12.8)')
     &    'Summed: ', (HECOMP(14,IS),IS=1,NSYM+1)
        WRITE(6,*)
      END IF

      end subroutine MS_STRANS
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MS_STRANS_BLK(ICASE,ISYM,NAS,IISTA,IIEND,V1,V2,OVL,
     &                         TG1,TG2,TG3,SCAL)
      USE SUPERINDEX, only: MTU, MTUV, MTGEU, MTGTU
      use caspt2_module, only: NAES, NASHT, NTUVES, NTUES, NTGEUES,
     &                         NTGTUES
      use definitions, only: wp, iwp
      use Constants, only: Zero, Two, Four, Eight
C Compute a contribution to the coupling Hamiltonian element (HEL)
C defined as HEL = < ROOT1 | H * OMEGA | ROOT2 >. The contribution
C arises from the block V_(A,I), with A=1,NAS and I=IISTA,IIEND,
C with A the active superindex and I the inactive superindex. Since
C the inactive superindex is partitioned over processes, each process
C only computes part of the HEL value, which is then sum reduced in the
C calling subroutine.
      implicit none

      integer(kind=iwp), intent(in) :: ICASE, ISYM, NAS, IISTA, IIEND
      real(kind=wp), intent(in) :: V1(*), OVL, SCAL
      real(kind=wp), intent(inout) :: V2(*), TG1(NASHT,NASHT),
     &  TG2(NASHT,NASHT,NASHT,NASHT), TG3(*)
C The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

      integer(kind=iwp) :: NISBLK, IAS, IASABS, ITABS, IUABS, IVABS,
     &  JAS, JASABS, IXABS, IYABS, IZABS, IND1, IND2, IND3, JND1,
     &  JND2, JND3, ITG3, NAS1, IAS1, IAS2, JAS1, JAS2
      real(kind=wp) :: TMP, SA, SC, SBtuxy, SBtuyx, SBP, SBM, GUTXY,
     &  SD11, SD12, SD21, SD22, GUY, SE, SFtuxy, SFtuyx, SFP, SFM, SG

      IF (IISTA <= 0) RETURN

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
            IF(IND2 > IND3) THEN
              IF(IND1 > IND2) THEN
                JND1=IND1
                JND2=IND2
                JND3=IND3
              ELSE IF(IND1 > IND3) THEN
                JND1=IND2
                JND2=IND1
                JND3=IND3
              ELSE
                JND1=IND2
                JND2=IND3
                JND3=IND1
              END IF
            ELSE
              IF(IND1 > IND3) THEN
                JND1=IND1
                JND2=IND3
                JND3=IND2
              ELSE IF(IND1 > IND2) THEN
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
            IF(IYABS == IUABS) THEN
              TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
            END IF
            IF(IYABS == ITABS) THEN
              TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
              IF(IXABS == IUABS) THEN
                TMP=TMP+TG1(IVABS,IZABS)
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              TMP=TMP+TG2(IVABS,ITABS,IYABS,IZABS)
            END IF
C SA is the negative of this, and then some correction:
            SA=-TMP
            IF(IXABS == ITABS) THEN
              SA=SA+Two*TG2(IVABS,IUABS,IYABS,IZABS)
              IF(IYABS == IUABS) THEN
                SA=SA+Two*TG1(IVABS,IZABS)
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
            IF(IND2 > IND3) THEN
              IF(IND1 > IND2) THEN
                JND1=IND1
                JND2=IND2
                JND3=IND3
              ELSE IF(IND1 > IND3) THEN
                JND1=IND2
                JND2=IND1
                JND3=IND3
              ELSE
                JND1=IND2
                JND2=IND3
                JND3=IND1
              END IF
            ELSE
              IF(IND1 > IND3) THEN
                JND1=IND1
                JND2=IND3
                JND3=IND2
              ELSE IF(IND1 > IND2) THEN
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
            IF(IYABS == IUABS) THEN
              TMP=TMP+TG2(IVABS,IZABS,IXABS,ITABS)
            END IF
            IF(IYABS == ITABS) THEN
              TMP=TMP+TG2(IVABS,IUABS,IXABS,IZABS)
              IF(IXABS == IUABS) THEN
                TMP=TMP+TG1(IVABS,IZABS)
              END IF
            END IF
            IF(IXABS == IUABS) THEN
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
            SBtuxy=Two*TG2(IXABS,ITABS,IYABS,IUABS)
            SBtuyx=Two*TG2(IYABS,ITABS,IXABS,IUABS)
            IF(IXABS == ITABS) THEN
              SBtuxy=SBtuxy-Four*TG1(IYABS,IUABS)
              SBtuyx=SBtuyx+Two*TG1(IYABS,IUABS)
              IF(IYABS == IUABS) THEN
                SBtuxy=SBtuxy+Eight*OVL
                SBtuyx=SBtuyx-Four*OVL
              END IF
            END IF
            IF(IYABS == IUABS) THEN
              SBtuxy=SBtuxy-Four*TG1(IXABS,ITABS)
              SBtuyx=SBtuyx+Two*TG1(IXABS,ITABS)
            END IF
            IF(IYABS == ITABS) THEN
              SBtuxy=SBtuxy+Two*TG1(IXABS,IUABS)
              SBtuyx=SBtuyx-Four*TG1(IXABS,IUABS)
              IF(IXABS == IUABS) THEN
                SBtuxy=SBtuxy-Four*OVL
                SBtuyx=SBtuyx+Eight*OVL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              SBtuxy=SBtuxy+Two*TG1(IYABS,ITABS)
              SBtuyx=SBtuyx-Four*TG1(IYABS,ITABS)
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
            SBtuxy=Two*TG2(IXABS,ITABS,IYABS,IUABS)
            SBtuyx=Two*TG2(IYABS,ITABS,IXABS,IUABS)
            IF(IXABS == ITABS) THEN
              SBtuxy=SBtuxy-Four*TG1(IYABS,IUABS)
              SBtuyx=SBtuyx+Two*TG1(IYABS,IUABS)
              IF(IYABS == IUABS) THEN
                SBtuxy=SBtuxy+Eight*OVL
                SBtuyx=SBtuyx-Four*OVL
              END IF
            END IF
            IF(IYABS == IUABS) THEN
              SBtuxy=SBtuxy-Four*TG1(IXABS,ITABS)
              SBtuyx=SBtuyx+Two*TG1(IXABS,ITABS)
            END IF
            IF(IYABS == ITABS) THEN
              SBtuxy=SBtuxy+Two*TG1(IXABS,IUABS)
              SBtuyx=SBtuyx-Four*TG1(IXABS,IUABS)
              IF(IXABS == IUABS) THEN
                SBtuxy=SBtuxy-Four*OVL
                SBtuyx=SBtuyx+Eight*OVL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              SBtuxy=SBtuxy+Two*TG1(IYABS,ITABS)
              SBtuyx=SBtuyx-Four*TG1(IYABS,ITABS)
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
            SD11=Two*GUTXY
            SD12= -GUTXY
            SD21= -GUTXY
            SD22= -TG2(IXABS,ITABS,IUABS,IYABS)
            IF(ITABS == IXABS) THEN
              GUY=TG1(IUABS,IYABS)
              SD11=SD11+Two*GUY
              SD12=SD12 -GUY
              SD21=SD21 -GUY
              SD22=SD22+Two*GUY
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
            IF(IXABS == ITABS) SE=SE+Two*OVL
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
            IF(IXABS == ITABS) SE=SE+Two*OVL
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
            SFtuxy=Two*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=Two*TG2(ITABS,IYABS,IUABS,IXABS)

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
            SFtuxy=Two*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=Two*TG2(ITABS,IYABS,IUABS,IXABS)

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
            Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
C           Do i = 1, NISBLK
C             V2(JAS+NAS*(i-1)) = SG*V1(IAS+NAS*(i-1))
C           End Do
          END DO
        END DO
************************************************************************
      CASE(12)
        IF(ABS(OVL) >= 1.0e-12_wp) THEN
C         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      CASE(13)
        IF(ABS(OVL) >= 1.0e-12_wp) THEN
C         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
************************************************************************
      END SELECT
      Return
      end subroutine MS_STRANS_BLK
C
C-----------------------------------------------------------------------
C
      Subroutine LoadCI_XMS(Bas,Mode,CI,iState,U0)

      use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: IFXMS, IFRMS, NCONF, NSTATE
      use Constants, only: Zero

      implicit none

      character(len=1), intent(in) :: Bas
      integer(kind=iwp), intent(in) :: Mode, iState
      real(kind=wp), intent(inout) :: CI(Nconf)
      real(kind=wp), intent(in) :: U0(nState,nState)

      integer(kind=iwp) :: ID, I
      real(kind=wp),allocatable :: WRK(:)
C
C     MODE=0 is equivalent to LoadCI (XMS basis)
C     MODE=1 constructs the CI vector in CASSCF basis (back-transformed)
C     CSF in natural (Bas=N) or quasi-canonical (Bas=C) orbital basis
C
      If (Bas == 'N' .or. Bas == 'n') Then
        ID = IDCIEX !! natural
      Else If (Bas == 'C' .or. Bas == 'c') Then
        ID = IDTCEX !! quasi-canonical
      ELse
        write (u6,*)"the first argument in LoadCI_XMS should be either",
     *              "N (natural) or C (quasi-canonical)"
        call abend
      End If
C
      If (Mode == 0 .or. (.not.IFXMS .and. .not.IFRMS)) Then
        do I = 1, iState-1
          call ddafile(LUCIEX,0,CI,Nconf,ID)
        end do
        call ddafile(LUCIEX,2,CI,Nconf,ID)
      Else If (Mode == 1) Then
        CI(1:nconf) = Zero
        call mma_allocate(WRK,nconf,Label='WRK')
        do I = 1, nState
          call ddafile(LUCIEX,2,WRK,Nconf,ID)
          CI(1:nconf) = CI(1:nconf) + U0(iState,I)*WRK(1:nconf)
        end do
        call mma_deallocate(WRK)
      End If

      End Subroutine LoadCI_XMS
C
C-----------------------------------------------------------------------
C
      Subroutine XMS_Grad(H0,U0,UEFF,OMGDER)

      use caspt2_global, only: do_nac, do_csf, iRoot1, iRoot2,
     *                           CLag,CLagFull,OLag,DPT2_tot,
     *                           FIFA_all,FIFASA_all
      use caspt2_global, only: FIFA, TORB, NDREF
      use caspt2_global, only: CMOPT2, if_equalW, weight
      use gugx, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: ENERGY, IFXMS, IFRMS, IFDW, STSYM, NCONF,
     &                         NFRO, NISH, NASHT, NDEL, NBAS, NBAST,
     &                         NBSQT, NSTATE, ZETA, ORBIN
      use Constants, only: Zero, One, Half, Two, Quart
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      real(kind=wp), intent(in) :: H0(nState,nState), U0(nState,nState),
     &                             OMGDER(nState,nState)
      real(kind=wp), intent(inout) :: UEFF(nState,nState)

      real(kind=wp),allocatable :: CI1(:),CI2(:),SGM1(:),SGM2(:),TG1(:),
     *                             TG2(:),DG1(:),DG2(:),DG3(:),G1(:,:),
     *                             RDMEIG(:,:),DPT2(:),Trf(:),
     *                             RDMSA(:,:),WRK1(:),WRK2(:),DPT2_AO(:)

      real(kind=wp) :: SLag(nState*nState), Scal, OVL, EEI, EEJ, fact,
     &  Wgt, TRC, EINACT, EDIFF
      integer(kind=iwp) :: nLev, iStat, jStat, NTG1, NTG3, iState,
     &  kStat, lStat, nOrbI, nBasI, iAsh, jAsh, nCor, I, II
      real(kind=wp), external :: ddot_

      nLev = SGS%nLev
C
C     The XMS rotation applies to any variants: XMS-CASPT2, XDW-CASPT2,
C     and RMS-CASPT2.
C
      If (IFXMS .or. IFRMS) Then
        call mma_allocate(CI1,nConf,Label='CI1')
        call mma_allocate(CI2,nConf,Label='CI2')
        call mma_allocate(SGM1,nConf,Label='SGM1')
        call mma_allocate(SGM2,nConf,Label='SGM2')
        call mma_allocate(TG1,nAshT**2,Label='TG1')
        call mma_allocate(TG2,nAshT**4,Label='TG2')
C
        call mma_allocate(DG1,nAshT**2,Label='DG1')
        call mma_allocate(DG2,nAshT**4,Label='DG2')
C
!       ----- Construct pseudo-density matrix -----

        !! First, we need to consider the derivative of the rotation
        !! of the CASSCF energy (Heff[1]).

        !! Forward transformation of UEFF (CASSCF basis to XMS basis)
        !! SLag is used as a working array
        Call DGEMM_('T','N',nState,nState,nState,
     &              One,U0,nState,UEFF,nState,
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C
        !! Then actual calculation
        CLag(1:nconf,1:nstate) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,CI1,iStat,U0)
          Do jStat = 1, nState
            Call LoadCI_XMS('C',0,CI2,jStat,U0)
C
            Call Dens2T_RPT2(CI1,CI2,SGM1,SGM2,
     *                       TG1,TG2,nAshT)
            TG1(:) = TG1(:)*Half
            TG2(:) = TG2(:)*Half
            DG1(:) = Zero
            DG2(:) = Zero
            Call CnstInt(1,DG1,DG2)
C
            If (do_nac) Then
              Scal =(UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
     *             + UEFF(jStat,iRoot1)*UEFF(iStat,iRoot2))*Half
            Else
              Scal = UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
            End If
C       write (*,*) " scal in xms"
C       write (*,*) istat,jstat,scal
            if (IFDW .and. zeta >= Zero) then
              scal = scal + OMGDER(iStat,jStat)
            end if
            DG1(:) = DG1(:)*Scal
            DG2(:) = DG2(:)*Scal
C
            call mma_allocate(DG3,nAshT**6,Label='DG3')
            DG3(:) = Zero
C
            NTG1 = NASHT**2
            NTG3 = (NTG1*(NTG1+1)*(NTG1+2))/6
            OVL = Zero
            CALL DERTG3(.False.,STSYM,STSYM,CI1,CI2,OVL,
     &                  DG1,DG2,NTG3,DG3,CLag(1,iStat),CLag(1,jStat))
            call mma_deallocate(DG3)
          End Do
        End Do
        call mma_deallocate(SGM2)
        call mma_deallocate(TG2)
C
        !! Back transformation of UEFF (XMS basis to CASSCF basis)
        !! SLag is used as a working array
        Call DGEMM_('N','N',nState,nState,nState,
     &              One,U0,nState,UEFF,nState,
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C     write (6,*) "ueff in casscf basis"
C     call sqprt(ueff,nstate)
C
        call mma_deallocate(DG1)
        call mma_deallocate(DG2)
C
      !! Add to the full CI Lagrangian
      !! CLag is in quasi-canonical basis, so transformation
      !! to natural basis is required.
        IF(ORBIN == 'TRANSFOR') Then
          Do iState = 1, nState
            Call CLagX_TrfCI(CLag(1,iState))
          End Do
        End If
        CLagFull(:,:) = CLagFull(:,:) + CLag(:,:)
C
        !! Compute the Lagrange multiplier for XMS
        !! The diagonal element is always zero.
        !! The code has an additional scaling with 0.5,
        !! because some contributions are doubled.
C     write (*,*) "istat,jstat,scal"
        SLag(:) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,CI1,iStat,U0)
          EEI = H0(iStat,iStat)
          Do jStat = 1, iStat
            If (iStat == jStat) Then
              SLag(iStat+nState*(jStat-1)) = Zero
            Else
              Call LoadCI_XMS('N',0,CI2,jStat,U0)
              EEJ = H0(jStat,jStat)
C
              Scal = DDOT_(nConf,CI1,1,CLagFull(1,jStat),1)
     *             - DDOT_(nConf,CI2,1,CLagFull(1,iStat),1)
C      write (*,*) "original scal = ", scal
              If (do_csf) Then
                !! JCTC 2017, 13, 2561: eq.(66)
                !! iStat and jStat: XMS
                !! kStat and lStat: CASSCF
                fact = Zero
                Do kStat = 1, nState
                  Do lStat = 1, nState
                    fact = fact
     *                   + (UEFF(kStat,iRoot1)*UEFF(lStat,iRoot2)
     *                   -  UEFF(kStat,iRoot2)*UEFF(lStat,iRoot1))*Half
     *                   * U0(kStat,iStat)*U0(lStat,jStat)
                  End Do
                End Do
                Scal = Scal+fact*(ENERGY(iRoot2)-ENERGY(iRoot1))*Two
C      write (*,*) "scal after= ", scal
C      write (*,*) fact,energy(iroot2)-energy(iroot1)
              End If
              Scal = Quart*Scal/(EEJ-EEI)
C           write (*,*) istat,jstat,scal
              SLag(iStat+nState*(jStat-1)) = Scal
              SLag(jStat+nState*(iStat-1)) = Scal
            End If
          End Do
        End Do
C
      !! Either subtract the CI derivative of Heff(1)
      !! or add off-diagonal couplings in rhs_sa (Z-vector)
      !CLagFull = CLagFull - CLag
C
        !! Finally, construct the pseudo-density matrix
        !! d = \sum_{ST} w_{ST} * d_{ST}
        call mma_allocate(G1,nAshT,nAshT,Label='G1')
        G1(:,:) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,CI1,iStat,U0)
          Do jStat = 1, nState
            If (ABS(SLag(iStat+nState*(jStat-1))) <= 1.0e-10_wp) Cycle
            Call LoadCI_XMS('C',0,CI2,jStat,U0)
C
            Call Dens1T_RPT2(CI1,CI2,
     *                       SGM1,TG1,nLev)
            Scal = SLag(iStat+nState*(jStat-1))*Two
C         write (*,*) "istat,jstat=",istat,jstat
C         write (*,*) "scal = ", scal
            Call DaXpY_(nAshT**2,Scal,TG1,1,G1,1)
          End Do
        End Do
C     write (*,*) "G1"
C     call sqprt(G1,nasht)
C
        call mma_deallocate(CI1)
        call mma_deallocate(CI2)
        call mma_deallocate(SGM1)
        call mma_deallocate(TG1)
C
C       ----- Calculate orbital derivatives -----
C
        call mma_allocate(RDMEIG,nAshT,nAshT,Label='RDMEIG')
        call mma_allocate(DPT2,NBSQT,Label='DPT2')
        DPT2(:) = Zero
C
        call mma_allocate(Trf,NBSQT,Label='TRFMAT')
        call mma_allocate(RDMSA,nAshT,nAshT,Label='RDMSA')
        call mma_allocate(WRK1,NBSQT,Label='WRK1')
        call mma_allocate(WRK2,NBSQT,Label='WRK2')
C
        !! Construct always state-averaged density; XMS basis is always
        !! generated with the equally-averaged density.
        WRK1(1:nDRef) = Zero
        call mma_allocate(CI1,nConf,Label='CI1')
        Wgt  = One/nState
        Do iState = 1, nState
          Call LoadCI(CI1,iState)
          call POLY1(CI1,nConf)
          call GETDREF(WRK2,nDRef)
          WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*WRK2(1:nDRef)
        End Do
        call mma_deallocate(CI1)
        Call SQUARE(WRK1,RDMSA,1,nAshT,nAshT)
C
        nOrbI = nBas(1) - nDel(1) !! nOrb(1)
        nBasI = nBas(1)
C       Call SQUARE(FIFA,FIFA_all,1,nOrbI,nOrbI)
        !! FIFASA_all is in natural orbital basis
        Trf(1:NBSQT) = Zero
        Call CnstTrf(TOrb,Trf)
C
        !! FIFA: natural -> quasi-canonical
        If (IFDW .or. IFRMS) Then
          Call DGemm_('T','N',nOrbI,nOrbI,nOrbI,
     *                One,Trf,nBasI,FIFASA_all,nOrbI,
     *                Zero,WRK1,nOrbI)
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,
     *                One,WRK1,nOrbI,Trf,nBasI,
     *                Zero,FIFASA_all,nOrbI)
        End If
!       Call DCopy_(NBSQT,FIFASA_all,1,FIFA_all,1)
        FIFA_all(1:NBSQT) = FIFASA_all(1:NBSQT)
C
        !! Orbital derivatives of FIFA
        !! Both explicit and implicit orbital derivatives are computed
        !! Also, compute G(D) and put the active contribution to RDMEIG
        OLag(:) = Zero
        Call EigDer2(RDMEIG,Trf,FIFA_all,RDMSA,G1,WRK1,WRK2)
C
        !! Add to PT2 density
        !! No inactive contributions. Correct as long as CASSCF CI
        !! vector are orthogonal.
        Call AddDEPSA(DPT2,G1)
        Call DPT2_TrfStore(1.0D+00,DPT2,DPT2_tot,Trf,WRK1)
C
        !! Finalize OLag (anti-symetrize) and construct WLag
        Call OLagFinal(OLag,Trf)
C
        If (.not.if_equalW) then
          call mma_allocate(DPT2_AO,NBSQT,Label='DPT2_AO')
C
          !! AddDEPSA considers the frozen orbital, whereas DPT2_Trf
          !! does not. In any case, construct DPT2 again.
          DPT2(:) = Zero
          CALL DPT2_Trf(DPT2,DPT2_AO,CMOPT2,G1,WRK1)
          !! Construct the SCF density
          WRK1(1:nDRef) = Zero
          call mma_allocate(CI1,nConf,Label='CI1')
          Do iState = 1, nState
            Call LoadCI_XMS('N',1,CI1,iState,U0)
            call POLY1(CI1,nConf)
            call GETDREF(WRK2,nDRef)
            wgt = Weight(iState)
            WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*WRK2(1:nDRef)
          End Do
          call mma_deallocate(CI1)
          !! WRK2 is the SCF density (for nstate=nroots)
          Call SQUARE(WRK1,WRK2,1,nAshT,nAshT)
          Call DaXpY_(nAshT**2,-One,WRK2,1,RDMSA,1)
          !! Construct the SS minus SA density matrix in WRK1
          Call OLagFroD(WRK1,WRK2,RDMSA,Trf)
          !! Subtract the inactive part
          Call DaXpY_(nBasT**2,-One,WRK2,1,WRK1,1)
          !! Save
          Call CnstAB_SSDM(DPT2_AO,WRK1)
          call mma_deallocate(DPT2_AO)
        End If
C
        call mma_deallocate(RDMSA)
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)
C
        call mma_deallocate(DPT2)
C
C       ----- Calculate CI derivatives -----
C
        !! use quasi-canonical CSF rather than natural CSF
C     ISAV = IDCIEX
C     IDCIEX = IDTCEX
        CLag(:,:) = Zero
C
        !! 1) Explicit CI derivative
        !! a: Extract FIFA in the AS for explicit CI derivative
        !!    Note that this FIFA uses state-averaged density matrix
C       call sqprt(fifa_all,nbast)
        Do iAsh = 1, nAshT
          Do jAsh = 1, nAshT
            G1(iAsh,jAsh) = FIFA_all( nFro(1)
     &          + nIsh(1) + iAsh + nBas(1)*(nFro(1)+nIsh(1)+jAsh-1))
          End Do
        End Do
C       call sqprt(g1),nasht)
        !! Transform quasi-canonical to natural
        nCor = nFro(1)+nIsh(1)
C     write (*,*) nfro(1),nish(1)
C     call sqprt(trf,nbast)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              One,Trf(1+nBasT*nCor+nCor),nBasT,G1,nAshT,
     *              Zero,FIFA_all,nAshT)
C     call sqprt(FIFA_all,nasht)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,
     *              One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,
     *              Zero,G1,nAshT)
C     call sqprt(g1,nasht)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,
     *              One,Trf(1+nBasT*nCor+nCor),nBasT,RDMEIG,nAshT,
     *              Zero,FIFA_all,nAshT)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,
     *              One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,
     *              Zero,RDMEIG,nAshT)
        call mma_deallocate(Trf)
C
        !! b: FIFA in the inactive space
        TRC=Zero
C     DO ISYM=1,NSYM
        DO I=1,NISH(1) ! ISYM)
C         II=IOFF(ISYM)+(I*(I+1))/2
          II=(I*(I+1))/2
          TRC=TRC+FIFA(II)
        END DO
C     END DO
* Contribution from inactive orbitals:
        EINACT=Two*TRC
        !! This EINACT may be wrong. Perhaps, FIFA has to be
        !! back-transformed to natural orbital basis. However, this does
        !! not contribute to the final gradient as long as all the
        !! (internal) CI vectors are orthogonal.
C
        !! c: Finally, compute explicit CI derivative
        !! y_{I,T} = w_{ST} <I|f|S>
        !! Here, G1 is FIFA = ftu
        Call CLagEigT(CLag,G1,SLag,EINACT)
C
        !! 2) Implicit CI derivative
        Call CLagEig(.False.,.True.,CLag,RDMEIG,nLev)
C
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          call GADSUM(CLag,nConf*nState)
        end if
#endif
C
        call mma_deallocate(RDMEIG)
        call mma_deallocate(G1)
C
      End If
C
      If (do_csf) Then
        !! Eq (68)
        !! I think this contributes even in XMS-CASPT2.
        !! However, maybe I implemented other terms in a different way
        !! from BAGEL.
        call mma_allocate(CI1,nConf,Label='CI1')
        !! Do in the natural orbital in the XMS basis
        !! because CLag is like that
        !! CASSCF -> XMS
        Call DGEMM_('T','N',nState,nState,nState,
     *              One,U0,nState,UEFF,nState,
     *              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C
        EDIFF = ENERGY(iRoot2)-ENERGY(iRoot1)
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,CI1,iStat,U0)
          Do jStat = 1, nState
            Scal = (UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
     *           -  UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1))*Half
            Scal = Scal * EDIFF
            Call DaXpY_(nConf,Scal,CI1,1,CLag(1,jStat),1)
          End Do
        End Do
C
        !! XMS -> CASSCF
        Call DGEMM_('N','N',nState,nState,nState,
     *              One,U0,nState,UEFF,nState,
     *              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
C
        call mma_deallocate(CI1)
      End If
C
      !! Finally, add to the total CI derivative array
      CLagFull(:,:) = CLagFull(:,:) + CLag(:,:)
C     IDCIEX = ISAV
C
      End Subroutine XMS_Grad
C
C-----------------------------------------------------------------------
C
      !! From poly3
      SUBROUTINE CLagEigT(CLag,RDMEIG,SLag,EINACT)

      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NCONF, ISCF, NSTATE
      use Constants, only: One, Two
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: nProcs, Is_Real_Par
#endif

      implicit none

      real(kind=wp), intent(inout) :: CLag(nConf,nState)
      real(kind=wp), intent(in) :: RDMEIG(*), SLag(*), EINACT

      real(kind=wp),allocatable :: CI1(:),CI2(:)
      real(kind=wp) :: Scal
      integer(kind=iwp) :: iStat, jStat

      !! RDMEIG
      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')
      Do iStat = 1, nState
        If (ISCF == 0) Then
          Call LoadCI(CI1,iStat)
        Else
          CI1(1) = One
        End If
        !! Skip iStat = jStat because SLag is zero
        Do jStat = 1, nState !! iStat-1
          If (ISCF == 0) Then
            Call LoadCI(CI2,jStat)
          Else
            CI2(1) = One
          End If
          !! One of doubling is due to the scaling factor
          !! One of doubling is due to the symmetry of iStat and jStat
C         Scal = SLag(iStat+nState*(jStat-1))*4.0d+00
          Scal = SLag(iStat+nState*(jStat-1))*Two
          If (ABS(Scal) <= 1.0e-09_wp) Cycle
C
          Call Poly1_CLagT(CI1,CI2,
     *                     CLag(1,iStat),CLag(1,jStat),RDMEIG,Scal)
          !! Inactive terms
#ifdef _MOLCAS_MPP_
          !! The inactive contributions are computed in all processes,
          !! whereas GADSUM will be done later, so divide
          if (is_real_par()) Scal = Scal/real(nProcs,kind=wp)
#endif
          CLag(1:nconf,jStat) = CLag(1:nconf,jStat)
     &      + Scal*EINACT*CI1(1:nconf)
          CLag(1:nconf,iStat) = CLag(1:nconf,iStat)
     &      + Scal*EINACT*CI2(1:nconf)
        End Do
      End Do
C
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
C
      Return
C
      End Subroutine CLagEigT
C
C-----------------------------------------------------------------------
C
      SUBROUTINE POLY1_CLagT(CI1,CI2,CLag1,CLag2,RDMEIG,Scal)
      use gugx, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: nConf
      use pt2_guga, only: MxCI, iAdr10, cLab10

      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES THE 1-EL DENSITY
* MATRIX FOR A CASSCF WAVE FUNCTION.

      real(kind=wp), intent(in) :: CI1(NCONF), CI2(NCONF), RDMEIG(*),
     &                             Scal
      real(kind=wp), intent(inout) :: CLag1(*), CLag2(*)

      real(kind=wp),allocatable :: SGM1(:)
      integer(kind=iwp) :: nLev, I

      nLev = SGS%nLev

      IF(NLEV > 0) THEN
        call mma_allocate(SGM1,MXCI,Label='SGM1')
        CALL DENS1T_RPT2_CLag(CI1,CI2,SGM1,
     *                        CLag1,CLag2,RDMEIG,Scal,nLev)
      END IF

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in pt2_guga.F90
* CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV > 0) THEN
        call mma_deallocate(SGM1)
      END IF

      RETURN
      end subroutine POLY1_CLagT
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DENS1T_RPT2_CLag(CI1,CI2,SGM1,CLag1,CLag2,RDMEIG,SCAL,
     &                            nLev)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use gugx, only: SGS, L2ACT, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: nConf, STSym, Mul
      use pt2_guga, only: MxCI

      IMPLICIT NONE

      real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI),
     &                             RDMEIG(NLEV,NLEV), SCAL
      real(kind=wp), intent(inout) :: SGM1(MXCI), CLag1(nConf),
     &                                CLag2(nConf)
      integer(kind=iwp), intent(in) :: nLev

      integer(kind=iwp),allocatable :: TASK(:,:)

      logical(kind=iwp), external :: RSV_TSK
      integer(kind=iwp) :: ID, IST, ISU, ISTU, IT, IU, LT, LU, ITASK,
     &                     NTASKS, ISSG, NSGM

* Purpose: Compute the 1-electron density matrix array G1.

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.F90

* SVC20100311: set up a task table with LT,LU
* SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks=(nLev**2+nLev)/2
      nTasks = nLev**2
      CALL mma_allocate (Task,nTasks,2,Label='TASK')

      iTask=0
      ! First, IL < JL pairs.
      Do LT = 1, nLev-1
        Do LU = LT+1, nLev
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      ! Then, IL = JL pairs.
      Do LT = 1, nLev
        iTask = iTask + 1
        TASK(iTask,1)=LT
        TASK(iTask,2)=LT
      End Do
      ! Last, IL > JL pairs.
      Do LT = 2, nLev
        Do LU = 1, LT-1
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      IF (iTask /= nTasks) WRITE(u6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

* SVC20100311: BEGIN SEPARATE TASK EXECUTION
      do while (Rsv_Tsk(ID,iTask))
* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
        LT=TASK(iTask,1)
        IST=SGS%ISM(LT)
        IT=L2ACT(LT)
        LU=Task(iTask,2)
        ! LTU=iTask
        ISU=SGS%ISM(LU)
        IU=L2ACT(LU)
        ISTU=MUL(IST,ISU)
        ISSG=MUL(ISTU,STSYM)
        NSGM=CIS%NCSF(ISSG)
        IF(NSGM == 0) cycle
* GETSGM2 computes E_UT acting on CI and saves it on SGM1
        IF(ISTU == 1) THEN
          CALL GETSGM2(LU,LT,STSYM,CI1,SGM1)
          CLag2(1:NSGM) = CLag2(1:NSGM)
     &      + SCAL*RDMEIG(IT,IU)*SGM1(1:NSGM)
          CALL GETSGM2(LU,LT,STSYM,CI2,SGM1)
          CLag1(1:NSGM) = CLag1(1:NSGM)
     &      + SCAL*RDMEIG(IT,IU)*SGM1(1:NSGM)
        END IF
* SVC: The master node now continues to only handle task scheduling,
*      needed to achieve better load balancing. So it exits from the task
*      list. It has to do it here since each process gets at least one
*      task.
      end do
      CALL Free_Tsk(ID)
      CALL mma_deallocate(Task)

      end subroutine DENS1T_RPT2_CLag
C
C-----------------------------------------------------------------------
C
      Subroutine CnstInt(Mode,INT1,INT2)
C
      Use CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only: FIMO_all
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp,wp
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NASH, NASHT,
     &                         NSSH, NBAS, NBAST, NBSQT, NBTCH, NBTCHES
      use Constants, only: Zero, One, Half
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      integer(kind=iwp), intent(in) :: Mode
      real(kind=wp), intent(out) :: INT1(nAshT,nAshT),
     &                              INT2(nAshT,nAshT,nAshT,nAshT)

      integer(kind=iwp),allocatable :: BGRP(:,:)
      real(kind=wp),allocatable :: WRK1(:),WRK2(:),KET(:)
C
      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) ::  nSh(8,3), iSym, nFroI, nIshI, nCorI, nBasI,
     &  iAshI, jAshI, iSymA, iSymI, iSymB, iSymJ, JSYM, IB, IB1, IB2,
     &  MXBGRP, IBGRP, NBGRP, NCHOBUF, MXPIQK, NADDBUF, IBSTA, IBEND,
     &  NV, nKET, kAshI, lAshI, iT, iU, iTU, iV, iX, iVX, iOrb, jOrb
      real(kind=wp) :: Val, SCAL
C
      INT1(:,:) = Zero
      Int2(:,:,:,:) = Zero
C
      iSym=1
      nFroI = nFro(iSym)
      nIshI = nIsh(iSym)
      nCorI = nFroI+nIshI
      nBasI = nBas(iSym)
      ! nOrbI = nOrb(iSym)

      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')
C
C     --- One-Electron Integral
C
      !! Read H_{\mu \nu}
      Do iAshI = 1, nAsh(iSym)
        Do jAshI = 1, nAsh(iSym)
          Val = FIMO_all(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
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
C
      If (IfChol) Then
        nSh(1:nSym,Inactive) = NISH(1:nSym)
        nSh(1:nSym,Active  ) = NASH(1:nSym)
        nSh(1:nSym,Virtual ) = NSSH(1:nSym)
        DO JSYM=1,NSYM
          IB1=NBTCHES(JSYM)+1
          IB2=NBTCHES(JSYM)+NBTCH(JSYM)
C
          MXBGRP=IB2-IB1+1
          IF (MXBGRP <= 0) CYCLE
          call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
          IBGRP=1
          DO IB=IB1,IB2
           BGRP(1,IBGRP) = IB
           BGRP(2,IBGRP) = IB
           IBGRP=IBGRP+1
          END DO
          NBGRP=MXBGRP
C
          CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,
     &                         NCHOBUF,MXPIQK,NADDBUF)
          call mma_allocate(KET,NCHOBUF,Label='KETBUF')
          Do IBGRP=1,NBGRP
C
            IBSTA=BGRP(1,IBGRP)
            IBEND=BGRP(2,IBGRP)
C
            NV=0
            DO IB=IBSTA,IBEND
              NV=NV+NVLOC_CHOBATCH(IB)
            END DO
C
            !! int2(tuvx) = (tu|vx)/2
            !! This can be computed without frozen orbitals
            Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                                KET,nKet,
     &                                IBSTA,IBEND)
C
            If (IBGRP == 1) SCAL = Zero
            If (IBGRP /= 1) SCAL = One
            Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,
     *                  Half,KET,NASH(JSYM)**2,KET,NASH(JSYM)**2,
     *                  SCAL   ,INT2,NASH(JSYM)**2)
          End Do
          call mma_deallocate(KET)
          call mma_deallocate(BGRP)
        End Do
      Else
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT2
            Do kAshI = 1, nAsh(iSym)
              Do lAshI = 1, nAsh(iSym)
                INT2(iAshI,jAshI,kAshI,lAshI)
     *        = INT2(iAshI,jAshI,kAshI,lAshI)
     *        + WRK1(nCorI+kAshI+nBasT*(nCorI+lAshI-1))*Half
              End Do
            End Do
          End Do
        End Do
      End If
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
      If (Mode == 0) Then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          iTU = iT + nAshT*(iU-1)
          Do iV = 1, nAshT
            Do iX = 1, nAshT
              iVX = iV + nAshT*(iX-1)
              If (iVX > iTU) Then
               INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX) + INT2(iV,iX,iT,iU)
               INT2(iV,iX,iT,iU) = Zero
              End If
            End Do
          End Do
        End Do
      End Do
      End If
C
#ifdef _MOLCAS_MPP_
      if (is_real_par()) CALL GADSUM (INT2,nAshT**4)
#endif
C
      Return
C
      End Subroutine CnstInt
C
C-----------------------------------------------------------------------
C
      Subroutine CnstAntiC(DPT2Canti,UEFF,U0)

      use caspt2_global, only: iRoot1, iRoot2, OLagFull
      use gugx, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: ENERGY, NSYM, NCONF, NFRO, NISH, NASH,
     &                         NASHT, NDEL, NORB, NBAS, NBAST, NBSQT,
     &                         NSTATE
      use Constants, only: Zero, One, Half

      implicit none

      real(kind=wp), intent(inout) :: DPT2Canti(*)
      real(kind=wp), intent(in) :: UEFF(nState,nState), U0(*)

      real(kind=wp),allocatable :: CI1(:),CI2(:),SGM1(:),SGM2(:),TG1(:),
     *                             G1(:,:),WRK1(:),WRK2(:)
      integer(kind=iwp) :: nLev, iStat, jStat, iMO1, iMO2, iSym, nOrbI1,
     &  nOrbI2, iOrb0, iOrb2, jOrb0, jOrb2, i, j
      real(kind=wp) :: Scal

      nLev = SGS%nLev
C
      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')
      call mma_allocate(SGM1,nConf,Label='SGM1')
      call mma_allocate(SGM2,nConf,Label='SGM2')
      call mma_allocate(TG1,nAshT**2,Label='TG1')
      call mma_allocate(G1,nAshT,nAshT,Label='G1')
C
      G1(:,:) = Zero
      Do iStat = 1, nState
        !! UEFF is in the CASSCF basis, so the CI coefficient
        !! has to be transformed(-back) accordingly
        Call LoadCI_XMS('N',1,CI1,iStat,U0)
        Do jStat = 1, nState
          If (iStat == jStat) Cycle
          Call LoadCI_XMS('N',1,CI2,jStat,U0)
C
          Call Dens1T_RPT2(CI1,CI2,SGM1,TG1,nLev)
          Scal = UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1)
     *         - UEFF(jStat,iRoot2)*UEFF(iStat,iRoot1)
          Scal = Scal*Half
C         write (*,*) istat,jstat,scal
C         call sqprt(tg1,5)
          Call DaXpY_(nAshT**2,Scal,TG1,1,G1,1)
        End Do
      End Do
C     call sqprt(g1,5)
C
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
      call mma_deallocate(SGM1)
      call mma_deallocate(SGM2)
      call mma_deallocate(TG1)
C
      !! The DPT2Canti computed so far has been doubled
      DPT2Canti(1:nBasT**2) = DPT2Canti(1:nBasT**2)*Half
C
      iMO1 = 1
      iMO2 = 1
      DO iSym = 1, nSym
        nOrbI1 = nOrb(iSym)
        nOrbI2 = nBas(iSym)-nDel(iSym)
        If (nOrbI2 > 0) Then
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
     *          + G1(iOrb0,jOrb0)
            End Do
          End Do
        END IF
        iMO1 = iMO1 + nOrbI1*nOrbI1
        iMO2 = iMO2 + nOrbI2*nOrbI2
      End Do
C     write (*,*) "dpt2anti after"
C     call sqprt(dpt2canti,nbast)
C
      call mma_deallocate(G1)
C
      !! Add orbital response
      call mma_allocate(WRK1,NBSQT,Label='WRK1')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')
C
      !! Scale with the CASPT2 energy difference
      Scal = ENERGY(iRoot2)-ENERGY(iRoot1)
      WRK1(1:nBasT**2) = Scal*DPT2Canti(1:nBasT**2)
      !! anti-symmetrize the orbital response
      Call DGeSub(WRK1,nBas(1),'N',
     &            WRK1,nBas(1),'T',
     &            WRK2,nBas(1),
     &            nBas(1),nBas(1))
C      write (*,*) "wrk1"
C      call sqprt(wrk),nbast)
C
      !! Probably, the way CSF term is computed in MOLCAS is different
      !! from in BAGEL; see Eqs.(51)--(53). In the active space, i.e.
      !! for SA-CASSCF, the orbital response in the active space cancel
      !! (see JCP 2004, 120, 7322), but not in off-diagonal blocks.
      !! The way MOLCAS computes adds more than the one BAGEL does, so
      !! the orbital response has to be subtracted?
      OLagFull(1:nBasT**2) = OLagFull(1:nBasT**2) - WRK1(1:nBasT**2)
C
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
C
      !! S[x]^A:D = S[x]:D^A, according to alaska_util/csfgrad.f
      !! so antisymmetrize D
      Do i = 1, nBasT
        Do j = 1, i-1
          Scal = DPT2Canti(i+nBasT*(j-1))
     *         - DPT2Canti(j+nBasT*(i-1))
          DPT2Canti(i+nBasT*(j-1)) =  Scal*Half
          DPT2Canti(j+nBasT*(i-1)) = -Scal*Half
        End Do
      End Do
C     write (*,*) "dpt2anti sym"
C     call sqprt(dpt2canti,nbast)
C     write (*,*) "dpt2c"
C     call sqprt(dpt2c_tot,nbast)
C
      Return
C
      End Subroutine CnstAntiC
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DENS1T_RPT2 (CI1,CI2,SGM1,G1,NLEV)

      use caspt2_global, only:iPrGlb
      use gugx, only: SGS, L2ACT, CIS
      use PrintLevel, only: debug
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: iSCF, nActEl, nAshT, STSym, Mul
      use pt2_guga, only: MxCI, nG1
      use Constants, only: Zero, One, Two

      IMPLICIT NONE

      real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI)
      real(kind=wp), intent(out) :: SGM1(MXCI), G1(NLEV,NLEV)
      integer(kind=iwp), intent(in):: nLev

      integer(kind=iwp), allocatable :: TASK(:,:)

      integer(kind=iwp) :: ID, IST, ISU, ISTU, IT, IU, LT, LU, ITASK,
     &  NTASKS, ISSG, NSGM
      real(kind=wp) :: GTU

      logical(kind=iwp), external :: RSV_TSK
      real(kind=wp), external :: ddot_, dnrm2_

c Purpose: Compute the 1- and 2-electron density matrix
c arrays G1 and G2.
      G1(:,:) = Zero

C For the special cases, there is no actual CI-routines involved:
c Special code for hi-spin case:
      IF(ISCF == 2) THEN
        DO IT = 1, NASHT
          G1(IT,IT) = One
        END DO
      else if (ISCF == 1 .and. NACTEL > 0) then
c Special code for closed-shell:
        DO IT = 1, NASHT
          G1(IT,IT) = Two
        END DO
      else

* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.F90

C-SVC20100311: set up a task table with LT,LU
        nTasks=(nLev**2+nLev)/2
        nTasks= nLev**2
        CALL mma_allocate (Task,nTasks,2,Label='TASK')

        iTask=0
        DO LT=1,nLev
          DO LU=1,nLev!LT
            iTask=iTask+1
            TASK(iTask,1)=LT
            TASK(iTask,2)=LU
          ENDDO
        ENDDO
        IF (iTask /= nTasks) WRITE(u6,*) "ERROR nTasks"

        Call Init_Tsk(ID, nTasks)

C-SVC20100311: BEGIN SEPARATE TASK EXECUTION
        do while (Rsv_Tsk(ID,iTask))
* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
C         LTU=0
C         DO 140 LT=1,NLEV
          LT=TASK(iTask,1)
          IST=SGS%ISM(LT)
          IT=L2ACT(LT)
C         DO 130 LU=1,LT
          LU=Task(iTask,2)
C         LTU=LTU+1
          ! LTU=iTask
          ISU=SGS%ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,STSYM)
          NSGM=CIS%NCSF(ISSG)
          IF(NSGM == 0) cycle
          CALL GETSGM2(LU,LT,STSYM,CI1,SGM1)
          IF(ISTU == 1) THEN
            GTU=DDOT_(NSGM,CI2,1,SGM1,1)
            G1(IT,IU)=G1(IT,IU)+GTU
          END IF

CSVC: The master node now continues to only handle task scheduling,
C     needed to achieve better load balancing. So it exits from the task
C     list.  It has to do it here since each process gets at least one
C     task.
        end do
        CALL Free_Tsk(ID)
        CALL mma_deallocate(Task)
        CALL GAdSUM (G1,NG1)
      end if

      IF(iPrGlb >= DEBUG) THEN
        WRITE(u6,'("DEBUG> ",A)')
     &   "DENS1_RPT2: norms of the density matrices:"
        WRITE(u6,'("DEBUG> ",A,1X,ES21.14)') "G1:", DNRM2_(NG1,G1,1)
      ENDIF

      end subroutine DENS1T_RPT2
C
C-----------------------------------------------------------------------
C
      Subroutine DWDER(OMGDER,HEFF,SLag)
C
      use definitions, only: wp, iwp
      use caspt2_module, only: NSTATE, DWTYPE, ZETA
      use Constants, only: Zero, Two, Half
C
      implicit none
C
      real(kind=wp), intent(in) :: OMGDER(nState,nState),
     &  HEFF(nState,nState)
      real(kind=wp), intent(inout) :: SLag(nState,nState)

      integer(kind=iwp) :: ilStat, jlStat, klStat
      real(kind=wp) :: Ebeta, Ealpha, Factor, Egamma, Dag, Hag, DEROMG,
     &  Scal, DERAB, DERAC, Dbg, Hbg
C
C     Computes the derivative of weight factor in XDW-CASPT2
C
      Do ilStat = 1, nState
        Ebeta = HEFF(ilStat,ilStat)
        Do jlStat = 1, nState
          Ealpha = HEFF(jlStat,jlStat)
C
          Factor = Zero
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWType == 1) Then
              Factor = Factor + exp(-zeta*(Ealpha - Egamma)**2)
            Else If (DWType == 2) Then
              Factor = Factor
     *          + exp(-zeta*(Ealpha/HEFF(jlStat,klStat))**2)
            Else If (DWType == 3) Then
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
          If (DWType == 1) Then
            DERAB = EXP(-ZETA*(Ealpha-Ebeta)**2)/Factor
            Scal = -Two*ZETA*DERAB*(Ealpha-Ebeta)*DEROMG
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
          Else If (DWType == 2) Then
            DERAB = EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
            DERAB = DERAB*Two*DEROMG*ZETA
     *            *(Ealpha/HEFF(jlStat,ilStat))**2
            Scal  = -DERAB/Ealpha
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            Scal  = +DERAB/HEFF(jlStat,ilStat)
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          Else If (DWType == 3) Then
            Dag = abs(Ealpha - Ebeta) + 1.0e-9_wp
            Hag = abs(HEFF(jlStat,ilStat))
            DERAB = EXP(-ZETA*Dag/(sqrt(Hag)+Tiny(Hag)))/Factor
            DERAB = -ZETA*DERAB*Dag/(sqrt(Hag)+tiny(Hag))
     *              *DEROMG
            Scal = DERAB/Dag
            If (Ealpha-Ebeta <= Zero) Scal = -Scal
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
            Scal  = -DERAB/(sqrt(Hag)+tiny(Hag))/sqrt(Hag)*Half
            If (HEFF(jlStat,ilStat) <= Zero) Scal = -Scal
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          End If
C
          !! derivative of alpha-gamma
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWtype == 1) Then
              DERAC = EXP(-ZETA*(Ealpha-Ebeta)**2)/(Factor*Factor)
     *               *EXP(-ZETA*(Ealpha-Egamma)**2)
              Scal = Two*ZETA*DERAC*(Ealpha-Egamma)*DEROMG
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
            Else If (DWType == 2) Then
              DERAC = EXP(-ZETA*(Ealpha/HEFF(jlStat,klStat))**2)/Factor
     *              * EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
              DERAC =-DERAC*Two*DEROMG*ZETA
     *              *(Ealpha/HEFF(jlStat,klStat))**2
              Scal  = -DERAC/Ealpha
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              Scal  = +DERAC/HEFF(jlStat,klStat)
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            Else IF (DWType == 3) Then
              Dbg = abs(Ealpha - Egamma) + 1.0e-9_wp
              Hbg = abs(HEFF(jlStat,klStat))
              DERAC =EXP(-ZETA*Dag/(sqrt(Hag)+tiny(Hag)))/Factor
     *              *EXP(-ZETA*Dbg/(sqrt(Hbg)+tiny(Hbg)))/Factor
              DERAC = ZETA*DERAC*Dbg/(sqrt(Hbg)+tiny(Hbg))
     *                *DEROMG
              Scal = DERAC/Dbg
              If (Ealpha-Egamma <= Zero) Scal = -Scal
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
              Scal  = -DERAC/(sqrt(Hbg)+tiny(Hbg))/sqrt(Hbg)*Half
              If (HEFF(jlStat,klStat) <= Zero) Scal = -Scal
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            End If
          End Do
        End Do
      End Do
C
      Return
C
      End Subroutine DWDER
