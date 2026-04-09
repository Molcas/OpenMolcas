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

      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV, only: IVECW, IVECC
      use stdalloc, only: mma_allocate,mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: STSYM, NCONF, NASHT, ISCF, NSTATE, JSTATE
      use pt2_guga, only: MXCI
      use Constants, only: Zero

      implicit none

      integer(kind=iwp) ::  IST, JST, I, NTG1, NTG2, NTG3, IDCI
      real(kind=wp) :: OVL, DUMMY(1)

      real(kind=wp), intent(inout) :: CLag(nConf,nState)
      real(kind=wp), intent(in) :: VECROT(*)

      real(kind=wp),allocatable :: DTG1(:),DTG2(:),DTG3(:),CI1(:),
     &  CI2(:),CI3(:)
!     return

! We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
! Note: Need proper allocation even if unused, sinced allocated
! arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      call mma_allocate(DTG1,NTG1,Label='DTG1')
      call mma_allocate(DTG2,NTG2,Label='DTG2')
      call mma_allocate(DTG3,NTG3,Label='DTG3')
      DTG1(:) = Zero
      DTG2(:) = Zero
      DTG3(:) = Zero

      !! OVL will contain the derivative contribution?
      !! It should be ignored
      OVL = Zero
      CALL DerHeffX(IVECW,IVECC,OVL,DTG1,DTG2,DTG3)

      call mma_allocate(CI1,MXCI,Label='MCCI1')
      call mma_allocate(CI2,MXCI,Label='MCCI2')
      call mma_allocate(CI3,MXCI,Label='MCCI3')

      IF(ISCF == 0) THEN
        IDCI=IDTCEX
        JST = jState
        CI1(1:NCONF) = Zero
        DO I=1,NSTATE
          IST = I
          IF (IST == JST) THEN
            CALL DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
          Else If (ABS(VECROT(IST)) <= 1.0e-12_wp) Then
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,2,CI3,NCONF,IDCI)
            CI1(1:NCONF) = CI1(1:NCONF) + VECROT(IST)*CI3(1:NCONF)
          END IF
        END DO

        CI3(1:NCONF) = Zero
        CALL DERTG3(.TRUE.,STSYM,STSYM,CI1,CI2,OVL,
     &              DTG1,DTG2,NTG3,DTG3,CI3,CLag(1,JST))

        DO I=1,NSTATE
          IST = I
          IF (IST == JST) THEN
            CYCLE
          Else If (ABS(VECROT(IST)) <= 1.0e-12_wp) Then
            CYCLE
          ELSE
            CLag(1:NCONF,IST) = CLag(1:NCONF,IST)
     &        + VECROT(IST)*CI3(1:NCONF)
          END IF
        END DO
      END IF

      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
      call mma_deallocate(CI3)

      call mma_deallocate(DTG1)
      call mma_deallocate(DTG2)
      call mma_deallocate(DTG3)

      End Subroutine DerHEff
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DerHeffX(IVEC,JVEC,OVL,DTG1,DTG2,DTG3)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NSYM, NASHT, NASUP, NISUP, NINDEP
      use definitions, only: wp, iwp, u6

      implicit none
! Compute the coupling Hamiltonian element defined as
!     HEL = < ROOT1 | H * OMEGA | ROOT2 >
! assuming that IVEC contains a contravariant representation of
! H|ROOT1>, JVEC contains a contravariant representation of
! OMEGA|ROOT2>, and OVL, TG1, TG2, TG3 contain the overlap (normally
! expected to be 0 or 1) and active transition density matrices of ROOT1
! and ROOT2. See also subroutine TSVEC for explanations.

! SVC (March 2014): modification of original code to handle distributed
! RHS arrays. There is now a main HCOUP subroutine that loops over cases
! and irreps and gets access to the process-specific block of the RHS.
! The coupling for that block is computed by the subroutine HCOUP_BLK.

      integer(kind=iwp), intent(in) :: IVEC, JVEC
      real(kind=wp), intent(out) :: OVL
      real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT),
     &  DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(*)

      integer(kind=iwp) :: ICASE, ISYM, NAS, NIN, NIS, lg_V1, lg_V2,
     &  iLo1, iHi1, jLo1, jHi1, MV1, iLo2, iHi2, jLo2, jHi2, MV2
! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif


! Sketch of procedure:
!  Loop over every (case/symmetry)-block.
!           If (No such vector block) Skip to end of loop
!           Allocate two places for this block, VEC1 and VEC2
!           Read VEC1 as IVEC component from file.
!           Read VEC2 as JVEC component from file.
!           Loop nest, computing
!              HEL := HEL + VEC1*GOM*VEC2
!           End of loop nest
!           Deallocate VEC1 and VEC2
!  End of loop.

      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIN=NINDEP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)

          IF(NAS*NIS == 0) cycle
          IF(NIN == 0) cycle

          CALL RHS_ALLO (NAS,NIS,lg_V1)
          CALL RHS_ALLO (NAS,NIS,lg_V2)
          CALL RHS_READ (NAS,NIS,lg_V1,ICASE,ISYM,IVEC)
          CALL RHS_READ (NAS,NIS,lg_V2,ICASE,ISYM,JVEC)
          CALL RHS_ACCESS(NAS,NIS,lg_V1,iLo1,iHi1,jLo1,jHi1,MV1)
          CALL RHS_ACCESS(NAS,NIS,lg_V2,iLo2,iHi2,jLo2,jHi2,MV2)

          IF ((iLo1 /= iLo2) .OR. (iHi1 /= iHi2) .OR.
     &        (jLo1 /= jLo2) .OR. (jHi1 /= jHi2)) THEN
            WRITE(u6,'(1X,A)') 'HCOUP: Error: block mismatch, abort...'
            CALL ABEND()
          END IF

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                      DBL_MB(MV1),DBL_MB(MV2),OVL,
     &                      DTG1,DTG2,DTG3)
          ELSE
#endif
            CALL DerHEffX_BLK(ICASE,ISYM,NAS,jLo1,jHi1,
     &                        GA_Arrays(MV1)%A,
     &                        GA_Arrays(MV2)%A,OVL,
     &                        DTG1,DTG2,DTG3)
#ifdef _MOLCAS_MPP_
          END IF
#endif

          CALL RHS_RELEASE (lg_V1,iLo1,iHi1,jLo1,jHi1)
          CALL RHS_RELEASE (lg_V2,iLo2,iHi2,jLo2,jHi2)
          CALL RHS_FREE (lg_V1)
          CALL RHS_FREE (lg_V2)
        END DO
      END DO

      end subroutine DerHeffX
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DerHEffX_BLK(ICASE,ISYM,NAS,IISTA,IIEND,V1,V2,OVL,
     &                        DTG1,DTG2,DTG3)

      USE SUPERINDEX, only: MTU, MTUV, MTGEU, MTGTU
      use caspt2_module, only: NAES, NASHT, NTUVES, NTUES, NTGEUES,
     &                         NTGTUES
      use definitions, only: wp, iwp
      use Constants, only: Zero, Two, Four, Eight
! Compute a contribution to the coupling Hamiltonian element (HEL)
! defined as HEL = < ROOT1 | H * OMEGA | ROOT2 >. The contribution
! arises from the block V_(A,I), with A=1,NAS and I=IISTA,IIEND,
! with A the active superindex and I the inactive superindex. Since
! the inactive superindex is partitioned over processes, each process
! only computes part of the HEL value, which is then sum reduced in the
! calling subroutine.
      implicit none

      integer(kind=iwp), intent(in) :: ICASE, ISYM, NAS, IISTA, IIEND
      real(kind=wp), intent(in) :: V1(*), V2(*)
      real(kind=wp), intent(out) :: OVL
      real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT),
     &  DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(*)
! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

      integer(kind=iwp) :: NISBLK, IAS, IASABS, ITABS, IUABS, IVABS,
     &  JAS, JASABS, IXABS, IYABS, IZABS, IND1, IND2, IND3, JND1,
     &  JND2, JND3, ITG3, NAS1, IAS1, IAS2, JAS1, JAS2
      real(kind=wp) :: VAL, VAL11, VAL12, VAL21, VAL22
      real(kind=wp), external :: ddot_

      IF (IISTA <= 0) RETURN

      OVL = Zero
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
! Compute and use SA(ITABS IUABS IVABS, IXABS IYABS IZABS)
! Formulae used:
!  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz -
!         - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
! Gvuxtyz is stored using full permutation symmetry of three pairs
! (vu),(xt), and (yz):
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
!  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz -
!         - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
! Compute TMP=Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            IF(IXABS == ITABS) THEN
              DTG2(IVABS,IUABS,IYABS,IZABS)
     &          = DTG2(IVABS,IUABS,IYABS,IZABS) + Two*VAL
              IF(IYABS == IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + Two*VAL
              END IF
            END IF
            VAL = -VAL
            DTG3(ITG3) = DTG3(ITG3) + VAL
            IF(IYABS == IUABS) THEN
              DTG2(IVABS,IZABS,IXABS,ITABS)
     &          = DTG2(IVABS,IZABS,IXABS,ITABS) + VAL
            END IF
            IF(IYABS == ITABS) THEN
              DTG2(IVABS,IUABS,IXABS,IZABS)
     &          = DTG2(IVABS,IUABS,IXABS,IZABS) + VAL
              IF(IXABS == IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + VAL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              DTG2(IVABS,ITABS,IYABS,IZABS)
     &          = DTG2(IVABS,ITABS,IYABS,IZABS) + VAL
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
! Compute and use SC(IXABS IUABS IVABS, ITABS IYABS IZABS)
! In SBMAT, the formula is written as SC(tuv,xyz)
!    = Gvutxyz +dyu Gvztx + dyx Gvutz + dtu Gvxyz + dtu dyx Gvz
! Rewritten, in order to reuse same quantities as in SA:
!  SC(xuv,tyz)
!    = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
! Gvuxtyz is stored using full permutation symmetry of three pairs
! (vu),(xt), and (yz):
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
!  SC(xuv,tyz) (rewritten, swapping x and t)
!    = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG3(ITG3) = DTG3(ITG3) + VAL
            IF(IYABS == IUABS) THEN
              DTG2(IVABS,IZABS,IXABS,ITABS)
     &          = DTG2(IVABS,IZABS,IXABS,ITABS) + VAL
            END IF
            IF(IYABS == ITABS) THEN
              DTG2(IVABS,IUABS,IXABS,IZABS)
     &          = DTG2(IVABS,IUABS,IXABS,IZABS) + VAL
              IF(IXABS == IUABS) THEN
                DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS) + VAL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              DTG2(IVABS,ITABS,IYABS,IZABS)
     &          = DTG2(IVABS,ITABS,IYABS,IZABS) + VAL
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
! Formulae used:
!    SB(tu,xy)=
!    = 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu
!      -4dxu dyt + 2dxu Gyt
!    SB(tu,yx)=
!    = 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu
!      -4dyu dxt + 2dyu Gxt
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(IXABS,ITABS,IYABS,IUABS)
     &        = DTG2(IXABS,ITABS,IYABS,IUABS) + Two*VAL
            DTG2(IYABS,ITABS,IXABS,IUABS)
     &        = DTG2(IYABS,ITABS,IXABS,IUABS) + Two*VAL
            IF(IXABS == ITABS) THEN
              DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS) - Four*VAL + Two*VAL
              IF(IYABS == IUABS) THEN
                OVL = OVL + Eight*VAL - Four*VAL
              END IF
            END IF
            IF(IYABS == IUABS) THEN
              DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - Four*VAL + Two*VAL
            END IF
            IF(IYABS == ITABS) THEN
              DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS) + Two*VAL - Four*VAL
              IF(IXABS == IUABS) THEN
                OVL = OVL - Four*VAL + Eight*VAL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS) + Two*VAL - Four*VAL
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
! Formulae used:
!    SB(tu,xy)=
!    = 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu
!      -4dxu dyt + 2dxu Gyt
!    SB(tu,yx)=
!    = 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu
!      -4dyu dxt + 2dyu Gxt
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(IXABS,ITABS,IYABS,IUABS)
     &        = DTG2(IXABS,ITABS,IYABS,IUABS) + Two*VAL
            DTG2(IYABS,ITABS,IXABS,IUABS)
     &        = DTG2(IYABS,ITABS,IXABS,IUABS) - Two*VAL
            IF(IXABS == ITABS) THEN
              DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS) - Four*VAL - Two*VAL
              IF(IYABS == IUABS) THEN
                OVL = OVL + Eight*VAL + Four*VAL
              END IF
            END IF
            IF(IYABS == IUABS) THEN
              DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - Four*VAL - Two*VAL
            END IF
            IF(IYABS == ITABS) THEN
              DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS) + Two*VAL + Four*VAL
              IF(IXABS == IUABS) THEN
                OVL = OVL - Four*VAL - Eight*VAL
              END IF
            END IF
            IF(IXABS == IUABS) THEN
              DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS) + Two*VAL + Four*VAL
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
! Formulae used:
!    SD11(tu1,xy1)=2*(Gutxy + dtx Guy)
!    SD12(tu2,xy1)= -(Gutxy + dtx Guy)
!    SD21(tu2,xy1)= -(Gutxy + dtx Guy)
!    SD22(tu2,xy2)= -Gxtuy +2*dtx Guy
            VAL11 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS1),NAS)
            VAL12 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS2),NAS)
            VAL21 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS1),NAS)
            VAL22 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS2),NAS)
            DTG2(IUABS,ITABS,IXABS,IYABS)
     &        = DTG2(IUABS,ITABS,IXABS,IYABS)
     &        + Two*VAL11 - VAL12 - VAL21
            DTG2(IXABS,ITABS,IUABS,IYABS)
     &        = DTG2(IXABS,ITABS,IUABS,IYABS) - VAL22
            IF(ITABS == IXABS) THEN
              DTG1(IUABS,IYABS) = DTG1(IUABS,IYABS)
     &          + Two*VAL11 - VAL12 - VAL21 + Two*VAL22
            END IF
          END DO
        END DO
************************************************************************
      CASE(6)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SE(t,x)=2*dxt - Dxt
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - VAL
            IF(IXABS == ITABS) OVL=OVL+Two*VAL
          END DO
        END DO
************************************************************************
      CASE(7)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SE(t,x)=2*dxt - Dxt
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS) - VAL
            IF(IXABS == ITABS) OVL=OVL+Two*VAL
          END DO
        END DO
************************************************************************
      CASE(8)
! ========================================================
! Compute and use SFP(ITABS IUABS , IXABS IYABS)
! and (later, similar) SFM(ITABS IUABS , IXABS IYABS)
        DO IAS=1,NAS
          IASABS=NTGEUES(ISYM)+IAS
          ITABS=MTGEU(1,IASABS)
          IUABS=MTGEU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGEUES(ISYM)+JAS
            IXABS=MTGEU(1,JASABS)
            IYABS=MTGEU(2,JASABS)
! Formulae used:
!    SF(tu,xy)= 2 Gtxuy
!    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
!    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(ITABS,IXABS,IUABS,IYABS)
     &        = DTG2(ITABS,IXABS,IUABS,IYABS) + Two*VAL
            DTG2(ITABS,IYABS,IUABS,IXABS)
     &        = DTG2(ITABS,IYABS,IUABS,IXABS) + Two*VAL
          END DO
        END DO
************************************************************************
      CASE(9)
! ========================================================
! Compute and use SFM(ITABS IUABS, IXABS ,IYABS)
        DO IAS=1,NAS
          IASABS=NTGTUES(ISYM)+IAS
          ITABS=MTGTU(1,IASABS)
          IUABS=MTGTU(2,IASABS)
          DO JAS=1,NAS
            JASABS=NTGTUES(ISYM)+JAS
            IXABS=MTGTU(1,JASABS)
            IYABS=MTGTU(2,JASABS)
! Formulae used:
!    SF(tu,xy)= 4 Ptxuy
!    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
!    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
            VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
            DTG2(ITABS,IXABS,IUABS,IYABS)
     &        = DTG2(ITABS,IXABS,IUABS,IYABS) + Two*VAL
            DTG2(ITABS,IYABS,IUABS,IXABS)
     &        = DTG2(ITABS,IYABS,IUABS,IXABS) - Two*VAL
          END DO
        END DO
************************************************************************
! CASES GP, GM
! Compute and use SG(ITABS , IXABS) (Same for cases GP and GM)
************************************************************************
      CASE(10)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SG(t,x)= Gtx
!           SG= TG1(ITABS,IXABS)
            DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)
     &        + DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(11)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SG(t,x)= Gtx
!           SG= TG1(ITABS,IXABS)
            DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)
     &        + DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
************************************************************************
      CASE(12)
        OVL = OVL + DDOT_(NAS*NISBLK,V2,1,V1,1)
!       IF(ABS(OVL) >= 1.0e-12_wp) THEN
!         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
!       END IF
************************************************************************
      CASE(13)
        OVL = OVL + DDOT_(NAS*NISBLK,V2,1,V1,1)
!       IF(ABS(OVL) >= 1.0e-12_wp) THEN
!         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
!       END IF
************************************************************************
      END SELECT
      Return
      end subroutine DerHEffX_BLK
!
!-----------------------------------------------------------------------
!
      SUBROUTINE DERTG3(DOG3,LSYM1,LSYM2,CI1,CI2,OVL,DTG1,DTG2,NTG3,
     &                  DTG3,CLAG1,CLAG2)
      use Symmetry_Info, only: Mul
      use gugx, only: SGS, L2ACT, CIS, EXS
      use stdalloc, only: mma_MaxDBLE, mma_allocate, mma_deallocate
      use definitions, only: iwp,wp,u6
      use caspt2_module, only: NACTEL, NASHT, ISCF, IASYM
      use pt2_guga, only: MXCI
      use Constants, only: Zero, One

      implicit none

      logical(kind=iwp), intent(in) :: DOG3
      integer(kind=iwp), intent(in) :: LSYM1, LSYM2, NTG3
      real(kind=wp), intent(in) :: CI1(MXCI), CI2(MXCI), OVL
      real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT),
     &  DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(NTG3), CLAG1(MXCI),
     &  CLAG2(MXCI)

      real(kind=wp), allocatable :: TG3WRK(:),BUF1(:),DTU(:,:),DYZ(:,:)
      integer(kind=iwp), allocatable :: P2LEV(:)
      integer(kind=iwp) :: nLev, LP2LEV1, LP2LEV2, IP, IL, JL, IP1,
     &  IT, IU, ITU, ITS, IUS, IS1, IP2, IV, IX, IVX, IVS, IXS, IS2,
     &  IP3, IY, IZ, IYS, IZS, IS3, IYZ, JTU, JVX, JYZ, JTUVXYZ, NCI1,
     &  NVECS, NTG3WRK, NYZBUF, NTUBUF, LSGM1, LTAU, LSGM2, IP3STA,
     &  IP3END, LTO, ISSG2, IP1STA, IP1END, ISSG1, LFROM, LFROMD, IM,
     &  JM, ISTAU, NTAU, ibuf
      real(kind=wp) :: VAL

      nLev = SGS%nLev
! Procedure for computing 1-body, 2-body, and 3-body transition
! density elements with active indices only.

! In: Wave functions CI1, with symmetry LSYM1, and CI2, with
!  symmetry LSYM2.
!
! Out: Transition density matrices, denoted here TG1, TG2 and TG3.
! Storage: TG1 and TG2 are simple two- and four-index arrays, and
! includes also such zeroes that are implied by symmetry.
! But TG3 is quite large, and while it is stored with zeroes, it
! is made more compact by the following addressing:

! <Psi1|E_tuvxyz|Psi2> is stored in TG3(ITG3) where
!    ITG3= ((i+1)*i*(i-1))/6 + (j*(j-1))/2 + k
!     i  = max(tu,vx,yz)
!     j  = mid(tu,vx,yz)
!     k  = min(tu,vx,yz)
! tu stands for the pair index tu= t + NASHT*(u-1), etc., and t is
! the usual active orbital number, when they are enumerated across
! all the symmetries (The ''absolute'' active index).


! Put in zeroes. Recognize special cases:
!     OVL=1.0D0
      IF(NASHT == 0) return
      IF(NACTEL == 0) return
!     IF(LSYM1 /= LSYM2) OVL=Zero
!     TG1(:,:) = Zero
!     TG2(:,:,:,:) = Zero
!     TG3(1:NTG3) = Zero

      IF(ISCF /= 0) then
! -Special code for the closed-shell or hi-spin cases:
! ISCF=1 for closed-shell, =2 for hispin
        write (u6,*) 'Here is the special case'
        write (u6,*) 'not yet'
        call abend()
      end if
! Here, for regular CAS or RAS cases.

! Special pair index allows true RAS cases to be handled:
      CALL mma_allocate(P2LEV,2*NASHT**2,Label='P2LEV')
      LP2LEV1=1
      LP2LEV2=1+NASHT**2
      IP=0
! First, IL < JL pairs.
      DO IL=1,NLEV-1
       DO JL=IL+1,NLEV
        IP=IP+1
        P2LEV(LP2LEV1-1+IP)=IL
        P2LEV(LP2LEV2-1+IP)=JL
       END DO
      END DO
! Then, IL = JL pairs.
      DO IL=1,NLEV
        IP=IP+1
        P2LEV(LP2LEV1-1+IP)=IL
        P2LEV(LP2LEV2-1+IP)=IL
      END DO
! Last, IL > JL pairs.
      DO IL=2,NLEV
       DO JL=1,IL-1
        IP=IP+1
        P2LEV(LP2LEV1-1+IP)=IL
        P2LEV(LP2LEV2-1+IP)=JL
       END DO
      END DO

! First, the 3-particle density matrix:
! <PSI1|E(T,U,V,X,Y,Z)|PSI2>  = <PSI1|E(TU)E(VX)E(YZ)|PSI2>
! -D(Y,X)*(TG2(T,U,V,Z)+D(V,U)*TG1(T,Z))
! -D(V,U)*TG2(T,X,Y,Z) C -D(Y,U)*TG2(V,X,T,Z)
      IF (DOG3) THEN
       DO IP1=1,NASHT**2
        IT=L2ACT(P2LEV(LP2LEV1-1+IP1))
        IU=L2ACT(P2LEV(LP2LEV2-1+IP1))
        ITU=IT+NASHT*(IU-1)
        ITS=IASYM(IT)
        IUS=IASYM(IU)
        IS1=Mul(Mul(ITS,IUS),LSYM1)
        DO IP2=1,IP1
         IV=L2ACT(P2LEV(LP2LEV1-1+IP2))
         IX=L2ACT(P2LEV(LP2LEV2-1+IP2))
         IVX=IV+NASHT*(IX-1)
         IVS=IASYM(IV)
         IXS=IASYM(IX)
         IS2=Mul(Mul(IVS,IXS),IS1)
         DO IP3=1,IP2
          IY=L2ACT(P2LEV(LP2LEV1-1+IP3))
          IZ=L2ACT(P2LEV(LP2LEV2-1+IP3))
          IYS=IASYM(IY)
          IZS=IASYM(IZ)
          IS3=Mul(Mul(IYS,IZS),IS2)
          IF(IS3 == LSYM2) THEN
           IYZ=IY+NASHT*(IZ-1)
           IF(ITU < IVX) THEN
             IF(ITU >= IYZ) THEN
               JTU=IVX
               JVX=ITU
               JYZ=IYZ
             ELSE IF(IVX < IYZ) THEN
                 JTU=IYZ
                 JVX=IVX
                 JYZ=ITU
             ELSE
                 JTU=IVX
                 JVX=IYZ
                 JYZ=ITU
             END IF
           ELSE
             IF(ITU < IYZ) THEN
               JTU=IYZ
               JVX=ITU
               JYZ=IVX
             ELSE IF (IVX >= IYZ) THEN
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
           IF(IY == IX) THEN
            DTG2(IT,IU,IV,IZ)=DTG2(IT,IU,IV,IZ)-VAL
            IF(IV == IU) THEN
             DTG1(IT,IZ)=DTG1(IT,IZ)-VAL
            END IF
           END IF
           IF(IV == IU) THEN
            DTG2(IT,IX,IY,IZ)=DTG2(IT,IX,IY,IZ)-VAL
           END IF
           IF(IY == IU) THEN
            DTG2(IV,IX,IT,IZ)=DTG2(IV,IX,IT,IZ)-VAL
           END IF
          END IF
         END DO
        END DO
       END DO
      END IF
!
! Then, the 2-particle density matrix:
! <PSI1|E(T,U,V,X)|PSI2>  = <PSI1|E(TU)E(VX)|PSI2> - D(V,U)*TG2(T,U,V,X)
      DO IP1=1,NASHT**2
       IT=L2ACT(P2LEV(LP2LEV1-1+IP1))
       IU=L2ACT(P2LEV(LP2LEV2-1+IP1))
       DO IP2=1,IP1
        IV=L2ACT(P2LEV(LP2LEV1-1+IP2))
        IX=L2ACT(P2LEV(LP2LEV2-1+IP2))
        IF (IP1 /= IP2) Then
          DTG2(IT,IU,IV,IX)=DTG2(IT,IU,IV,IX)+DTG2(IV,IX,IT,IU)
          DTG2(IV,IX,IT,IU) = Zero
        End If
        IF(IV == IU) DTG1(IT,IX)=DTG1(IT,IX)-DTG2(IT,IU,IV,IX)
       END DO
      END DO
!
! If now any matrix element E(t1u1)E(t2u2)..E(tnun) is arranged
! such that the pair indices are non-decreasing, then the matrix
! element can be correctly computed by performing explicit
! excitations within the RAS space.
! But we also need the 'usual' pair index in order to use the
! packed addressing.

      NCI1=CIS%NCSF(LSYM1)
! Overlap:
!     IF(LSYM1 == LSYM2) OVL=DDOT_(NCI1,CI1,1,CI2,1)
      IF(LSYM1 == LSYM2) THEN
        CLAG2(1:NCI1) = CLAG2(1:NCI1) + OVL*CI1(1:NCI1)
        CLAG1(1:NCI1) = CLAG1(1:NCI1) + OVL*CI2(1:NCI1)
      END IF
!     write (*,*) 'overlap = ',DDOT_(NCI1,CI1,1,CI2,1)
! Allocate as many vectors as possible:
! Wishful thinking:
      NVECS=2*NASHT**2+1
! But what is really available?
      CALL mma_MaxDBLE(NTG3WRK)
      NTG3WRK=NTG3WRK-3*MXCI ! for BUF1, DTU, and DYZ, allocated later
      NTG3WRK=NTG3WRK/2
      NTG3WRK=MIN(MXCI*NVECS,NTG3WRK)
      NVECS=NTG3WRK/MXCI
      NTG3WRK=NVECS*MXCI
! Find optimal subdivision of available vectors:
      NYZBUF=NINT(real(NVECS-1,kind=wp)/real(NASHT,kind=wp),kind=iwp)
      NYZBUF=MAX(1,NYZBUF)
      NTUBUF=MIN(NASHT**2,NVECS-1-NYZBUF)
      NYZBUF=NVECS-1-NTUBUF
! Insufficient memory?
      IF(NTUBUF <= 0) THEN
        WRITE(u6,*)' Too little memory left for MKTG3.'
        WRITE(u6,*)' Need at least 6 vectors of length MXCI=',MXCI
        CALL ABEND()
      END IF
      IF(NTUBUF <= (NASHT**2)/5) THEN
        WRITE(u6,*)' WARNING: MKTG3 will be inefficient owing to'
        WRITE(u6,*)' small memory.'
      END IF
      CALL mma_allocate(TG3WRK,NTG3WRK,Label='TG3WRK')
      CALL mma_allocate(BUF1,MXCI,Label='BUF1')

      call mma_allocate(DTU,MXCI,NTUBUF,Label='DTU')
      call mma_allocate(DYZ,MXCI,NYZBUF,Label='DYZ')

! And divide it up:
      !! LSGM1: NTUBUF vectors
      !! LTAU : 1 vector
      !! LSGM2: NYZBUF vectors
      LSGM1=1
      LTAU=LSGM1+NTUBUF*MXCI
      LSGM2=LTAU+MXCI

! Sectioning loops over pair indices IP3 (ket side):
      DO IP3STA=1,NASHT**2,NYZBUF
       IP3END=MIN(NASHT**2,IP3STA-1+NYZBUF)
! Compute a section of sigma vectors E(YZ)*PSI2 to memory:
       LTO=LSGM2
       DO IP3=IP3STA,IP3END
! Translate to levels in the SGUGA coupling order:
        IL=P2LEV(LP2LEV1-1+IP3)
        JL=P2LEV(LP2LEV2-1+IP3)
        IY=L2ACT(IL)
        IZ=L2ACT(JL)
        IYS=IASYM(IY)
        IZS=IASYM(IZ)
        ISSG2=Mul(Mul(IYS,IZS),LSYM2)
        TG3WRK(LTO:LTO+MXCI-1) = Zero
! LTO is first element of Sigma2 = E(YZ) Psi2
        CALL SIGMA1(SGS,CIS,EXS,
     &              IL,JL,One,LSYM2,CI2,TG3WRK(LTO))
        IF(ISSG2 == LSYM1 .AND. DTG1(IY,IZ) /= Zero) THEN
          !! It is possible to calculate the contribution using
          !! DGEMV, but DAXPY seems to be faster than DGEMV
          CLAG1(1:NCI1) = CLAG1(1:NCI1)
     &      + DTG1(IY,IZ)*TG3WRK(LTO:LTO+NCI1-1)
        END IF
        LTO=LTO+MXCI
       END DO

       DYZ(1:MXCI,1:NYZBUF) = Zero
! Sectioning loops over pair indices IP1 (bra side):
       DO IP1STA=IP3STA,NASHT**2,NTUBUF
        IP1END=MIN(NASHT**2,IP1STA-1+NTUBUF)
! Compute a section of sigma vectors E(UT)*PSI1 to memory:
        LTO=LSGM1
        !! <Psi1|E(TU)
        DO IP1=IP1STA,IP1END
! Translate to levels:
         JL=P2LEV(LP2LEV1-1+IP1)
         IL=P2LEV(LP2LEV2-1+IP1)
         IT=L2ACT(IL)
         IU=L2ACT(JL)
         ITS=IASYM(IT)
         IUS=IASYM(IU)
         ISSG1=Mul(Mul(ITS,IUS),LSYM1)
         TG3WRK(LTO:LTO+MXCI-1) = Zero
         CALL SIGMA1(SGS,CIS,EXS,
     &               IL,JL,One,LSYM1,CI1,TG3WRK(LTO))
         IF (ISSG1 == LSYM1 .AND. DTG1(IU,IT) /= Zero
     &       .AND. IP3STA == 1) THEN
          CLAG2(1:NCI1) = CLAG2(1:NCI1)
     &      + DTG1(IU,IT)*TG3WRK(LTO:LTO+NCI1-1)
         END IF
         LTO=LTO+MXCI
        END DO

        DTU(1:MXCI,1:NTUBUF) = Zero
! Now compute as many elements as possible:
        LFROM=LSGM2
        LFROMD=1
        DO IP3=IP3STA,IP3END
         IY=L2ACT(P2LEV(LP2LEV1-1+IP3))
         IZ=L2ACT(P2LEV(LP2LEV2-1+IP3))
! LFROM will be start element of Sigma2=E(YZ) Psi2
         IYZ=IY+NASHT*(IZ-1)
         IYS=IASYM(IY)
         IZS=IASYM(IZ)
         ISSG2=Mul(Mul(IYS,IZS),LSYM2)
         IM=P2LEV(LP2LEV1-1+IP3)
         JM=P2LEV(LP2LEV2-1+IP3)
         DO IP2=IP3,IP1END
          IL=P2LEV(LP2LEV1-1+IP2)
          JL=P2LEV(LP2LEV2-1+IP2)
          IV=L2ACT(IL)
          IX=L2ACT(JL)
          IVX=IV+NASHT*(IX-1)
          IVS=IASYM(IV)
          IXS=IASYM(IX)
          ISTAU=Mul(Mul(IVS,IXS),ISSG2)
          NTAU=CIS%NCSF(ISTAU)
          TG3WRK(LTAU:LTAU+MXCI-1) = Zero
! LTAU  will be start element of Tau=E(VX) Sigma2=E(VX) E(YZ) Psi2
          !! LTAU = EvxEyz|Psi2>
          CALL SIGMA1(SGS,CIS,EXS,
     &                IL,JL,One,ISSG2,TG3WRK(LFROM),TG3WRK(LTAU))
          IF(ISTAU == LSYM1 .AND. DTG2(IV,IX,IY,IZ) /= Zero) THEN
!          DTG2(IV,IX,IY,IZ)=DDOT_(NTAU,TG3WRK(LTAU),1,CI1,1)
           !! For left derivative: <I|Evx Eyz|Psi2>
           CLAG1(1:NTAU) = CLAG1(1:NTAU)
     &       + DTG2(IV,IX,IY,IZ)*TG3WRK(LTAU:LTAU+NTAU-1)
           !! For right derivative: <Psi1|Evx Eyz|I>
           IF (IP2 >= IP1STA.AND.IP2 <= IP1END) THEN
              ibuf = lsgm1+mxci*(ip2-ip1sta)
              DYZ(1:MXCI,LFROMD) = DYZ(1:MXCI,LFROMD)
     &          + DTG2(IV,IX,IY,IZ)*TG3WRK(IBUF:IBUF+MXCI-1)
           ELSE
         CALL SIGMA1(SGS,CIS,EXS,
     &               JL,IL,DTG2(IV,IX,IY,IZ),ISSG2,CI1,DYZ(1,LFROMD))
           END IF
           DTG2(IV,IX,IY,IZ) = Zero
          END IF
          IF (DOG3) THEN
            BUF1(1:MXCI) = Zero
            DO IP1=MAX(IP2,IP1STA),IP1END
             IT=L2ACT(P2LEV(LP2LEV1-1+IP1))
             IU=L2ACT(P2LEV(LP2LEV2-1+IP1))
             ITS=IASYM(IT)
             IUS=IASYM(IU)
             ISSG1=Mul(Mul(ITS,IUS),LSYM1)
             IF(ISSG1 == ISTAU) THEN
!             L=LSGM1+MXCI*(IP1-IP1STA)
!             VAL=DDOT_(NTAU,TG3WRK(LTAU),1,TG3WRK(L),1)
              ITU=IT+NASHT*(IU-1)
! Here VAL is the value <PSI1|E(IT1,IU1)E(IT2,IU2)E(IT3,IU3)|PSI2>
! Code to put it in correct place:
              IF(ITU < IVX) THEN
                IF(ITU >= IYZ) THEN
                  JTU=IVX
                  JVX=ITU
                  JYZ=IYZ
                ELSE IF(IVX < IYZ) THEN
                  JTU=IYZ
                  JVX=IVX
                  JYZ=ITU
                ELSE
                  JTU=IVX
                  JVX=IYZ
                  JYZ=ITU
                END IF
              ELSE
                IF(ITU < IYZ) THEN
                  JTU=IYZ
                  JVX=ITU
                  JYZ=IVX
                ELSE IF (IVX >= IYZ) THEN
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
              IF (DTG3(JTUVXYZ) /= Zero) then
                !! For left derivative: <I|Evx Eyz|Psi2> * Dtuvxyz
                !! I don't understand, but this is much faster than
                !! processing all possible vectors at once with DGEMM
                !! (and DGER) after finishing the IP1 loop below
                DTU(1:MXCI,1+IP1-IP1STA) = DTU(1:MXCI,1+IP1-IP1STA)
     &            + DTG3(JTUVXYZ)*TG3WRK(LTAU:LTAU+MXCI-1)
                !! For right derivative: <Psi1|Etu|I> * Dtuvxyz
                !! This is also (slightly) faster than DGEMV, apparently
                Call DaXpY_(MXCI,DTG3(JTUVXYZ),
     &                      TG3WRK(LSGM1+MXCI*(IP1-IP1STA)),1,BUF1,1)
               END IF
! End of symmetry requirement IF-clause:
             END IF
! End of IP1 loop.
            END DO
            !! Second operator for the right derivative:
            !! <Psi1|Etu Evx|I> * Dtuvxyz
            CALL SIGMA1(SGS,CIS,EXS,
     &                  JL,IL,One,ISTAU,BUF1,DYZ(1,LFROMD))
          END IF !! End of DOG3 clause
! End of IP2 loop.
         END DO
         LFROM=LFROM+MXCI
         LFROMD=LFROMD+1
! End of IP3 loop.
        END DO

        LTO=1
        !! <I|Etu Evx Eyz|Psi2> * Dtuvxyz
        DO IP1=IP1STA,IP1END
! Translate to levels:
         IL=P2LEV(LP2LEV1-1+IP1)
         JL=P2LEV(LP2LEV2-1+IP1)
         IT=L2ACT(IL)
         IU=L2ACT(JL)
         ITS=IASYM(IT)
         IUS=IASYM(IU)
         ISSG1=Mul(Mul(ITS,IUS),LSYM1)
         CALL SIGMA1(SGS,CIS,EXS,IL,JL,One,LSYM1,DTU(1,LTO),CLAG1)
         LTO=LTO+1
        END DO
! End of IP1STA sectioning loop
       END DO

       LTO=1
       DO IP3=IP3STA,IP3END
        IY=L2ACT(P2LEV(LP2LEV1-1+IP3))
        IZ=L2ACT(P2LEV(LP2LEV2-1+IP3))
! LFROM will be start element of Sigma2=E(YZ) Psi2
        IYZ=IY+NASHT*(IZ-1)
        IYS=IASYM(IY)
        IZS=IASYM(IZ)
        ISSG2=Mul(Mul(IYS,IZS),LSYM2)
        IM=P2LEV(LP2LEV1-1+IP3)
        JM=P2LEV(LP2LEV2-1+IP3)
! LTO is first element of Sigma2 = E(YZ) Psi2
        CALL SIGMA1(SGS,CIS,EXS,JM,IM,One,LSYM2,DYZ(1,LTO),CLAG2)
        LTO=LTO+1
       END DO
! End of IP3STA sectioning loop
      END DO

      call mma_deallocate(TG3WRK)
      call mma_deallocate(BUF1)
      call mma_deallocate(DTU)
      call mma_deallocate(DYZ)
      call mma_deallocate(P2LEV)

      end subroutine DERTG3
