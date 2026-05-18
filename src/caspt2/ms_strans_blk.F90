!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      SUBROUTINE MS_STRANS_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,IISTA,   &
     &                         IIEND,V1,V2,OVL,TG1,TG2,TG3,SCAL)
      USE SUPERINDEX, only: MTU, MTUV, MTGEU, MTGTU
      use caspt2_module, only: NAES, NTUVES, NTUES, NTGEUES, NTGTUES
      use definitions, only: wp, iwp
      use Constants, only: Two, Four, Eight
! Compute a contribution to the coupling Hamiltonian element (HEL)
! defined as HEL = < ROOT1 | H * OMEGA | ROOT2 >. The contribution
! arises from the block V_(A,I), with A=1,NAS and I=IISTA,IIEND,
! with A the active superindex and I the inactive superindex. Since
! the inactive superindex is partitioned over processes, each process
! only computes part of the HEL value, which is then sum reduced in the
! calling subroutine.
      implicit none

      integer(kind=iwp), intent(in) :: ICASE, ISYM, NAS, nvlen, NASHT,  &
     &                                 NTG3, IISTA, IIEND
      real(kind=wp), intent(in) :: V1(nvlen), OVL, SCAL
      real(kind=wp), intent(inout) :: V2(nvlen), TG1(NASHT,NASHT),      &
     &  TG2(NASHT,NASHT,NASHT,NASHT), TG3(NTG3)
! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

      integer(kind=iwp) :: NISBLK, IAS, IASABS, ITABS, IUABS, IVABS,    &
     &  JAS, JASABS, IXABS, IYABS, IZABS, IND1, IND2, IND3, JND1,       &
     &  JND2, JND3, ITG3, NAS1, IAS1, IAS2, JAS1, JAS2
      real(kind=wp) :: TMP, SA, SC, SBtuxy, SBtuyx, SBP, SBM, GUTXY,    &
     &  SD11, SD12, SD21, SD22, GUY, SE, SFtuxy, SFtuyx, SFP, SFM, SG

      IF (IISTA <= 0) RETURN

      NISBLK=IIEND-IISTA+1
      SELECT CASE (ICASE)
!***********************************************************************
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
! SA is the negative of this, and then some correction:
            SA=-TMP
            IF(IXABS == ITABS) THEN
              SA=SA+Two*TG2(IVABS,IUABS,IYABS,IZABS)
              IF(IYABS == IUABS) THEN
                SA=SA+Two*TG1(IVABS,IZABS)
              END IF
            END IF
! SA has been computed.

!           HEBLK=HEBLK+SA*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SA*SCAL,V1(IAS),NAS,V2(JAS),NAS)
!           Do i = 1, NAS
!             V2(JAS+i-1) = SA*V1(iAS-1)
!           End Do
          END DO
        END DO
!***********************************************************************
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

!           HEBLK=HEBLK+SC*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SC*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
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
!    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
!    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)
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

!           HEBLK=HEBLK+SBP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SBP*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
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
!    SBP(tu,xy)=SB(tu,xy)+SB(tu,yx)
!    SBM(tu,xy)=SB(tu,xy)-SB(tu,yx)
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

!           HEBLK=HEBLK+SBM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SBM*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
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

!           HEBLK=HEBLK+SD11*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS1),NAS)
!           HEBLK=HEBLK+SD12*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS1),NAS)
!           HEBLK=HEBLK+SD21*DDOT_(NISBLK,V2(JAS1),NAS,V1(IAS2),NAS)
!           HEBLK=HEBLK+SD22*DDOT_(NISBLK,V2(JAS2),NAS,V1(IAS2),NAS)
            Call DaXpY_(NISBLK,SD11*SCAL,V1(IAS1),NAS,V2(JAS1),NAS)
            Call DaXpY_(NISBLK,SD12*SCAL,V1(IAS1),NAS,V2(JAS2),NAS)
            Call DaXpY_(NISBLK,SD21*SCAL,V1(IAS2),NAS,V2(JAS1),NAS)
            Call DaXpY_(NISBLK,SD22*SCAL,V1(IAS2),NAS,V2(JAS2),NAS)
          END DO
        END DO
!***********************************************************************
      CASE(6)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SE(t,x)=2*dxt - Dxt
            SE=-TG1(IXABS,ITABS)
            IF(IXABS == ITABS) SE=SE+Two*OVL
!           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SE*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
      CASE(7)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SE(t,x)=2*dxt - Dxt
            SE=-TG1(IXABS,ITABS)
            IF(IXABS == ITABS) SE=SE+Two*OVL
!           HEBLK=HEBLK+SE*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SE*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
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
            SFtuxy=Two*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=Two*TG2(ITABS,IYABS,IUABS,IXABS)

            SFP=SFtuxy + SFtuyx
!           HEBLK=HEBLK+SFP*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SFP*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
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
            SFtuxy=Two*TG2(ITABS,IXABS,IUABS,IYABS)
            SFtuyx=Two*TG2(ITABS,IYABS,IUABS,IXABS)

            SFM=SFtuxy - SFtuyx
!           HEBLK=HEBLK+SFM*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SFM*SCAL,V1(IAS),NAS,V2(JAS),NAS)
          END DO
        END DO
!***********************************************************************
! CASES GP, GM
! Compute and use SG(ITABS , IXABS) (Same for cases GP and GM)
!***********************************************************************
      CASE(10)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SG(t,x)= Gtx
            SG= TG1(ITABS,IXABS)

!           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
!           Do i = 1, NISBLK
!             V2(JAS+NAS*(i-1)) = SG*V1(IAS+NAS*(i-1))
!           End Do
          END DO
        END DO
!***********************************************************************
      CASE(11)
        DO IAS=1,NAS
          ITABS=IAS+NAES(ISYM)
          DO JAS=1,NAS
            IXABS=JAS+NAES(ISYM)
! Formula used: SG(t,x)= Gtx
            SG= TG1(ITABS,IXABS)

!           HEBLK=HEBLK+SG*DDOT_(NISBLK,V2(JAS),NAS,V1(IAS),NAS)
            Call DaXpY_(NISBLK,SG*SCAL,V1(IAS),NAS,V2(JAS),NAS)
!           Do i = 1, NISBLK
!             V2(JAS+NAS*(i-1)) = SG*V1(IAS+NAS*(i-1))
!           End Do
          END DO
        END DO
!***********************************************************************
      CASE(12)
        IF(ABS(OVL) >= 1.0e-12_wp) THEN
!         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
!***********************************************************************
      CASE(13)
        IF(ABS(OVL) >= 1.0e-12_wp) THEN
!         HEBLK=HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
        END IF
!***********************************************************************
      END SELECT
      Return
      end subroutine MS_STRANS_BLK
