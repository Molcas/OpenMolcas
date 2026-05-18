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

subroutine DerHEffX_BLK(ICASE,ISYM,NAS,nvlen,NASHT,NTG3,IISTA,IIEND,V1,V2,OVL,DTG1,DTG2,DTG3)
! Compute a contribution to the coupling Hamiltonian element (HEL)
! defined as HEL = < ROOT1 | H * OMEGA | ROOT2 >. The contribution
! arises from the block V_(A,I), with A=1,NAS and I=IISTA,IIEND,
! with A the active superindex and I the inactive superindex. Since
! the inactive superindex is partitioned over processes, each process
! only computes part of the HEL value, which is then sum reduced in the
! calling subroutine.

use SUPERINDEX, only: MTGEU, MTGTU, MTU, MTUV
use caspt2_module, only: NAES, NTGEUES, NTGTUES, NTUES, NTUVES
use Constants, only: Zero, Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ICASE, ISYM, NAS, nvlen, NASHT, NTG3, IISTA, IIEND
real(kind=wp), intent(in) :: V1(nvlen), V2(nvlen)
real(kind=wp), intent(out) :: OVL
real(kind=wp), intent(inout) :: DTG1(NASHT,NASHT), DTG2(NASHT,NASHT,NASHT,NASHT), DTG3(NTG3)
integer(kind=iwp) :: IAS, IAS1, IAS2, IASABS, IND1, IND2, IND3, ITABS, ITG3, IUABS, IVABS, IXABS, IYABS, IZABS, JAS, JAS1, JAS2, &
                     JASABS, JND1, JND2, JND3, NAS1, NISBLK
real(kind=wp) :: VAL, VAL11, VAL12, VAL21, VAL22
real(kind=wp), external :: ddot_

! The dimension of TG3 is NTG3=(NASHT**2+2 over 3)

if (IISTA <= 0) return

OVL = Zero
NISBLK = IIEND-IISTA+1
select case (ICASE)
  !*********************************************************************
  case (1)
    do IAS=1,NAS
      IASABS = NTUVES(ISYM)+IAS
      ITABS = MTUV(1,IASABS)
      IUABS = MTUV(2,IASABS)
      IVABS = MTUV(3,IASABS)
      do JAS=1,NAS
        JASABS = NTUVES(ISYM)+JAS
        IXABS = MTUV(1,JASABS)
        IYABS = MTUV(2,JASABS)
        IZABS = MTUV(3,JASABS)
        ! Compute and use SA(ITABS IUABS IVABS, IXABS IYABS IZABS)
        ! Formulae used:
        !  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
        ! Gvuxtyz is stored using full permutation symmetry of three pairs
        ! (vu),(xt), and (yz):
        IND1 = IVABS+NASHT*(IUABS-1)
        IND2 = IXABS+NASHT*(ITABS-1)
        IND3 = IYABS+NASHT*(IZABS-1)
        if (IND2 > IND3) then
          if (IND1 > IND2) then
            JND1 = IND1
            JND2 = IND2
            JND3 = IND3
          else if (IND1 > IND3) then
            JND1 = IND2
            JND2 = IND1
            JND3 = IND3
          else
            JND1 = IND2
            JND2 = IND3
            JND3 = IND1
          end if
        else
          if (IND1 > IND3) then
            JND1 = IND1
            JND2 = IND3
            JND3 = IND2
          else if (IND1 > IND2) then
            JND1 = IND3
            JND2 = IND1
            JND3 = IND2
          else
            JND1 = IND3
            JND2 = IND2
            JND3 = IND1
          end if
        end if
        ITG3 = ((JND1+1)*JND1*(JND1-1))/6+(JND2*(JND2-1))/2+JND3
        !  SA(tuv,xyz) =  -Gvuxtyz -dyu Gvzxt - dyt Gvuxz - dxu Gvtyz - dxu dyt Gvz +2 dtx Gvuyz + 2 dtx dyu Gvz
        ! Compute TMP=Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        if (IXABS == ITABS) then
          DTG2(IVABS,IUABS,IYABS,IZABS) = DTG2(IVABS,IUABS,IYABS,IZABS)+Two*VAL
          if (IYABS == IUABS) DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS)+Two*VAL
        end if
        VAL = -VAL
        DTG3(ITG3) = DTG3(ITG3)+VAL
        if (IYABS == IUABS) DTG2(IVABS,IZABS,IXABS,ITABS) = DTG2(IVABS,IZABS,IXABS,ITABS)+VAL
        if (IYABS == ITABS) then
          DTG2(IVABS,IUABS,IXABS,IZABS) = DTG2(IVABS,IUABS,IXABS,IZABS)+VAL
          if (IXABS == IUABS) DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS)+VAL
        end if
        if (IXABS == IUABS) DTG2(IVABS,ITABS,IYABS,IZABS) = DTG2(IVABS,ITABS,IYABS,IZABS)+VAL
      end do
    end do
    !*******************************************************************
  case (4)
    do IAS=1,NAS
      IASABS = NTUVES(ISYM)+IAS
      IXABS = MTUV(1,IASABS)
      IUABS = MTUV(2,IASABS)
      IVABS = MTUV(3,IASABS)
      do JAS=1,NAS
        JASABS = NTUVES(ISYM)+JAS
        ITABS = MTUV(1,JASABS)
        IYABS = MTUV(2,JASABS)
        IZABS = MTUV(3,JASABS)
        ! Compute and use SC(IXABS IUABS IVABS, ITABS IYABS IZABS)
        ! In SBMAT, the formula is written as
        !  SC(tuv,xyz) = Gvutxyz +dyu Gvztx + dyx Gvutz + dtu Gvxyz + dtu dyx Gvz
        ! Rewritten, in order to reuse same quantities as in SA:
        !  SC(xuv,tyz) = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
        ! Gvuxtyz is stored using full permutation symmetry of three pairs
        ! (vu),(xt), and (yz):
        IND1 = IVABS+NASHT*(IUABS-1)
        IND2 = IXABS+NASHT*(ITABS-1)
        IND3 = IYABS+NASHT*(IZABS-1)
        if (IND2 > IND3) then
          if (IND1 > IND2) then
            JND1 = IND1
            JND2 = IND2
            JND3 = IND3
          else if (IND1 > IND3) then
            JND1 = IND2
            JND2 = IND1
            JND3 = IND3
          else
            JND1 = IND2
            JND2 = IND3
            JND3 = IND1
          end if
        else
          if (IND1 > IND3) then
            JND1 = IND1
            JND2 = IND3
            JND3 = IND2
          else if (IND1 > IND2) then
            JND1 = IND3
            JND2 = IND1
            JND3 = IND2
          else
            JND1 = IND3
            JND2 = IND2
            JND3 = IND1
          end if
        end if
        ITG3 = ((JND1+1)*JND1*(JND1-1))/6+(JND2*(JND2-1))/2+JND3
        !  SC(xuv,tyz) (rewritten, swapping x and t)
        !    = Gvuxtyz +dyu Gvzxt + dyt Gvuxz + dxu Gvtyz + dxu dyt Gvz
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG3(ITG3) = DTG3(ITG3)+VAL
        if (IYABS == IUABS) DTG2(IVABS,IZABS,IXABS,ITABS) = DTG2(IVABS,IZABS,IXABS,ITABS)+VAL
        if (IYABS == ITABS) then
          DTG2(IVABS,IUABS,IXABS,IZABS) = DTG2(IVABS,IUABS,IXABS,IZABS)+VAL
          if (IXABS == IUABS) DTG1(IVABS,IZABS) = DTG1(IVABS,IZABS)+VAL
        end if
        if (IXABS == IUABS) DTG2(IVABS,ITABS,IYABS,IZABS) = DTG2(IVABS,ITABS,IYABS,IZABS)+VAL
      end do
    end do
    !*******************************************************************
  case (2)
    do IAS=1,NAS
      IASABS = NTGEUES(ISYM)+IAS
      ITABS = MTGEU(1,IASABS)
      IUABS = MTGEU(2,IASABS)
      do JAS=1,NAS
        JASABS = NTGEUES(ISYM)+JAS
        IXABS = MTGEU(1,JASABS)
        IYABS = MTGEU(2,JASABS)
        ! Formulae used:
        !    SB(tu,xy)= 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu -4dxu dyt + 2dxu Gyt
        !    SB(tu,yx)= 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu -4dyu dxt + 2dyu Gxt
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG2(IXABS,ITABS,IYABS,IUABS) = DTG2(IXABS,ITABS,IYABS,IUABS)+Two*VAL
        DTG2(IYABS,ITABS,IXABS,IUABS) = DTG2(IYABS,ITABS,IXABS,IUABS)+Two*VAL
        if (IXABS == ITABS) then
          DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS)-Four*VAL+Two*VAL
          if (IYABS == IUABS) OVL = OVL+Eight*VAL-Four*VAL
        end if
        if (IYABS == IUABS) DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)-Four*VAL+Two*VAL
        if (IYABS == ITABS) then
          DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS)+Two*VAL-Four*VAL
          if (IXABS == IUABS) OVL = OVL-Four*VAL+Eight*VAL
        end if
        if (IXABS == IUABS) DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS)+Two*VAL-Four*VAL
      end do
    end do
    !*******************************************************************
  case (3)
    do IAS=1,NAS
      IASABS = NTGTUES(ISYM)+IAS
      ITABS = MTGTU(1,IASABS)
      IUABS = MTGTU(2,IASABS)
      do JAS=1,NAS
        JASABS = NTGTUES(ISYM)+JAS
        IXABS = MTGTU(1,JASABS)
        IYABS = MTGTU(2,JASABS)
        ! Formulae used:
        !    SB(tu,xy)= 2 Gxtyu -4dxt Gyu -4dyu Gxt +2dyt Gxu + 8 dxt dyu -4dxu dyt + 2dxu Gyt
        !    SB(tu,yx)= 2 Gytxu -4dyt Gxu -4dxu Gyt +2dxt Gyu + 8 dyt dxu -4dyu dxt + 2dyu Gxt
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG2(IXABS,ITABS,IYABS,IUABS) = DTG2(IXABS,ITABS,IYABS,IUABS)+Two*VAL
        DTG2(IYABS,ITABS,IXABS,IUABS) = DTG2(IYABS,ITABS,IXABS,IUABS)-Two*VAL
        if (IXABS == ITABS) then
          DTG1(IYABS,IUABS) = DTG1(IYABS,IUABS)-Four*VAL-Two*VAL
          if (IYABS == IUABS) OVL = OVL+Eight*VAL+Four*VAL
        end if
        if (IYABS == IUABS) DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)-Four*VAL-Two*VAL
        if (IYABS == ITABS) then
          DTG1(IXABS,IUABS) = DTG1(IXABS,IUABS)+Two*VAL+Four*VAL
          if (IXABS == IUABS) OVL = OVL-Four*VAL-Eight*VAL
        end if
        if (IXABS == IUABS) DTG1(IYABS,ITABS) = DTG1(IYABS,ITABS)+Two*VAL+Four*VAL
      end do
    end do
    !*******************************************************************
  case (5)
    NAS1 = NAS/2
    do IAS1=1,NAS1
      IAS2 = IAS1+NAS1
      IASABS = NTUES(ISYM)+IAS1
      ITABS = MTU(1,IASABS)
      IUABS = MTU(2,IASABS)
      do JAS1=1,NAS1
        JAS2 = JAS1+NAS1
        JASABS = NTUES(ISYM)+JAS1
        IXABS = MTU(1,JASABS)
        IYABS = MTU(2,JASABS)
        ! Formulae used:
        !    SD11(tu1,xy1)=2*(Gutxy + dtx Guy)
        !    SD12(tu2,xy1)= -(Gutxy + dtx Guy)
        !    SD21(tu2,xy1)= -(Gutxy + dtx Guy)
        !    SD22(tu2,xy2)= -Gxtuy +2*dtx Guy
        VAL11 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS1),NAS)
        VAL12 = DDOT_(NISBLK,V1(IAS1),NAS,V2(JAS2),NAS)
        VAL21 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS1),NAS)
        VAL22 = DDOT_(NISBLK,V1(IAS2),NAS,V2(JAS2),NAS)
        DTG2(IUABS,ITABS,IXABS,IYABS) = DTG2(IUABS,ITABS,IXABS,IYABS)+Two*VAL11-VAL12-VAL21
        DTG2(IXABS,ITABS,IUABS,IYABS) = DTG2(IXABS,ITABS,IUABS,IYABS)-VAL22
        if (ITABS == IXABS) DTG1(IUABS,IYABS) = DTG1(IUABS,IYABS)+Two*VAL11-VAL12-VAL21+Two*VAL22
      end do
    end do
    !*******************************************************************
  case (6)
    do IAS=1,NAS
      ITABS = IAS+NAES(ISYM)
      do JAS=1,NAS
        IXABS = JAS+NAES(ISYM)
        ! Formula used: SE(t,x)=2*dxt - Dxt
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)-VAL
        if (IXABS == ITABS) OVL = OVL+Two*VAL
      end do
    end do
    !*******************************************************************
  case (7)
    do IAS=1,NAS
      ITABS = IAS+NAES(ISYM)
      do JAS=1,NAS
        IXABS = JAS+NAES(ISYM)
        ! Formula used: SE(t,x)=2*dxt - Dxt
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG1(IXABS,ITABS) = DTG1(IXABS,ITABS)-VAL
        if (IXABS == ITABS) OVL = OVL+Two*VAL
      end do
    end do
    !*******************************************************************
  case (8)
    ! ========================================================
    ! Compute and use SFP(ITABS IUABS , IXABS IYABS)
    ! and (later, similar) SFM(ITABS IUABS , IXABS IYABS)
    do IAS=1,NAS
      IASABS = NTGEUES(ISYM)+IAS
      ITABS = MTGEU(1,IASABS)
      IUABS = MTGEU(2,IASABS)
      do JAS=1,NAS
        JASABS = NTGEUES(ISYM)+JAS
        IXABS = MTGEU(1,JASABS)
        IYABS = MTGEU(2,JASABS)
        ! Formulae used:
        !    SF(tu,xy)= 2 Gtxuy
        !    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
        !    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG2(ITABS,IXABS,IUABS,IYABS) = DTG2(ITABS,IXABS,IUABS,IYABS)+Two*VAL
        DTG2(ITABS,IYABS,IUABS,IXABS) = DTG2(ITABS,IYABS,IUABS,IXABS)+Two*VAL
      end do
    end do
    !*******************************************************************
  case (9)
    ! ========================================================
    ! Compute and use SFM(ITABS IUABS, IXABS ,IYABS)
    do IAS=1,NAS
      IASABS = NTGTUES(ISYM)+IAS
      ITABS = MTGTU(1,IASABS)
      IUABS = MTGTU(2,IASABS)
      do JAS=1,NAS
        JASABS = NTGTUES(ISYM)+JAS
        IXABS = MTGTU(1,JASABS)
        IYABS = MTGTU(2,JASABS)
        ! Formulae used:
        !    SF(tu,xy)= 4 Ptxuy
        !    SFP(tu,xy)=SF(tu,xy)+SF(tu,yx)
        !    SFM(tu,xy)=SF(tu,xy)-SF(tu,yx)
        VAL = DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
        DTG2(ITABS,IXABS,IUABS,IYABS) = DTG2(ITABS,IXABS,IUABS,IYABS)+Two*VAL
        DTG2(ITABS,IYABS,IUABS,IXABS) = DTG2(ITABS,IYABS,IUABS,IXABS)-Two*VAL
      end do
    end do
    !*******************************************************************
    ! CASES GP, GM
    ! Compute and use SG(ITABS , IXABS) (Same for cases GP and GM)
    !*******************************************************************
  case (10)
    do IAS=1,NAS
      ITABS = IAS+NAES(ISYM)
      do JAS=1,NAS
        IXABS = JAS+NAES(ISYM)
        ! Formula used: SG(t,x)= Gtx
        !           SG= TG1(ITABS,IXABS)
        DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)+DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
      end do
    end do
    !*******************************************************************
  case (11)
    do IAS=1,NAS
      ITABS = IAS+NAES(ISYM)
      do JAS=1,NAS
        IXABS = JAS+NAES(ISYM)
        ! Formula used: SG(t,x)= Gtx
        !           SG= TG1(ITABS,IXABS)
        DTG1(ITABS,IXABS) = DTG1(ITABS,IXABS)+DDOT_(NISBLK,V1(IAS),NAS,V2(JAS),NAS)
      end do
    end do
    !*******************************************************************
  case (12)
    OVL = OVL+DDOT_(NAS*NISBLK,V2,1,V1,1)
    !if (abs(OVL) >= 1.0e-12_wp) HEBLK = HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
    !*******************************************************************
  case (13)
    OVL = OVL+DDOT_(NAS*NISBLK,V2,1,V1,1)
    !if (abs(OVL) >= 1.0e-12_wp) HEBLK = HEBLK+OVL*DDOT_(NAS*NISBLK,V2,1,V1,1)
  !*********************************************************************
end select

return

end subroutine DerHEffX_BLK
