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

      Subroutine CLagDXC_DP(iSym,nAS,nAshT,BDER,SDER,DG1,DG2,DF1,DF2,   &
     &                      DEPSA,DEASUM,iLo,iHi,jLo,jHi,LDC,           &
     &                      G1,G2,SC,SC2,lg_S)

      USE SUPERINDEX, only: MTUV
      use caspt2_global, only:ipea_shift
      use caspt2_module, only: EASUM, EPSA, NTUVES
      use Constants, only: Zero, Half, Four
      use definitions, only: wp, iwp
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, nProcs
#endif

      implicit none

#ifdef _MOLCAS_MPP_
#include "global.fh"
#else
#include "macros.fh"
#endif

      integer(kind=iwp), intent(in) :: iSym, nAS, nAshT, iLo, iHi, jLo, &
     &                                 jHi, LDC, lg_S
      real(kind=wp), intent(in) :: BDER((iHi-iLo+1)*(jHi-jLo+1)),       &
     &  SC((iHi-iLo+1)*(jHi-jLo+1)), SC2((iHi-iLo+1)*(jHi-jLo+1))
      real(kind=wp), intent(inout) :: SDER((iHi-iLo+1)*(jHi-jLo+1)),    &
     &  DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT),                 &
     &  DF1(nAshT,nAshT), DF2(nAshT,nAshT,nAshT,nAshT),                 &
     &  DEPSA(nAshT,nAshT), DEASUM, G1(nAshT,nAshT),                    &
     &  G2(nAshT,nAshT,nAshT,nAshT)

      integer(kind=iwp) :: ISADR, NROW, iLoS, jLoS, IXYZ, IXYZABS,      &
     &  IXABS, IYABS, IZABS, ITUV, ITUVABS, ITABS, IUABS, IVABS,        &
     &  ISADR2, iWabs, iTWV, iXWZ
      real(kind=wp) :: EU, EY, EYU, FACT, ValB, bsBDER, ValS
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: irank, iHiS, jHiS
#endif
!
!     LDC == 0, if not parallel; SC is triangular
!     LDC /= 0, if parallel    ; SC is square
!     In both cases, BDER and SDER are square
!
      ISADR=0
      NROW = 0
      iLoS = 0
      jLoS = 0
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        irank = 0
        CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
        NROW = jHiS-jLoS+1 !! = NAS
        CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
      end if
#else
      unused_var(lg_S)
#endif
      DO IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EY=EPSA(IYABS)
        DO ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          EU=EPSA(IUABS)
          EYU=EY + EU
          FACT=EYU-EASUM
          IF (LDC == 0) THEN
            ISADR=1+iTUV-iLo+NAS*(iXYZ-jLo)
          ELSE
            ISADR=1+iTUV-iLo+LDC*(iXYZ-jLo)
          END IF
          ValB=BDER(ISADR)

          If (iTUV == iXYZ .and. ipea_shift /= Zero) Then
!           !! BC in the next equation refers to the active overlap
!    !! ipea_shift*Half*BC(ISADR)*(Four-DREF(IDT)-DREF(IDV)+DREF(IDU))
            bsBDER = ipea_shift*Half*ValB
            SDER(iSAdr) = SDER(iSAdr) + bsBDER*(Four                    &
     &       -G1(iTabs,iTabs)+G1(iUabs,iUabs)-G1(iVabs,iVabs))
            IF (LDC == 0) THEN
              iSAdr2 = iTUV*(iTUV+1)/2
            ELSE
              ISADR2 = 1+iTUV-iLo+LDC*(iTUV-jLo)
            END IF
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) - bsBDER*SC(iSAdr2)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) + bsBDER*SC(iSAdr2)
            DG1(iVabs,iVabs) = DG1(iVabs,iVabs) - bsBDER*SC(iSAdr2)
          End If

          !! First VALUE contribution in MKBC_DP (FACT)
!         IF (LDC == 0) ISADR=ITUV*(ITUV-1)/2+IXYZ
          SDER(ISADR) = SDER(ISADR) + FACT*ValB
          ValS=SDER(ISADR)

          !! DEPSA can be computed simultaneously with the F3
          !! contributions, but remember that SC here is the overlap,
          !! whereas SC in FG3 subroutines are just G3 matrix
          if (ldc == 0) then
            Do iWabs = 1, nAshT
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              iSAdr2 = Max(iTWV,iXYZ)*(Max(iTWV,iXYZ)-1)/2              &
     &               + Min(iTWV,iXYZ)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)                   &
     &          + ValB*SC(iSAdr2)

              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              iSAdr2 = Max(iTUV,iXWZ)*(Max(iTUV,iXWZ)-1)/2              &
     &               + Min(iTUV,iXWZ)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)                   &
     &          + ValB*SC(iSAdr2)
            End Do
          else
            !! assume SC has all columns (which is reasonable)
!           iTWV = iTabs+nAshT**2*(iVabs-1)
!           ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
!           Call DaXpY_(nAshT,ValB,SC2(iSAdr2),nAshT,DEPSA(1,iUabs),1)
!           iXWZ = iXabs+nAshT**2*(iZabs-1)
!           ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
!           Call DaXpY_(nAshT,ValB,SC(iSAdr2),LDC*nAshT,
!    *                             DEPSA(1,iYabs),1)
            Do iWabs = 1, nAshT
              !! we do not have all elements, so use distributed memory
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)                   &
     &          + ValB*SC2(iSAdr2)

              !! we have all elements, so local memory
              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)                   &
     &          + ValB*SC(iSAdr2)
            End Do
          end if

          !! If non-parallel, overlap is triangular
          !! If parallel, the index is the same to that of SDER
          IF (LDC == 0) THEN
            iSAdr = Max(iTUV,iXYZ)*(Max(iTUV,iXYZ)-1)/2                 &
     &            + Min(iTUV,iXYZ)
          END IF
          DEASUM = DEASUM - ValB*SC(iSAdr)

!         dyu ( Fvztx - EPSA(u)*Gvztx )
!         dyu Gvztx
          IF(IYABS == IUABS) THEN
            !! VALUE=VALUE+Two*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iZabs,iTabs,iXabs)                                &
     &        = DF2(iVabs,iZabs,iTabs,iXabs) + ValB
            DG2(iVabs,iZabs,iTabs,iXabs)                                &
     &        = DG2(iVabs,iZabs,iTabs,iXabs) - EU*ValB

            !! VALUE=VALUE+Two*PREF(IP)
            DG2(iVabs,iZabs,iTabs,iXabs)                                &
     &        = DG2(iVabs,iZabs,iTabs,iXabs) + ValS
          END IF
          DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs)                       &
     &      - ValB*G2(iVabs,iZabs,iTabs,iXabs)

!         dyx ( Fvutz - EPSA(y)*Gvutz )
!         dyx Gvutz -> dut Gzyxv
          IF(IYABS == IXABS) THEN
            !! VALUE=VALUE+Two*(FP(IP)-EY*PREF(IP))
            DF2(iVabs,iUabs,iTabs,iZabs)                                &
     &        = DF2(iVabs,iUabs,iTabs,iZabs) + ValB
            DG2(iVabs,iUabs,iTabs,iZabs)                                &
     &        = DG2(iVabs,iUabs,iTabs,iZabs) - EY*ValB

            !! VALUE=VALUE+Two*PREF(IP)
            DG2(iVabs,iUabs,iTabs,iZabs)                                &
     &        = DG2(iVabs,iUabs,iTabs,iZabs) + ValS
          END IF
          DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs)                       &
     &      - ValB*G2(iVabs,iUabs,iTabs,iZabs)

!         dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz -
!                (EPSA(u)+EPSA(y)*dyz Gvz)
!         dtu Gvxyz + dtu dyx Gvz
          IF(ITABS == IUABS) THEN
            !! VALUE=VALUE+Two*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iXabs,iYabs,iZabs)                                &
     &        = DF2(iVabs,iXabs,iYabs,iZabs) + ValB
            DG2(iVabs,iXabs,iYabs,iZabs)                                &
     &        = DG2(iVabs,iXabs,iYabs,iZabs) - EU*ValB

            !! VALUE=VALUE+Two*PREF(IP)
            DG2(iVabs,iXabs,iYabs,iZabs)                                &
     &        = DG2(iVabs,iXabs,iYabs,iZabs) + ValS
            IF(IYABS == IXABS) THEN
              !! VALUE=VALUE+FD(ID)-EYU*DREF(ID)
              DF1(iVabs,iZabs) = DF1(iVabs,iZabs) + ValB
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) - EYU*ValB

              !! VALUE=VALUE+DREF((ID1*(ID1-1))/2+ID2)
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) + ValS
            END IF
          END IF
          DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs)                       &
     &      - ValB*G2(iVabs,iXabs,iYabs,iZabs)
          If (iYabs == iXabs)                                           &
     &    DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs) - ValB*G1(iVabs,iZabs)
          If (iTabs == iUabs)                                           &
     &    DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs) - ValB*G1(iVabs,iZabs)
        end do
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. IXYZ == iHiS .and. iRank/=NPROCS-1) then
          irank = irank + 1
          CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
          CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
        end if
#endif
      end do

      Return

      End Subroutine CLagDXC_DP
