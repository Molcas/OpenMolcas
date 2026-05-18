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

subroutine CLagDXA_DP(iSym,nAS,nAshT,BDER,SDER,DG1,DG2,DF1,DF2,DEPSA,DEASUM,iLo,iHi,jLo,jHi,LDA,G1,G2,SA,SA2,lg_S)

use SUPERINDEX, only: MTUV
use caspt2_global, only: ipea_shift
use caspt2_module, only: EASUM, EPSA, NTUVES
use Constants, only: Zero, Half, Two
use definitions, only: wp, iwp
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif

implicit none
#ifdef _MOLCAS_MPP_
#include "global.fh"
#else
#include "macros.fh"
#endif
integer(kind=iwp), intent(in) :: iSym, nAS, nAshT, iLo, iHi, jLo, jHi, LDA, lg_S
real(kind=wp), intent(in) :: BDER((iHi-iLo+1)*(jHi-jLo+1)), SA((iHi-iLo+1)*(jHi-jLo+1)), SA2((iHi-iLo+1)*(jHi-jLo+1))
real(kind=wp), intent(inout) :: SDER((iHi-iLo+1)*(jHi-jLo+1)), DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DF1(nAshT,nAshT), &
                                DF2(nAshT,nAshT,nAshT,nAshT), DEPSA(nAshT,nAshT), DEASUM, G1(nAshT,nAshT), &
                                G2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: ISADR, NROW, iLoS, jLoS, IXYZ, IXYZABS, IXABS, IYABS, IZABS, ITUV, ITUVABS, ITABS, IUABS, IVABS, ISADR2, &
                     iWabs, iTWV, iXWZ, iWYZ, iWUV
real(kind=wp) :: ET, EU, ETU, EX, EY, FACT, ValB, bsBDER, ValS
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: irank, iHiS, jHiS
#endif

! LDA == 0, if not parallel; SA is triangular
! LDA /= 0, if parallel    ; SA is square
! In both cases, BDER and SDER are square

ISADR = 0
NROW = 0
iLoS = 0
jLoS = 0
#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  irank = 0
  call GA_Distribution(lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
  NROW = jHiS-jLoS+1 !! = NAS
  call GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SA2,NROW)
end if
#else
unused_var(lg_S)
#endif
do IXYZ=jLo,jHi
  IXYZABS = IXYZ+NTUVES(ISYM)
  IXABS = MTUV(1,IXYZABS)
  IYABS = MTUV(2,IXYZABS)
  IZABS = MTUV(3,IXYZABS)
  EX = EPSA(IXABS)
  EY = EPSA(IYABS)
  do ITUV=iLo,iHi
    ITUVABS = ITUV+NTUVES(ISYM)
    ITABS = MTUV(1,ITUVABS)
    IUABS = MTUV(2,ITUVABS)
    IVABS = MTUV(3,ITUVABS)
    ET = EPSA(ITABS)
    EU = EPSA(IUABS)
    ETU = ET+EU
    FACT = EY+EU+EX+ET-EASUM
    if (LDA == 0) then
      ISADR = 1+iTUV-iLo+NAS*(iXYZ-jLo)
    else
      ISADR = 1+iTUV-iLo+LDA*(iXYZ-jLo)
    end if
    ValB = BDER(ISADR)

    if ((iTUV == iXYZ) .and. (ipea_shift /= Zero)) then
      !! BA in the next equation refers to the active overlap
      ! ipea_shift*Half*BA(ISADR)*(Two-DREF(IDV)+DREF(IDT)+DREF(IDU))
      bsBDER = ipea_shift*Half*ValB
      SDER(iSAdr) = SDER(iSAdr)+bsBDER*(Two+G1(iTabs,iTabs)+G1(iUabs,iUabs)-G1(iVabs,iVabs))
      if (LDA == 0) then
        iSAdr2 = iTUV*(iTUV+1)/2
      else
        ISADR2 = 1+iTUV-iLo+LDA*(iTUV-jLo)
      end if
      DG1(iTabs,iTabs) = DG1(iTabs,iTabs)+bsBDER*SA(iSAdr2)
      DG1(iUabs,iUabs) = DG1(iUabs,iUabs)+bsBDER*SA(iSAdr2)
      DG1(iVabs,iVabs) = DG1(iVabs,iVabs)-bsBDER*SA(iSAdr2)
    end if

    !! First VALUE contribution in MKBC_DP (FACT)
    SDER(ISADR) = SDER(ISADR)+FACT*ValB
    ValS = SDER(ISADR)

    if (lda == 0) then
      do iWabs=1,nAshT
        !! EU derivative
        iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
        iSAdr2 = max(iTWV,iXYZ)*(max(iTWV,iXYZ)-1)/2+min(iTWV,iXYZ)
        DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)+ValB*SA(iSAdr2)

        !! EY derivative
        iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
        iSAdr2 = max(iTUV,iXWZ)*(max(iTUV,iXWZ)-1)/2+min(iTUV,iXWZ)
        DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)+ValB*SA(iSAdr2)

        !! EX derivative
        iWYZ = iWabs+nAshT*(iYabs-1)+nAshT**2*(iZabs-1)
        iSAdr2 = max(iTUV,iWYZ)*(max(iTUV,iWYZ)-1)/2+min(iTUV,iWYZ)
        DEPSA(iWabs,iXabs) = DEPSA(iWabs,iXabs)+ValB*SA(iSAdr2)

        !! ET derivative
        iWUV = iWabs+nAshT*(iUabs-1)+nAshT**2*(iVabs-1)
        iSAdr2 = max(iWUV,iXYZ)*(max(iWUV,iXYZ)-1)/2+min(iWUV,iXYZ)
        DEPSA(iWabs,iTabs) = DEPSA(iWabs,iTabs)+ValB*SA(iSAdr2)
      end do
    else
      do iWabs=1,nAshT
        !! EU derivative
        iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
        ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
        DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)+ValB*SA2(iSAdr2)

        !! EY derivative
        iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
        ISADR2 = 1+iTUV-iLo+LDA*(iXWZ-jLo)
        DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)+ValB*SA(iSAdr2)

        !! EX derivative
        iWYZ = iWabs+nAshT*(iYabs-1)+nAshT**2*(iZabs-1)
        ISADR2 = 1+iTUV-iLo+LDA*(iWYZ-jLo)
        DEPSA(iWabs,iXabs) = DEPSA(iWabs,iXabs)+ValB*SA(iSAdr2)

        !! ET derivative
        iWUV = iWabs+nAshT*(iUabs-1)+nAshT**2*(iVabs-1)
        ISADR2 = 1+iWUV-jLoS+NROW*(iXYZ-iLoS)
        DEPSA(iWabs,iTabs) = DEPSA(iWabs,iTabs)+ValB*SA2(iSAdr2)
      end do
    end if

    if (LDA == 0) iSAdr = max(iTUV,iXYZ)*(max(iTUV,iXYZ)-1)/2+min(iTUV,iXYZ)
    DEASUM = DEASUM-ValB*SA(iSAdr)

    ! 2dtx ( Fvuyz-Et*Gvuyz )
    ! 2 dtx Gvuyz + 2 dtx dyu Gvz
    if (iTabs == iXabs) then
      !! VALUE=VALUE+Four*(FP(IP)-ET*PREF(IP))
      DF2(iVabs,iUabs,iYabs,iZabs) = DF2(iVabs,iUabs,iYabs,iZabs)+Two*ValB
      DG2(iVabs,iUabs,iYabs,iZabs) = DG2(iVabs,iUabs,iYabs,iZabs)-Two*ET*ValB

      !! VALUE=VALUE+Four*PREF(IP)
      DG2(iVabs,iUabs,iYabs,iZabs) = DG2(iVabs,iUabs,iYabs,iZabs)+Two*ValS
      !! VALUE=VALUE+Two*DREF(ID)
      if (iYabs == iUabs) DG1(iVabs,iZabs) = DG1(iVabs,iZabs)+Two*ValS
    end if
    DEPSA(iTabs,iXabs) = DEPSA(iTabs,iXabs) -Two*ValB*G2(iVabs,iUabs,iYabs,iZabs)

    ! dxu ( -Fvtyz + Eu*Gvtyz )
    ! -dxu Gvtyz -dxu dyt Gvz
    if (iXabs == iUabs) then
      !! VALUE=VALUE-Two*(FP(IP)-EU*PREF(IP))
      DF2(iVabs,iTabs,iYabs,iZabs) = DF2(iVabs,iTabs,iYabs,iZabs)-ValB
      DG2(iVabs,iTabs,iYabs,iZabs) = DG2(iVabs,iTabs,iYabs,iZabs)+EU*ValB
      !! VALUE=VALUE - Two*PREF(IP)
      DG2(iVabs,iTabs,iYabs,iZabs) = DG2(iVabs,iTabs,iYabs,iZabs)-ValS
      !! VALUE=VALUE - DREF(ID)
      if (iYabs == iTabs) DG1(iVabs,iZabs) = DG1(iVabs,iZabs)-ValS
    end if
    DEPSA(iXabs,iUabs) = DEPSA(iXabs,iUabs) +ValB*G2(iVabs,iTabs,iYabs,iZabs)

    ! dyt ( -Fvuxz + Et*Gvuxz +dxu (-Fvz+(Et+Eu)*Gvz))
    ! -dyt Gvuxz
    if (iYabs == iTabs) then
      !! VALUE=VALUE-Two*(FP(IP)-ET*PREF(IP))
      DF2(iVabs,iUabs,iXabs,iZabs) = DF2(iVabs,iUabs,iXabs,iZabs)-ValB
      DG2(iVabs,iUabs,iXabs,iZabs) = DG2(iVabs,iUabs,iXabs,iZabs)+ET*ValB
      if (iXabs == iUabs) then
        !! VALUE=VALUE - (FD(ID)-ETU*DREF(ID))
        DF1(iVabs,iZabs) = DF1(iVabs,iZabs)-ValB
        DG1(iVabs,iZabs) = DG1(iVabs,iZabs)+ETU*ValB
      end if

      !! VALUE=VALUE - Two*PREF(IP)
      DG2(iVabs,iUabs,iXabs,iZabs) = DG2(iVabs,iUabs,iXabs,iZabs)-ValS
    end if
    DEPSA(iYabs,iTabs) = DEPSA(iYabs,iTabs) +ValB*G2(iVabs,iUabs,iXabs,iZabs)
    if (iYabs == iTabs) DEPSA(iXabs,iUabs) = DEPSA(iXabs,iUabs)+ValB*G1(iVabs,iZabs)
    if (iXabs == iUabs) DEPSA(iYabs,iTabs) = DEPSA(iYabs,iTabs)+ValB*G1(iVabs,iZabs)

    ! -dyu Gvzxt
    if (iYabs == iUabs) then
      !! VALUE=VALUE-Two*(FP(IP)-EU*PREF(IP))
      DF2(iVabs,iZabs,iXabs,iTabs) = DF2(iVabs,iZabs,iXabs,iTabs)-ValB
      DG2(iVabs,iZabs,iXabs,iTabs) = DG2(iVabs,iZabs,iXabs,iTabs)+EU*ValB
      if (iXabs == iTabs) then
        !! VALUE=VALUE+Two*(FD(ID)-ETU*DREF(ID))
        DF1(iVabs,iZabs) = DF1(iVabs,iZabs)+Two*ValB
        DG1(iVabs,iZabs) = DG1(iVabs,iZabs)-Two*ETU*ValB
      end if

      !! VALUE=VALUE - Two*PREF(IP)
      DG2(iVabs,iZabs,iXabs,iTabs) = DG2(iVabs,iZabs,iXabs,iTabs)-ValS
    end if
    DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs) +ValB*G2(iVabs,iZabs,iXabs,iTabs)
    if (iYabs == iUabs) DEPSA(iXabs,iTabs) = DEPSA(iXabs,iTabs)-Two*ValB*G1(iVabs,iZabs)
    if (iXabs == iTabs) DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs)-Two*ValB*G1(iVabs,iZabs)
  end do
# ifdef _MOLCAS_MPP_
  if (is_real_par() .and. (IXYZ == iHiS) .and. (iRank /= NPROCS-1)) then
    irank = irank+1
    call GA_Distribution(lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
    call GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SA2,NROW)
  end if
# endif
end do

return

end subroutine CLagDXA_DP
