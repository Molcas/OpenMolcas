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

subroutine CLagDXC_DP(iSym,nAS,nAshT,BDER,SDER,DG1,DG2,DF1,DF2,DEPSA,DEASUM,iLo,iHi,jLo,jHi,LDC,G1,G2,SC,SC2,lg_S)

use SUPERINDEX, only: MTUV
use caspt2_global, only: ipea_shift
use caspt2_module, only: EASUM, EPSA, NTUVES
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, nProcs
#endif
use Constants, only: Zero, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, nAS, nAshT, iLo, iHi, jLo, jHi, LDC, lg_S
real(kind=wp), intent(in) :: BDER((iHi-iLo+1)*(jHi-jLo+1)), SC((iHi-iLo+1)*(jHi-jLo+1)), SC2((iHi-iLo+1)*(jHi-jLo+1))
real(kind=wp), intent(inout) :: SDER((iHi-iLo+1)*(jHi-jLo+1)), DG1(nAshT,nAshT), DG2(nAshT,nAshT,nAshT,nAshT), DF1(nAshT,nAshT), &
                                DF2(nAshT,nAshT,nAshT,nAshT), DEPSA(nAshT,nAshT), DEASUM, G1(nAshT,nAshT), &
                                G2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: iLoS, ISADR, ISADR2, ITABS, ITUV, ITUVABS, iTWV, IUABS, IVABS, iWabs, IXABS, iXWZ, IXYZ, IXYZABS, IYABS, &
                     IZABS, jLoS, NROW
real(kind=wp) :: bsBDER, EU, EY, EYU, FACT, ValB, ValS
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: iHiS, irank, jHiS
#include "global.fh"
#else
#include "macros.fh"
unused_var(lg_S)
#endif

! LDC == 0, if not parallel; SC is triangular
! LDC /= 0, if parallel    ; SC is square
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
  call GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
end if
#endif
do IXYZ=jLo,jHi
  IXYZABS = IXYZ+NTUVES(ISYM)
  IXABS = MTUV(1,IXYZABS)
  IYABS = MTUV(2,IXYZABS)
  IZABS = MTUV(3,IXYZABS)
  EY = EPSA(IYABS)
  do ITUV=iLo,iHi
    ITUVABS = ITUV+NTUVES(ISYM)
    ITABS = MTUV(1,ITUVABS)
    IUABS = MTUV(2,ITUVABS)
    IVABS = MTUV(3,ITUVABS)
    EU = EPSA(IUABS)
    EYU = EY+EU
    FACT = EYU-EASUM
    if (LDC == 0) then
      ISADR = 1+iTUV-iLo+NAS*(iXYZ-jLo)
    else
      ISADR = 1+iTUV-iLo+LDC*(iXYZ-jLo)
    end if
    ValB = BDER(ISADR)

    if ((iTUV == iXYZ) .and. (ipea_shift /= Zero)) then
      !! BC in the next equation refers to the active overlap
      !! ipea_shift*Half*BC(ISADR)*(Four-DREF(IDT)-DREF(IDV)+DREF(IDU))
      bsBDER = ipea_shift*Half*ValB
      SDER(iSAdr) = SDER(iSAdr)+bsBDER*(Four-G1(iTabs,iTabs)+G1(iUabs,iUabs)-G1(iVabs,iVabs))
      if (LDC == 0) then
        iSAdr2 = iTUV*(iTUV+1)/2
      else
        ISADR2 = 1+iTUV-iLo+LDC*(iTUV-jLo)
      end if
      DG1(iTabs,iTabs) = DG1(iTabs,iTabs)-bsBDER*SC(iSAdr2)
      DG1(iUabs,iUabs) = DG1(iUabs,iUabs)+bsBDER*SC(iSAdr2)
      DG1(iVabs,iVabs) = DG1(iVabs,iVabs)-bsBDER*SC(iSAdr2)
    end if

    !! First VALUE contribution in MKBC_DP (FACT)
    !if (LDC == 0) ISADR = ITUV*(ITUV-1)/2+IXYZ
    SDER(ISADR) = SDER(ISADR)+FACT*ValB
    ValS = SDER(ISADR)

    !! DEPSA can be computed simultaneously with the F3
    !! contributions, but remember that SC here is the overlap,
    !! whereas SC in FG3 subroutines are just G3 matrix
    if (ldc == 0) then
      do iWabs=1,nAshT
        iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
        iSAdr2 = max(iTWV,iXYZ)*(max(iTWV,iXYZ)-1)/2+min(iTWV,iXYZ)
        DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)+ValB*SC(iSAdr2)

        iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
        iSAdr2 = max(iTUV,iXWZ)*(max(iTUV,iXWZ)-1)/2+min(iTUV,iXWZ)
        DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)+ValB*SC(iSAdr2)
      end do
    else
    !! assume SC has all columns (which is reasonable)
    !iTWV = iTabs+nAshT**2*(iVabs-1)
    !ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
    !call DaXpY_(nAshT,ValB,SC2(iSAdr2),nAshT,DEPSA(1,iUabs),1)
    !iXWZ = iXabs+nAshT**2*(iZabs-1)
    !ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
    !call DaXpY_(nAshT,ValB,SC(iSAdr2),LDC*nAshT,DEPSA(1,iYabs),1)
      do iWabs=1,nAshT
        !! we do not have all elements, so use distributed memory
        iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
        ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
        DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)+ValB*SC2(iSAdr2)

        !! we have all elements, so local memory
        iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
        ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
        DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)+ValB*SC(iSAdr2)
      end do
    end if

    !! If non-parallel, overlap is triangular
    !! If parallel, the index is the same to that of SDER
    if (LDC == 0) iSAdr = max(iTUV,iXYZ)*(max(iTUV,iXYZ)-1)/2+min(iTUV,iXYZ)
    DEASUM = DEASUM-ValB*SC(iSAdr)

    ! dyu ( Fvztx - EPSA(u)*Gvztx )
    ! dyu Gvztx
    if (IYABS == IUABS) then
      !VALUE = VALUE+Two*(FP(IP)-EU*PREF(IP))
      DF2(iVabs,iZabs,iTabs,iXabs) = DF2(iVabs,iZabs,iTabs,iXabs)+ValB
      DG2(iVabs,iZabs,iTabs,iXabs) = DG2(iVabs,iZabs,iTabs,iXabs)-EU*ValB

      !VALUE = VALUE+Two*PREF(IP)
      DG2(iVabs,iZabs,iTabs,iXabs) = DG2(iVabs,iZabs,iTabs,iXabs)+ValS
    end if
    DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs) -ValB*G2(iVabs,iZabs,iTabs,iXabs)

    ! dyx ( Fvutz - EPSA(y)*Gvutz )
    ! dyx Gvutz -> dut Gzyxv
    if (IYABS == IXABS) then
      !VALUE = VALUE+Two*(FP(IP)-EY*PREF(IP))
      DF2(iVabs,iUabs,iTabs,iZabs) = DF2(iVabs,iUabs,iTabs,iZabs)+ValB
      DG2(iVabs,iUabs,iTabs,iZabs) = DG2(iVabs,iUabs,iTabs,iZabs)-EY*ValB

      !VALUE = VALUE+Two*PREF(IP)
      DG2(iVabs,iUabs,iTabs,iZabs) = DG2(iVabs,iUabs,iTabs,iZabs)+ValS
    end if
    DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs) -ValB*G2(iVabs,iUabs,iTabs,iZabs)

    ! dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz -
    !        (EPSA(u)+EPSA(y)*dyz Gvz)
    ! dtu Gvxyz + dtu dyx Gvz
    if (ITABS == IUABS) then
      !VALUE = VALUE+Two*(FP(IP)-EU*PREF(IP))
      DF2(iVabs,iXabs,iYabs,iZabs) = DF2(iVabs,iXabs,iYabs,iZabs)+ValB
      DG2(iVabs,iXabs,iYabs,iZabs) = DG2(iVabs,iXabs,iYabs,iZabs)-EU*ValB

      !VALUE = VALUE+Two*PREF(IP)
      DG2(iVabs,iXabs,iYabs,iZabs) = DG2(iVabs,iXabs,iYabs,iZabs)+ValS
      if (IYABS == IXABS) then
        !VALUE = VALUE+FD(ID)-EYU*DREF(ID)
        DF1(iVabs,iZabs) = DF1(iVabs,iZabs)+ValB
        DG1(iVabs,iZabs) = DG1(iVabs,iZabs)-EYU*ValB

        !VALUE = VALUE+DREF((ID1*(ID1-1))/2+ID2)
        DG1(iVabs,iZabs) = DG1(iVabs,iZabs)+ValS
      end if
    end if
    DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs) -ValB*G2(iVabs,iXabs,iYabs,iZabs)
    if (iYabs == iXabs) DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs)-ValB*G1(iVabs,iZabs)
    if (iTabs == iUabs) DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs)-ValB*G1(iVabs,iZabs)
  end do
# ifdef _MOLCAS_MPP_
  if (is_real_par() .and. (IXYZ == iHiS) .and. (iRank /= NPROCS-1)) then
    irank = irank+1
    call GA_Distribution(lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
    call GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
  end if
# endif
end do

return

end subroutine CLagDXC_DP
