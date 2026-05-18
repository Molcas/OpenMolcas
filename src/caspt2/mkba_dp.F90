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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MKBA_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,BA,MBA,iLo,iHi,jLo,jHi,LDA)

use SUPERINDEX, only: MTUV
use caspt2_global, only: ipea_shift
use caspt2_module, only: EASUM, EPSA, NASHT, NTUVES, NTUVES
use Constants, only: Two, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, iSYM, MBA, iLo, iHi, jLo, jHi, LDA
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), FD(NDREF), FP(NPREF)
real(kind=wp), intent(inout) :: BA(MBA)
integer(kind=iwp) :: ID, ID1, ID2, IDT, IDU, IDV, IP, IP1, IP2, ISADR, ITABS, ITUV, ITUVABS, IUABS, IVABS, IVT, IVU, IVZ, IXABS, &
                     IXT, IXYZ, IXYZABS, IXZ, IYABS, IYZ, IZABS
real(kind=wp) :: ET, ETU, EU, EX, EY, FACT, Val

!SV.20100831: fill in the F2 and F1 corrections for this BA block
! on entry, BA should contain SA!!
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
    ISADR = 1+iTUV-iLo+LDA*(iXYZ-jLo)
    if (LDA == 0) then
      if (iXYZ <= iTUV) then
        ISADR = (ITUV*(ITUV-1))/2+IXYZ
      else
        cycle
      end if
    end if
    Val = FACT*BA(ISADR)
    if (IYABS == IUABS) then
      IVZ = IVABS+NASHT*(IZABS-1)
      IXT = IXABS+NASHT*(ITABS-1)
      IP1 = max(IVZ,IXT)
      IP2 = min(IVZ,IXT)
      IP = (IP1*(IP1-1))/2+IP2
      Val = Val-Two*(FP(IP)-EU*PREF(IP))
      if (IXABS == ITABS) then
        ID1 = max(IVABS,IZABS)
        ID2 = min(IVABS,IZABS)
        ID = (ID1*(ID1-1))/2+ID2
        Val = Val+Two*(FD(ID)-ETU*DREF(ID))
      end if
    end if
    ! Add  dyt ( -Fvuxz + Et*Gvuxz +dxu (-Fvz+(Et+Eu)*Gvz))
    if (IYABS == ITABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      IXZ = IXABS+NASHT*(IZABS-1)
      IP1 = max(IVU,IXZ)
      IP2 = min(IVU,IXZ)
      IP = (IP1*(IP1-1))/2+IP2
      Val = Val-Two*(FP(IP)-ET*PREF(IP))
      if (IXABS == IUABS) then
        ID1 = max(IVABS,IZABS)
        ID2 = min(IVABS,IZABS)
        ID = (ID1*(ID1-1))/2+ID2
        Val = Val-(FD(ID)-ETU*DREF(ID))
      end if
    end if
    ! Add  dxu ( -Fvtyz + Eu*Gvtyz )
    if (IXABS == IUABS) then
      IVT = IVABS+NASHT*(ITABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP1 = max(IVT,IYZ)
      IP2 = min(IVT,IYZ)
      IP = (IP1*(IP1-1))/2+IP2
      Val = Val-Two*(FP(IP)-EU*PREF(IP))
    end if
    ! Add  2dtx ( Fvuyz-Et*Gvuyz )
    if (ITABS == IXABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP1 = max(IVU,IYZ)
      IP2 = min(IVU,IYZ)
      IP = (IP1*(IP1-1))/2+IP2
      Val = Val+Four*(FP(IP)-ET*PREF(IP))
    end if
    !GG.Nov03
    if (ITUV == IXYZ) then
      IDT = (ITABS*(ITABS+1))/2
      IDU = (IUABS*(IUABS+1))/2
      IDV = (IVABS*(IVABS+1))/2
      Val = Val+ipea_shift*Half*BA(ISADR)*(Two-DREF(IDV)+DREF(IDT)+DREF(IDU))
    end if
    !GG End
    BA(ISADR) = Val
  end do
end do

end subroutine MKBA_DP
