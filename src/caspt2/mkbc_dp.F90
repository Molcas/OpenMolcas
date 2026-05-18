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

subroutine MKBC_DP(DREF,NDREF,PREF,NPREF,FD,FP,iSYM,BC,MBC,iLo,iHi,jLo,jHi,LDC)

use definitions, only: iwp, wp
use constants, only: Half, Two, Four
use SUPERINDEX, only: MTUV
use caspt2_global, only: ipea_shift
use caspt2_module, only: EASUM, NASHT, NTUVES, EPSA

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, iSYM, MBC, iLo, iHi, jLo, jHi, LDC
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(in) :: FD(NDREF), FP(NPREF)
real(kind=wp), intent(inout) :: BC(MBC)
integer(kind=iwp) IXYZ, IXYZABS, IXABS, IYABS, IZABS, ITUV, ITUVABS, ITABS, IUABS, IVABS, ISADR, IVZ, ITX, IVU, ITZ, IP1, IP2, IP, &
                  ID1, ID2, ID, IDT, IDU, IDV, IVX, IYZ
real(kind=wp) EY, EU, EYU, FACT, value

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
    ISADR = 1+iTUV-iLo+LDC*(iXYZ-jLo)
    if (LDC == 0) then
      if (iXYZ <= iTUV) then
        ISADR = (ITUV*(ITUV-1))/2+IXYZ
      else
        cycle
      end if
    end if
    value = FACT*BC(ISADR)
    !value = Fvutxyz+(EPSA(y)+EPSA(u))*SC(tuv,xyz)
    ! Add  dyu ( Fvztx - EPSA(u)*Gvztx )
    if (IYABS == IUABS) then
      IVZ = IVABS+NASHT*(IZABS-1)
      ITX = ITABS+NASHT*(IXABS-1)
      IP1 = max(IVZ,ITX)
      IP2 = min(IVZ,ITX)
      IP = (IP1*(IP1-1))/2+IP2
      value = value+Two*(FP(IP)-EU*PREF(IP))
    end if
    ! Add  dyx ( Fvutz - EPSA(y)*Gvutz )
    if (IYABS == IXABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      ITZ = ITABS+NASHT*(IZABS-1)
      IP1 = max(IVU,ITZ)
      IP2 = min(IVU,ITZ)
      IP = (IP1*(IP1-1))/2+IP2
      value = value+Two*(FP(IP)-EY*PREF(IP))
    end if
    ! Add  dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz - (EPSA(u)+EPSA(y)*dyz Gvz)
    if (ITABS == IUABS) then
      IVX = IVABS+NASHT*(IXABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP1 = max(IVX,IYZ)
      IP2 = min(IVX,IYZ)
      IP = (IP1*(IP1-1))/2+IP2
      value = value+Two*(FP(IP)-EU*PREF(IP))
      if (IYABS == IXABS) then
        ID1 = max(IVABS,IZABS)
        ID2 = min(IVABS,IZABS)
        ID = (ID1*(ID1-1))/2+ID2
        value = value+FD(ID)-EYU*DREF(ID)
      end if
    end if
    !GG.Nov03
    if (ITUV == IXYZ) then
      IDT = (ITABS*(ITABS+1))/2
      IDU = (IUABS*(IUABS+1))/2
      IDV = (IVABS*(IVABS+1))/2
      value = value+ipea_shift*Half*BC(ISADR)*(Four-DREF(IDT)-DREF(IDV)+DREF(IDU))
    end if
    !GG End
    BC(ISADR) = value
  end do
end do

end subroutine MKBC_DP
