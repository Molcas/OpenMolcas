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

use Index_Functions, only: iTri, nTri_Elem
use SUPERINDEX, only: MTUV
use caspt2_global, only: ipea_shift
use caspt2_module, only: EASUM, EPSA, NASHT, NTUVES
use Constants, only: Two, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, iSYM, MBC, iLo, iHi, jLo, jHi, LDC
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF), FD(NDREF), FP(NPREF)
real(kind=wp), intent(inout) :: BC(MBC)
integer(kind=iwp) :: ID, IDT, IDU, IDV, IP, ISADR, ITABS, ITUV, ITUVABS, ITX, ITZ, IUABS, IVABS, IVU, IVX, IVZ, IXABS, IXYZ, &
                     IXYZABS, IYABS, IYZ, IZABS
real(kind=wp) :: EU, EY, EYU, FACT, Val

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
      if (iXYZ > iTUV) cycle
      ISADR = iTri(ITUV,IXYZ)
    end if
    Val = FACT*BC(ISADR)
    !Val = Fvutxyz+(EPSA(y)+EPSA(u))*SC(tuv,xyz)
    ! Add  dyu ( Fvztx - EPSA(u)*Gvztx )
    if (IYABS == IUABS) then
      IVZ = IVABS+NASHT*(IZABS-1)
      ITX = ITABS+NASHT*(IXABS-1)
      IP = iTri(IVZ,ITX)
      Val = Val+Two*(FP(IP)-EU*PREF(IP))
    end if
    ! Add  dyx ( Fvutz - EPSA(y)*Gvutz )
    if (IYABS == IXABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      ITZ = ITABS+NASHT*(IZABS-1)
      IP = iTri(IVU,ITZ)
      Val = Val+Two*(FP(IP)-EY*PREF(IP))
    end if
    ! Add  dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz - (EPSA(u)+EPSA(y)*dyz Gvz)
    if (ITABS == IUABS) then
      IVX = IVABS+NASHT*(IXABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP = iTri(IVX,IYZ)
      Val = Val+Two*(FP(IP)-EU*PREF(IP))
      if (IYABS == IXABS) then
        ID = iTri(IVABS,IZABS)
        Val = Val+FD(ID)-EYU*DREF(ID)
      end if
    end if
    !GG.Nov03
    if (ITUV == IXYZ) then
      IDT = nTri_Elem(ITABS)
      IDU = nTri_Elem(IUABS)
      IDV = nTri_Elem(IVABS)
      Val = Val+ipea_shift*Half*BC(ISADR)*(Four-DREF(IDT)-DREF(IDV)+DREF(IDU))
    end if
    !GG End
    BC(ISADR) = Val
  end do
end do

end subroutine MKBC_DP
