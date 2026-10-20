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

subroutine MKSC_DP(DREF,NDREF,PREF,NPREF,iSYM,SC,NSC,iLo,iHi,jLo,jHi,LDC)
! In parallel, this subroutine is called on a local chunk of memory
! and LDC is set. In serial, the whole array is passed but then the
! storage uses a triangular scheme, and the LDC passed is zero.

use Index_Functions, only: iTri
use SUPERINDEX, only: MTUV
use caspt2_module, only: NASHT, nTUVES
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, iSYM, NSC, iLo, iHi, jLo, jHi, LDC
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(inout) :: SC(NSC)
integer(kind=iwp) :: ID, IP, ISADR, ITABS, ITUV, ITUVABS, ITX, ITZ, IUABS, IVABS, IVU, IVX, IVZ, IXABS, IXYZ, IXYZABS, IYABS, IYZ, &
                     IZABS
real(kind=wp) :: Val

ISADR = 0
!-SVC20100831: fill in the G2 and G1 corrections for this SC block
do IXYZ=jLo,jHi
  IXYZABS = IXYZ+NTUVES(ISYM)
  IXABS = MTUV(1,IXYZABS)
  IYABS = MTUV(2,IXYZABS)
  IZABS = MTUV(3,IXYZABS)
  do ITUV=iLo,iHi
    ITUVABS = ITUV+NTUVES(ISYM)
    ITABS = MTUV(1,ITUVABS)
    IUABS = MTUV(2,ITUVABS)
    IVABS = MTUV(3,ITUVABS)
    if (LDC /= 0) then
      Val = SC(1+iTUV-iLo+LDC*(iXYZ-jLo))
    else
      if (IXYZ <= ITUV) then
        ISADR = iTri(ITUV,IXYZ)
        Val = SC(ISADR)
      else
        cycle
      end if
    end if
    ! Add  dyu Gvztx
    if (IYABS == IUABS) then
      IVZ = IVABS+NASHT*(IZABS-1)
      ITX = ITABS+NASHT*(IXABS-1)
      IP = iTri(IVZ,ITX)
      Val = Val+Two*PREF(IP)
    end if
    ! Add  dyx Gvutz
    if (IYABS == IXABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      ITZ = ITABS+NASHT*(IZABS-1)
      IP = iTri(IVU,ITZ)
      Val = Val+Two*PREF(IP)
    end if
    ! Add  dtu Gvxyz + dtu dyx Gvz
    if (ITABS == IUABS) then
      IVX = IVABS+NASHT*(IXABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP = iTri(IVX,IYZ)
      Val = Val+Two*PREF(IP)
      if (IYABS == IXABS) then
        ID = iTri(IVABS,IZABS)
        Val = Val+DREF(ID)
      end if
    end if
    if (LDC /= 0) then
      SC(1+iTUV-iLo+LDC*(iXYZ-jLo)) = Val
    else
      SC(ISADR) = Val
    end if
  end do
end do

end subroutine MKSC_DP
