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

subroutine MKSA_DP(DREF,NDREF,PREF,NPREF,iSYM,SA,NSA,iLo,iHi,jLo,jHi,LDA)
! In parallel, this subroutine is called on a local chunk of memory
! and LDA is set. In serial, the whole array is passed but then the
! storage uses a triangular scheme, and the LDA passed is zero.

use Index_Functions, only: iTri
use SUPERINDEX, only: MTUV
use caspt2_module, only: NASHT, nTUVES
use Constants, only: Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NDREF, NPREF, iSYM, NSA, iLo, iHi, jLo, jHi, LDA
real(kind=wp), intent(in) :: DREF(NDREF), PREF(NPREF)
real(kind=wp), intent(inout) :: SA(NSA)
integer(kind=iwp) :: ID, IP, ISADR, ITABS, ITUV, ITUVABS, IUABS, IVABS, IVT, IVU, IVZ, IXABS, IXT, IXYZ, IXYZABS, IXZ, IYABS, IYZ, &
                     IZABS
real(kind=wp) :: Val

ISADR = 0
!-SVC20100831: fill in the G2 and G1 corrections for SA
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
    ! Add  2 dtx Gvuyz + 2 dtx dyu Gvz
    if (LDA /= 0) then
      Val = SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo))
    else if (IXYZ > ITUV) then
      cycle
    else
      ISADR = iTri(ITUV,IXYZ)
      Val = SA(ISADR)
    end if
    if (ITABS == IXABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP = iTri(IVU,IYZ)
      Val = Val+Four*PREF(IP)
      if (IYABS == IUABS) then
        ID = iTri(IVABS,IZABS)
        Val = Val+Two*DREF(ID)
      end if
    end if
    ! Add  -dxu Gvtyz -dxu dyt Gvz
    if (IXABS == IUABS) then
      IVT = IVABS+NASHT*(ITABS-1)
      IYZ = IYABS+NASHT*(IZABS-1)
      IP = iTri(IVT,IYZ)
      Val = Val-Two*PREF(IP)
      if (IYABS == ITABS) then
        ID = iTri(IVABS,IZABS)
        Val = Val-DREF(ID)
      end if
    end if
    ! Add  -dyt Gvuxz
    if (IYABS == ITABS) then
      IVU = IVABS+NASHT*(IUABS-1)
      IXZ = IXABS+NASHT*(IZABS-1)
      IP = iTri(IVU,IXZ)
      Val = Val-Two*PREF(IP)
    end if
    ! Add -dyu Gvzxt
    if (IYABS == IUABS) then
      IVZ = IVABS+NASHT*(IZABS-1)
      IXT = IXABS+NASHT*(ITABS-1)
      IP = iTri(IVZ,IXT)
      Val = Val-Two*PREF(IP)
    end if
    if (LDA /= 0) then
      SA(1+(iTUV-iLo)+LDA*(iXYZ-jLo)) = Val
    else
      SA(ISADR) = Val
    end if
  end do
end do

end subroutine MKSA_DP
