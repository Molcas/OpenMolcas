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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MLTSCA_DH(IMLTOP,LST1,LST2,X,NXI,NXA,F,NFI,NFA,Y,NAS2,jYLo,jYHi)
! Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L11,L12,L13,L14 for short, for a given
! item, and with V1=VAL1(L14), and similar for the other list,
! compute, for IMLTOP=0 or 1 respectively,
!     X(L11,L21) := Add V1*V2*F(L12,L22)*Y(L13,L23)
!  or Y(L13,L23) := Add V1*V2*F(L12,L22)*X(L11,L21)
! or for IMLTOP=2, compute
!     F(L12,L22) := Add V1*V2*X(L11,L21)*Y(L13,L23)

use Sigma_data, only: NLST1, NLST2, VAL1, VAL2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, LST1(4,NLST1), LST2(4,NLST2), NXI, NXA, NFI, NFA, NAS2, jYLo, jYHi
real(kind=wp), intent(inout) :: X(NXI,NXA), F(NFI,NFA), Y(NAS2,jYHi-jYLo+1)
integer(kind=iwp) :: ILST1, ILST2, JY, L11, L12, L13, L14, L21, L22, L23, L24
real(kind=wp) :: V1, V2

if (IMLTOP == 0) then
  do ILST1=1,NLST1
    L11 = LST1(1,ILST1)
    L12 = LST1(2,ILST1)
    L13 = LST1(3,ILST1)
    L14 = LST1(4,ILST1)
    V1 = VAL1(L14)
    if ((L13 >= jYLo) .and. (L13 <= jYHi)) then
      JY = L13-jYLo+1
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        X(L11,L21) = X(L11,L21)+V1*V2*F(L12,L22)*Y(L23,JY)
      end do
    end if
  end do
else if (IMLTOP == 1) then
  do ILST1=1,NLST1
    L11 = LST1(1,ILST1)
    L12 = LST1(2,ILST1)
    L13 = LST1(3,ILST1)
    L14 = LST1(4,ILST1)
    V1 = VAL1(L14)
    if ((L13 >= jYLo) .and. (L13 <= jYHi)) then
      JY = L13-jYLo+1
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        Y(L23,JY) = Y(L23,JY)+V1*V2*F(L12,L22)*X(L11,L21)
      end do
    end if
  end do
else
  do ILST1=1,NLST1
    L11 = LST1(1,ILST1)
    L12 = LST1(2,ILST1)
    L13 = LST1(3,ILST1)
    L14 = LST1(4,ILST1)
    V1 = VAL1(L14)
    if ((L13 >= jYLo) .and. (L13 <= jYHi)) then
      JY = L13-jYLo+1
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        F(L12,L22) = F(L12,L22)+V1*V2*X(L11,L21)*Y(L23,JY)
      end do
    end if
  end do
end if

!NFSCA = NFSCA+4*NLST1*NLST2

end subroutine MLTSCA_DH
