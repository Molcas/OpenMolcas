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

subroutine DIELMV(ICASE,nICASE,JCASE,nJCASE,NUP,NDWN,EMU)

use sguga, only: CIS, SGS
use caspt2_module, only: ETA
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nICASE, ICASE(nICASE), nJCASE, JCASE(nJCASE), NUP, NDWN
real(kind=wp), intent(inout) :: EMU(NUP,NDWN)
integer(kind=iwp) :: I, IC, IC1, II, IOC, ISTEP, LEV, LV1, nIpWlk, nLev
real(kind=wp) :: rSUM

nLev = SGS%nLev
nIpWlk = CIS%nIpWlk

do I=1,NUP
  II = NIPWLK*(I-1)
  rSUM = Zero
  do LV1=SGS%MIDLEV+1,NLEV,15
    II = II+1
    IC = ICASE(II)
    do LEV=LV1,min(LV1+14,NLEV)
      IC1 = IC/4
      ISTEP = IC-4*IC1
      IOC = (ISTEP+1)/2
      rSUM = rSUM+real(IOC,kind=wp)*ETA(LEV)
      IC = IC1
    end do
  end do
  EMU(I,:) = EMU(I,:)+rSUM
end do
! THEN THE LOWER HALF:
do I=1,NDWN
  II = NIPWLK*(I-1)
  rSUM = Zero
  do LV1=1,SGS%MIDLEV,15
    II = II+1
    IC = JCASE(II)
    do LEV=LV1,min(LV1+14,SGS%MIDLEV)
      IC1 = IC/4
      ISTEP = IC-4*IC1
      IOC = (ISTEP+1)/2
      rSUM = rSUM+real(IOC,kind=wp)*ETA(LEV)
      IC = IC1
    end do
  end do
  EMU(:,I) = EMU(:,I)+rSUM
end do

end subroutine DIELMV
