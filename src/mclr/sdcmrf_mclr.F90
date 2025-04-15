!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine SDCMRF_MCLR(CSD,CCM,IWAY,IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB)
! Change a block of coefficients bwtween combination format and
! Slater determinant format
!
!     IWAY = 1 : SD => Combinations
!     IWAY = 2 : Combinations => SD
!
! Input
! =====
! CSD       : Block in determinant form
! CCM       : Block in combination  form
! IWAY      : as above
! IATP,IBTP : type of alpha- and beta- string
! NA,NB     : Number of alpha- and beta- strings
! IDC       : Combination type
! PS        : Spin combination sign
! PL        : Ml   combination sign
! ISGVST    : Ml reflection of strings

use Index_Functions, only: nTri_Elem
use Constants, only: Two, One, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: CSD(*), CCM(*), PS, PL
integer(kind=iwp) :: IWAY, IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST(*), LDET, LCOMB
integer(kind=iwp) :: IPACK
real(kind=wp) :: FACTOR, SGN
real(kind=wp), parameter :: SQRT2 = sqrt(Two), SQRT2I = sqrt(Half)

! Is combination array packed ?

IPACK = 0
FACTOR = One

if (((IDC == 2) .or. (IDC == 4)) .and. (IASM == IBSM)) then
  SGN = PS
  FACTOR = SQRT2
  if (IATP == IBTP) IPACK = 1
else if ((IDC == 4) .and. (IASM == ISGVST(IBSM))) then
  if (IATP == IBTP) IPACK = 1
  SGN = PS*PL
  FACTOR = Two
end if

LDET = NA*NB
if (IPACK == 0) then
  LCOMB = LDET
else
  LCOMB = nTri_Elem(NA)
end if
if ((IDC == 4) .and. (IPACK == 0)) FACTOR = SQRT2
if (IWAY == 2) FACTOR = One/FACTOR

! SD => combination transformation

if (IWAY == 1) then
  if (IPACK == 1) then
    ! Pack to triangular form
    call TRIPK2(CSD,CCM,1,NA,NA,SGN)
    !    TRIPK2(AUTPAK,APAK,IWAY,MATDIM,NDIM,SGN)
  else
    CCM(1:NA*NB) = CSD(1:NA*NB)
  end if
  ! Scale
  if (FACTOR /= One) then
    CCM(1:LCOMB) = FACTOR*CCM(1:LCOMB)
    if (IPACK == 1) call SCLDIA(CCM,SQRT2I,NA,1)
  end if
end if

! Combination => SD transformation

if (IWAY == 2) then
  if (IPACK == 1) then
    ! Unpack from triangular form
    call TRIPK2(CSD,CCM,2,NA,NA,SGN)
  else
    CSD(1:NA*NB) = CCM(1:NA*NB)
  end if
  ! Scale
  if (FACTOR /= One) then
    CSD(1:LDET) = FACTOR*CSD(1:LDET)
    if (IPACK == 1) call SCLDIA(CSD,SQRT2,NA,0)
  end if
end if

return

end subroutine SDCMRF_MCLR
