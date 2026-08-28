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

subroutine SDCMRF_MCLR_1(CSD,CCM,IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB)
! Change a block of coefficients bwtween combination format and
! Slater determinant format
!
!     SD => Combinations
!
! Input
! =====
! CSD       : Block in determinant form
! CCM       : Block in combination form
! IATP,IBTP : type of alpha- and beta- string
! NA,NB     : Number of alpha- and beta- strings
! IDC       : Combination type
! PS        : Spin combination sign
! PL        : Ml   combination sign
! ISGVST    : Ml reflection of strings

use Constants, only: One, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CSD(*), PS, PL
real(kind=wp), intent(_OUT_) :: CCM(*)
integer(kind=iwp), intent(in) :: IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST
integer(kind=iwp), intent(out) :: LDET, LCOMB
integer(kind=iwp) :: IPACK
real(kind=wp) :: FACTOR, SGN
real(kind=wp), parameter :: SQRT2I = sqrt(Half)

call SDCMRF_PRE(IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB,IPACK,FACTOR,SGN)

! SD => combination transformation

if (IPACK == 1) then
  ! Pack to triangular form
  call TRIPK2_1(CSD,CCM,NA,NA)
  !    TRIPK2_1(AUTPAK,APAK,MATDIM,NDIM)
else
  CCM(1:NA*NB) = CSD(1:NA*NB)
end if
! Scale
if (FACTOR /= One) then
  CCM(1:LCOMB) = FACTOR*CCM(1:LCOMB)
  if (IPACK == 1) call SCLDIA(CCM,SQRT2I,NA,1)
end if

end subroutine SDCMRF_MCLR_1
