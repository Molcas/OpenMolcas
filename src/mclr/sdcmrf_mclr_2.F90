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

subroutine SDCMRF_MCLR_2(CSD,CCM,IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB)
! Change a block of coefficients bwtween combination format and
! Slater determinant format
!
!     Combinations => SD
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

use Constants, only: One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CSD(*)
real(kind=wp), intent(in) :: CCM(*), PS, PL
integer(kind=iwp), intent(in) :: IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST
integer(kind=iwp), intent(out) :: LDET, LCOMB
integer(kind=iwp) :: IPACK
real(kind=wp) :: FACTOR, SGN
real(kind=wp), parameter :: SQRT2 = sqrt(Two)

call SDCMRF_PRE(IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB,IPACK,FACTOR,SGN)

! Combination => SD transformation

if (IPACK == 1) then
  ! Unpack from triangular form
  call TRIPK2_2(CSD,CCM,NA,NA,SGN)
else
  CSD(1:NA*NB) = CCM(1:NA*NB)
end if
! Scale
if (FACTOR /= One) then
  CSD(1:LDET) = CSD(1:LDET)/FACTOR
  if (IPACK == 1) call SCLDIA(CSD,SQRT2,NA,0)
end if

end subroutine SDCMRF_MCLR_2
