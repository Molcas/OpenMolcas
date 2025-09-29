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

!#define _DEBUGPRINT_
subroutine SDCMRF(CSD,CCM,IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
! Change a block of coefficients bwtween combination format and
! Slater determinant format
!
! Combinations => SD
!
! Input
! =====
! CSD : Block in determinant form
! CCM : Block in combination form
! IATP,IBTP : type of alpha- and beta- string
! NA,NB : Number of alpha- and beta- strings
! IDC  : Combination type
! PS   : Spin combination sign
! PL   : Ml   combination sign
! ISGVST : Ml reflection of strings
!
! If ISCALE == 0, no overall scaling is performed,
!                 the overall scale factor is returned as SCLFAC

use Index_Functions, only: nTri_Elem
use Constants, only: One, Two, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CSD(*)
real(kind=wp), intent(in) :: CCM(*), PS, PL
integer(kind=iwp), intent(in) :: IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST(IBSM), ISCALE
integer(kind=iwp), intent(out) :: LDET, LCOMB
real(kind=wp), intent(out) :: SCLFAC
integer(kind=iwp) :: IPACK
real(kind=wp) :: FACTOR, SGN
logical(kind=iwp) :: Test1, Test2
real(kind=wp), parameter :: SQRT2 = sqrt(Two), SQRT2I = sqrt(Half)

! Is combination array packed ?

SCLFAC = One
IPACK = 0
FACTOR = One

Test1 = (IDC == 2) .or. (IDC == 4)
if (.not. Test1) then
  Test2 = IDC == 4
  if (Test2) Test2 = Test2 .and. (IASM == ISGVST(IBSM))
else
  Test2 = .false.
end if

if (Test1) then
  SGN = PS
  FACTOR = SQRT2
  if ((IASM == IBSM) .and. (IATP == IBTP)) IPACK = 1
else if (Test2) then
  if (IATP == IBTP) IPACK = 1
  SGN = PS*PL
  FACTOR = Two
end if

LDET = NA*NB
#ifdef _DEBUGPRINT_
write(u6,*) ' SDCMRF : NA, NB =',NA,NB
#endif
if (IPACK == 0) then
  LCOMB = LDET
else
  LCOMB = nTri_Elem(NA)
end if
if ((IDC == 4) .and. (IPACK == 0)) FACTOR = SQRT2
FACTOR = One/FACTOR

! Combination => SD transformation

if (IPACK == 1) then
  ! Unpack from triangular form
  call TRIPK32(CSD,CCM,NA,NA,SGN)
else
  CSD(1:NA*NB) = CCM(1:NA*NB)
end if
! Scale
if (FACTOR /= One) then
  if (ISCALE == 1) then
    SCLFAC = One
    CSD(1:LDET) = FACTOR*CSD(1:LDET)
  else
    SCLFAC = FACTOR
  end if
  if (IPACK == 1) call SCLDIA(CSD,SQRT2,NA,0)
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Information from SDCMRF'

write(u6,'(A,6I4)') ' IATP IBTP IASM IBSM IDC ',IATP,IBTP,IASM,IBSM,IDC
write(u6,'(A,I4,3X,2ES15.8)') ' IPACK FACTOR SIGN',IPACK,FACTOR,SGN
write(u6,*) ' Slater determinant block'
call WRTMAT(CSD,NA,NB,NA,NB)
write(u6,*)
write(u6,*) ' Combination block'
if (IPACK == 1) then
  call PRSM2(CCM,NA)
else
  call WRTMAT(CCM,NA,NB,NA,NB)
end if
#endif

end subroutine SDCMRF
