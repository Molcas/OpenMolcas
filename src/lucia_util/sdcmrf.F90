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

subroutine SDCMRF(CSD,CCM,IWAY,IATP,IBTP,IASM,IBSM,NA,NB,IDC,PS,PL,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
! Change a block of coefficients bwtween combination format and
! Slater determinant format
!
! IWAY = 1 : SD => Combinations
! IWAY = 2 : Combinations => SD
!
! Input
! =====
! CSD : Block in determinant form
! CCM : Block in combination  form
! IWAY : as above
! IATP,IBTP : type of alpha- and beta- string
! NA,NB : Number of alpha- and beta- strings
! IDC  : Combination type
! PS   : Spin combination sign
! PL   : Ml   combination sign
! ISGVST : Ml reflection of strings
!
! If ISCALE == 0, no overall scaling is performed,
!                 the overall scale factor is returned as SCLFAC

use Constants, only: One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CSD(*), CCM(*), PS, PL, SCLFAC
integer(kind=iwp) :: IWAY, IATP, IBTP, IASM, IBSM, NA, NB, IDC, ISGVST(*), LDET, LCOMB, ISCALE
integer(kind=iwp) :: IPACK, NTEST
real(kind=wp) :: FACTOR, SGN
logical(kind=iwp) :: Test1, Test2
real(kind=wp), parameter :: SQRT2 = sqrt(Two), SQRT2I = sqrt(Half)

NTEST = 0

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
if (NTEST >= 100) write(u6,*) ' SDCMRF : NA, NB =',NA,NB
if (IPACK == 0) then
  LCOMB = LDET
else
  LCOMB = NA*(NA+1)/2
end if
if ((IDC == 4) .and. (IPACK == 0)) FACTOR = SQRT2
if (IWAY == 2) FACTOR = One/FACTOR

! SD => combination transformation

if (IWAY == 1) then
  if (IPACK == 1) then
    ! Pack to triangular form
    call TRIPK3(CSD,CCM,1,NA,NA,SGN)
    !    TRIPK3(AUTPAK,APAK,IWAY,MATDIM,NDIM,SGN)
  else
    CCM(1:NA*NB) = CSD(1:NA*NB)
  end if
  ! Scale
  if (FACTOR /= One) then
    if (ISCALE == 1) then
      SCLFAC = One
      CCM(1:LCOMB) = FACTOR*CCM(1:LCOMB)
    else
      SCLFAC = FACTOR
    end if
    if (IPACK == 1) call SCLDIA(CCM,SQRT2I,NA,1)
  end if
end if

! Combination => SD transformation

if (IWAY == 2) then
  if (IPACK == 1) then
    ! Unpack from triangular form
    call TRIPK3(CSD,CCM,2,NA,NA,SGN)
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
end if

NTEST = 0
!if ((NTEST /= 0) .and. (IWAY == 1)) then
if (NTEST /= 0) then
  write(u6,*) ' Information from SDCMRF'

  write(u6,'(A,6I4)') ' IWAY IATP IBTP IASM IBSM IDC ',IWAY,IATP,IBTP,IASM,IBSM,IDC
  write(u6,'(A,I4,3X,2ES15.8)') ' IPACK FACTOR SIGN',IPACK,FACTOR,SGN
  if (NTEST >= 100) then
    write(u6,*) ' Slater determinant block'
    call WRTMAT(CSD,NA,NB,NA,NB)
    write(u6,*)
    write(u6,*) ' Combination block'
    if (IPACK == 1) then
      call PRSM2(CCM,NA)
    else
      call WRTMAT(CCM,NA,NB,NA,NB)
    end if
  end if
end if

end subroutine SDCMRF
