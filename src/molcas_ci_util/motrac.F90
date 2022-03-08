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

subroutine MOTRAC(CMO,F,X1,X2)
! RASSCF PROGRAM IBM-3090 VERSION: CI SECTION
! PURPOSE: TRANSFORMS A MATRIX F WITH THE TRANSFORMATION MATRIX CMO
!          NORMALLY FROM AO TO MO BASIS. THE RESULT OVERWRITES THE
!          INITIAL MATRIX. SYMMETRY BLOCKED.
!          THE MATRIX F IS TRANSFORMED FROM FULL AO BASIS TO
!          ACTIVE ORBITAL MO BASIS.
!          INPUT AND OUTPUT MATRICES IN LOWER TRIANGULAR FORM.
!          CALLED FROM FCIN
!
! ********** IBM-3090 RELEASE 86 12 05 **********

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: CMO(*)
real(kind=wp), intent(inout) :: F(*)
real(kind=wp), intent(out) :: X1(*), X2(*)
integer(kind=iwp) :: ISTFA, ISTFP, ISYM, LMOP, LMOP1, NA, NB
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

LMOP = 1
ISTFA = 1
ISTFP = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  NA = NASH(ISYM)
  LMOP1 = LMOP+NB*(NISH(ISYM)+NFRO(ISYM))
  if (NA /= 0) then
    call SQUARE(F(ISTFP),X1,1,NB,NB)
    call DGEMM_('N','N',NB,NA,NB,One,X1,NB,CMO(LMOP1),NB,Zero,X2,NB)
    call DGEMM_Tri('T','N',NA,NA,NB,One,X2,NB,CMO(LMOP1),NB,Zero,F(ISTFA),NA)

    ISTFA = ISTFA+ITRI(NA+1)
  end if
  LMOP = LMOP+NB**2
  ISTFP = ISTFP+ITRI(NB+1)
end do

return

end subroutine MOTRAC
