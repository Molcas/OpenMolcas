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
function IBASSPC_FOR_CLS(ICLS)
! Obtain base space for occupation class ICLS

use lucia_data, only: ICMBSPC, IGSOCCX, LCMBSPC, NCMBSPC, NCMBSPC, NGAS
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IBASSPC_FOR_CLS
integer(kind=iwp), intent(in) :: ICLS(NGAS)
integer(kind=iwp) :: I_AM_OKAY, IBASE, IGAS, ISPC, JJSPC, JSPC, NEL

! Some dummy initializations

NEL = 0 ! jwk-cleanup

IBASE = 0
do ISPC=1,NCMBSPC
  do JJSPC=1,LCMBSPC(ISPC)
    JSPC = ICMBSPC(JJSPC,ISPC)
    ! Test for occupation constraints in CI space JSPC
    I_AM_OKAY = 1
    do IGAS=1,NGAS
      if (IGAS == 1) then
        NEL = ICLS(IGAS)
      else
        NEL = NEL+ICLS(IGAS)
      end if

      if ((NEL < IGSOCCX(IGAS,1,JSPC)) .or. (NEL > IGSOCCX(IGAS,2,JSPC))) I_AM_OKAY = 0
    end do
    ! End of loop over gasspaces for given cispace

    if ((I_AM_OKAY == 1) .and. (IBASE == 0)) IBASE = ISPC

  end do
  ! End of loop over cisspaces for given combination space
end do
! End of loop over combinations apaces

IBASSPC_FOR_CLS = IBASE

#ifdef _DEBUGPRINT_
write(u6,*) ' Occupation class and its basespace'
call IWRTMA(ICLS,1,NGAS,1,NGAS)
write(u6,*) IBASE
#endif

end function IBASSPC_FOR_CLS
