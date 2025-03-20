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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB)
! Reform between two ways of writing occupations
!
! IOCC_EXP : Occupation in expanded form, i.e. the orbital for each
!            electron is given
!
! IOCC_PCK  : Occupation is given in packed form, i.e. each occupied
!             orbitals is given once, and a negative index indicates
!             a double occupation
!
! Expanded to Packed form
!
! Jeppe Olsen, Nov. 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NEL, IOCC_EXP(NEL), NOCOB
integer(kind=iwp), intent(out) :: IOCC_PCK(NOCOB)
integer(kind=iwp) :: IEL, IOCC

! Expanded => Packed form

! Loop over electrons
IEL = 1
IOCC = 0
do
  !IEL = IEL+1
  if (IEL < NEL) then
    if (IOCC_EXP(IEL) == IOCC_EXP(IEL+1)) then
      IOCC = IOCC+1
      IOCC_PCK(IOCC) = -IOCC_EXP(IEL)
      IEL = IEL+2
    else
      IOCC = IOCC+1
      IOCC_PCK(IOCC) = IOCC_EXP(IEL)
      IEL = IEL+1
    end if
  else
    ! Last occupation was not identical to previous, so single occupied
    IOCC = IOCC+1
    IOCC_PCK(IOCC) = IOCC_EXP(IEL)
    IEL = IEL+1
  end if
  if (IEL > NEL) exit
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Reforming form of configuration'
write(u6,*) ' Expanded to packed form'
write(u6,*) ' IOCC_EXP :'
call IWRTMA(IOCC_EXP,1,NEL,1,NEL)
write(u6,*) ' IOCC_PCK :'
call IWRTMA(IOCC_PCK,1,NOCOB,1,NOCOB)
#endif

end subroutine REFORM_CONF_OCC
