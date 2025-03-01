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
subroutine NXT_CONF(ICONF,NEL,NORB,INI,NONEW)
! Next configuration of NEL electrons distributed in NORB orbitals
!
! A configuration is stored as the occupied orbitals
! in nonstrict ascending order - two consecutive orbitals are allowed
! to be identical
! allowing two
!
! IF INI = 1 : Generate initial configuration
!    NONOEW = 1 : No new configuration could be generated
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NEL, NORB, INI
integer(kind=iwp), intent(inout) :: ICONF(NEL)
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: I, IADD, IEL, INCREASE, JORB, N_DOUBLE, NDOUBLE

#ifdef _DEBUGPRINT_
write(u6,*) ' Input configuration to NXT_CONF'
call IWRTMA(ICONF,1,NEL,1,NEL)
write(u6,*) ' NEL, NORB = ',NEL,NORB
#endif
if (INI == 1) then
  ! Check that NEL electrons can be distributed in NORB orbitals
  if (NEL <= 2*NORB) then
    NONEW = 0
    N_DOUBLE = NEL/2
    do I=1,N_DOUBLE
      ICONF(2*I-1) = I
      ICONF(2*I) = I
    end do
    if (2*N_DOUBLE /= NEL) ICONF(2*N_DOUBLE+1) = N_DOUBLE+1
  else
    NONEW = 1
  end if

else if (INI == 0) then

  IADD = 1
  IEL = 0
  ! Increase orbital number of next electron
  do while (IADD == 1)
    IEL = IEL+1
    ! Can orbital number be increased for electron IEL ?
    INCREASE = 0
    if (IEL < NEL) then
      if (ICONF(IEL) < ICONF(IEL+1)-1) then
        INCREASE = 1
      else if (ICONF(IEL) == ICONF(IEL+1)-1) then
        ! If ICONF(IEL) is increased, ICONF(IEL) = ICONF(IEL+1), check if this is ok
        if (IEL == NEL-1) then
          INCREASE = 1
        else if (ICONF(IEL+1) /= ICONF(IEL+2)) then
          INCREASE = 1
        end if
      end if
    else if (IEL == NEL) then
      if (ICONF(IEL) < NORB) then
        INCREASE = 1
      else
        ! Nothing more to do
        NONEW = 1
        exit
      end if
    else
      ! Nothing more to do
      NONEW = 1
      exit
    end if

    if (INCREASE == 1) then
      ! Increase orbital for elec IEL
      NONEW = 0
      ICONF(IEL) = ICONF(IEL)+1
      ! Minimize orbital occupations
      NDOUBLE = (IEL-1)/2
      do JORB=1,NDOUBLE
        ICONF(2*JORB-1) = JORB
        ICONF(2*JORB) = JORB
      end do
      if (2*NDOUBLE < IEL-1) ICONF(IEL-1) = NDOUBLE+1
      IADD = 0
    end if
  end do
end if
! End if INI = 0

#ifdef _DEBUGPRINT_
if (NONEW == 1) then
  write(u6,*) ' No new configurations'
  write(u6,*) ' Input configuration'
  call IWRTMA(ICONF,1,NEL,1,NEL)
else
  write(u6,*) ' Next configurations'
  call IWRTMA(ICONF,1,NEL,1,NEL)
end if
#endif

end subroutine NXT_CONF
