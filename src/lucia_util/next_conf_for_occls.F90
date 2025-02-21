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
subroutine NEXT_CONF_FOR_OCCLS(ICONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
! Obtain next configuration for occupation class
!
! Jeppe Olsen, Nov. 2001

use lucia_data, only: MXPNGAS, MXPORB
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(inout) :: ICONF(*)
integer(kind=iwp), intent(in) :: NGAS, IOCCLS(NGAS), NOBPT(NGAS), INI
integer(kind=iwp), intent(out) :: NONEW
integer(kind=iwp) :: IBEL(MXPNGAS), IBORB(MXPNGAS), ICONF_GAS(MXPORB), IGAS, INI_L, JBEL, JBORB, JEL, JGAS, JORB, NEL_GAS, &
                     NONEW_L, NORB_GAS
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NEL

! Total number of electrons
NEL = sum(IOCCLS)
#endif
!write(u6,*) ' NEXT_CONF ... NEL, NGAS = ',NEL,NGAS
! Offset for orbitals and electrons
do IGAS=1,NGAS
  if (IGAS == 1) then
    IBORB(IGAS) = 1
    IBEL(IGAS) = 1
  else
    IBORB(IGAS) = IBORB(IGAS-1)+NOBPT(IGAS-1)
    IBEL(IGAS) = IBEL(IGAS-1)+IOCCLS(IGAS-1)
  end if
end do

NONEW = 1

if (INI == 1) then

  ! Initial configuration

  NONEW = 0
  INI_L = 1
  NONEW_L = 0
  do IGAS=1,NGAS
    ! Initial configuration for this GASSPACE
    NEL_GAS = IOCCLS(IGAS)
    NORB_GAS = NOBPT(IGAS)
    !write(u6,*) ' IGAS, NEL_GAS, NORB_GAS = ',IGAS,NEL_GAS,NORB_GAS
    call NXT_CONF(ICONF_GAS,NEL_GAS,NORB_GAS,INI_L,NONEW_L)
    if (NONEW_L == 1) then
      NONEW = 1
      exit
    else
      JBEL = IBEL(IGAS)
      JBORB = IBORB(IGAS)
      JEL = IOCCLS(IGAS)
      JORB = NOBPT(IGAS)
      call REFORM_CONF_FOR_GAS2(ICONF_GAS,ICONF,JBORB,JBEL,MXPORB,JEL)
    end if
  end do

# ifdef _DEBUGPRINT_
  if (IGAS > NGAS) then
    write(u6,*) ' Initial configuration'
    call IWRTMA(ICONF,1,NEL,1,NEL)
  end if
# endif

else

  ! Next configuration

  ! Loop over GAS spaces and find first GASspace where a new configuration
  ! could be obtained
  do IGAS=1,NGAS
    !write(u6,*) ' IGAS = ',IGAS
    ! Remove the offsets for this space
    JBEL = IBEL(IGAS)
    JBORB = IBORB(IGAS)
    JEL = IOCCLS(IGAS)
    JORB = NOBPT(IGAS)
    !write(u6,*) ' JBEL, JBORB, JEL, JORB = ',JBEL,JBORB,JEL,JORB
    call REFORM_CONF_FOR_GAS1(ICONF_GAS,ICONF,JBORB,JBEL,MXPORB,JEL)
    ! Generate next configuration for this space
    INI_L = 0
    call NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
    if (NONEW_L == 0) then
      NONEW = 0
      ! Configuration in space IGAS, was increased. Copy this and reset configurations
      ! in previous gasspaces to initial form
      call REFORM_CONF_FOR_GAS2(ICONF_GAS,ICONF,JBORB,JBEL,MXPORB,JEL)

      do JGAS=1,IGAS-1
        JBEL = IBEL(JGAS)
        JBORB = IBORB(JGAS)
        JEL = IOCCLS(JGAS)
        JORB = NOBPT(JGAS)
        call REFORM_CONF_FOR_GAS1(ICONF_GAS,ICONF,JBORB,JBEL,MXPORB,JEL)
        INI_L = 1
        call NXT_CONF(ICONF_GAS,JEL,JORB,INI_L,NONEW_L)
        call REFORM_CONF_FOR_GAS2(ICONF_GAS,ICONF,JBORB,JBEL,MXPORB,JEL)
      end do
      ! Get out of the loop
      exit
    end if
  end do
  ! End of loop over gasspaces
end if
! End if switch between initialization/next

#ifdef _DEBUGPRINT_
if (NONEW == 1) then
  write(u6,*) ' No new configuration'
else
  write(u6,*) ' New configuration'
  call IWRTMA(ICONF,1,NEL,1,NEL)
end if
#endif

end subroutine NEXT_CONF_FOR_OCCLS
