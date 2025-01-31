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
!               2024, Giovanni Li Manni                                *
!***********************************************************************

subroutine GEN_CONF_FOR_OCCLS(IOCCLS,IB_OCCLS,INITIALIZE_CONF_COUNTERS,NGAS,ISYM,MINOP,MAXOP,NSMST,IONLY_NCONF,NTORB,NOBPT, &
                              NCONF_OP,NCONF,IBCONF_REO,IBCONF_OCC,ICONF,IDOREO,IZ_CONF,NCONF_ALL_SYM,ireo,nconf_tot)
! IONLY_NCONF = 1 :
!
! Generate number of configurations of occclass IOCCLS and sym ISYM
!
! IONLY_NCONF = 0 :
!
! Generate number and actual configurations of occclass IOCCLS
! and sym ISYM
!
! Jeppe Olsen, Nov. 2001
!
! G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations

use lucia_data, only: MXPORB

implicit none
integer IB_OCCLS, INITIALIZE_CONF_COUNTERS, NGAS, ISYM, MINOP, MAXOP, NSMST, IONLY_NCONF, NTORB, NCONF, IDOREO, NCONF_ALL_SYM, &
        nconf_tot
! Input
! Number of electrons per gas space
integer IOCCLS(NGAS)
! Number of orbitals per gasspace
integer NOBPT(NGAS)
! Arc weights for configurations
integer IZ_CONF(*)
! Offset for reordering array and occupation array
integer IBCONF_REO(*), IBCONF_OCC(*)
! Output
! Number of configurations per number of open shells, all symmetries
integer NCONF_OP(MAXOP+1)
! And the actual configurations
integer ICONF(*)
! Reorder array : Lex number => Actual number
!integer IREO(*)
integer ireo(nconf_tot)
! Local scratch
integer JCONF(2*MXPORB)
integer IDUM_ARR(1)
integer, external :: IELSUM
integer NTEST, NEL, IZERO, INI, ISUM, I, NONEW, ISYM_CONF, NOPEN, NOCOB, IB_OCC, ILEXNUM, JREO, IDUM, ILEX_FOR_CONF_NEW, ISYMST, &
        NOP_FOR_CONF

NTEST = 0
! Total number of electrons
NEL = IELSUM(IOCCLS,NGAS)
if (INITIALIZE_CONF_COUNTERS == 1) then
  IZERO = 0
  call ISETVC(NCONF_OP,IZERO,MAXOP+1)
  NCONF_ALL_SYM = 0
end if
! Loop over configurations
INI = 1
NCONF = 0
ISUM = 0
call ISETVC(JCONF,IZERO,2*MXPORB)

1000 continue
! Generate an array of integers from 1 to NEL.
! It is the only CSF that matter for HS calculations.
! Skip any loop below... It is an overwhelmingly long loop.
if ((NEL == MINOP) .and. (NEL == NOBPT(2))) then
  do i=1,NEL
    JCONF(i) = i
  end do
  NONEW = 0
else
  call NEXT_CONF_FOR_OCCLS(JCONF,IOCCLS,NGAS,NOBPT,INI,NONEW)
end if
ISUM = ISUM+1
INI = 0
if (NONEW == 0) then
  ! Check symmetry and number of open orbitals for this space
  ISYM_CONF = ISYMST(JCONF,NEL)
  NOPEN = NOP_FOR_CONF(JCONF,NEL)
  NOCOB = NOPEN+(NEL-NOPEN)/2
  if ((NOPEN >= MINOP) .or. (IONLY_NCONF /= 0)) NCONF_ALL_SYM = NCONF_ALL_SYM+1
  if ((ISYM_CONF == ISYM) .and. (NOPEN >= MINOP)) then
    ! A new configuration to be included, reform and save in packed form
    NCONF = NCONF+1
    NCONF_OP(NOPEN+1) = NCONF_OP(NOPEN+1)+1
    if (IONLY_NCONF == 0) then
      ! Lexical number of this configuration
      IB_OCC = IBCONF_OCC(NOPEN+1)+(NCONF_OP(NOPEN+1)-1)*NOCOB
      call REFORM_CONF_OCC(JCONF,ICONF(IB_OCC),NEL,NOCOB,1)
      if (IDOREO /= 0) then
        ! Giovanni and Dongxia 2011.1.31
        ilexnum = ilex_for_conf_new(iconf(ib_occ),nocob,ntorb,nel,iz_conf,0,idum_arr,idum,idum)
        JREO = IBCONF_REO(NOPEN+1)-1+NCONF_OP(NOPEN+1)
        ! Giovanni and Dongxia 2011
        ireo(jreo) = ib_occls-1+ilexnum
      end if
    end if
  end if
  if ((NEL /= MINOP) .or. (NEL /= NOBPT(2))) goto 1000
end if
! End if nonew = 0

if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' ======================================'
  write(6,*) ' Results from configuration generator :'
  write(6,*) ' ======================================'
  write(6,*)
  write(6,*) ' Occupation class in action :'
  call IWRTMA(IOCCLS,1,NGAS,1,NGAS)
  write(6,*) ' Number of configurations of correct symmetry ',NCONF
  write(6,*) ' Number of configurations of all symmetries   ',NCONF_ALL_SYM
  write(6,*) ' Number of configurations for various number of open orbs'
  call IWRTMA(NCONF_OP,1,MAXOP+1,1,MAXOP+1)
  if (IONLY_NCONF == 0) then
    write(6,*) ' Updated list of configurations (may not be the final...)'
    call WRT_CONF_LIST(ICONF,NCONF_OP,MAXOP,NCONF,NEL)
    write(6,*) ' Updated reordering of conf, Lex=>Act (may not be the final'
    call IWRTMA(IREO,1,NCONF_tot,1,NCONF_tot)
  end if
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(NSMST)

end subroutine GEN_CONF_FOR_OCCLS
