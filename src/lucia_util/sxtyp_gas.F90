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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine SXTYP_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP)
! Two supergroups are given. Find single excitations that connects
! these supergroups
!
! Jeppe Olsen, July 1995

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: NSXTYP
integer(kind=iwp), intent(_OUT_) :: ITP(*), JTP(*)
integer(kind=iwp), intent(in) :: NGAS, ILTP(NGAS), IRTP(NGAS)
integer(kind=iwp) :: IANNI, IAS, ICREA, NANNI, NCREA
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ISX
#endif

! Some dummy initializations
ICREA = 0 ! jwk-cleanup
IANNI = 0 ! jwk-cleanup
! Differences between occupations :
NCREA = 0
NANNI = 0
do IAS=1,NGAS
  if (ILTP(IAS) > IRTP(IAS)) then
    NCREA = NCREA+ILTP(IAS)-IRTP(IAS)
    ICREA = IAS
  else if (IRTP(IAS) > ILTP(IAS)) then
    NANNI = NANNI+IRTP(IAS)-ILTP(IAS)
    IANNI = IAS
  end if
end do

if (NCREA > 1) then
  ! Sorry : No single excitation connects
  NSXTYP = 0
else if (NCREA == 1) then
  ! Supergroups differ by one sngle excitation.
  NSXTYP = 1
  ITP(1) = ICREA
  JTP(1) = IANNI
else if (NCREA == 0) then
  ! Supergroups are identical, connects with all diagonal excitations.
  NSXTYP = 0
  do IAS=1,NGAS
    if (IRTP(IAS) /= 0) then
      NSXTYP = NSXTYP+1
      ITP(NSXTYP) = IAS
      JTP(NSXTYP) = IAS
    end if
  end do
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from SXTYP_GAS :'
write(u6,*) ' ======================='
write(u6,*) ' Input left  supergroup'
call IWRTMA(ILTP,1,NGAS,1,NGAS)
write(u6,*) ' Input right supergroup'
call IWRTMA(IRTP,1,NGAS,1,NGAS)
write(u6,*) ' Number of connecting single excitations ',NSXTYP
if (NSXTYP /= 0) then
  write(u6,*) ' Connecting single excitations'
  do ISX=1,NSXTYP
    write(u6,*) ITP(ISX),JTP(ISX)
  end do
end if
#endif

end subroutine SXTYP_GAS
