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
! Copyright (C) 1988, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine TRACI_LUCIA(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
! A rotation matrix X is defining expansion from
! old to new orbitals
!        PHI(NEW) = PHI(OLD) * X
!
! change CI coefficients (sym IXSM, space IXSPC)
! so they corresponds to PHI(NEW) basis
!
! The input CI vector is on LUCIN and the transformed CI vector
! will be delivered on LUCOUT.
!
! Transformation as conceived by Per-AAke Malmquist
! (I.J.Q.C. vol XXX, p479, 1986 (OCTOBER ISSUE))
!
!  Jeppe Olsen 1988
!
! New LUCIA version of Jan 1998
!
! note The transformation matrix X is supposed to be in complete form
! as a matrix over NTOOB orbitals.

use Index_Functions, only: nTri_Elem
use CandS, only: ICSM, ICSPC, ISSM, ISSPC
use lucia_data, only: LUSC1, LUSC2, LUSC3, NSMOB, NTOOB, NTOOBS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: X(*)
integer(kind=iwp), intent(in) :: LUCIN, LUCOUT, IXSPC, IXSM
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)
real(kind=wp), allocatable :: SCR(:), LT(:)
integer(kind=iwp) :: IOFF, ISM

! Some dummy initializations
IOFF = 0 ! jwk-cleanup

#ifdef _DEBUGPRINT_
write(u6,*) ' ================'
write(u6,*) ' Welcome to TRACI'
write(u6,*) ' ================'
write(u6,*)
write(u6,*) ' IXSPC,IXSM = ',IXSPC,IXSM
#endif
! Memory allocation
! for a matrix T
call mma_allocate(LT,NTOOB**2,Label='LT')
! Scratch in PAMTMT
call mma_allocate(SCR,NTOOB**2+nTri_Elem(NTOOB),Label='SCR')
! Obtain T matrix used for transformation, for each symmetry separately
do ISM=1,NSMOB
  if (ISM == 1) then
    IOFF = 1
  else
    IOFF = IOFF+NTOOBS(ISM-1)**2
  end if
  if (NTOOBS(ISM) > 0) call PAMTMT(X(IOFF),LT(IOFF),SCR,NTOOBS(ISM))
end do
! Transform CI-vector
ICSPC = IXSPC
ICSM = ICSM
ISSPC = IXSPC
ISSM = IXSM

call TRACID(LT,LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,VEC1,VEC2)

call mma_deallocate(SCR)
call mma_deallocate(LT)

end subroutine TRACI_LUCIA
