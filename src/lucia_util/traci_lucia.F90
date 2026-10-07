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

subroutine TRACI_LUCIA(LT,LUCIN,LUCOUT,VEC1,VEC2)
! A rotation matrix X is defining expansion from
! old to new orbitals
!        PHI(NEW) = PHI(OLD) * X
!
! change CI coefficients (sym ISSM, space ISSPC)
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

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUCIN, LUCOUT
real(kind=wp), intent(in) :: LT(*)
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)

! Transform CI-vector

call TRACID(LT,LUCIN,LUCOUT,VEC1,VEC2)

end subroutine TRACI_LUCIA
