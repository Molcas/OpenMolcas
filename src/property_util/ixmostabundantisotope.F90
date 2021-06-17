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
! Copyright (C) 2000, Per-Olof Widmark                                 *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

function ixMostAbundantIsotope(Z,Rc,Opt)

!***********************************************************************
!                                                                      *
! Routine: ixMostAbundantIsotope                                       *
! Purpose: The mass number for the most abundant isotope of the atom   *
!          with charge Z is returned. For radioactive atoms, the       *
!          most stable isotope is returned.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden.                                    *
! Written: March 2000                                                  *
! History: March 2017, use isotopes module, Ignacio Fdez. Galvan       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Algorithm: The nuclear mass number is tabulated for most charges.    *
!            For Z=0, A=1 is returned, corresponding to a neutron.     *
!            For untabulated Z, A=176+Z.                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
! Z   - The nuclear charge for which the nuclear mass number is        *
!       returned.                                                      *
! Rc  - Return code.                                                   *
! Opt - Options.                                                       *
!                                                                      *
!***********************************************************************

use isotopes, only: Initialize_Isotopes, ElementList, MaxAtomNum

implicit none
#include "proputil.fh"
integer ixMostAbundantIsotope
!----------------------------------------------------------------------*
! Parameters.                                                          *
!----------------------------------------------------------------------*
integer StopOnError
parameter(StopOnError=_OPT_STOP_ON_ERROR_)
!----------------------------------------------------------------------*
! Dummy parameters.                                                    *
!----------------------------------------------------------------------*
integer Z
integer Rc
integer Opt
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
integer A

!----------------------------------------------------------------------*
! Compute A.                                                           *
!----------------------------------------------------------------------*
call Initialize_Isotopes()
if (Z < 0) then
  write(6,'(a)') '***'
  write(6,'(a)') '*** ixMostAbundantIsotope: error'
  write(6,'(a)') '***    Charge less than zero!'
  write(6,'(a)') '***'
  if (iand(Opt,StopOnError) /= 0) call Quit_OnUserError()
  A = 1
else if (Z == 0) then
  A = 1
else if (Z > MaxAtomNum) then
  A = 176+Z
else
  A = ElementList(Z)%Isotopes(1)%A
end if
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
ixMostAbundantIsotope = A

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(Rc)

end function ixMostAbundantIsotope
