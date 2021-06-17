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

real*8 function dxNuclearMass(Z,A,Rc,Opt)
!***********************************************************************
!                                                                      *
! Routine: dxNuclearMass                                               *
! Purpose: To compute mass for isotope with Z protons and A-Z neutrons.*
!          Nontabulated isotopes are computed by the semiempirical     *
!          mass formula.                                               *
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
! Parameters:                                                          *
! Z   - The nuclear charge for the isotope. Input.                     *
! A   - The number of nucleons for the isotope. Input.                 *
! Rc  - Return code. Output                                            *
! Opt - Options. Input.                                                *
!                                                                      *
!***********************************************************************

use Isotopes

implicit none
#include "proputil.fh"
#include "constants2.fh"
!----------------------------------------------------------------------*
! Parameters.                                                          *
!----------------------------------------------------------------------*
integer StopOnWarning
parameter(StopOnWarning=_OPT_STOP_ON_WARNING_)
integer StopOnError
parameter(StopOnError=_OPT_STOP_ON_ERROR_)
!----------------------------------------------------------------------*
! Dummy parameters.                                                    *
!----------------------------------------------------------------------*
integer Z
integer A
integer Rc
integer Opt
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
real*8 Coef(7)
data Coef/1.00781360d0,1.00866184d0,0.01685183d0,0.01928950d0,0.00075636d0,0.10146129d0,0.02449108d0/
save Coef
real*8 t

!----------------------------------------------------------------------*
! Search table.                                                        *
!----------------------------------------------------------------------*
dxNuclearMass = NuclideMass(Z,A)
!----------------------------------------------------------------------*
! Optionally use the semi-empirical mass formula.                      *
!----------------------------------------------------------------------*
if (dxNuclearMass < 0.0d0) then
  write(6,'(a)') '***'
  write(6,'(a)') '*** dxNuclearMass: warning'
  write(6,'(a,2i6)') '*** semi empirical mass formula used for nuclei (Z,A)=',Z,A
  write(6,'(a)') '***'
  if (iand(StopOnWarning,Opt) /= 0) call Quit_OnUserError()
  t = 0.0d0
  t = t+Coef(1)*Z
  t = t+Coef(2)*(A-Z)
  t = t-Coef(3)*A
  t = t+Coef(4)*A**(2.0d0/3.0d0)
  t = t+Coef(5)*Z*(Z-1)/dble(A)**(1.0d0/3.0d0)
  t = t+Coef(6)*(Z-0.5d0*A)**2/dble(A)
  if ((mod(Z,2) == 0) .and. (mod(A,2) == 0)) then
    t = t-Coef(7)/dble(A)**(0.75)
  end if
  if ((mod(Z,2) == 1) .and. (mod(A,2) == 0)) then
    t = t+Coef(7)/dble(A)**(0.75)
  end if
  dxNuclearMass = uToau*t
end if
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
return
! Avoid unused argument warnings
if (.false.) call Unused_integer(Rc)

end function dxNuclearMass
