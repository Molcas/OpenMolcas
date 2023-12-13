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

function dNuclearMass(Z,A)
!***********************************************************************
!                                                                      *
! Routine: dNuclearMass                                                *
! Purpose: Return the mass for isotope with Z protons and A-Z neutrons.*
!          Nontabulated isotopes are computed by the semiempirical     *
!          mass formula. The mass is in atomic units, not in dalton.   *
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
!                                                                      *
!***********************************************************************

use Isotopes, only: NuclideMass
use Constants, only: Zero, One, Two, Three, Half, uToau
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: dNuclearMass
integer(kind=iwp), intent(in) :: Z, A
real(kind=wp) :: t
real(kind=wp), parameter :: Coef(7) = [1.00781360_wp,1.00866184_wp,0.01685183_wp,0.01928950_wp,0.00075636_wp,0.10146129_wp, &
                                       0.02449108_wp]

!----------------------------------------------------------------------*
! Search table.                                                        *
!----------------------------------------------------------------------*
dNuclearMass = NuclideMass(Z,A)
!----------------------------------------------------------------------*
! Optionally use the semi-empirical mass formula.                      *
!----------------------------------------------------------------------*
if (dNuclearMass < Zero) then
  write(u6,'(a)') '***'
  write(u6,'(a)') '*** dNuclearMass: warning'
  write(u6,'(a,2i6)') '*** semi empirical mass formula used for nuclei (Z,A)=',Z,A
  write(u6,'(a)') '***'
  t = Zero
  t = t+Coef(1)*Z
  t = t+Coef(2)*(A-Z)
  t = t-Coef(3)*A
  t = t+Coef(4)*A**(Two/Three)
  t = t+Coef(5)*Z*(Z-1)/real(A,kind=wp)**(One/Three)
  t = t+Coef(6)*(Z-Half*A)**2/real(A,kind=wp)
  if ((mod(Z,2) == 0) .and. (mod(A,2) == 0)) then
    t = t-Coef(7)/real(A,kind=wp)**(0.75_wp)
  end if
  if ((mod(Z,2) == 1) .and. (mod(A,2) == 0)) then
    t = t+Coef(7)/real(A,kind=wp)**(0.75_wp)
  end if
  dNuclearMass = uToau*t
end if
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
return

end function dNuclearMass
