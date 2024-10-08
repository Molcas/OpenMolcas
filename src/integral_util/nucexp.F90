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
! Copyright (C) 1999, Per-Olof Widmark                                 *
!***********************************************************************

function NucExp(A)
!***********************************************************************
!                                                                      *
! Routine: NucExp.                                                     *
! Purpose: This routine computes a nuclear radius in the form of a     *
!          Gaussian exponent. The exponent is a function of the        *
!          nuclear mass number (A).                                    *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden.                                    *
! Written: September 1999                                              *
! History: none.                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
! A  - The nuclear mass number for the nucleus.                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Algorithm: Empirical relation                                        *
!               rms(r)/[fm] = 0.836*A^(1/3) + 0.570                    *
!            Equate empirical rms with Gaussian rms                    *
!               rms(r) = 3/(2*xi)                                      *
!                                                                      *
! References: W. R. Johnson, G. Soff. At. Data Nucl. Data Tables, 33,  *
!             (1985) 405-446. doi:10.1016/0092-640X(85)90010-5         *
!                                                                      *
!             L. Visscher, K. G. Dyall. At. Data Nucl. Data Tables,    *
!             67 (1997) 207-224. doi:10.1006/adnd.1997.0751            *
!                                                                      *
!***********************************************************************

use Constants, only: One, Three, OneHalf, rBohr
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: NucExp
integer(kind=iwp), intent(in) :: A
real(kind=wp) :: A3, R, Xi

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
A3 = real(A,kind=wp)**(One/Three)

R = 0.836_wp*A3+0.570_wp ! fm
R = R*1.0e-15_wp         ! m
R = R/rBohr              ! bohr

Xi = OneHalf/R**2

NucExp = Xi
!----------------------------------------------------------------------*
! Done.                                                                *
!----------------------------------------------------------------------*
return

end function NucExp
