!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module RMat

use Definitions, only: wp, iwp

implicit none
private

!                                                                      *
!***********************************************************************
!                                                                      *
! RmatR    : radius of the R-matrix sphere (bohr)
! Epsabs   : absolute precision of numerical radial integration
! Epsrel   : relative precision of numerical radial integration
! qCoul    : effective charge of the target molecule
! dipol(3) : effective dipole moment of the target molecule
! dipol1   : abs sum of dipol(i)
! epsq     : minimal value of qCoul and/or dipol1 to be considered
! bParm    : Bloch term parameter
!
! keyr: option to the dqag quadpack routine
!
! Quadpack: logical flag to indicate use of QUADPACK routines
! NagInt  : logical flag to indicate use of NAG routines
! testint : logical flag to indicate use of both QUADPACK and NAG routines
! RMat_On : Logical flag to signal that R-matrix type integrals are to be computed.
!
! K.P.
! Dieses File enthaelt Feld fuer die theta/phi  Integration der
! Ein-Elektronen-Integrale fuer die R-Matrix.
! lgamma, n_gam, m_gam, lgammax, gammath, gammaph
!
! Dieses File enthaelt Inputdaten fuer die Modifikation der
! Ein-Elektronen-Integrale fuer die R-Matrix.
! expsum, l, lcost, lsint, lcosf, lsinf
!                                                                      *
!***********************************************************************
!                                                                      *
! The following section should be activated if a NAG library is
! present. At that point I suggest that it gets its own include
! file. -RL-
!
!integer(kind=iwp), parameter :: INTparm = 200
!integer(kind=iwp) :: Iwrk1(INTparm)
!real(kind=wp) :: wrk1(4*INTparm)

integer(kind=iwp), parameter :: lgammax = 15
integer(kind=iwp) :: keyr, l, lgamma
real(kind=wp) :: bParm, dipol(3), Dipol1, Epsabs, Epsq, Epsrel, expsum, gammaph(-2:2*lgammax+3,-2:2*lgammax+4), &
                 gammath(-2:2*lgammax+3,-2:2*lgammax+4), qCoul, RmatR
logical(kind=iwp) :: nagint, Quadpack, RMat_On, RMat_Type_Integrals, testint

public :: bParm, dipol, Dipol1, Epsabs, Epsq, Epsrel, expsum, gammaph, gammath, keyr, l, lgamma, nagint, qCoul, Quadpack, RMat_On, &
          RMat_Type_Integrals, RmatR, testint

end module RMat
