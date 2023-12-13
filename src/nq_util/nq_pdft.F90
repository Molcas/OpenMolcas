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

module nq_pdft

use Definitions, only: wp, iwp

implicit none
private

! ThrsRho: threshold of total density
! ThrsOMR: threshold of (1 - R)
! ThrsFT : threshold for doing full translation in ft functionals
!          a.k.a. R0 in the ft paper
! ThrsNT : threshold for not doing any translation in ft functionals
!          a.k.a. R1 in the ft paper

real(kind=wp), parameter :: fta = -475.6065601_wp, ftb = -379.4733192_wp, ftc = -85.3814968_wp, ThrsFT = 0.9_wp, ThrsNT = 1.15_wp, &
                            ThrsOMR = 1.0e-15_wp, ThrsRho = 1.0e-15_wp
logical(kind=iwp) :: lft = .false., lGGA = .false.
real(kind=wp), allocatable :: d2RdRho2(:), d2RdRhodPi(:), d2ZdR2(:), dEdPi(:), dEdPiMO(:,:), dEdPix(:), dEdPiy(:), dEdPiz(:), &
                              dEdRho(:), dEdRhox(:), dEdRhoy(:), dEdRhoz(:), dF_dRhoamb(:), dF_dRhoapb(:), dF_dRhoxamb(:), &
                              dF_dRhoxapb(:), dF_dRhoyamb(:), dF_dRhoyapb(:), dF_dRhozamb(:), dF_dRhozapb(:), dRdPi(:), dRdRho(:), &
                              dRhodX(:), dRhodY(:), dRhodZ(:), dZdR(:), dZdRho(:), GdEdPiMO(:,:), GradPidFdRho(:), GradRdFdRho(:), &
                              GradRhodFdRho(:), MOas(:,:), MOax(:,:), MOay(:,:), MOaz(:,:), OneMZ(:), OnePZ(:), RatioA(:), &
                              RhoAB(:), ZetaA(:)
logical(kind=iwp), allocatable :: Pass1(:), Pass2(:), Pass3(:)

public :: d2RdRho2, d2RdRhodPi, d2ZdR2, dEdPi, dEdPiMO, dEdPix, dEdPiy, dEdPiz, dEdRho, dEdRhox, dEdRhoy, dEdRhoz, dF_dRhoamb, &
          dF_dRhoapb, dF_dRhoxamb, dF_dRhoxapb, dF_dRhoyamb, dF_dRhoyapb, dF_dRhozamb, dF_dRhozapb, dRdPi, dRdRho, dRhodX, dRhodY, &
          dRhodZ, dZdR, dZdRho, fta, ftb, ftc, GdEdPiMO, GradPidFdRho, GradRdFdRho, GradRhodFdRho, lft, lGGA, MOas, MOax, MOay, &
          MOaz, OneMZ, OnePZ, Pass1, Pass2, Pass3, RatioA, RhoAB, ThrsFT, ThrsNT, ThrsOMR, ThrsRho, ZetaA

end module nq_pdft
