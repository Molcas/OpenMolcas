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

! ThrsRho: threshold of total density
! ThrsOMR: threshold of (1 - R)
! ThrsFT : threshold for doing full translation in ft functionals
!          a.k.a R0 in the ft paper
! ThrsNT : threshold for not doing any translation in ft functionals
!          a.k.a R1 in the ft paper

real*8 :: ThrsRho = 1.00d-15
real*8 :: ThrsOMR = 1.00d-15
real*8 :: ThrsFT = 0.90d0
real*8 :: ThrsNT = 1.15d0
real*8 :: fta = -4.756065601d2
real*8 :: ftb = -3.794733192d2
real*8 :: ftc = -8.538149682d1

logical :: lGGA = .false., lft = .false.
logical, dimension(:), allocatable :: Pass1, Pass2, Pass3
real*8, dimension(:), allocatable :: RhoAB, OnePZ, OneMZ, RatioA, ZetaA
real*8, dimension(:), allocatable :: dZdR, dRdRho, dZdRho, dRdPi
real*8, dimension(:), allocatable :: dRhoadZ, dRhoaxdZ, dRhoaydZ, dRhoazdZ
real*8, dimension(:), allocatable :: dRhodX, dRhodY, dRhodZ
real*8, dimension(:), allocatable :: dF_dRhoapb, dF_dRhoamb
real*8, dimension(:), allocatable :: dF_dRhoxapb, dF_dRhoxamb
real*8, dimension(:), allocatable :: dF_dRhoyapb, dF_dRhoyamb
real*8, dimension(:), allocatable :: dF_dRhozapb, dF_dRhozamb
real*8, dimension(:), allocatable :: GradRhodFdRho, GradRdFdRho, GradPidFdRho
real*8, dimension(:), allocatable :: dEdRho, dEdRhox, dEdRhoy, dEdRhoz
real*8, dimension(:), allocatable :: dEdPi, dEdPix, dEdPiy, dEdPiz
real*8, dimension(:), allocatable :: dEdPiMO, GdEdPiMO
real*8, dimension(:), allocatable :: d2RdRho2, d2RdRhodPi, d2ZdR2
real*8, dimension(:), allocatable :: MOas, MOax, MOay, MOaz

end module nq_pdft
