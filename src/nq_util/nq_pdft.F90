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
Module nq_pdft

! ThrsRho: threshold of total density
! ThrsOMR: threshold of (1 - R)
! ThrsFT : threshold for doing full translation in ft functionals
!          a.k.a R0 in the ft paper
! ThrsNT : threshold for not doing any translation in ft functionals
!          a.k.a R1 in the ft paper

Real*8 :: ThrsRho  =  1.00d-15
Real*8 :: ThrsOMR  =  1.00d-15
Real*8 :: ThrsFT   =  0.90d0
Real*8 :: ThrsNT   =  1.15d0
Real*8 :: fta      = -4.756065601d2
Real*8 :: ftb      = -3.794733192d2
Real*8 :: ftc      = -8.538149682d1

Logical :: lGGA=.False., lft=.False.
Logical,DIMENSION(:),Allocatable::Pass1,Pass2,Pass3
Real*8 ,DIMENSION(:),Allocatable::RhoAB,OnePZ,OneMZ,RatioA,ZetaA
Real*8 ,DIMENSION(:),Allocatable::dZdR,dRdRho,dZdRho,dRdPi
Real*8 ,DIMENSION(:),Allocatable::dRhoadZ,dRhoaxdZ,dRhoaydZ,dRhoazdZ
Real*8 ,DIMENSION(:),Allocatable::dRhodX,dRhodY,dRhodZ
Real*8 ,DIMENSION(:),Allocatable::dF_dRhoapb,dF_dRhoamb
Real*8 ,DIMENSION(:),Allocatable::dF_dRhoxapb,dF_dRhoxamb
Real*8 ,DIMENSION(:),Allocatable::dF_dRhoyapb,dF_dRhoyamb
Real*8 ,DIMENSION(:),Allocatable::dF_dRhozapb,dF_dRhozamb
Real*8 ,DIMENSION(:),Allocatable::GradRhodFdRho,GradRdFdRho,GradPidFdRho
Real*8 ,DIMENSION(:),Allocatable::dEdRho,dEdRhox,dEdRhoy,dEdRhoz
Real*8 ,DIMENSION(:),Allocatable::dEdPi,dEdPix,dEdPiy,dEdPiz
Real*8 ,DIMENSION(:),Allocatable::dEdPiMO,GdEdPiMO
Real*8 ,DIMENSION(:),Allocatable::d2RdRho2,d2RdRhodPi,d2ZdR2
Real*8 ,DIMENSION(:),Allocatable::MOas,MOax,MOay,MOaz
End Module nq_pdft
