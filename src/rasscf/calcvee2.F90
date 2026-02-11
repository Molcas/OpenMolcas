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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine CalcVee2(Vee,GD,Gtuvx)
      use rasscf_global, only: lRoots, NAC
      Implicit None


#include "warnings.h"
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(lRoots)::Vee
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx

      INTEGER I,t,u,v,x,III
      DO I=1,lRoots
       Vee(I)=0.0d0
       III=I*(I+1)/2
       Do t=1,nac
       Do u=1,nac
       Do v=1,nac
       Do x=1,nac
        Vee(I)=Vee(I)+GD(III,t,u)*GD(III,v,x)*Gtuvx(t,u,v,x)
       End Do
       End Do
       End Do
       End Do
       Vee(I)=Vee(I)/2.0d0
      END DO
      END SUBROUTINE CalcVee2
