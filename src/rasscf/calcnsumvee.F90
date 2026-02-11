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
      Function CalcNSumVee(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
      use rasscf_global, only: lRoots
      Implicit None

#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RotMat
      Real*8,DIMENSION(:),Allocatable::Vee
      Real*8 CalcNSumVee
      INTEGER IState

      CALL mma_allocate(Vee,lRoots)
      CalcNSumVee=0.0d0
      CALL CalcVee(Vee,RotMat,DDg)
      DO IState=1,lRoots
       CalcNSumVee=CalcNSumVee+Vee(IState)
      END DO
      CALL mma_deallocate(Vee)
      END Function CalcNSumVee
