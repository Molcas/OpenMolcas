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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************
      Subroutine NewCar_Kriging(kIter,SaveBMx,Error)
      use Slapaf_Info, only: Cx, BMx, BMx_kriging
      use Slapaf_Parameters, only: PrQ, Numerical, mTtAtm, Force_dB
      Implicit None
#include "stdalloc.fh"
      Integer :: kIter
      Logical :: SaveBMx, Error

      Real*8, Allocatable :: Coor(:,:), BMx_Tmp(:,:)
      Logical :: Numerical_Save,PrQ_Save
!
      Call mma_allocate(Coor,3,SIZE(Cx,2),Label='Coor')
      Coor(:,:) = Cx(:,:,kIter)

      Call mma_allocate(BMx_tmp,SIZE(BMx,1),SIZE(BMx,2),Label='BMx_tmp')
      BMx_tmp(:,:) = BMx(:,:)
!
      Numerical_Save=Numerical
      PrQ_Save=PrQ
      Numerical=.False.
      PrQ_Save=PrQ
      PrQ=.False.
!
      Force_dB=SaveBMx
!
      Call NewCar(kIter,SIZE(Coor,2),Coor,mTtAtm,Error)
!
      Numerical=Numerical_Save
      PrQ=PrQ_Save
      Force_dB=.False.
      Call mma_deallocate(Coor)

!     Stash away this B-matrix for later us with EDiff constraints

      If (Allocated(BMx_kriging)) Call mma_deallocate(BMx_kriging)
      Call mma_allocate(BMx_kriging,SIZE(BMx,1),SIZE(BMx,2),            &
     &                  Label='BMx_kriging')
      BMx_kriging(:,:) = BMx(:,:)

      If (.NOT.SaveBMx) BMx(:,:) = BMx_tmp(:,:)

      Call mma_deallocate(BMx_tmp)
!
      End Subroutine NewCar_Kriging
