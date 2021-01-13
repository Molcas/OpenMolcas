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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
Subroutine Energy_Kriging(x0_,y_,ndimx)
  use kriging_mod
  Implicit None
  Integer nDimx
  Real*8 x0_(ndimx),y_
!
!x0 is the n-dimensional vector of the coordinates at which the energy is evaluated
! subroutine
  x0(:) = x0_(:)
!
  call covarvector(0) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(0)
  y_=pred
!
  return
End Subroutine Energy_Kriging
