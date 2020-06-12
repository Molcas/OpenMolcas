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
Subroutine Energy_Kriging(x_,y_,ndimx)
  use kriging_mod
  Implicit None
  Integer nInter,nPoints,nDimx
  Real*8 x_(ndimx,1),y_
!
  nPoints=nPoints_save
  nInter=nInter_save
!
  npx = npxAI
!nx is the n-dimensional vector of the last iteration computed in update_sl
! subroutine
  nx(:,:) = x_
!
  call covarvector(0,nPoints,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
  call predict(0,nPoints,nInter)
  y_=pred(npx)
!
  return
End Subroutine Energy_Kriging
