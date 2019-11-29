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
Subroutine Init_Kriging()
  use kriging_mod
!
! Initiate Kriging parameters.
!
  Kriging = .False.
  nspAI = 1
  anMd = .True.
  pAI = 2
  npxAI = 1
  lb(1) = 20.0D0
  lb(2) = 20.0D0
  lb(3) = 1
  miAI = 50
  meAI = 1.0D-8
  blAI = .False.
  mblAI = .False.
  blaAI = .True.
  blavAI=0.50D0
  set_l=.False.
!
  Return
End Subroutine Init_Kriging
