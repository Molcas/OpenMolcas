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

module output

  use definitions, only: wp,iwp

  implicit none

  ! amount of output written
  Integer(kind=iwp),parameter :: silent  = 0
  Integer(kind=iwp),parameter :: terse   = 1
  Integer(kind=iwp),parameter :: usual   = 2
  Integer(kind=iwp),parameter :: verbose = 3
  Integer(kind=iwp),parameter :: debug   = 4
  Integer(kind=iwp),parameter :: insane  = 5

  Integer(kind=iwp) :: iPrGlb

  ! thresholds for printing
  Real(kind=wp) :: dnmThr,cmpThr,cntThr

  Real(kind=wp) :: EMP2
  Real(kind=wp) :: STrA, STrF, STrX

end module
