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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
<<<<<<< HEAD
Module Fockof
use definitions, only: iwp, wp
Private
Public:: IOFFIT,IOFFIA,IOFFTA
Public:: FIT, FIT_Full
Public:: FTI, FTI_Full
Public:: FIA, FIA_Full
Public:: FAI, FAI_Full
Public:: FTA, FTA_Full
Public:: FAT, FAT_Full

Type rPointers
     Real(kind=wp), Pointer:: A(:)=>Null()
End Type rPointers

Type (rPointers):: FIT(8)
Real(kind=wp), Allocatable, Target:: FIT_Full(:)
Type (rPointers):: FTI(8)
Real(kind=wp), Allocatable, Target:: FTI_Full(:)
Type (rPointers):: FIA(8)
Real(kind=wp), Allocatable, Target:: FIA_Full(:)
Type (rPointers):: FAI(8)
Real(kind=wp), Allocatable, Target:: FAI_Full(:)
Type (rPointers):: FTA(8)
Real(kind=wp), Allocatable, Target:: FTA_Full(:)
Type (rPointers):: FAT(8)
Real(kind=wp), Allocatable, Target:: FAT_Full(:)

Integer(kind=iwp):: IOFFIT(8),IOFFIA(8),IOFFTA(8)
End Module Fockof
=======

module Fockof

use Data_Structures, only: V1
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: IOFFIT(8), IOFFIA(8), IOFFTA(8)
type(V1) :: FAI(8), FAT(8), FIA(8), FIT(8), FTA(8), FTI(8)
real(kind=wp), allocatable, target :: FAI_Full(:), FAT_Full(:), FIA_Full(:), FIT_Full(:), FTA_Full(:), FTI_Full(:)

public :: FAI, FAI_Full, FAT, FAT_Full, FIA, FIA_Full, FIT, FIT_Full, FTA, FTA_Full, FTI, FTI_Full, IOFFIA, IOFFIT, IOFFTA

end module Fockof
>>>>>>> upstream-openmolcas/master
