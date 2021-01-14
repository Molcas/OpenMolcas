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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_

Module Phase_Info
Implicit None
Integer ::iPhase(3,0:7)=Reshape ([ 1, 1, 1,-1, 1, 1, 1,-1, &
                                   1,-1,-1, 1, 1, 1,-1,-1, &
                                   1,-1, 1,-1,-1,-1,-1,-1], [3,8])
End Module Phase_Info
