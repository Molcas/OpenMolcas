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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

real*8 function Get_LDFAccuracy()

implicit none
#include "localdf.fh"
logical LDF_X_IsSet
external LDF_X_IsSet

if (.not. LDF_X_IsSet()) call Get_dScalar('LDF Accuracy',Thr_Accuracy)
Get_LDFAccuracy = Thr_Accuracy

end function Get_LDFAccuracy
