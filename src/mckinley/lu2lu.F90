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
      Subroutine Lu2Lu(Filename,LuInput)
      Character FileName*(*), Line*180
      Logical Exist
#include "warnings.h"
!
      Call f_inquire(Filename,Exist)
      If (.Not.Exist) Then
         Write (6,*) 'SuperMac: Missing ',Filename
         Call Finish(_RC_INTERNAL_ERROR_)
      End If
!
      LuSpool2 = 77
      LuSpool2 = IsFreeUnit(LuSpool2)
      Call Molcas_Open(LuSpool2, Filename)
!
 100  Continue
         Read (LuSpool2,'(A)',End=900) Line
         Write(LuInput,'(A)') Line
         Go To 100
 900  Continue
!
      Close(LuSpool2)
!
      Return
      End
