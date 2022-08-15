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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************
      Subroutine Request_MCLR_Run(Run_MCLR,ireturn,iPrint)
      Logical Run_MCLR
      Character*16 StdIn
#include "warnings.h"
!
      If (Run_MCLR) Then
!
!        McKinley will automatically generate the input for MCLR
!        and signal to AUTO (iRC=2) to run the input file Stdin.x.
!
         If (iPrint.ge.6) Then
         Write (6,*)
         Write (6,*)                                                    &
     &     ' McKinley requests the MCLR module to be executed!'
         Write (6,*)
         End If
!
         LuInput=11
         LuInput=IsFreeUnit(LuInput)
         Call StdIn_Name(StdIn)
         Call Molcas_Open(LuInput,StdIn)
         Write (LuInput,'(A)') ' &MCLR &End'
         Write (LuInput,'(A)') 'End of Input'
         Close(LuInput)
         ireturn=_RC_INVOKED_OTHER_MODULE_
      Else
         ireturn=_RC_ALL_IS_WELL_
      End if
!
      Return
!
      Write (6,*)
      Write (6,*) ' Error opening Stdin.x'
      Write (6,*)
      Call Quit(_RC_INPUT_EMIL_ERROR_)
!
      End
