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
      Integer Function RPA_iUHF()
      Implicit None
#include "rpa_config.fh"
      Integer iUHF
      If (Reference(1:1).eq.'R') Then
         iUHF=1
      Else If (Reference(1:1).eq.'U') Then
         iUHF=2
      Else
         Write(6,'(A,A)') 'Reference=',Reference
         Call RPA_Warn(3,'Unable to determine iUHF in RPA')
         iUHF=-1
      End If
      RPA_iUHF=iUHF
      End
