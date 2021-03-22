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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine RPA_RdOrb()
!
!     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
!     Read orbitals and orbital energies from InpOrb or from Runfile.
!
      Implicit None
#include "rpa_config.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_RdOrb')


      If (LumOrb) Then
         ! read from InpOrb
         Call RPA_RdOrb_FromInpOrb()
      Else
         ! read from Runfile
         Call RPA_RdOrb_FromRunfile()
      End If


      End
