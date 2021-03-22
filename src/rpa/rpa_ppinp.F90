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
      Subroutine RPA_PPInp()
!
!     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
!     Input postprocessing.
!
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"

      Character*9 SecNam
      Parameter (SecNam='RPA_PPInp')

      ! set RPAModel
      If (dRPA) Then
         If (SOSEX) Then
            RPAModel='SOSX@'//Reference(1:3)
         Else
            RPAModel='dRPA@'//Reference(1:3)
         End If
      Else
         ! this should never happen
         Call RPA_Warn(3,SecNam//': internal error [RPAModel]')
         RPAModel='None@Non'
      End If

      ! freeze orbitals
      Call RPA_Freezer()

      ! print config after input processing
      Call RPA_PrInp()

      End
