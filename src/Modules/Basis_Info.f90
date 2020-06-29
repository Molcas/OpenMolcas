!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
      Module Basis_Info
      Integer, Private, Parameter :: Mxdbsc=2000
!     Work in progress
      Type Distinct_Basis_set_centers
          Integer:: ipCntr
          Integer:: nCntr
      End Type Distinct_Basis_set_centers
!
      Type (Distinct_Basis_set_centers) :: dbsc(Mxdbsc)
      End Module Basis_Info
