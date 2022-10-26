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
! Copyright (C) 1999, Thorstein Thorsteinsson                          *
!***********************************************************************
      Subroutine Real2Int(RealArg,IntArg)
!***********************************************************************
!                                                                      *
!    Purpose: Convert Real argument to integer                         *
!                                                                      *
!    Calling parameters:                                               *
!    RealArg: contains real to be converted                            *
!    IntArg : contains on return the corresponding integer array       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     T. Thorsteinsson                                                 *
!     University of Lund, Sweden, 1999                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
!
#ifndef _I8_
      Integer IntArg(2),iTemp(2)
#else
      Integer IntArg(1),iTemp(1)
#endif
      Real*8 RealArg,Temp
      Equivalence (iTemp,Temp)

      Temp=RealArg
      IntArg(1)=iTemp(1)
#ifndef _I8_
      IntArg(2)=iTemp(2)
#endif
      Return
      End
