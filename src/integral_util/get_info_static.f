************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine Get_Info_Static()
************************************************************************
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
      use Symmetry_Info, only: Symmetry_Info_Get
      use Sizes_of_Seward, only: Size_Get
      use DKH_Info, only: DKH_Info_Get
      use Real_Info, only: Real_Info_Get
      use RICD_Info, only: RICD_Info_Get
      use Logical_Info, only: Logical_Info_Get

      Call Symmetry_Info_Get()
      Call Size_Get()
      Call DKH_Info_Get()
      Call Real_Info_Get()
      Call RICD_Info_Get()
      Call Logical_Info_Get()

      Return
      End
