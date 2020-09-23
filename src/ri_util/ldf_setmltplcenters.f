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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetMltplCenters(MltplOrder)
C
C     Thomas Bondo Pedersen, July 2012.
C
C     Set multipole centers.
C
      use Sizes_of_Seward, only: S
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"

      If (MltplOrder.ge.0) Then
         S%nMltpl=MltplOrder
         Call Get_dArray('Center of Mass',CoM,3)
         Call SetMltplCenters()
      End If

      End
      Subroutine LDF_UnsetMltplCenters(MltplOrder)
C
C     Thomas Bondo Pedersen, July 2012.
C
C     Unset multipole centers.
C
      use MpmC
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"

      If (MltplOrder.ge.0) Call mma_deallocate(Coor_MPM)

      End
