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
      Subroutine SetMltplCenters()
C
C     Thomas Bondo Pedersen, July 2012.
C
C     Set multipole centers.
C
      use MpmC
      use Sizes_of_Seward, only: S
      use Real_Info, only: CoM
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"

      ! Check
      If (S%nMltpl.lt.0) Then
         Call WarningMessage(2,'SetMltplCenters: illegal input')
         Write(6,'(A,I10)') 'S%nMltpl=',S%nMltpl
         Call Abend()
      End If

      ! Allocate center array
      Call mma_allocate(Coor_MPM,3,S%nMltpl+1,label='Coor_MPM')

      ! Set origin as center for overlap (0th order)
      Call fZero(Coor_MPM(1,1),3)

      ! Set origin as center for dipole (1st order) and
      ! center of mass as center for higher-order multipoles
      If (S%nMltpl.gt.0) Then
         Call fZero(Coor_MPM(1,2),3)
         Do i=2,S%nMltpl
            Call dCopy_(3,CoM,1,Coor_MPM(1,i+1),1)
         End Do
      End If

      End
