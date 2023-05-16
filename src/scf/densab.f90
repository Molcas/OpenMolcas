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
      Subroutine DensAB(nBT,nDens,nD,Dens)
      use stdalloc, only: mma_allocate, mma_deallocate
!***********************************************************************
!                                                                      *
!     purpose: calculate Density alpha+beta, alpha-beta and dump them  *
!              to runfile                                              *
!                                                                      *
!***********************************************************************
      Implicit None
      Integer nDens,nD,nBT
      Real*8 Dens(nBT,nD,nDens)

      Real*8, Dimension(:), Allocatable:: Dtemp
!
      If (nD.eq.1) Then
         Call Put_dArray('D1ao',Dens(1,1,nDens),nBT)
      Else
         Call mma_allocate(DTemp,nBT,Label='DTemp')
!
         DTemp(:)=Dens(:,1,nDens)+Dens(:,2,nDens)
         Call Put_dArray('D1ao',DTemp,nBT)
!
         DTemp(:)=Dens(:,1,nDens)-Dens(:,2,nDens)
         Call Put_dArray('D1Sao',DTemp,nBT)
!
         Call mma_deallocate(DTemp)
      End If
!
      Return
      End Subroutine DensAB
