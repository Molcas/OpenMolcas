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
      Subroutine AppFld_NonEq_1(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer lMax
      Real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
      Real*8, Allocatable:: CavSph(:,:)
      Real*8 Radius, Eps, EpsInf
      Logical NonEq
!
      Call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
      Call AppFld_1(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
      Call mma_deallocate(CavSph)
!
      Return
      End Subroutine AppFld_NonEq_1
