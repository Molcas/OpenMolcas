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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

      Subroutine CalcDg(Dgorbit,GDorbit,Gtuvx,nGD,nTUVX,NAC,lRoots)
      Implicit None
      INTEGER nGD,nTUVX,NAC,lRoots
      Real*8 Dgorbit(nGD),GDorbit(nGD),Gtuvx(nTUVX)

      INTEGER NAC2,lRoots2

      NAC2=NAC**2
      lRoots2=lRoots**2

      CALL DGEMM_('T','N',NAC2,lRoots2,NAC2,1.0d0,                      &
     &              Gtuvx,NAC2,GDorbit,NAC2,0.0d0,                      &
     &                         Dgorbit,NAC2)

      RETURN
      End Subroutine
