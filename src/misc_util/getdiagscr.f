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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Jul 01, 2022, created this file.               *
! ****************************************************************


      Subroutine GetDiagScr(nScr,Mat,EigVal,nDim)
      INTEGER nScr,nDim,INFO
      Real*8 Mat(nDim**2)
      Real*8 EigVal(nDim)
      Real*8 Scr(2)

      CALL DSYEV_('V','U',nDim,Mat,nDim,EigVal,Scr,-1,INFO)
      NScr=INT(Scr(1))
      RETURN
      End Subroutine


