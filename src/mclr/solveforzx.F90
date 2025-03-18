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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      subroutine SolveforzX(zX,AXX,bX)
      use stdalloc, only : mma_allocate, mma_deallocate
      use cmslag,   only : ResQaaLag2
      use Constants, only: Pi
      use input_mclr, only: nRoots,Eps
      Implicit None
#include "warnings.h"
!***** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::zX
!***** Input
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
      Real*8,DIMENSION(((nRoots-1)*nRoots/2)**2)::AXX
!***** Assistants
      Real*8,DIMENSION(:),Allocatable::EigVal,bxscr,zXscr,Scr
      INTEGER NDim,nSPair,iPair,nScr,INFO


      NDim=((nRoots-1)*nRoots/2)
      nSPair=nDim
      ResQaaLag2=0.0d0
      CALL mma_allocate(EigVal,nDim)
      CALL mma_allocate(bxScr ,nDim)
      CALL mma_allocate(zXScr ,nDim)

      CALL GetDiagScr(nScr,AXX,EigVal,nDim)
      CALL mma_allocate(Scr   ,nScr)

      CALL DSYEV_('V','U',nDim,AXX,nDim,EigVal,Scr,nScr,INFO)

      CALL DGEMM_('n','n',1,nDim,nDim,1.0d0,bx,1,AXX,nDim,              &
     &                                0.0d0,bxScr,1)


      DO iPair=1,nDim
       zxScr(iPair)=-bxScr(iPair)/EigVal(iPair)
       IF(Abs(zxScr(iPair)).gt.2.0d0*Pi) THEN
        zxScr(iPair)=0.0d0
        ResQaaLag2=ResQaaLag2+bxScr(iPair)**2
       END IF
      END DO

      write(6,'(6X,A37,2X,ES17.9)')                                     &
     & 'Residual in Qaa Lagrange Multipliers:',SQRT(ResQaaLag2)
      IF(ResQaaLag2.gt.Eps**2) THEN
        write(6,*)
        write(6,'(6X,A)')                                               &
     &    'ERROR: RESIDUAL(S) FOR INTERMEDIATE STATE TOO BIG!'
        write(6,*)
        write(6,'(6X,A)')                                               &
     &    'This may come from a linear molecular or a linear'
        write(6,'(6X,A)')                                               &
     &    'fragment.'
        write(6,'(6X,A)')                                               &
     &    'CMS-PDFT Lagrange multipliers are not solved.'
        CALL WarningMessage(2,                                          &
     &    'Residual in Lagrange Multipliers for Qaa Too Big')
        CALL Quit(_RC_EXIT_EXPECTED_)
      END IF

      CALL DGEMM_('n','t',    1,nSPair,nSPair,                          &
     &            1.0d0,zXScr,1,AXX,nSPair,                             &
     &            0.0d0,zx   ,1)

      CALL mma_deallocate(EigVal)
      CALL mma_deallocate(bxScr )
      CALL mma_deallocate(zXScr )
      CALL mma_deallocate(Scr   )
      END SUBROUTINE SolveforzX
