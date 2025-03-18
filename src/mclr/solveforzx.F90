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

use stdalloc, only: mma_allocate, mma_deallocate
use cmslag, only: ResQaaLag2
use Constants, only: Pi
use input_mclr, only: nRoots, Eps

implicit none
#include "warnings.h"
! Output
real*8, dimension((nRoots-1)*nRoots/2) :: zX
! Input
real*8, dimension((nRoots-1)*nRoots/2) :: bX
real*8, dimension(((nRoots-1)*nRoots/2)**2) :: AXX
! Assistants
real*8, dimension(:), allocatable :: EigVal, bxscr, zXscr, Scr
integer NDim, nSPair, iPair, nScr, INFO

NDim = ((nRoots-1)*nRoots/2)
nSPair = nDim
ResQaaLag2 = 0.0d0
call mma_allocate(EigVal,nDim)
call mma_allocate(bxScr,nDim)
call mma_allocate(zXScr,nDim)

call GetDiagScr(nScr,AXX,EigVal,nDim)
call mma_allocate(Scr,nScr)

call DSYEV_('V','U',nDim,AXX,nDim,EigVal,Scr,nScr,INFO)

call DGEMM_('n','n',1,nDim,nDim,1.0d0,bx,1,AXX,nDim,0.0d0,bxScr,1)

do iPair=1,nDim
  zxScr(iPair) = -bxScr(iPair)/EigVal(iPair)
  if (abs(zxScr(iPair)) > 2.0d0*Pi) then
    zxScr(iPair) = 0.0d0
    ResQaaLag2 = ResQaaLag2+bxScr(iPair)**2
  end if
end do

write(6,'(6X,A37,2X,ES17.9)') 'Residual in Qaa Lagrange Multipliers:',sqrt(ResQaaLag2)
if (ResQaaLag2 > Eps**2) then
  write(6,*)
  write(6,'(6X,A)') 'ERROR: RESIDUAL(S) FOR INTERMEDIATE STATE TOO BIG!'
  write(6,*)
  write(6,'(6X,A)') 'This may come from a linear molecular or a linear'
  write(6,'(6X,A)') 'fragment.'
  write(6,'(6X,A)') 'CMS-PDFT Lagrange multipliers are not solved.'
  call WarningMessage(2,'Residual in Lagrange Multipliers for Qaa Too Big')
  call Quit(_RC_EXIT_EXPECTED_)
end if

call DGEMM_('n','t',1,nSPair,nSPair,1.0d0,zXScr,1,AXX,nSPair,0.0d0,zx,1)

call mma_deallocate(EigVal)
call mma_deallocate(bxScr)
call mma_deallocate(zXScr)
call mma_deallocate(Scr)

end subroutine SolveforzX
