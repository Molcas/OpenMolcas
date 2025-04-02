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

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: ResQaaLag2
use input_mclr, only: nRoots, Eps
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Pi
use Definitions, only: u6

implicit none
#include "warnings.h"
! Output
real*8, dimension(nTri_Elem(nRoots-1)) :: zX
! Input
real*8, dimension(nTri_Elem(nRoots-1)) :: bX
real*8, dimension(nTri_Elem(nRoots-1)**2) :: AXX
! Assistants
real*8, dimension(:), allocatable :: EigVal, bxscr, zXscr, Scr
integer NDim, nSPair, iPair, nScr, INFO

NDim = nTri_Elem(nRoots-1)
nSPair = nDim
ResQaaLag2 = Zero
call mma_allocate(EigVal,nDim)
call mma_allocate(bxScr,nDim)
call mma_allocate(zXScr,nDim)

call GetDiagScr(nScr,AXX,EigVal,nDim)
call mma_allocate(Scr,nScr)

call DSYEV_('V','U',nDim,AXX,nDim,EigVal,Scr,nScr,INFO)

call DGEMM_('n','n',1,nDim,nDim,One,bx,1,AXX,nDim,Zero,bxScr,1)

do iPair=1,nDim
  zxScr(iPair) = -bxScr(iPair)/EigVal(iPair)
  if (abs(zxScr(iPair)) > Two*Pi) then
    zxScr(iPair) = Zero
    ResQaaLag2 = ResQaaLag2+bxScr(iPair)**2
  end if
end do

write(u6,'(6X,A37,2X,ES17.9)') 'Residual in Qaa Lagrange Multipliers:',sqrt(ResQaaLag2)
if (ResQaaLag2 > Eps**2) then
  write(u6,*)
  write(u6,'(6X,A)') 'ERROR: RESIDUAL(S) FOR INTERMEDIATE STATE TOO BIG!'
  write(u6,*)
  write(u6,'(6X,A)') 'This may come from a linear molecular or a linear'
  write(u6,'(6X,A)') 'fragment.'
  write(u6,'(6X,A)') 'CMS-PDFT Lagrange multipliers are not solved.'
  call WarningMessage(2,'Residual in Lagrange Multipliers for Qaa Too Big')
  call Quit(_RC_EXIT_EXPECTED_)
end if

call DGEMM_('n','t',1,nSPair,nSPair,One,zXScr,1,AXX,nSPair,Zero,zx,1)

call mma_deallocate(EigVal)
call mma_deallocate(bxScr)
call mma_deallocate(zXScr)
call mma_deallocate(Scr)

end subroutine SolveforzX
