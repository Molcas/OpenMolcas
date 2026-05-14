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
use input_mclr, only: Eps, nRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: zX(nTri_Elem(nRoots-1))
real(kind=wp), intent(inout) :: AXX(nTri_Elem(nRoots-1),nTri_Elem(nRoots-1))
real(kind=wp), intent(in) :: bX(nTri_Elem(nRoots-1))
#include "warnings.h"
integer(kind=iwp) :: INFO, iPair, NDim, nScr, nSPair
real(kind=wp), allocatable :: bxscr(:), EigVal(:), Scr(:), zXscr(:)

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
