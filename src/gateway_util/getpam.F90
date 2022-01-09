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
! Copyright (C) 2001,2020,  Roland Lindh                               *
!               Sergey Gusarov                                         *
!***********************************************************************

subroutine GetPAM(lUnit,iCnttp)
!***********************************************************************
!                                                                      *
!    Objective: To read potential information, for DMFT calculations   *
!               This means that we read (from input stream)            *
!               the potential terms (PAM)                              *
!                                                                      *
! Called from: GetBs                                                   *
!                                                                      *
! Calling    : RecPrt, Rdbsl                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!                                                                      *
!     Modified: Sergey Gusarov SPb State Ubiv, Russia.                 *
!***********************************************************************

use Basis_Info

implicit real*8(A-H,O-Z)
#include "itmax.fh"
#include "real.fh"
#include "stdalloc.fh"
real*8, allocatable :: Array(:)
character*180 Line, Get_Ln
!external Get_Ln
integer nPAM2
logical test
#ifdef _DEBUGPRINT_
data test/.true./
#else
data test/.false./
#endif

if (test) write(6,*) ' Reading PAM potencials'
nArray = 10000
call mma_Allocate(Array,nArray,Label='Array')

iStrt = 1
Line = Get_Ln(lUnit)
if (index(Line,'PAM') == 0) then
  call WarningMessage(2,'ERROR: Keyword PAM expected, offending line : '//Line)
  call Quit_OnUserError()
end if
Line = Get_Ln(lUnit)
call Get_i1(1,nPAM2)
dbsc(iCnttp)%nPAM2 = nPAM2
do iPAM_Ang=0,nPAM2
  Line = Get_Ln(lUnit)
  call Get_i1(1,nPAM2Prim)
  call Get_i1(2,nPAM2Bas)
  Array(iStrt) = dble(nPAM2Prim)
  Array(iStrt+1) = dble(nPAM2Bas)
  iStrt = iStrt+2
  iEnd = iStrt+nPAM2Prim-1

  ! Read exponents

  if (nPAM2Prim > 0) then
    call read_v(lUnit,Array,iStrt,iEnd,1,Ierr)
    if (Ierr /= 0) then
      call WarningMessage(2,'GetPAM: Error reading GPA exponents')
      call Abend()
    end if
  end if

  ! Read coefficents

  iStrt = iEnd+1
  iEnd = iStrt+nPAM2Prim*nPAM2Bas-1
  do iPrim=0,nPAM2Prim-1
    call Read_v(lUnit,Array,iStrt+iPrim,iEnd,nPAM2Prim,ierr)
    if (ierr /= 0) then
      call WarningMessage(2,'GetPAM: Error in reading GPA!!!')
      call Abend()
    end if
  end do
  iStrt = iEnd+1
end do

call mma_allocate(dbsc(iCnttp)%PAM2,iEnd,Label='PAM2')
dbsc(iCnttp)%PAM2(:) = Array(:)

call mma_deAllocate(Array)

return

end subroutine GetPAM
