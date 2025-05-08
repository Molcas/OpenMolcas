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
!  IniMem
!
!> @brief
!>   Initialize memory for Molcas
!>
!> @details
!> Initialize memory for Molcas.
!***********************************************************************

subroutine IniMem()

use stdalloc, only: MxMem
use Definitions, only: iwp, u6

implicit none
#include "warnings.h"
#include "mama.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iRc
interface
  function allocmem(ref,intof,dblof,chrof,size_) bind(C,name='allocmem_')
    use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
    integer(kind=MOLCAS_C_INT) :: allocmem
    real(kind=MOLCAS_C_REAL), intent(in) :: ref(*)
    integer(kind=MOLCAS_C_INT), intent(out) :: intof, dblof, chrof, size_
  end function allocmem
end interface

!----------------------------------------------------------------------*
!     Initialize the manager the first time it is referenced           *
!----------------------------------------------------------------------*
MemStat = .true.

!----------------------------------------------------------------------*
!     Grab from the system a pointer to the dynamic work area          *
!----------------------------------------------------------------------*
iRc = allocmem(Work,iofint,iofdbl,iofchr,MxMem)
if (iRc /= 0) then
  write(u6,'(A,I3,A)') 'The initialization of the memory manager failed ( iRc=',iRc,' ).'
  call Quit(_RC_MEMORY_ERROR_)
end if
!----------------------------------------------------------------------*
!     exit                                                             *
!----------------------------------------------------------------------*

end subroutine IniMem
