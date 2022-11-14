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

subroutine check_use(nToc,i_run_used,Label)

use RunFile_data, only: lw
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

integer(kind=iwp), intent(in) :: nToc, i_run_used(nToc)
character(len=*), intent(in) :: Label
integer(kind=iwp) :: i, nData, RecTyp
character(len=60) :: Line
character(len=lw), allocatable :: RecLab(:)
integer(kind=iwp), parameter :: MakeWarn = 40 !, MakeErr = 100
#ifdef _HAVE_EXTRA_
! FIXME: This include file should be created, as it contains a shared common block
!# include "lfalcon.fh"
#else
logical(kind=iwp), parameter :: isFalcon = .false.
#endif

do i=1,nToc
  if ((i_run_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    if (.not. allocated(RecLab)) then
      call mma_allocate(RecLab,nToc,label='RecLab')
      call ffRun(Label//' labels',nData,RecTyp)
      call cRdRun(Label//' labels',RecLab,lw*nToc)
    end if
    !write(u6,*) Label//' label ',i,' used ',i_run_used(i), ' times'
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',RecLab(i),';was used ',i_run_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
if (allocated(RecLab)) call mma_deallocate(RecLab)

end subroutine check_use
