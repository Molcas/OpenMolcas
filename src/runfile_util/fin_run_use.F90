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

subroutine fin_run_use()

use RunFile_data, only: i_run_CA_used, i_run_DA_used, i_run_DS_used, i_run_IA_used, i_run_IS_used, nTocCA, nTocDA, nTocDS, nTocIA, &
                        nTocIS
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i
character(len=60) :: Line
character(len=16) :: Label
integer(kind=iwp), parameter :: MakeWarn = 40 !, MakeErr = 100
#ifdef _HAVE_EXTRA_
! FIXME: This include file should be created, as it contains a shared common block
!#include "lfalcon.fh"
#else
logical(kind=iwp), parameter :: isFalcon = .false.
#endif
logical(kind=iwp), external :: Reduce_Prt

if (Reduce_Prt()) return
!need_abend = 0
do i=1,nTocCA
  if ((i_run_CA_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    !write(u6,*) 'cArray label ',i,' used ',i_run_CA_used(i),' times'
    call lookup_label(i,'cArray labels',Label)
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,';was used ',i_run_CA_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
do i=1,nTocDA
  if ((i_run_DA_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    !write(u6,*) 'dArray label ',i,' used ',i_run_DA_used(i),' times'
    call lookup_label(i,'dArray labels',Label)
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,';was used ',i_run_DA_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
do i=1,nTocDS
  if ((i_run_DS_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    !write(u6,*) 'dScalar label ',i,' used ',i_run_DS_used(i),' times'
    call lookup_label(i,'dScalar labels',Label)
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,';was used ',i_run_DS_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
do i=1,nTocIA
  if ((i_run_IA_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    !write(u6,*) 'iArray label ',i,' used ',i_run_IA_used(i),' times'
    call lookup_label(i,'iArray labels',Label)
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,';was used ',i_run_IA_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
do i=1,nTocIS
  if ((i_run_IS_used(i) > MakeWarn) .and. (.not. isFalcon)) then
    !write(u6,*) 'iScalar label ',i,' used ',i_run_IS_used(i), ' times'
    call lookup_label(i,'iScalar labels',Label)
    write(Line,'(A,A,A,I8,A)') 'RunFile label ',Label,';was used ',i_run_IS_used(i),' times'
    call WarningMessage(1,Line)
  end if
end do
!if (need_abend == 1) call abend()

return

end subroutine fin_run_use
