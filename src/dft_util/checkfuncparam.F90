!**********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2023, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Oct. 20, 2023, created this file.               *
!*****************************************************************

subroutine CheckFuncParam(FileExtParam)

use Functionals, only: check_n_ext_params
use libxc_parameters, only: FuncExtParams, lExtParams
use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
character(len=128), intent(in) :: FileExtParam
integer(kind=iwp) :: iFunc, istat, LUFile, MaxParam, NFuncs
integer(kind=iwp), allocatable :: nParam(:)
integer(kind=iwp), external :: IsFreeUnit

LUFile = IsFreeUnit(100)
call Molcas_Open(LUFile,FileExtParam)

read(LUFile,*,iostat=istat) NFuncs
if (istat /= 0) then
  call WarningMessage(2,'Error Reading NFuncs in External Parameter File!')
  call Quit_OnUserError()
end if

call mma_allocate(nParam,nFuncs)
read(LUFile,*,iostat=istat) NParam(:)
if (istat /= 0) then
  call WarningMessage(2,'Error Reading NParam in External Parameter File!')
  call Quit_OnUserError()
end if

MaxParam = maxval(NParam(:))
call mma_allocate(FuncExtParams,MaxParam,nFuncs)
FuncExtParams(:,:) = Zero
do iFunc=1,nFuncs
  read(LUFile,*,iostat=istat) FuncExtParams(1:nParam(iFunc),iFunc)
  if (istat /= 0) then
    call WarningMessage(2,'Error Reading Parameters in Ext Param File!')
    call Quit_OnUserError()
  end if
end do
close(LUFile)

call check_n_ext_params(nFuncs,nParam)

lExtParams = .true.
!call mma_deallocate(FuncExtParams) !It will be freed in Set_External_Params.
call mma_deallocate(nParam)

end subroutine CheckFuncParam
