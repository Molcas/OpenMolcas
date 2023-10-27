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

  use Definitions,    only: wp, iwp, u6
  use Constants,      only: Zero
  use stdalloc,       only: mma_allocate, mma_deallocate
  use Functionals,    only: check_n_ext_params
  use libxc_parameters, only: FuncExtParams, lExtParams

  CHARACTER(len=128),intent(in) :: FileExtParam
  INTEGER(kind=iwp) :: NFuncs, MaxParam, istat, iFunc, I
  INTEGER(kind=iwp),DIMENSION(:),Allocatable :: nParam

LUFile = IsFreeUnit(100)
CALL Molcas_Open(LUFile, FileExtParam)

READ(LUFile, *, iostat=istat) NFuncs
IF (istat /= 0) THEN
  CALL WarningMessage(2, 'Error Reading NFuncs in External Parameter File!')
  CALL Quit_OnUserError()
END IF


CALL mma_allocate(nParam, nFuncs)
READ(LUFile, *, iostat=istat) (NParam(I), I = 1, nFuncs)
IF (istat /= 0) THEN
  CALL WarningMessage(2, 'Error Reading NParam in External Parameter File!')
  CALL Quit_OnUserError()
END IF
write(u6, *) 'Number of External Parameters in each functional'
write(u6, *) nParam

MaxParam = 0
DO iFunc = 1,  nFuncs
  If(nParam(iFunc) .gt. MaxParam) Then
    MaxParam = nParam(iFunc)
  End If
END DO

CALL mma_allocate(FuncExtParams, MaxParam, nFuncs)
FuncExtParams(:,:) = Zero
DO iFunc = 1,  nFuncs
  READ(LUFile, *, iostat=istat) (FuncExtParams(I, iFunc), I=1, nParam(iFunc))
  IF (istat /= 0) THEN
    CALL WarningMessage(2, 'Error Reading Parameters in Ext Param File!')
    CALL Quit_OnUserError()
  END IF
END DO
CLOSE(LUFile)

write(u6,'(A12,10(F9.6,1A))') 'Parameters:'
CALL RecPrt(' ','(10(F9.6))', FuncExtParams, MaxParam, nFuncs)

CALL check_n_ext_params(nFuncs, nParam)

lExtParams = .true.
!CALL mma_deallocate(FuncExtParams) !It will be freed in Set_External_Params.
CALL mma_deallocate(nParam)
end subroutine CheckFuncParam

