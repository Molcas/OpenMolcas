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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Oct. 20, 2023, created this file.               *
!*****************************************************************

subroutine FuncParam_PDFT(Functional,FileExtParam)

  use xc_f03_lib_m,   only: xc_f03_func_t, xc_f03_func_info_t, xc_f03_func_get_info, xc_f03_func_info_get_n_ext_params, &
  xc_f03_func_init_, XC_UNPOLARIZED
  use Definitions,    only: wp, iwp, u6
  use stdalloc,       only: mma_allocate, mma_deallocate

  CHARACTER(len=*),   intent(in) :: Functional
  CHARACTER(len=128), intent(in) :: FileExtParam
  integer(kind=iwp)              :: func_id, LUFile, IReadStat, K, N_Ext_Param
  integer(kind=iwp), external    :: isFreeUnit
  type(xc_f03_func_t)            :: func
  type(xc_f03_func_info_t)       :: info
  real(kind=wp), dimension(:), allocatable :: Ext_Params
  Logical(kind=iwp)               :: Test
 INTERFACE
     Function Get_Func(XcLabel, test)
       integer(kind=LibxcInt) :: get_func
       character(len=*), intent(in) :: xcLabel
       logical(kind=iwp), intent(in), optional :: test
     End Function Get_Func
 End INTERFACE

write(u6,*) 'getting functional id for ', Functional
func_id = get_func(Functional, test=.true.)
write(u6,*) 'functional ', adjustl(trim(Functional)), ' is ', func_id
write(u6,*) 'calling xc_f03_func_init'
call xc_f03_func_init(func, func_id, XC_UNPOLARIZED)
write(u6,*) 'calling xc_f03_func_get_info'
info = xc_f03_func_get_info(func)
write(u6,*) 'calling xc_f03_func_get_n_ext_params'
n_ext_param = xc_f03_func_info_get_n_ext_params(info)

CALL mma_allocate(Ext_Params, N_Ext_Param)
LUFile = IsFreeUnit(100)
CALL Molcas_Open(LUFile, FileExtParam)

READ(LUFile, *, iostat=IReadStat) (Ext_Params(K), K=1, N_Ext_Param)

IF (IReadStat /= 0) THEN
  CALL WarningMessage(2, ' MC-PDFT Error Reading External Parameter File! ')
  write(u6,'(A12,10(F9.6,1A))') 'Parameters: ', Ext_Params
  CALL Quit_OnUserError()
END IF
write(u6,'(A12,10(F9.6,1A))') 'Parameters: ', Ext_Params
CLOSE(LUFile)
CALL mma_deallocate(Ext_Params)

end subroutine FuncParam_PDFT
