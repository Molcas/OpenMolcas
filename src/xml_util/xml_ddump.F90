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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine dumps doubles in xml format.                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_dDump(TagName,Appear,Units,Level,Content,nx,ny)

use Definitions, only: wp, iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName, Appear, Units
real(kind=wp), intent(in) :: Content(*)
integer(kind=iwp), intent(in) :: nx, ny, Level
interface
  subroutine xml_ddumpc(name_,nx_name,appear,nx_appear,units,nx_units,Level,data_,nxx,nyx) bind(C,name='xml_ddumpc_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT, MOLCAS_C_REAL
    character(kind=c_char) :: name_(*), appear(*), units(*)
    integer(kind=MOLCAS_C_INT) :: nx_name, nx_appear, nx_units, Level, nxx, nyx
    real(kind=MOLCAS_C_REAL) :: data_(*)
  end subroutine xml_ddumpc
end interface
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call xml_ddumpc(TagName,len(TagName),Appear,len(Appear),Units,len(Units),Level,Content,nx,ny)
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_dDump
