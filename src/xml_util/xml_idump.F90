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
! This routine dumps integers in xml format.                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_iDump(TagName,Appear,Units,Level,Content,nx,ny)

use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName, Appear, Units
integer(kind=iwp), intent(in) :: Content(*), nx, ny, Level
interface
  subroutine xml_iDumpc(name_,nx_name,appear,nx_appear,units,nx_units,Level,data_,nxx,nyx) bind(C,name='xml_idumpc_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: name_(*), appear(*), units(*)
    integer(kind=MOLCAS_C_INT) :: nx_name, nx_appear, nx_units, Level, data_(*), nxx, nyx
  end subroutine xml_iDumpc
end interface
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
call xml_iDumpc(TagName,len(TagName),Appear,len(Appear),Units,len(Units),Level,Content,nx,ny)
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_iDump
