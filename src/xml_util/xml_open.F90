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
! This routine opens an xml container.                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_Open(TagName,Appear,Units,Level,Content)

use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName, Appear, Units, Content
integer(kind=iwp), intent(in) :: Level
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
character(len=16) :: myName
interface
  subroutine xml_Openc(name_,nx_name,appear,nx_appear,units,nx_units,Level,value_,nx_value) bind(C,name='xml_openc_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: name_(*), appear(*), units(*), value_(*)
    integer(kind=MOLCAS_C_INT) :: nx_name, nx_appear, nx_units, Level, nx_value
  end subroutine xml_Openc
end interface
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
myName = TagName
call Upcase(myName)
if (myName == 'MODULE') then
  call poke_iScalar("xml opened",1)
end if
call xml_Openc(TagName,len(TagName),Appear,len(Appear),Units,len(Units),Level,Content,len(Content))
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_Open
