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
! This routine closes an xml container.                                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

subroutine xml_Close(TagName)

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
character(len=*), intent(in) :: TagName
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
character(len=16) :: myName
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
myName = TagName
call Upcase(myName)
if (myName == 'MODULE') then
  call poke_iScalar('xml opened',0)
end if
call xml_Closec(TagName,len(TagName))
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine xml_Close
