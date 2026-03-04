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

subroutine Ext_Rank(FileName)

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank
use Definitions, only: iwp
#endif

implicit none
character(len=*), intent(inout) :: FileName
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: Length, NameLength
integer(kind=iwp), external :: StrnLn
#endif

#ifdef _MOLCAS_MPP_
Length = len(FileName)
NameLength = StrnLn(FileName)
Filename(NameLength+1:NameLength+1) = '_'
call WrNumber(FileName(NameLength+2:Length),MyRank)
#else
#include "macros.fh"
unused_var(FileName)
#endif

end subroutine Ext_Rank
