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

subroutine RdTraOne()
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     read the header of the transformed one-electron integral file    *
!                                                                      *
!***********************************************************************

#include "Molcas.fh"
#include "motra.fh"
#include "SysDef.fh"
integer iDisk, LuTraOne
integer TocTraOne(64)
character*(LENIN8) BsLbl(MxOrb)

LuTraOne = 3

call DaName(LuTraOne,'TRAONE')

!ulf
iDisk = 0

call WR_MOTRA_Info(LuTraOne,2,iDisk,TocTraOne,64,EcorX,nSymX,nBasX,nOrbX,nFroX,nDelX,8,BsLbl,LENIN8*MxOrb)

call Daclos(LuTraOne)

return

end subroutine RdTraOne
