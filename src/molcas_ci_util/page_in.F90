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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine page_in(KeyWord,nConf,Vector,LuDavid)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Load any vector for further use by the Davidson diagonalization  *
!     which has been saved by the write through cache mechanism        *
!                                                                      *
!     calling arguments:                                               *
!     KeyWord : character*16                                           *
!               record identifier                                      *
!     nConf   : integer                                                *
!               length of the vector H_diag                            *
!     Vector  : array of real*8                                        *
!               any vector of length nConf                             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
real*8 Vector(nConf)
character*16 KeyWord
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"

! check input arguments
if (nConf < 0) then
  write(6,*) 'page_in: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if

! search for a metching record identifier
nStk = 0
do iStk=1,(mxMemStk+mxDiskStk)
  if (LblStk(iStk) == KeyWord) nStk = iStk
end do
if (nStk == 0) then
  write(6,*) 'page_in: nStk equal 0'
  write(6,*) 'nStk = ',nStk
  call Abend()
end if

if (nStk <= mxMemStk) then
  iMem = memory_address(nStk)
  call dCopy_(nConf,Work(iMem),1,Vector,1)
else
  iDisk = disk_address(nStk-mxMemStk)
  call DDaFile(LuDavid,2,Vector,nConf,iDisk)
end if

return

end subroutine page_in
