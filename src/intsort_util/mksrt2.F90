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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1991, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MkSrt2()
!***********************************************************************
!                                                                      *
!     Purpose: Initialize counters and offsets required                *
!              for bin sorting algorithm                               *
!                                                                      *
!     Called from: Sort1                                               *
!                                                                      *
!     Calls to : none                                                  *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use sort_data, only: iDIBin, iDVBin, mSyBlk, n_Int, nRec, nSln
use Definitions, only: iwp, u6

implicit none
#include "print.fh"
integer(kind=iwp) :: iBin, iPrint, iRout, iSlice, iSyBlk, nSlice

iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 10) write(u6,*) ' >>> Enter MKSRT2 <<<'

!----------------------------------------------------------------------*
!     initialize various pointers, counters and disk adresses          *
!----------------------------------------------------------------------*

iBin = 0
do iSyBlk=1,mSyBlk
  nSlice = nSln(iSyBlk)
  if (nSlice /= 0) then
    do iSlice=1,nSlice
      iBin = iBin+1
      iDIBin(2,iBin) = -1
      iDVBin(2,iBin) = -1
      iDVBin(3,iBin) = -1
      iDVBin(4,iBin) = -1
      n_Int(iBin) = 0
      nRec(iBin) = 0
    end do
  end if
end do

return

end subroutine MkSrt2
