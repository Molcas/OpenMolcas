************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
*----------------------------------------------------------------------*
* This subroutine reads from the formatted output of mpprop the        *
* coordinates of the expansion centers.                                *
*----------------------------------------------------------------------*
      Subroutine Get_Centers(nAt,xyz)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "warnings.fh"

      Dimension xyz(MxAt,MxAt,3)
      Character*13 TheLine
      Logical Exist

*Open the file
      Lu=40
      Lu=IsFreeUnit(40)
      Call Opnfl('MPPROP',Lu,Exist)
      If(.not.Exist) then
        Write(6,*)
        Write(6,*)' Can not locate output file from MpProp. '
        Call Quit(_RC_IO_ERROR_READ_)
      Endif
      Rewind(Lu)

*Read until you get standard line
10    Continue
        Read(Lu,'(A)') TheLine
      If(TheLine.ne.'* All centers') Go To 10
      Read(Lu,*) nCent

*Read atom centers.
      Do 15, i=1,nAt
        Read(Lu,'(A)') TheLine
        Read(Lu,*)(xyz(i,i,k),k=1,3)
        Do 25, j=1,10
          Read(Lu,'(A)') TheLine
25      Continue
15    Continue

*Read bond centers.
      Do 30, i=2,nAt
        Do 32, j=1,i-1
          Read(Lu,'(A)') TheLine
          Read(Lu,*)(xyz(i,j,k),k=1,3)
          Do 35, jj=1,10
            Read(Lu,'(A)') TheLine
35        Continue
32      Continue
30    Continue

*Square xyz for later convinience
      Do 40, i=2,nAt
        Do 42, j=1,i-1
          xyz(j,i,1)=xyz(i,j,1)
          xyz(j,i,2)=xyz(i,j,2)
          xyz(j,i,3)=xyz(i,j,3)
42      Continue
40    Continue

      Close(Lu)
      Return
      End
