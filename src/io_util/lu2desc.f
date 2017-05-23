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
      Subroutine LU2DESC(Lu,Desc)
      Implicit Integer (A-Z)
#include "ctl.fh"

      Integer handle,n

      handle=lu2handle(Lu)

      n=1
10    If(CtlBlk(pHndle,n).ne.handle) Then
        n=n+1
        If(n.gt.MxFile) Then
           Return
        End If
        Go To 10
      End If
      nFile=n
      Desc=CtlBlk(pDesc,nFile)

      Return
      End


      integer function lu2handle(lu)

#include "fio.fh"

      integer lu

      lu2handle=FSCB(lu)

      Return
      End
