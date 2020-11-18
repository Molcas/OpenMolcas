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
      Subroutine cp_SpcInt
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Character*14  qLbl, filnam*16
      Logical Exist
      Real*8, Allocatable:: rK(:)

*     Write (*,*) ' Copy files'
*
*.... Copy file SPCINX to SPCIN1
*
      filnam='SPCINX'
      Call f_Inquire(filnam,Exist)
      If (Exist) Then
      LuTmp1=11
      LuTmp2=12
      call molcas_binaryopen_vanilla(luTmp1,filnam)
c      Open(luTmp1,File=filnam,Form='unformatted',Status='unknown')
      filnam='SPCIN1'
      call molcas_binaryopen_vanilla(luTmp2,filnam)
c      Open(luTmp2,File=filnam,Form='unformatted',Status='unknown')
      ReWind (LuTmp1)
      ReWind (LuTmp2)
*
      Read (LuTmp1) nq,nQQ
      Write(LuTmp2) nq,nQQ
      Call mma_allocate(rK,nQQ,Label='rK')
      Do iq = 1, nq
         Read (LuTmp1) qLbl,(rK(i),i=1,nQQ)
         Write(LuTmp2) qLbl,(rK(i),i=1,nQQ)
      End Do
      Call mma_deallocate(rK)
*
      Close  (LuTmp1)
      Close  (LuTmp2)
      End If
*
      Return
      End
