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
      Subroutine RBuf_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS)
      Implicit real*8(a-h,o-z)
*     Integer  MEMX,KBUF,BLKSZ,NBLCK,NPASS,DLNGTH

#include "SysDef.fh"
      Dimension W(*)
*      Call GetMem('MaxMem','MAX','REAL',KBUF,MEMX)
*      BLKSZ=(NOTU-1)*IADXS
*      NBLCK=(MEMX)/BLKSZ

      IST=1
*      if(NBLCK.gt.1.and.LBuf.lt.256) Then
*      Call RBufF_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS,MEMX)
*      Else
      Call RBufO_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS)
*      End If
      IST=1
      Return
      End
