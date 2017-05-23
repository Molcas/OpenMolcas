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
      Subroutine RBufO_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS)
      Implicit real*8(a-h,o-z)

#include "SysDef.fh"
      Dimension W(*)
      Call GetMem('MaxMem','MAX','REAL',KDUM,MEMX)
      IADX=(KKTU-1)*IADXS
      IST=1
      Length=LBuf
      IEnd=LBuf
   52 CALL dDAFILE(LUHLFX,2,W(IST),Length,IADX)
      IST=IST+LBuf
      IEnd=IEnd+LBuf
      If(IEnd.gt.LL)Length=mod(LL,LBuf)
      IADX=IADX+(NOTU-1)*IADXS
      IF(IST.LE.LL) GO TO 52
      IST=1
      Return
      End
