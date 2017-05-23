************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************

      SubRoutine ChoMP2g_MOReOrd(CMO,COrb1,COrb2,iMoType1,iMOType2)
C
C     Jonas Bostrom, Jan. 2010. (modified from ChoMP2_MOReOrd)
C
C     Purpose: Make CMO:s of appropriate length, transpose COrb2,
C
C              CMO(alpha,p) -> COrb1(p,alpha)
C              CMO(alpha,q) -> COrb2(alpha,q)
C
#include "implicit.fh"
      Real*8 COrb1(*), COrb2(*), CMO(*)
      Integer nOrb1(8), nOrb2(8)
      Integer nOffOrb1(8), nOffOrb2(8)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2g.fh"
#include "choorb.fh"

      Do iSym = 1, nSym
         nOffOrb1(iSym) = 0
         nOffOrb2(iSym) = 0
         Do i = 1, iMOType1-1
            nOffOrb1(iSym) = nOffOrb1(iSym) + nMo(iSym,i)
         End Do
         Do i = 1, iMOType2-1
            nOffOrb2(iSym) = nOffOrb2(iSym) + nMo(iSym,i)
         End Do
         nOrb1(iSym) = nMo(iSym,iMOType1)
         nOrb2(iSym) = nMo(iSym,iMOType2)
      End Do
*
      iCount = 0
      Do iSym = 1,nSym
*
         jCount = iCount + nOffOrb1(iSym)*nBas(iSym)
*
         Do i = 1,nOrb1(iSym)
            kOff1 = jCount + nBas(iSym)*(i-1) + 1
            kOff2 = iMoAo(iSym,iSym,iMoType1) + i
           Call dCopy_(nBas(iSym),CMO(kOff1),1,COrb1(kOff2),nOrb1(iSym))
         End Do

         kOff1 = iCount + nOffOrb2(iSym)*nBas(iSym) + 1
         kOff2 = iAoMo(iSym,iSym,iMoType2) + 1
         Call dCopy_(nBas(isym)*nOrb2(iSym),CMO(kOff1),1,COrb2(kOff2),1)

         iCount = iCount + nBas(iSym)*nBas(iSym)

      End Do

      End
