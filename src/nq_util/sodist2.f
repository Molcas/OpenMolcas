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
      Subroutine SODist2(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
     &                  nSOs,iAO,TmpCMOs,nCMO,TmpDoIt)
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      Implicit None
      Integer mAO, nCoor, mBas, nCmp, nDeg, nSOs, iAO, nCMO
      Real*8  SOValue(mAO*nCoor,mBas,nCmp*nDeg),
     &        SO(mAO*nCoor,nSOs),
     &        TmpCMOs(nCMO)
      Integer TmpDoIt(nSOs)
*     Local
      Integer k, iOff, i, j, iBas, ii
*
      Do k=1,nSOs
         TmpDoIt(k) = 1
      End Do
*
      TmpCMOs(:)=0.0D0
*
      iOff=0
      Do i=0,nIrrep-1
         iBas=nBas(i)
         Do j=1,iBas
            ii = iOff+(j-1)*iBas + j
            TmpCMOs(ii) = 1.0d0
         End Do
         iOff=iOff+nBas(i)**2
      End Do
*
      Call SODist(SOValue,mAO,nCoor,mBas,nCmp,nDeg,SO,
     &            nSOs,iAO,TmpCMOs,nCMO,TmpDoIt)
*
      Return
      End

