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
      Subroutine WR_MOTRA_Info(Lu,iOpt,iDisk,
     &                         TCONEMO,nTCONEMO,ECOR,NSYM,
     &                         NBAS,NORB,NFRO,NDEL,MxSym,BSLBL,nBSLBL)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
      Integer TCONEMO(nTCONEMO), nBas(MxSym), nOrb(MxSym), nFro(MxSym),
     &        nDel(MxSym)
      Character BSLBL(nBSLBL)*1
*
      Call iDafile(Lu,iOpt,TCONEMO,nTCONEMO,iDisk)
      Call dDafile(Lu,iOpt,ECor,   1,       iDisk)
      Call iDafile(Lu,iOpt,nSym,   1,       iDisk)
      Call iDafile(Lu,iOpt,nBas,   MxSym,   iDisk)
      Call iDafile(Lu,iOpt,nOrb,   MxSym,   iDisk)
      Call iDafile(Lu,iOpt,nFro,   MxSym,   iDisk)
      Call iDafile(Lu,iOpt,nDel,   MxSym,   iDisk)
      Call cDafile(Lu,iOpt,BSLBL,  nBSLBL,  iDisk)
*
      Return
      End
