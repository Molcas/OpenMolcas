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
      FUNCTION nPre(ks)

#include "Input.fh"
      iOut=0
      Do is=1,nSym
       jS=iEOR(iS-1,kS-1)+1
       nRest=nOrb(js)-nIsh(js)
       iOut=iOut+nIsh(is)*nRest*(nRest+1)
       nRest=nOrb(js)-nRs1(js)
       iOut=iOut+nRs1(is)*nRest*(nRest+1)
       nRest=nOrb(js)-nRs2(js)
       iOut=iOut+nRs2(is)*nRest*(nRest+1)
       nRest=nOrb(js)-nRs3(js)
       iOut=iOut+nRs3(is)*nRest*(nRest+1)
      End Do
      nPre=iOut
      Return
      End
