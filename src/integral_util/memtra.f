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
c -------------------------------------------------------------------
c The following function returns the size needed by ptrans for
c temporaries. It also puts into common some offsets and stuff.
c -------------------------------------------------------------------
      integer function memtra(npam)
      use pso_stuff
#include "etwas.fh"
      Integer nPam(4,0:7)
      intrinsic max
*
*     iQ = 1
*     Call qEnter('MemTra')
      mxact=0
      mxS1=0
      mxS2=0
      mxS3=0
      mxS4=0
      do 10 isym=0,mirrep-1
        na=nash(isym)
        if(na.eq.0) goto 10
        mxact=max(mxact,na)
        mxS1=max(mxS1,npam(1,isym))
        mxS2=max(mxS2,npam(2,isym))
        mxS3=max(mxS3,npam(3,isym))
        mxS4=max(mxS4,npam(4,isym))
  10  continue
      mxS=max(mxS1,mxS2,mxS3,mxS4)
      mxa2=mxact*mxact
      mxa3=mxa2*mxact
      mxa4=mxa3*mxact
      mxS34=mxS3*mxS4
      mxS234=mxS2*mxS34
c Max sizes, in common, needed for certain temporaries:
      ncred=Max(1,mxS*mxact)

      nscr1=mxa4
      nscr2=mxa3*mxS4
      nscr3=mxa2*mxS34
      nscr4=mxact*mxS234
      nscr5=mxS1*mxS234
      nScr1=max(1,nscr1,nscr3,nscr5)
      nScr2=max(1,nscr2,nscr4)
      memtra=nCred+2*nScr1+nScr2+3
*
*     Call qExit('MemTra')
      Return
      End
