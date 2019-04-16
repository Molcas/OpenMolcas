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
      Subroutine Memory_PolPrep(Que,ixx,iyy,izz,irr3,ixxi,iyyi
     &                         ,izzi,iGri,nPol,nPart)
      Implicit Real*8 (a-h,o-z)

#include "numbers.fh"
#include "WrkSpc.fh"

      Character*(*) Que

      nSize=nPart*nPol
      Call GetMem('xx',Que,'Real',ixx,nSize**2)
      Call GetMem('yy',Que,'Real',iyy,nSize**2)
      Call GetMem('zz',Que,'Real',izz,nSize**2)
      Call GetMem('ixx',Que,'Real',ixxi,nSize**2)
      Call GetMem('iyy',Que,'Real',iyyi,nSize**2)
      Call GetMem('izz',Que,'Real',izzi,nSize**2)
      Call GetMem('irr3',Que,'Real',irr3,nSize**2)
      Call GetMem('iGri',Que,'Real',iGri,nSize**2)

      If(Que(1:4).eq.'Allo') then
        call dcopy_(nSize**2,[ZERO],iZERO,Work(ixx),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(iyy),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(izz),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(ixxi),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(iyyi),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(izzi),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(irr3),iONE)
        call dcopy_(nSize**2,[ZERO],iZERO,Work(iGri),iONE)
      Endif

      Return
      End
