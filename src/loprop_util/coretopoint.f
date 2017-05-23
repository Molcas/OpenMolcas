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
      Subroutine CoreToPoint(nAt,ipMP,iTP)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension ByggareBas(6)

      Logical GoHere

      Data ByggareBas/2.0d0,8.0d0,8.0d0,18.0d0,18.0d0,32.0d0/

*
*-- A crude but to the point algorithm to put core electrons and
*   nuclei together, and separating the presumably more diffuse
*   part of the charge distribution in a separate chunk.
*
      Call GetMem('NucC','Allo','Real',iNucCh,nAt)
      Call Get_dArray('Nuclear charge',Work(iNucCh),nAt)
      kAt=0
      dToPoint=0.0d0
      Do iAt=1,nAt
        ByggMeraHus=0.0d0
        GoHere=.true.
        dScaleOffSave=Work(iNucCh+iAt-1)
        Do i=1,6
          dScaleOff=dScaleOffSave-ByggareBas(i)
          If(dScaleOff.le.0.0D0.and.GoHere) then
            dToPoint=ByggMeraHus
            GoHere=.false.
          Endif
          dScaleOffSave=dScaleOff
          ByggMeraHus=ByggMeraHus+ByggareBas(i)
        Enddo
        kAt=kAt+iAt
        Work(ipMP+kAt-1)=Work(ipMP+kAt-1)+dToPoint
        Work(iTP+iAt-1)=Work(iNucCh+iAt-1)-dToPoint
      Enddo
      Call GetMem('NucC','Free','Real',iNucCh,nAt)

      Return
      End
