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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      Subroutine GenerateTab_ptr(nAtoms,nOrb2Loc,nBas_Start,Name,
     &                           iTab_ptr,Debug,PA)
c
c     Author: Y. Carissan.
c
      Implicit Real*8 (a-h,o-z)
      Logical Debug
#include "WrkSpc.fh"
#include "Molcas.fh"
      Character*(LENIN8) Name(*),PALbl
      Integer iTab_ptr(*),nBas_Start(*)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)

      nSize=nOrb2Loc**2
      Do iAt=1,nAtoms
        If (Debug) Write(6,*) 'Atom : ',iAt
        PALbl='PA__'//Name(nBas_Start(iAt))(1:LENIN)
        Call GetMem(PALbl,'Allo','Real',ip,nSize)
        Call dCopy_(nSize,[0.0d0],0,Work(ip),1)
        iTab_ptr(iAt)=ip
        If (Debug) Write(6,*) 'gen Atom ip',iAt,ip
      End Do

      Return
      End
