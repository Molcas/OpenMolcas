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
      Subroutine ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
c
c     Author: Y. Carissan
c
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Logical Debug
c
      Functional=Zero
      Do iAt=1,nAtoms
        Do iMO_s=1,nOrb2Loc
          Functional=Functional+PA(iMO_s,iMO_s,iAt)**2
        End Do
      End Do
c
      If (Debug) Then
         Write(6,*) 'ComputeFunc: Functional: ',Functional
      End If
c
      Return
      End
