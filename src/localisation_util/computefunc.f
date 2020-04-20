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
#define _TEST1_
#ifdef _TEST1_
      Subroutine ComputeFunc(nAtoms,nOrb2Loc,iTab_ptr,PA,Functional,
     &                       Debug)
#else
      Subroutine ComputeFunc(nAtoms,nOrb2Loc,iTab_ptr,Functional,
     &                       Debug)
#endif
c
c     Author: Y. Carissan
c
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "real.fh"
#ifdef _TEST1_
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
#endif
      Integer iTab_ptr(*)
      Logical Debug
c
#ifdef _TEST1_
      Functional=Zero
      Do iAt=1,nAtoms
        Do iMO_s=1,nOrb2Loc
          Functional=Functional+PA(iMO_s,iMO_s,iAt)**2
        End Do
      End Do
*        Write(6,*) 'ComputeFunc: Functional: ',Functional
#else
      Functional=Zero
      Do iAt=1,nAtoms
        ip=iTab_ptr(iAt)
        Do iMO_s=1,nOrb2Loc
          Functional=Functional+(Work(ip+(iMO_s-1)*nOrb2Loc+iMO_s-1))**2
        End Do
      End Do
*        Write(6,*) 'ComputeFunc: Functional: ',Functional
#endif
c
      If (Debug) Then
         Write(6,*) 'ComputeFunc: Functional: ',Functional
      End If
c
      Return
      End
