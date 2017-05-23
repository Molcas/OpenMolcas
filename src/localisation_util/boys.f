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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine Boys(Functional,CMO,Thrs,ThrRot,ThrGrad,
     &                nBas,nOrb2Loc,nFro,nSym,nMxIter,
     &                Maximisation,Converged,Debug,Silent)
C
C     Author: T.B. Pedersen
C
C     Purpose: Boys localisation of occupied orbitals.
C
      Implicit Real*8 (a-h,o-z)
      Real*8   CMO(*)
      Integer nBas(nSym), nOrb2Loc(nSym), nFro(nSym)
      Logical Maximisation, Converged, Debug, Silent
#include "WrkSpc.fh"

      Character*4 SecNam
      Parameter (SecNam = 'Boys')

      Parameter (nComp = 3) ! 3 components of dipole operator
      Integer ipLbl(nComp), ipLbl_MO(nComp)
      Character*8 Label, AlloLbl(nComp), AlloLbl_MO(nComp)

C     Symmetry is NOT allowed!!
C     -------------------------

      If (nSym .ne. 1) Then
         Call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
      End If

C     Initializations.
C     ----------------

      Functional = -9.9d9

      nBasT     = nBas(1)
      nOrb2LocT = nOrb2Loc(1)
      nFroT     = nFro(1)

      Converged = .False.

C     Read AO dipole moment integrals.
C     --------------------------------

      lLbl = nBasT*nBasT
      Do iComp = 1,nComp
         Write(AlloLbl(iComp),'(A,I2)') 'Dipole',iComp
         Call GetMem(AlloLbl(iComp),'Allo','Real',ipLbl(iComp),lLbl)
      End Do

      lAux = nBasT*(nBasT+1)/2 + 4
      Call GetMem('DipAux','Allo','Real',ipAux,lAux)
      Label = 'Mltpl  1'
      Do iComp = 1,nComp
         irc = -1
         iOpt = 2
         iSym = 1
         Call RdOne(irc,iOpt,Label,iComp,Work(ipAux),iSym)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': RdOne returned ',irc
            Write(6,*) 'Label = ',Label,'   Component = ',iComp
            Call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
         End If
         If (Debug) Then
            Write(6,*)
            Write(6,*) ' Triangular dipole matrix at start'
            Write(6,*) ' ---------------------------------'
            Write(6,*) ' Component: ',iComp
            Call TriPrt(' ',' ',Work(ipAux),nBasT)
         End If
         Call Tri2Rec(Work(ipAux),Work(ipLbl(iComp)),nBasT,Debug)
      End Do
      Call GetMem('DipAux','Free','Real',ipAux,lAux)

C     Allocate MO arrays.
C     -------------------

      lLbl_MO = nOrb2LocT*nOrb2LocT
      Do iComp = 1,nComp
         Write(AlloLbl_MO(iComp),'(A,I2)') 'MO dip',iComp
         Call GetMem(AlloLbl_MO(iComp),'Allo','Real',ipLbl_MO(iComp),
     &                                               lLbl_MO)
      End Do

C     Localise orbitals.
C     ------------------

      kOffC = nBasT*nFroT + 1
      Call Boys_Iter(Functional,CMO(kOffC),Thrs,ThrRot,ThrGrad,
     &               ipLbl,ipLbl_MO,
     &               nBasT,nOrb2LocT,nComp,nMxIter,
     &               Maximisation,Converged,Debug,Silent)

C     De-allocations.
C     ---------------

      Do iComp = nComp,1,-1
         Call GetMem(AlloLbl_MO(iComp),'Free','Real',ipLbl_MO(iComp),
     &                                               lLbl_MO)
      End Do
      Do iComp = nComp,1,-1
         Call GetMem(AlloLbl(iComp),'Free','Real',ipLbl(iComp),lLbl)
      End Do

      End
