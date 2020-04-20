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
*               Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine PipekMezey(Functional,CMO,Thrs,ThrRot,ThrGrad,
     &                      Name,
     &                      nBas,nOrb2Loc,nFro,
     &                      nSym,nAtoms,nMxIter,
     &                      Maximisation,Converged,Debug,Silent)
C
C     Author: Y. Carissan [modified by T.B. Pedersen].
C
C     Purpose: Pipek-Mezey localisation of occupied orbitals.
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Real*8  CMO(*)
      Integer nBas(nSym), nOrb2Loc(nSym), nFro(nSym)
      Logical Maximisation, Converged, Debug, Silent
      Character*(LENIN8) Name(*) ! dimension should be tot. #bf
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: PA(:,:,:)

      Character*10 SecNam
      Parameter (SecNam = 'PipekMezey')

      Character*8 Label

C     Symmetry is NOT allowed!!
C     -------------------------

      If (nSym .ne. 1) Then
         Call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
      End If

C     Initializations.
C     ----------------

      Functional = -9.9D9

      nBasT     = nBas(1)
      nOrb2LocT = nOrb2Loc(1)
      nFroT     = nFro(1)

      Converged = .False.

C     Read overlap matrix.
C     --------------------

      lOaux = nBasT*(nBasT+1)/2 + 4
      lOvlp = nBasT**2
      Call GetMem('Ovlp','Allo','Real',ipOvlp,lOvlp)
      Call GetMem('AuxOvlp','Allo','Real',ipOaux,lOaux)

      irc    = -1
      iOpt   = 2
      iComp  = 1
      iSyLbl = 1
      Label  = 'Mltpl  0'
      Call RdOne(irc,iOpt,Label,iComp,Work(ipOaux),iSyLbl)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': RdOne returned ',irc
         Write(6,*) 'Label = ',Label,'  iSyLbl = ',iSyLbl
         Call SysAbendMsg(SecNam,'I/O error in RdOne',' ')
      End If

      If (Debug) Then
         Write (6,*)
         Write (6,*) ' Triangular overlap matrix at start'
         Write (6,*) ' ----------------------------------'
         Call TriPrt('Overlap',' ',Work(ipOaux),nBasT)
      End If

      Call Tri2Rec(Work(ipOaux),Work(ipOvlp),nBasT,Debug)
      Call GetMem('AuxOvlp','Free','Real',ipOaux,lOaux)

C     Allocate and get index arrays for basis functions per atom.
C     -----------------------------------------------------------

      l_nBas_per_Atom = nAtoms
      l_nBas_Start    = l_nBas_per_Atom
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 Name,nBasT,nAtoms,Debug)

C     Allocate PA array.
C     ------------------
      Call mma_Allocate(PA,nOrb2LocT,nOrb2LocT,nAtoms,Label='PA')
      PA(:,:,:)=0.0D0

C     Localise orbitals.
C     ------------------

      kOffC = nBasT*nFroT + 1
      Call PipekMezey_Iter(Functional,CMO(kOffC),
     &                     Work(ipOvlp),Thrs,ThrRot,ThrGrad,PA,
     &                     iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                     Name,nBasT,nOrb2LocT,nAtoms,nMxIter,
     &                     Maximisation,Converged,Debug,Silent)

C     De-allocations.
C     ---------------

      Call mma_deallocate(PA)
      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

      End
