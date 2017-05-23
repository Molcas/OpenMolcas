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
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Domain_Localisation(irc)
C
C     Thomas Bondo Pedersen, January 2006.
C
C     Purpose: set up orbital domains and pair domains. Find number of
C              strong, weak, distant, and very distant pairs.
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "debug.fh"
#include "WrkSpc.fh"

      Character*19 SecNam
      Parameter (SecNam = 'Domain_Localisation')

      Integer iCount(0:3)
      Real*8  ThrPD(3)

C     Set return code.
C     ----------------

      irc = 0

C     Check for symmetry (not allowed).
C     ---------------------------------

      If (nSym .ne. 1) Then
         irc = -1
         Return
      End If

C     Initializations.
C     ----------------

      nBasT = nBas(1)
      nOcc  = nOrb2Loc(1)
      nnOcc = nOcc*(nOcc+1)/2
      nAtom = nAtoms

      l_nBas_per_Atom = 0
      l_nBas_Start = 0
      l_iDomain = 0
      l_QD = 0
      l_f = 0
      l_iPairDomain = 0
      l_iClass = 0
      l_Rmin = 0
      l_Coord = 0

C     There must be at least 2 atoms and 2 orbitals.
C     ----------------------------------------------

      If (nAtom.lt.2 .or. nOcc.lt.2) Then
         irc = -2
         Return
      End If

C     Allocate and get index arrays for indexation of basis functions on
C     each atom.
C     ------------------------------------------------------------------

      l_nBas_per_Atom = nAtom
      l_nBas_Start    = nAtom
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 Name,nBasT,nAtom,Debug)

C     Define domains.
C     ---------------

      l_iDomain = (nAtom+1)*nOcc
      l_QD = nOcc
      l_f = nOcc
      Call GetMem('iDomain','Allo','Inte',ip_iDomain,l_iDomain)
      Call GetMem('QD','Allo','Real',ip_QD,l_QD)
      Call GetMem('f','Allo','Real',ip_f,l_f)

      kC = ipCMO + nBasT*nFro(1)
      Call DefineDomain(irc,iWork(ip_iDomain),Work(ip_QD),Work(ip_f),
     &                  Work(kC),ThrDomain,
     &                  iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                  nAtom,nBasT,nOcc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ERROR: DefineDomain returned ',irc
         Go To 1 ! return after deallocations
      End If

      If (Debug) Then
         Write(6,*) SecNam,': checking domain definitions...'
         Call CheckDomain(irc,iWork(ip_iDomain),nAtom,nOcc)
         If (irc .eq. 0) Then
            Write(6,*) '....OK!'
         Else
            Write(6,*) '....Ooops. Buggy domain definition!'
            irc = 2
            Go To 1 ! return after deallocations
         End If
      End If

C     Define pair domains.
C     Make sure that ThrPairDomain is in ascending order.
C     ---------------------------------------------------

      Call dCopy_(3,ThrPairDomain,1,ThrPD,1)
      Call Cho_Order(ThrPD,3,1)
      iChange = 0
      i = 0
      Do While (i.lt.2 .and. iChange.eq.0)
         i = i + 1
         Tst = ThrPairDomain(i) - ThrPD(i)
         If (abs(Tst) .gt. 1.0d-15) Then
            iChange = 1
         End If
      End Do

      l_iPairDomain = (nAtom+1)*nnOcc
      l_iClass = nnOcc
      l_Rmin = nnOcc
      l_Coord = 3*nAtom
      Call GetMem('iPairDomain','Allo','Inte',ip_iPairDomain,
     &                                        l_iPairDomain)
      Call GetMem('iClass','Allo','Inte',ip_iClass,l_iClass)
      Call GetMem('Rmin','Allo','Real',ip_Rmin,l_Rmin)
      Call GetMem('Coord','Allo','Real',ip_Coord,l_Coord)

      Call Get_dArray('Unique Coordinates',Work(ip_Coord),l_Coord)
      Call DefinePairDomain(irc,iWork(ip_iPairDomain),iWork(ip_iClass),
     &                      Work(ip_Rmin),iWork(ip_iDomain),ThrPD,
     &                      Work(ip_Coord),nAtom,nOcc,3)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ERROR: DefinePairDomain returned ',irc
         Go To 1 ! return after deallocations
      End If

      If (Debug) Then
         Write(6,*) SecNam,': checking pair domain definitions...'
         Call CheckDomain(irc,iWork(ip_iPairDomain),nAtom,nnOcc)
         If (irc .eq. 0) Then
            Write(6,*) '....OK!'
         Else
            Write(6,*) '....Ooops. Buggy pair domain definition!'
            irc = 3
            Go To 1 ! return after deallocations
         End If
      End If

C     Print info.
C     -----------

      Call Domain_Histogram(iWork(ip_iDomain),nAtom,nOcc,
     &                      'Histogram of domain sizes')
      Call Domain_Histogram(iWork(ip_iPairDomain),nAtom,nnOcc,
     &                      'Histogram of pair domain sizes')

      Call Cho_Head('Pair domain classification','=',80,6)
      Do i = 0,3
         iCount(i) = 0
      End Do
      Do ij = 0,nnOcc-1
         iC = iWork(ip_iClass+ij)
         iCount(iC) = iCount(iC) + 1
      End Do
      Write(6,'(/,A)') 'Definition:'
      If (iChange .ne. 0) Then
         Write(6,'(A,A)')
     &   'Notice: the input thresholds were re-ordered to ascending ',
     &   'order'
         Write(6,'(A,1P,3(1X,D15.5))') 'Your input order was:',
     &   (ThrPairDomain(i),i=1,3)
      End If
      Write(6,'(A,1P,D15.5)')
     & 'Strong       pairs:                   R <= ',ThrPD(1)
      Write(6,'(A,1P,D15.5,A,D15.5)')
     & 'Weak         pairs: ',ThrPD(1),' < R <= ',ThrPD(2)
      Write(6,'(A,1P,D15.5,A,D15.5)')
     & 'Distant      pairs: ',ThrPD(2),' < R <= ',ThrPD(3)
      Write(6,'(A,1P,D15.5,A)')
     & 'Very distant pairs: ',ThrPD(3),' < R'
      Write(6,'(/,A)') 'Classification:'
      Fac = 1.0d2/dble(nnOcc)
      Write(6,'(A,I9,3X,F7.2,A)')
     & 'Number of strong       pairs: ',iCount(0),
     & Fac*dble(iCount(0)),'%'
      Write(6,'(A,I9,3X,F7.2,A)')
     & 'Number of weak         pairs: ',iCount(1),
     & Fac*dble(iCount(1)),'%'
      Write(6,'(A,I9,3X,F7.2,A)')
     & 'Number of distant      pairs: ',iCount(2),
     & Fac*dble(iCount(2)),'%'
      Write(6,'(A,I9,3X,F7.2,A,/)')
     & 'Number of very distant pairs: ',iCount(3),
     & Fac*dble(iCount(3)),'%'

C     Analysis of individual domains (if requested).
C     ----------------------------------------------

      If (AnaDomain) Then
         Call Analysis_Domain(iWork(ip_iDomain),Work(ip_QD),Work(ip_f),
     &                        Work(ip_Coord),Name,
     &                        iWork(ip_nBas_Start),nAtom,nBas,nOcc)
      End If

C     Deallocations.
C     --------------

    1 If (l_Coord .gt. 0) Then
         Call GetMem('Coord','Free','Real',ip_Coord,l_Coord)
      End If
      If (l_Rmin .gt. 0) Then
         Call GetMem('Rmin','Free','Real',ip_Rmin,l_Rmin)
      End If
      If (l_iClass .gt. 0) Then
         Call GetMem('iClass','Free','Inte',ip_iClass,l_iClass)
      End If
      If (l_iPairDomain .gt. 0) Then
         Call GetMem('iPairDomain','Free','Inte',ip_iPairDomain,
     &                                        l_iPairDomain)
      End If
      If (l_f .gt. 0) Then
         Call GetMem('f','Free','Real',ip_f,l_f)
      End If
      If (l_QD .gt. 0) Then
         Call GetMem('QD','Free','Real',ip_QD,l_QD)
      End If
      If (l_iDomain .gt. 0) Then
         Call GetMem('iDomain','Free','Inte',ip_iDomain,l_iDomain)
      End If
      If (l_nBas_Start .gt. 0) Then
         Call GetMem('nB_Start','Free','Inte',
     &               ip_nBas_Start,l_nBas_Start)
      End If
      If (l_nBas_per_Atom .gt. 0) Then
         Call GetMem('nB_per_Atom','Free','Inte',
     &               ip_nBas_per_Atom,l_nBas_per_Atom)
      End If

      End
