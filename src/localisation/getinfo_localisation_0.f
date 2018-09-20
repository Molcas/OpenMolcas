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
      SubRoutine GetInfo_Localisation_0()
C
C     Author: T.B. Pedersen
C
C     Purpose: read basic info from runfile and INPORB.
C
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "inflocal.fh"
#include "WrkSpc.fh"

      Character*22 SecNam
      Parameter (SecNam = 'GetInfo_Localisation_0')
      Character*80 Txt
      Character*512 FName

C     Read number of irreps.
C     ----------------------

      Call Get_iScalar('nSym',nSym)
      If (nSym.lt.1 .or. nSym.gt.MxSym) Then
         Write(Txt,'(A,I9)') 'nSym =',nSym
         Call SysAbendMsg(SecNam,'Number of irreps out of bounds!',Txt)
      End If

C     Read number of basis functions.
C     -------------------------------

      Call Get_iArray('nBas',nBas,nSym)
      nBasT = nBas(1)
      Do iSym = 2,nSym
         nBasT = nBasT + nBas(iSym)
      End Do
      If (nBasT.lt.1 .or. nBasT.gt.MxBas) Then
         Write(Txt,'(A,I9)') 'nBasT =',nBasT
         Call SysAbendMsg(SecNam,'Basis set limits exceeded!',Txt)
      End If

C     Set number of orbitals equal to nBas.
C     -------------------------------------

      Call Icopy(nSym,nBas,1,nOrb,1)
      nOrbT = nOrb(1)
      Do iSym = 2,nSym
         nOrbT = nOrbT + nOrb(iSym)
      End Do
      If (nOrbT.lt.1 .or. nOrbT.gt.MxBas) Then
         Write(Txt,'(A,I9)') 'nOrbT =',nOrbT
         Call SysAbendMsg(SecNam,'Orbital limits exceeded!',Txt)
      End If
      Do iSym=1,nSym
         If (nOrb(iSym) .gt. nBas(iSym)) Then
            Write(Txt,'(A,I2,2(1X,I9))')
     &      'iSym,nOrb,nBas:',iSym,nOrb(iSym),nBas(iSym)
            Call SysAbendMsg(SecNam,'#orb > #bas:',Txt)
         End If
      End Do

C     Read MO coefficients, orbital occupations, and orbital energies
C     from INPORB.
C     ---------------------------------------------------------------

      n2Bas = nBas(1)**2
      Do iSym = 2,nSym
         n2Bas = n2Bas + nBas(iSym)**2
      End Do

      nCMO = n2Bas
      lOcc = nBasT
      lEor = nBasT
      lInd = nBasT
      Call GetMem('CMO','Allo','Real',ipCMO,nCMO)
      Call GetMem('Occup','Allo','Real',ipOcc,lOcc)
      Call GetMem('OrbEn','Allo','Real',ipEor,lEor)
      Call GetMem('IndT','Allo','Inte',ipInd,lInd)
      FName=LC_FileOrb
      If (mylen(FName).eq.0) FName='INPORB' ! file name
      Call RdVec_Localisation(nSym,nBas,nOrb,iWork(ipInd),
     &                        Work(ipCMO),Work(ipOcc),Work(ipEor),
     &                        FName(:mylen(FName)))

C     Set number of occupied and virtual orbitals according to the
C     occupation numbers from INPORB (assuming that occupied orbitals
C     precede virtual ones on file).
C     ---------------------------------------------------------------

      kOff = ipOcc - 1
      Do iSym = 1,nSym
         nOccInp(iSym) = 0
         i = 0
         Do While (i .lt. nOrb(iSym))
            i = i + 1
            If (Work(kOff+i) .gt. 0.0d0) Then
               nOccInp(iSym) = nOccInp(iSym) + 1
            Else
               i = nOrb(iSym) ! break while loop
            End If
         End Do
         nVirInp(iSym) = nOrb(iSym) - nOccInp(iSym)
         If (nVirInp(iSym) .lt. 0) Then
            Write(Txt,'(3(A,I9))')
     &      'No. of occupied: ',nOccInp(iSym),
     &      ' No. of orbitals: ',nOrb(iSym),
     &      ' Symmetry: ',iSym
            Call SysAbendMsg(SecNam,'#occ > #orb:',Txt)
         End If
         kOff = kOff + nBas(iSym)
      End Do

C     Read number of atoms, atomic labels, and basis function labels
C     from runfile.
C     --------------------------------------------------------------

      Call Get_nAtoms_All(nAtoms)
      If (nAtoms.lt.1 .or. nAtoms.gt.MxAtom) Then
         Write(Txt,'(A,I9)') 'nAtoms =',nAtoms
         Call SysAbendMsg(SecNam,'Atom limit exceeded!',Txt)
      End If
c     Call Get_cArray('Unique Atom Names',AtomLbl,4*nAtoms)
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBasT)
      Do i = 1,nBasT
         Atom(i) = Name(i)(1:LENIN)
         Type(i) = Name(i)(LENIN1:LENIN8)
      End Do

      End
