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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Sort_Localisation(CMO,nBas,nOcc,nFro,nSym)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: sort CMOs according to Cholesky orbital ordering.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(*)
      Integer nBas(nSym), nOcc(nSym), nFro(nSym)
#include "WrkSpc.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Sort_Localisation')

      Character*8  Label
      Character*80 Txt

C     Static setting of decomposition threshold.
C     ------------------------------------------

      ThrCho = 1.0d-12

C     Get a copy of the occupied orbitals: X=CMO.
C     -------------------------------------------

      lX = nBas(1)*nOcc(1)
      Do iSym = 2,nSym
         lX = lX + nBas(iSym)*nOcc(iSym)
      End Do
      Call GetMem('XCho','Allo','Real',ipX,lX)
      k1 = 1
      kX = ipX
      Do iSym = 1,nSym
         kC = k1 + nBas(iSym)*nFro(iSym)
         Call dCopy_(nBas(1)*nOcc(1),CMO(kC),1,Work(kX),1)
         k1 = k1 + nBas(iSym)**2
         kX = kX + nBas(iSym)*nOcc(iSym)
      End Do

C     Get the overlap matrix.
C     -----------------------

      lOAux = nBas(1)*(nBas(1)+1)/2
      lOvlp = nBas(1)*nBas(1)
      Do iSym = 1,nSym
         lOaux = lOaux + nBas(iSym)*(nBas(iSym)+1)/2
         lOvlp = lOvlp + nBas(iSym)*nBas(iSym)
      End Do
      lOaux = lOaux + 4
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

      kTri = ipOaux
      kSqr = ipOvlp
      Do iSym = 1,nSym
         Call Tri2Rec(Work(kTri),Work(kSqr),nBas(iSym),.False.)
         kTri = kTri + nBas(iSym)*(nBas(iSym)+1)/2
         kSqr = kSqr + nBas(iSym)*nBas(iSym)
      End Do
      Call GetMem('AuxOvlp','Free','Real',ipOaux,lOaux)

C     Sort each symmetry block.
C     -------------------------

      kX = ipX
      kC = 1
      kS = ipOvlp
      Do iSym = 1,nSym

C        Cycle loop for empty symmetry blocks.
C        -------------------------------------

         If (nBas(iSym).lt.1 .or. nOcc(iSym).lt.1) Go To 100

C        Allocations.
C        ------------

         lDen = nBas(iSym)*nBas(iSym)
         lU = nOcc(iSym)*nOcc(iSym)
         lScr = nBas(iSym)*nOcc(iSym)
         Call GetMem('SrtDen','Allo','Real',ipDen,lDen)
         Call GetMem('SrtU','Allo','Real',ipU,lU)
         Call GetMem('SrtScr','Allo','Real',ipScr,lScr)

C        Cholesky decompose D=C^TC and thus define the ordering.
C        At this stage, X contains the original MOs (CMO).
C        After the decomposition, X contains the Cholesky MOs.
C        -------------------------------------------------------

         Call GetDens_Localisation(Work(ipDen),Work(kX),nBas(iSym),
     &                             nOcc(iSym))
         irc = -1
         Call ChoLoc(irc,Work(ipDen),Work(kX),ThrCho,xNrm,nBas(iSym),
     &               nOcc(iSym))
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoLoc returned ',irc
            Write(6,*) 'Symmetry block: ',iSym
            Write(6,*) 'Unable to continue...'
            Write(Txt,'(A,I6)') 'ChoLoc return code:',irc
            Call SysAbendMsg(SecNam,
     &                       'Density Cholesky decomposition failed!',
     &                       Txt)
         End If

C        Compute address in CMO, skipping frozen orbitals.
C        -------------------------------------------------

         k1 = kC + nBas(iSym)*nFro(iSym)

C        Compute U = X^TSC.
C        ------------------

         Call GetUmat_Localisation(Work(ipU),Work(kX),Work(kS),CMO(k1),
     &                             Work(ipScr),lScr,nBas(iSym),
     &                             nOcc(iSym))

C        Sort.
C        -----

         Call Sort_Localisation_1(CMO(k1),Work(ipU),nBas(iSym),
     &                            nOcc(iSym))

C        Update counters.
C        ----------------

         kX = kX + nBas(iSym)*nOcc(iSym)
         kC = kC + nBas(iSym)**2
         kS = kS + nBas(iSym)**2

C        De-allocations.
C        ---------------

         Call GetMem('SrtScr','Free','Real',ipScr,lScr)
         Call GetMem('SrtU','Free','Real',ipU,lU)
         Call GetMem('SrtDen','Free','Real',ipDen,lDen)

C        Loop cycling (empty symmetries jump here).
C        ------------------------------------------

  100    Continue

      End Do

C     De-allocations.
C     ---------------

      Call GetMem('XCho','Free','Real',ipX,lX)
      Call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

      End
      SubRoutine Sort_Localisation_1(CMO,U,nBas,nOcc)
C
C     Thomas Bondo Pedersen, November 2005.
C
C     Purpose: sort CMO columns according to U.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(nBas,nOcc), U(nOcc,nOcc)
#include "WrkSpc.fh"

      I1(i)=iWork(ipI1-1+i)
      I2(i)=iWork(ipI2-1+i)

C     Allocations.
C     ------------

      lI1 = nOcc
      lI2 = nOcc
      lC  = nBas*nOcc
      Call GetMem('Sr1I1','Allo','Inte',ipI1,lI1)
      Call GetMem('Sr1I2','Allo','Inte',ipI2,lI2)
      Call GetMem('Sr1C','Allo','Real',ipC,lC)

C     Find max U element in each row.
C     -------------------------------

      ip1 = ipI1 - 1
      Do i = 1,nOcc
         iWork(ip1+i) = i
      End Do

      ip2 = ipI2 - 1
      Do i = 1,nOcc
         jmax = 0
         Umax = -1.0d15
         Do j = 1,nOcc
            If (I1(j) .eq. j) Then
               Utst = abs(U(i,j))
               If (Utst .gt. Umax) Then
                  jmax = j
                  Umax = Utst
               End If
            End If
         End Do
         If (jmax .eq. 0) Then
            Call SysAbendMsg('Sort_Localisation_1','Error:','jmax=0')
         Else
            iWork(ip1+jmax) = 0
            iWork(ip2+i) = jmax
         End If
      End Do

C     Swap MOs according to I2.
C     -------------------------

      Call dCopy_(nBas*nOcc,CMO,1,Work(ipC),1)
      Do i = 1,nOcc
         kOff = ipC + nBas*(I2(i)-1)
         Call dCopy_(nBas,Work(kOff),1,CMO(1,i),1)
      End Do

C     De-allocate.
C     ------------

      Call GetMem('Sr1C','Free','Real',ipC,lC)
      Call GetMem('Sr1I2','Free','Inte',ipI2,lI2)
      Call GetMem('Sr1I1','Free','Inte',ipI1,lI1)

      End
