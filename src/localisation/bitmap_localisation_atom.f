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
      SubRoutine BitMap_Localisation_Atom(PreFix)
      Implicit Real*8 (a-h,o-z)
      Character*2 PreFix
#include "Molcas.fh"
#include "inflocal.fh"
#include "shinf.fh"
#include "WrkSpc.fh"

      Character*24 SecNam
      Parameter (SecNam = 'BitMap_Localisation_Atom')

      Character*4  Typ(2)
      Character*12 BasNam

      Integer kC0(2)

      Logical Debug

      Debug = .False.

C     Symmetry is not allowed!
C     ------------------------

      If (nSym .ne. 1) Then
         Call SysAbendMsg(SecNam,'Symmetry not implemented!','Sorry!')
      End If

C     Allocate max. sym. block of density matrix
C     and atom based density and CMO matrices.
C     ------------------------------------------

      lDen  = nBas(1)**2
      lDAt  = nAtoms**2
      lCAt  = nAtoms*nOrb2Loc(1)
      lXAt  = lCAt
      Call GetMem('BMpLoc','Allo','Real',ipDen,lDen)
      Call GetMem('DAt','Allo','Real',ipDAt,lDAt)
      Call GetMem('CAt','Allo','Real',ipCAt,lCAt)
      Call GetMem('XAt','Allo','Real',ipXAt,lXAt)

C     Allocate and get index arrays for basis functions per atom.
C     -----------------------------------------------------------

      l_nBas_per_Atom = nAtoms
      l_nBas_Start    = nAtoms
      Call GetMem('nB_per_Atom','Allo','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('nB_Start','Allo','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call BasFun_Atom(iWork(ip_nBas_per_Atom),iWork(ip_nBas_Start),
     &                 Name,nBas(1),nAtoms,Debug)

C     Compute density matrix, Den = CC^T, and set atom based matrices.
C     Generate bitmap and perform sparsity analysis.
C     ----------------------------------------------------------------

      kC1 = ipMOrig + nBas(1)*nFro(1)
      Call GetDens_Localisation(Work(ipDen),Work(kC1),nBas(1),
     &                          nOrb2Loc(1))
      Call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),
     &                        Work(ipDAt),nAtoms,2,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),
     &                        AnaNrm)
      Call GetAt_Localisation(Work(kC1),nBas(1),nOrb2Loc(1),
     &                        Work(ipCAt),nAtoms,1,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),
     &                        AnaNrm)
      kX1 = ipCMO + nBas(1)*nFro(1)
      Call GetAt_Localisation(Work(kX1),nBas(1),nOrb2Loc(1),
     &                        Work(ipXAt),nAtoms,1,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),
     &                        AnaNrm)
      Call GenBMp_Localisation(Work(ipDAt),Work(ipCAt),Work(ipXAt),
     &                         nAtoms,1,'r','r','r',PreFix)
      Call Anasize_Localisation(Work(ipDAt),Work(ipCAt),Work(ipXAt),
     &                          nAtoms,nOrb2Loc(1),1)
      Write(6,*) 'Bitmap files have been generated. Norm: ',AnaNrm

C     Allocate memory for nuclear coordinates.
C     Read nuclear coordinates from the runfile.
C     Generate gnuplot files.
C     ------------------------------------------

      lCoord = 3*nAtoms
      Call GetMem('NucCoord','Allo','Real',ipCoord,lCoord)
      Call Get_dArray('Unique Coordinates',Work(ipCoord),lCoord)

      Call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),
     &                        Work(ipDAt),nAtoms,2,
     &                        iWork(ip_nBas_per_Atom),
     &                        iWork(ip_nBas_Start),
     &                        AnaNrm)
      Write(BasNam,'(A2,A10)') PreFix,'TotDensity'
      Call GenGnu_Localisation(BasNam,Work(ipDAt),Work(ipCoord),nAtoms)
      Typ(1) = 'Dini'
      Typ(2) = 'Dloc'
      kC0(1) = ipMOrig
      kC0(2) = ipCMO
      Do iTyp = 1,2
         kC1 = kC0(iTyp) + nBas(1)*nFro(1)
         Do i = 1,nOrb2Loc(1)
            kC2 = kC1 + nBas(1)*(i-1)
            Call GetDens_Localisation(Work(ipDen),Work(kC2),nBas(1),1)
            Call GetAt_Localisation(Work(ipDen),nBas(1),nBas(1),
     &                              Work(ipDAt),nAtoms,2,
     &                              iWork(ip_nBas_per_Atom),
     &                              iWork(ip_nBas_Start),
     &                              AnaNrm)
            Write(BasNam,'(A2,A4,I6)') PreFix,Typ(iTyp),i
            Call GenGnu_Localisation(BasNam,Work(ipDAt),Work(ipCoord),
     &                               nAtoms)
         End Do
      End Do
      Write(6,*) 'Gnuplot files have been generated. Norm: ',AnaNrm

C     De-allocations.
C     ---------------

      Call GetMem('NucCoord','Free','Real',ipCoord,lCoord)
      Call GetMem('nB_Start','Free','Inte',
     &            ip_nBas_Start,l_nBas_Start)
      Call GetMem('nB_per_Atom','Free','Inte',
     &            ip_nBas_per_Atom,l_nBas_per_Atom)
      Call GetMem('XAt','Free','Real',ipXAt,lXAt)
      Call GetMem('CAt','Free','Real',ipCAt,lCAt)
      Call GetMem('DAt','Free','Real',ipDAt,lDAt)
      Call GetMem('BMpLoc','Free','Real',ipDen,lDen)

      End
