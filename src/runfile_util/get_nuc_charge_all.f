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
* Copyright (C) Luca De Vico                                           *
************************************************************************
      Subroutine Get_Nuc_Charge_All(Charges_All,nAtoms_All)
************************************************************
*
*   <DOC>
*     <Name>Get\_Nuc\_Charge\_All</Name>
*     <Syntax>Call Get\_Nuc\_Charge\_All(Charges\_All,nAtoms\_All)</Syntax>
*     <Arguments>
*       \Argument{Charges\_All}{Array of charges}{Real*8 (nAtoms\_All)}{out}
*       \Argument{nAtoms\_All}{Number of atoms}{Integer}{in}
*     </Arguments>
*     <Purpose>Get Nuclear charges from RUNFILE</Purpose>
*     <Dependencies></Dependencies>
*     <Author>L. De Vico</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Place nuclear charges (in a.u.) into array Charges\_All(*)
*       Based on Get\_Coord\_All
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Charges_All(nAtoms_All)
*
      Call Get_nAtoms_All(nAtoms_Allx)
      If (nAtoms_All.ne.nAtoms_Allx) Then
         Write (6,*) 'Get_Nuc_Charge_All: nAtoms_All.ne.nAtoms_Allx'
         Write (6,*) 'nAtoms_All=',nAtoms_All
         Write (6,*) 'nAtoms_Allx=',nAtoms_Allx
         Call QTrace
         Call Abend
      End If

      Call Get_iScalar('Unique atoms',nAtoms)

      Call Allocate_Work(ipCU,3*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ipCU),3*nAtoms)

      Call Allocate_Work(ipCMu, nAtoms)
      Call Get_dArray('Nuclear charge',Work(ipCMu),nAtoms)

      Call Get_Nuc_Charge_All_(Work(ipCU),Work(ipCMu),nAtoms,
     &                    Charges_All,nAtoms_All)

      Call Free_Work(ipCMu)
      Call Free_Work(ipCU)
*
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Get_Nuc_Charge_All_(Coord_Unique,Charges_Unique,
     &                       nUnique_Atoms,Charges_All,nAll_Atoms)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer iGen(3), iCoSet(0:7,0:7), iStab(0:7)
      Real*8  Coord_Unique(3,nUnique_Atoms),
     &        Charges_Unique(nUnique_Atoms),
     &        Charges_All(nAll_Atoms)
      integer is_nSym, nSym
      integer is_iOper, iOper(0:7)
      save is_nSym, is_iOper
      data is_nSym/0/, is_iOper/0/
      save nSym, iOper
*     Write (*,*) 'Enter Get_Nuc_Charge_All_'
*                                                                      *
************************************************************************
*                                                                      *
      if(is_nSym.eq.0) then
       Call Get_iScalar('nSym',nSym)
       is_nSym=1
      endif
      nIrrep=nSym
*     Write (*,*) 'Get_Nuc_Charge_All_: nSym=',nSym
      if(is_iOper.eq.0) then
       Call Get_iArray('Symmetry operations',iOper,nSym)
       is_iOper=1
      endif
*     Write (*,*) 'Get_Nuc_Charge_All_: iOper=',(iOper(i),i=0,nSym-1)
*                                                                      *
************************************************************************
*                                                                      *
      nGen=0
      If (nSym.eq.2) nGen=1
      If (nSym.eq.4) nGen=2
      If (nSym.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.eq.3) iGen(3)=iOper(4)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of all nuclear charges
*
      iAll_Atom=0
      MaxDCR=0
      Do iUnique_Atom = 1, nUnique_Atoms
         iChAtom=iChxyz(Coord_Unique(1,iUnique_Atom),iGen,nGen)
         Call Stblz(iChAtom,iOper,nIrrep,nStab,iStab,MaxDCR,iCoSet)
         nCoSet=nIrrep/nStab

         Charge_Old = Charges_Unique(iUnique_Atom)

         Do iCo = 0, nCoSet-1
            iAll_Atom = iAll_Atom + 1
            Charges_All(iAll_Atom) = Charge_Old
         End Do
*
      End Do
*
*     Call RecPrt('Charges_Unique',' ',Charges_Unique,1,nUnique_Atoms)
*     Call RecPrt('Charges_All',' ',Charges_All,1,nAll_Atoms)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
