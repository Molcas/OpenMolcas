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
      Subroutine LoProp_Print(rMP,nij,nElem,nAtoms,Q_Nuc,LblCnt,lSave)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 rMP(nij,nElem), Q_Nuc(nAtoms)
      Character*(LENIN4) LblCnt(nAtoms)
      Character*120 Banner_Line(3)
*
      Real*8 E_Charge(MxAtom), Q_Charge(MxAtom)
      Character*(LENIN) Lbl(MxAtom)
      Logical lSave, Reduce_Prt
      External Reduce_Prt
*                                                                      *
************************************************************************
*                                                                      *
*     Get the print level
*
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
      If (iPL.lt.2) Return
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Binom()
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over all domains
*
      Write (6,*)
      Banner_Line(1)='LoProp Charges per center'
#ifdef _DEBUG_
      Banner_Line(2)=' '
      Banner_Line(3)='Note that this charge analysis only makes '
     &             //'sense if the orbital basis is of true AO '
     &             //'type!'
      Call Banner (Banner_Line,3,120)
#else
      Write(6,'(6X,A)') Banner_Line(1)(:mylen(Banner_Line(1)))
#endif
*
*---- Collect data
*
      mAtoms=0
      ij = 0
      Do iAtom = 1, nAtoms
         ij = iAtom*(iAtom+1)/2
         If (LblCnt(iAtom)(LENIN1:LENIN4).eq.':E  ' .or.
     &       LblCnt(iAtom)(LENIN1:LENIN4).eq.'    ') Then
            mAtoms=mAtoms+1
            Q_Charge(mAtoms)=Q_nuc(iAtom)
            E_Charge(mAtoms)=rMP(ij,1)
            Lbl(mAtoms)=LblCnt(iAtom)(1:LENIN)
         End If
      End Do
      If (lSave) Then
         Call GetMem('LoProp Chg','Allo','Real',ipLPChg,mAtoms)
         Call dCopy_(mAtoms,Q_Charge,1,Work(ipLPChg),1)
         call daxpy_(mAtoms,One,E_Charge,1,Work(ipLPChg),1)
         Call Put_dArray('LoProp Charge',Work(ipLPChg),mAtoms)
         Call GetMem('LoProp Chg','Free','Real',ipLPChg,mAtoms)
      End If
*
*---- Print out the stuff!
*
      Inc=10
      Do iSt = 1, mAtoms, Inc
         iEnd=Min(mAtoms,iSt+Inc-1)
         Write(6,*)
         Write(6,'(/16X,10(3X,A))') (Lbl(i),i=iSt,iEnd)
         Write(6,'(6X,A,10F9.4)')'Nuclear   ',(Q_Charge(i),i=iSt,iEnd)
         Write(6,'(6X,A,10F9.4)')'Electronic',(E_Charge(i),i=iSt,iEnd)
         Write(6,*)
         Write(6,'(6X,A,10F9.4)')'Total     ',
     &        (Q_Charge(i)+E_Charge(i),i=iSt,iEnd)
      End Do
*
#ifdef _DEBUG_
      Write (6,*)
      Call Banner (Banner_Line,0,120)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
