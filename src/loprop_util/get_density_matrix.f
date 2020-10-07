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
      Subroutine Get_Density_Matrix(ip_D,nBas1,nBas2,nBasMax,nBas,nSym,
     &                      ipP,UserDen,PrintDen,SubtractDen,SubScale,
     &                      Q_Nuc,nAtoms,iPert,Restart,Utility,TDensity,
     &                      nStateI,nStateF)

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
* OBSERVE! MxState has to be the same as MxStat in cntrl.fh
*          in src/rassi.
      Parameter (MxState=200,MxStOT=MxState*(MxState+1)/2)
      Integer nBas(nSym)
      Character*16 Label, filename
      Dimension Q_Nuc(nAtoms)
      Dimension iToc(MxStOT)
      Logical UserDen, PrintDen, SubtractDen, Exist, Restart, Utility
      Logical TDensity, ok1, ok2, Found
      Character*8 Method
      Real*8, Allocatable:: DTmp(:), DSym(:)
*
*define _DEBUG_
*
      Write(Label,'(A,I1)') 'LoProp Dens ',iPert
      If (Restart) Then
*
* Retrieve density matrix from the runfile
*
         Call qpg_dArray(Label,Exist,nDens)
         If (.NOT. Exist .or. nDens .eq. 0) Then
            Call SysAbendMsg('get_density_matrix',
     &                       'Could not locate:',Label)
         End If
         Call Allocate_Work(ip_D,nDens)
         Call Get_dArray(Label,Work(ip_D),nDens)
*
         nSize=nDens
*
      Else If (nSym.eq.1) Then
*
         nSize=nBas(1)*(nBas(1)+1)/2
*
         If(UserDen) then
* Additions made by A.Ohrn to enable LoProp to read in a density matrix
* provided by the user. Generalized by P. Soderhjelm so that also
* perturbed density matrices may be read for polarizability calculation
           Lu_Ud=56
           Lu_Ud=IsFreeUnit(Lu_Ud)
           If(ipert.eq.0) Then
              filename='USERDEN'
           Else
              Write(filename,'(A7,I1)')'USERDEN',ipert
           EndIf
           Call OpnFl(filename,Lu_Ud,Exist)
           If(.not.Exist) then
             Write(6,*)
             Write(6,*)' Unable to locate user density matrix.'
             Call Abend()
           Endif
           Call GetMem('UserDen','Allo','Real',ipUser,nSize)
           Read(Lu_Ud,*)(Work(ipUser+k),k=0,nSize-1)
           Call Put_D1ao(Work(ipUser),nSize)
           Call GetMem('UserDen','Free','Real',ipUser,nSize)
           Close(Lu_Ud)
         Endif
* Addition of A.Ohrn, to collect transition density matrix from ToFile.
         If(TDensity) then
           LuIn=57
           LuIn=IsFreeUnit(LuIn)
           Call DaName(LuIn,'TOFILE')
           iDisk=0
*---  Table-of-contents
           Call iDaFile(LuIn,2,iToc,MxStOT,iDisk)
*---  Allocation of density
           Call GetMem('TDMden','Allo','Real',iTDMden,nSize)
*---  Loop to 'suck-out' the relevant matrix from the ToFile.
           nStateM=max(nStateI,nStateF)
           kaunter=0
           Do iS1=1,nStateM
             Do iS2=1,iS1
               kaunter=kaunter+1
               iDisk=iToc(kaunter)
               Call dDaFile(LuIn,2,Work(iTDMden),nSize,iDisk)
               ok1=iS1.eq.nStateI.and.iS2.eq.nStateF
               ok2=iS1.eq.nStateF.and.iS2.eq.nStateI
               If(ok1.or.ok2) then
                 Call Put_D1ao(Work(iTDMden),nSize)
               Endif
             Enddo
           Enddo
           Call GetMem('TDMden','Free','Real',iTDMden,nSize)
           Call DaClos(LuIn)
         Endif
* End addition A.Ohrn.

         If(SubtractDen) then
* Additions made by P.Soderhjelm to enable LoProp to read in a density matrix
* provided by the user and subtract that from the current one (which could
* in fact be provided by the Userdens keyword. After subtraction the matrix
* is scaled by a constant SubScale. Nuclear charges are set to zero.
           Lu_Ud=56
           Lu_Ud=IsFreeUnit(Lu_Ud)
           Call OpnFl('SUBDEN',Lu_Ud,Exist)
           If(.not.Exist) then
             Write(6,*)
             Write(6,*)' Unable to locate density matrix to subtract.'
             Call Abend()
           Endif
           Call GetMem('UserDen','Allo','Real',ipUser,nSize)
           Read(Lu_Ud,*)(Work(ipUser+k),k=0,nSize-1)
           Call mma_allocate(DTmp,nSize)
           Call Get_D1ao(Dtmp,nSize)
           Do k=0,nSize-1
              Dtmp(1+k)=SubScale*(Dtmp(1+k)-Work(ipUser+k))
           EndDo
           Call Put_D1ao(Dtmp,nSize)
           Call mma_deallocate(DTmp)
           Call GetMem('UserDen','Free','Real',ipUser,nSize)
           Close(Lu_Ud)
           Do k=1,nAtoms
              Q_Nuc(k)=0.0D0
           EndDo
         Endif
* End addition P.Soderhjelm.

* Addition by J.Bostrom to check if loprop should pick
* the variational density instead.
         Call Get_cArray('Relax Method',Method,8)
         iMp2Prpt = 0
         If(Method.eq.'MBPT2   ') Then
            Call Get_iScalar('mp2prpt',iMp2Prpt)
         End If
         Call getmem('D','Allo','Real',ip_D,nSize)
         If(iMp2Prpt.ne.0) Then
            Call Get_D1ao_var(Work(ip_D),nSize)
         else
* End Addition J.Bostrom (well, the End If too ofc.)
            Call Get_D1ao(Work(ip_D),nSize)
         End If
         nDens=nSize
#ifdef _DEBUG_
         Call RecPrt('D',' ',Work(ip_D),1,nDens)
#endif
         If(PrintDen) Then
* Addition made by P.Soderhjelm to enable LoProp to print the density matrix
* to a text file for later use.
           Lu_Ud=56
           Lu_Ud=IsFreeUnit(Lu_Ud)
           Call OpnFl('PRDEN',Lu_Ud,Exist)
           Write(Lu_Ud,'(10d25.16)')(Work(ip_D+k),k=0,nDens-1)
           Close(Lu_Ud)
         EndIf

*
      Else
*
         Call Allocate_Work(ip_D,nBas1*(nBas1+1)/2)
         Call Allocate_Work(ip_D_sq,nBas1**2)
         Call Allocate_Work(ip_Tmp,nBas2)
*
         Call Qpg_darray('D1ao',Found,nDens)
         If (Found .and. nDens/=0) Then
            Call mma_allocate(DSym,nDens,Label='DSym')
            Call Get_D1ao(DSym,nDens)
         Else
            Write (6,*) 'Get_density_matrix: not found.'
            Call Abend()
         End If
         iSyLbl=1
#ifdef _DEBUG_
         Call RecPrt('DSym',' ',DSym,1,nDens)
#endif
         iOfft = 1
         iOffs = ip_Tmp
         Do iSym = 1, nSym
            If (nBas(iSym).eq.0) Go To 99
#ifdef _DEBUG_
            Call TriPrt('DSym',' ',DSym(iOfft),nBas(iSym))
#endif
            Call Square(DSym(iOfft),Work(iOffs),1,nBas(iSym),nBas(iSym))
            Call DScal_(nBas(iSym)**2,Half,Work(iOffs),1)
            Call DScal_(nBas(iSym)   ,Two, Work(iOffs),nBas(iSym)+1)
#ifdef _DEBUG_
            Call RecPrt('DSym',' ',Work(iOffs),nBas(iSym),nBas(iSym))
#endif
            iOfft = iOfft + nBas(iSym)*(nBas(iSym)+1)/2
            iOffs = iOffs + nBas(iSym)**2
 99         Continue
         End Do
         Call mma_deallocate(DSym)
*
         nScr=nBasMax*nBas1
         Call Allocate_Work(ipScr,nScr)
         Call Desymmetrize(Work(ip_Tmp),nBas2,Work(ipScr),nScr,
     &                     Work(ip_D_sq),nBas,nBas1,Work(ipP),nSym,
     &                     iSyLbl)
         Call Free_Work(ipScr)
         Call Free_Work(ip_Tmp)
*
         Call Triangularize(Work(ip_D_sq),Work(ip_D),nBas1,.True.)
         Call Free_Work(ip_D_sq)
      End If
#ifdef _DEBUG_
      Call TriPrt('Density Matrix',' ',Work(ip_D),nBas1)
#endif
*
* Copy Density Matrix to the runfile
*
      If ((.Not. Restart) .AND. (.Not. Utility)) Then
*        Call put_dArray(Label,Work(ip_D),nDens)
         Call put_dArray(Label,Work(ip_D),nBas1*(nBas1+1)/2)
      End If
*
      Return
      End
      Subroutine Triangularize(A_sq, A_tr, n, Fold)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 A_sq(n,n), A_tr(n*(n+1)/2)
      Logical Fold
*
      Fact = One
      If (Fold) Fact = Two
      ij = 1
      Do i = 1, n
         Do j = 1, i-1
            A_tr(ij)=Fact*A_sq(i,j)
            ij = ij + 1
         End Do
         A_tr(ij)=A_sq(i,j)
         ij = ij + 1
      End Do
*
      Return
      End
*
