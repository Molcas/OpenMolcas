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
      SubRoutine PotGrd(Temp,nGrad)
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External PCMGrd1,PCMMmg
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"
      Character Method*8, Label*80
      Real*8 Temp(nGrad)
      Logical DiffOp
*
*-----Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
*...  Prologue
      iRout = 131
      iPrint = nPrint(iRout)
      Call CWTime(TCpu1,TWall1)
      Call qEnter('PotGrd')
*
*---- Allocate memory for density and Fock matrices
*
      nFock = 0
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nFock = nFock + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*...  Get the method label
*     print *,' Read Method label'
      Call Get_cArray('Relax Method',Method,8)
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*     print *,' Read density matrix'
      Call Get_D1ao_Var(ipD_var,Length)
      If ( length.ne.nDens ) Then
         Write (6,*) 'PotGrd: length.ne.nDens'
         Write (6,*) 'length,nDens=',length,nDens
         Call QTrace
         Call Abend()
      End If
*
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=ipD_Var
         Do iIrrep = 0, nIrrep - 1
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',Work(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*
*...  Read the generalized Fock matrix
*...  Fock matrix in AO/SO basis
*     print *,' Read Fock matrix'
      If (.Not.HF_Force) Then
         Call Get_Fock_Occ(ipFock,Length)
         If ( length.ne.nDens ) Then
            Write (6,*) 'PotGrd: length.ne.nDens'
            Write (6,*) 'length,nDens=',length,nDens
            Call QTrace
            Call Abend()
         End If
         If (iPrint.ge.99) then
            Write(6,*) 'generalized Fock matrix'
            ii=ipFock
            Do iIrrep = 0, nIrrep - 1
               Write(6,*) 'symmetry block',iIrrep
               Call TriPrt(' ',' ',Work(ii),nBas(iIrrep))
               ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
            End Do
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     Work(ip1): lOper of each component of the operator
*
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call GetMem('Coor','Allo','Real',ipC,3*nComp)
      Call GetMem('lOper','Allo','Inte',ip1,nComp)
      call dcopy_(nComp*3,[Zero],0,Work(ipC),1)
      iWork(ip1) = 1
      DiffOp = .True.
      Call dZero(Temp,nGrad)
      Call OneEl_g_mck(PCMGrd1,PCMMmG,Temp,nGrad,DiffOp,Work(ipC),
     &             Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,
     &             Label)
      Call PrGrad_mck(' TEST '
     &   //'(PCM) contribution',Temp,nGrad,ChDisp,5)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('lOper','Free','Inte',ip1,nComp)
      Call GetMem('Coor','Free','Real',ipC,3*nComp)
*                                                                      *
************************************************************************
*                                                                      *
*...  Epilogue, end
*
*
      If (.Not.HF_Force) Call GetMem('Fock','Free','Real',ipFock,nFock)
      Call GetMem('D0  ','Free','Real',ipD_Var,nDens)
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      Call qExit('PotGrd')
      Return
      End
