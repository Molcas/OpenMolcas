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
      SubRoutine Drvespf(Grad,Temp,nGrad,CCoor)
*
*     Driver to compute the ESPF B*dV contribution to the gradient
*     This is a hack of the alaska/drvh1 subroutine with a little
*     piece of (extinct) integral_util/drvprop subroutine
*
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "espf.fh"
*
      External BdVGrd
      External NAMmG
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
#include "wldata.fh"
#include "iavec.fh"
      Character Label*80
      Real*8 Grad(nGrad), Temp(nGrad)
      Real*8 Ccoor(*)
*     Logical DiffOp, DoRys
      Logical DiffOp
*
*-----Statement function
*
      nElem(i) = (i+1)*(i+2)/2
*
*...  Prologue
      iPrint = 1
      Call qEnter('Drvespf')
*
*     Set up the angular index vector
*
      i = 0
      Do 1000 iR = 0, iTabMx
         Do 2000 ix = iR, 0, -1
            Do 3000 iy = iR-ix, 0, -1
               iz = iR-ix-iy
               i = i + 1
               ixyz(1,i) = ix
               ixyz(2,i) = iy
               ixyz(3,i) = iz
 3000       Continue
 2000    Continue
 1000 Continue

*
*---- Allocate memory for density matrix
*
      nDens = 0
      Do iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
      End Do
*
*...  Read the variational 1st order density matrix
*...  density matrix in AO/SO basis
*
      Call Get_D1ao_Var(ipD_var,Length)
      If ( length.ne.nDens ) Then
         Write (6,*) 'Drvespf: length.ne.nDens'
         Write (6,*) 'length,nDens=',length,nDens
         Call Abend
      End If
      If (iPrint.ge.99) then
         Write(6,*) 'variational 1st order density matrix'
         ii=ipD_Var
         Do iIrrep = 0, nIrrep - 1
            Write(6,*) 'symmetry block',iIrrep
            Call TriPrt(' ',' ',Work(ii),nBas(iIrrep))
            ii = ii + nBas(iIrrep)*(nBas(iIrrep)+1)/2
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     nOrdOp: order/rank of the operator
*     Work(ip1): lOper of each component of the operator
*
      nPrint(112) = 5
      iPL = iPL_espf()
      If (iPL.ge.3) nPrint(112) = 15
      nOrdOp=0
      nComp = nElem(nOrdOp)
      Call GetMem('lOper','Allo','Inte',ip1,nComp)
      iWork(ip1) = 1
*
************************************************************************
*                                                                      *
*     Trace the "variational" first order density matrix with the      *
*     gradient of the external field integrals.                        *
*                                                                      *
************************************************************************
      DiffOp = .True.
      Label = ' The ESPF BdV contribution'
      Call OneEl_g(BdVGrd,NAMmG,Temp,nGrad,DiffOp,CCoor,
     &           Work(ipD_Var),nDens,iWork(ip1),nComp,nOrdOp,Label)
      Call DaXpY_(nGrad,One,Temp,1,Grad,1)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('lOper','Free','Inte',ip1,nComp)
      Call GetMem('D0  ','Free','Real',ipD_Var,nDens)
*
      Call qExit('Drvespf')
      Return
      End
