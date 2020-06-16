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
      Subroutine PCM_Hss(Hess,nHss)
      use PCM_arrays
      Implicit real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "angstr.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "rctfld.fh"
#include "periodic_table.fh"
      Dimension Hess(nHss)
      Character*2 Elements(MxAtom*8)
*
      iRout = 1
      iPrint = nPrint(iRout)
      Call QEnter('PCM_Hss')
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the PCM cavity contributions to Hessian in solution
*
*
*---- Retrieve atomic info
      Call Get_nAtoms_All(nAtoms)
      Call Allocate_Work(ipCoor,3*nAtoms)
      Call Get_Coord_All(Work(ipCoor),nAtoms)
      Call Get_Name_All(Elements)
      Call GetMem('ANr','Allo','Inte',ipANr,nAtoms)
      Do i = 1, nAtoms
         Do j = 0, Num_Elem
            If (PTab(j).eq.Elements(i)) iWork(ipANr+i-1)=j
         End Do
      End Do
      Call GetMem('Chrg','Allo','Real',ipChrg,nAtoms)
      Call Get_dArray('Nuclear charge',Work(ipChrg),nAtoms)
*
*----- Allocate space for electric field
      nComp=3
      Call GetMem('EF_n','Allo','Real',ip_EF_n,nComp*nTs)
      Call GetMem('EF_e','Allo','Real',ip_EF_e,nComp*nTs)
*
*----- Allocate space for the matrix derivative, total charges,
*      derivatives of charge, two temporary vectors with tesserae
*      dimensions, one vector for the derivative of the potential
*      on tesserae (with dim. nts*nat3),
*      then the PCM contribution to Hessian
      nAt3 = nAtoms * 3
      Call GetMem('DerMat','Allo','Real',ip_DerMat,nTs*nTs)
      Call GetMem('Qtot','Allo','Real',ip_Qtot,nTs)
      Call GetMem('Qder','Allo','Real',ip_Qder,nAt3*nTs)
      Call GetMem('Der1','Allo','Real',ip_Der1,nTs)
      Call GetMem('Der2','Allo','Real',ip_Der2,nTs)
      Call GetMem('VDer','Allo','Real',ip_VDer,nAt3*nTs)
      Call GetMem('HssPCM','Allo','Real',ip_HssPCM,nAt3*nAt3)

      Call PCM_Cav_Hss(Angstr,nAtoms,nAt3,nTs,nS,Eps,iWork(ipANr),
     &     Work(ipCoor),Work(ipChrg),Work(ip_EF_n),Work(ip_EF_e),
     &     PCMSph,iWork(ip_ISph),iWork(ip_N),PCMTess,
     &     PCM_SQ,Work(ip_Qtot),PCMDM,Work(ip_HssPCM),
     &     Work(ip_DerMat),dTes,dPnt,dRad,
     &     dCntr,Work(ip_QDer),Work(ip_Der1),Work(ip_Der2),
     &     Work(ip_VDer))
      Call GetMem('HssPCM','Free','Real',ip_HssPCM,nHss)
      Call GetMem('DerMat','Free','Real',ip_DerMat,nTs*nTs)
      Call GetMem('Qtot','Free','Real',ip_Qtot,nTs)
      Call GetMem('Qder','Free','Real',ip_Qder,nAt3*nTs)
      Call GetMem('Der1','Free','Real',ip_Der1,nTs)
      Call GetMem('Der2','Free','Real',ip_Der2,nTs)
      Call GetMem('VDer','Free','Real',ip_VDer,nAt3*nTs)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('PCM_Hss')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Hess)
      End
