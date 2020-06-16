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
      SubRoutine PCM_Cav_grd(Grad,nGrad)
      use PCM_arrays
      Implicit Real*8 (A-H,O-Z)
      Real*8 Grad(nGrad)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
*
      iRout = 1
      iPrint = nPrint(iRout)
      Call QEnter('PCM_Cav_grd')
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the geometric contributions to
*     derivatives in solution
*
      Call GetMem('DerDM','Allo','Real',ip_DerDM,nTs*nTs)
      Call Get_nAtoms_All(MaxAto)
      Call GetMem('PCMGrd','Allo','Real',ip_PCMGrd,3*MaxAto)
      LcNAtm = ISlPar(42)
      Call GeoDer(LcNAtm,Conductor,nTs,nS,Eps,Work(ip_Sph),
     &            iWork(ip_ISph),iWork(ip_N),Work(ip_Tess),
     &            PCM_SQ,Work(ip_DerDM),Work(ip_PCMGrd),
     &            dTes,dPnt,dRad,dCntr)
*     Call RecPrt('PCM_Cav_Grd','(5G20.10)',Work(ip_PCMGrd),3,MaxAto)
      Call GrdTr_Alaska(Work(ip_PCMGrd),MaxAto,Grad,nGrad)
      Call GetMem('PCMGrd','Free','Real',ip_PCMGrd,3*MaxAto)
      Call GetMem('DerDM','Free','Real',ip_DerDM,nTs*nTs)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('PCM_Cav_grd')
      Return
      End
