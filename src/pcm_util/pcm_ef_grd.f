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
      SubRoutine PCM_EF_grd(Grad,nGrad)
      use Basis_Info
      use Center_Info
      use PCM_arrays
      Implicit Real*8 (A-H,O-Z)
      Real*8 Grad(nGrad)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Logical Save_tmp
      Real*8 EF_Temp(3)
      Real*8, Allocatable :: Cord(:,:), Chrg(:), FactOp(:)
      Integer, Allocatable :: lOper(:)
*
      iRout = 1
      iPrint = nPrint(iRout)
      Call QEnter('PCM_EF_grd')
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Save_tmp=PrPrt
      PrPrt=.True.
      nOrdOp=1
      nComp=(nOrdOp+1)*(nOrdOp+2)/2
      Call GetMem('EF','Allo','Real',ip_EF,nComp*2*nTs)
      ip_EF_nuclear   =ip_EF
      ip_EF_electronic=ip_EF + nComp
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_nAtoms_All(MaxAto)
*
      Call mma_allocate(Cord,3,MaxAto)
      Call mma_allocate(Chrg,MaxAto)
*
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         If (dbsc(jCnttp)%Aux) Cycle
         Z = dbsc(jCnttp)%Charge
         mCnt = dbsc(jCnttp)%nCntr
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            Do i = 0, nIrrep/dc(ndc)%nStab - 1
               Call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),
     &                 Cord(1:3,nc))
               Chrg(nc) = Z
               nc = nc + 1
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute the the electric field on the tiles
*
*     1) The nuclear contribution
*
      ip_EF_n=ip_EF_nuclear
      ip_EF_e=ip_EF_electronic
      Do iTile = 1, nTs
         Call EFNuc(PCMTess(1,iTile),Chrg,Cord,MaxAto,EF_temp,nOrdOp)
         Work(ip_EF_n  )=EF_Temp(1)
         Work(ip_EF_n+1)=EF_Temp(2)
         Work(ip_EF_n+2)=EF_Temp(3)
         Work(ip_EF_e  )=Zero
         Work(ip_EF_e+1)=Zero
         Work(ip_EF_e+2)=Zero
         ip_EF_n = ip_EF_n + 2*nComp
         ip_EF_e = ip_EF_e + 2*nComp
      End Do
*
      Call mma_deallocate(Cord)
      Call mma_deallocate(Chrg)
*
*     2) The electronic contribution
*
*
*     Get the total 1st order AO density matrix
*
      Call Get_D1ao(ipD1ao,nDens)
*
      Call mma_allocate(FactOp,nTs)
      Call mma_allocate(lOper,nTs)
      Call DCopy_(nTs,[One],0,FactOp,1)
      Call ICopy(nTs,[255],0,lOper,1)
*
      Call Drv1_PCM(FactOp,nTs,Work(ipD1ao),nDens,
     &              PCMTess,lOper,Work(ip_EF),nOrdOp)
*
      Call mma_deallocate(lOper)
      Call mma_deallocate(FactOp)
      Call GetMem('D1ao','Free','Real',ipD1ao,nDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Now form the correct combinations
*
      Call Cmbn_EF_DPnt(Work(ip_EF),nTs,dPnt,MaxAto,
     &                  dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)
*
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('EF','Free','Real',ip_EF,nComp*2*nTs)
      PrPrt=Save_tmp
      Call Free_iSD()
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('PCM_EF_grd')
      Return
      End
