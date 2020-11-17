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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine NewCar_Kriging(kIter,nLines,nAtom,nDim,nInter,BMx,
     &                          dMass,Lbl,kShift,qInt,dqInt,Name,Cx,
     &                          SaveBMx,RefIter,Error)
      use Slapaf_Info, only: Gx
      Implicit None
#include "info_slapaf.fh"
#include "db.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Integer :: kIter,nLines,nAtom,nDim,nInter,RefIter
      Real*8 :: BMx(3*nAtom,nInter),dMass(nAtom),kShift(nInter,kIter),
     &          qInt(nInter,MaxItr),dqInt(nInter,MaxItr),
     &          Cx(3*nAtom,kIter+1)
      Character :: Lbl(nInter)*8,Name(nAtom)*(LENIN)

      Real*8, Allocatable :: DFC(:),dss(:),qTemp(:), Coor(:)
      Integer :: ipBMx,ip_qInt,ip_dqInt,ipShift
      Logical :: Numerical,PrQ,Error,SaveBMx
      Integer, External :: ip_of_Work
*
      Call mma_allocate(Coor,3*nAtom,Label='Coor')
      Call GetMem('BMx','ALLO','REAL',ipBMx,(3*nAtom)*nInter)
      Call dCopy_(3*nAtom,Cx(1,kIter),1,Coor,1)
      Call dCopy_((3*nAtom)*nInter,BMx,1,Work(ipBMx),1)
*
      ip_qInt=ip_of_Work(qInt(1,1))
      ip_dqInt=ip_of_Work(dqInt(1,1))
      ipShift=ip_of_Work(kShift(1,1))
*
      Numerical=.False.
      PrQ=.False.
*
      Call mma_allocate(DFC,3*nAtom,Label='DFC')
      Call mma_allocate(dss,nInter,Label='dss')
      Call mma_allocate(qTemp,nInter,Label='qTemp')
      Force_dB=SaveBMx
*
      Call NewCar(kIter,nBVec,nLines,nAtom,nDim,nInter,Coor,
     &            ipBMx,dMass,Lbl,kShift,ip_qInt,ip_dqInt,DFC,dss,
     &            qTemp,Name,iSym,Smmtrc,Degen,
     &            Gx,Cx,mTtAtm,iWork(ipANr),iOptH,
     &            User_Def,nStab,jStab,Curvilinear,Numerical,
     &            DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &            iCoSet,rHidden,ipRef,Redundant,nqInt,MaxItr,
     &            RefIter,Error)
*
      Force_dB=.False.
      Call mma_deallocate(DFC)
      Call mma_deallocate(dss)
      Call mma_deallocate(qTemp)
      If (SaveBMx) Call dCopy_((3*nAtom)*nInter,Work(ipBMx),1,BMx,1)
      Call mma_deallocate(Coor)
      Call GetMem('BMx','FREE','REAL',ipBMx,(3*nAtom)**2)
*
      End Subroutine NewCar_Kriging
