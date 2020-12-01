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
     &                          Lbl,Name,SaveBMx,RefIter,Error)
      use Slapaf_Info, only: Cx
      Implicit None
#include "info_slapaf.fh"
#include "db.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Integer :: kIter,nLines,nAtom,nDim,nInter,RefIter
      Real*8 :: BMx(3*nAtom,nInter)
      Character :: Lbl(nInter)*8,Name(nAtom)*(LENIN)

      Real*8, Allocatable :: DFC(:),dss(:),qTemp(:), Coor(:)
      Integer :: ipBMx
      Logical :: Numerical,PrQ,Error,SaveBMx
*
      Call mma_allocate(Coor,3*nAtom,Label='Coor')
      Call GetMem('BMx','ALLO','REAL',ipBMx,(3*nAtom)*nInter)
      Call dCopy_(3*nAtom,Cx(:,:,kIter),1,Coor,1)
      Call dCopy_((3*nAtom)*nInter,BMx,1,Work(ipBMx),1)
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
     &            ipBMx,Lbl,DFC,dss,
     &            qTemp,Name,iSym,Smmtrc,Degen,
     &            mTtAtm,iOptH,
     &            User_Def,nStab,jStab,Curvilinear,Numerical,
     &            DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &            iCoSet,rHidden,ipRef,Redundant,MaxItr,
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
