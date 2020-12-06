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
      Integer :: kIter,nLines,nAtom,nDim,nInter,RefIter
      Real*8 :: BMx(3*nAtom,nInter)
      Integer, External:: ip_of_Work
      Character :: Lbl(nInter)*8,Name(nAtom)*(LENIN)

      Real*8, Allocatable :: Coor(:,:), BMx_Tmp(:,:)
      Integer :: ipBMx
      Logical :: Numerical,PrQ,Error,SaveBMx
*
      Call mma_allocate(Coor,3,nAtom,Label='Coor')
      Coor(:,:) = Cx(:,:,kIter)

      Call mma_allocate(BMx_tmp,(3*nAtom),nInter,Label='BMx_tmp')
      BMx_tmp(:,:) = BMx(:,:)
      ipBMx = ip_of_Work(BMx(1,1))
*
      Numerical=.False.
      PrQ=.False.
*
      Force_dB=SaveBMx
*
      Call NewCar(kIter,nBVec,nLines,nAtom,nDim,nInter,Coor,
     &            ipBMx,Lbl,Name,iSym,Smmtrc,Degen,
     &            mTtAtm,iOptH,
     &            User_Def,nStab,jStab,Curvilinear,Numerical,
     &            DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &            iCoSet,rHidden,Redundant,MaxItr,
     &            RefIter,Error)
*
      Force_dB=.False.
      Call mma_deallocate(Coor)

      If (.NOT.SaveBMx) BMx(:,:) = BMx_tmp(:,:)

      Call mma_deallocate(BMx_tmp)
*
      End Subroutine NewCar_Kriging
