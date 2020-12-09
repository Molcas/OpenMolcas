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
      Subroutine NewCar(kIter,nBVec,nLines,nAtom,nDim,nInter,
     &                  Coor,Lbl,Name,iSym,Smmtrc,
     &                  Degen,mTtAtm,iOptH,User_Def,nStab,
     &                  jStab,Curvilinear,Numerical,DDV_Schlegel,HWRS,
     &                  Analytic_Hessian,iOptC,PrQ,mxdc,iCoSet,rHidden,
     &                  Redundant,MaxItr,iRef,Error)
      use Slapaf_Info, only: Shift, qInt
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8 Coor(3,nAtom), Degen(3*nAtom)
      Integer iSym(3), nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Character Lbl(nInter)*8, Name(nAtom)*(LENIN)
      Logical Smmtrc(3,nAtom), User_Def, Redundant,
     &        Curvilinear, Numerical, DDV_Schlegel, HWRS,
     &        Analytic_Hessian, PrQ, Error
      Real*8, Allocatable:: DFC(:), dss(:), Tmp(:)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(DFC,3*nAtom,Label='DFC')
      Call mma_allocate(dss,nInter,Label='dss')
      Call mma_allocate(Tmp,nInter,Label='Tmp')
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('NewCar: q',' ',qInt,nInter,kIter+1)
      Call RecPrt('NewCar: Shift',' ',Shift,nInter,kIter)
#endif
*
*-----Form the new set of symmetry distinct cartesian coordinates.
*
      rMax = Zero
      ThrR = 0.1d0
      Do i = 1, nInter
        If (Abs(Shift(i,kIter)).gt.rMax) rMax=Abs(Shift(i,kIter))
      End Do
*
      Tmp(:) = qInt(:,kIter)
      dss(:) = Shift(:,kIter)
      Tmp(:) = Tmp(:) + dss(:)
*                                                                      *
************************************************************************
*                                                                      *
      Call Int2Car(dss,Tmp,nInter,Coor,nAtom,nBVec,
     &             nLines,DFC,ndim,Lbl,Name,iSym,Smmtrc,
     &             Degen,kIter,mTtAtm,iOptH,
     &             User_Def,nStab,jStab,Curvilinear,Numerical,
     &             DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &             iCoSet,rHidden,Error,Redundant,MaxItr,
     &             iRef)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Tmp)
      Call mma_deallocate(dss)
      Call mma_deallocate(DFC)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
