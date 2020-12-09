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
*                                                                      *
************************************************************************
*                                                                      *
*-----Form the new set of symmetry distinct cartesian coordinates.
*
      Call Int2Car(nInter,Coor,nAtom,nBVec,
     &             nLines,ndim,Lbl,Name,iSym,Smmtrc,
     &             Degen,kIter,mTtAtm,iOptH,
     &             User_Def,nStab,jStab,Curvilinear,Numerical,
     &             DDV_Schlegel,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &             iCoSet,rHidden,Error,Redundant,MaxItr,
     &             iRef)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
