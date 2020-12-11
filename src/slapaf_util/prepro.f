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
      Subroutine PrePro(nLines,iInt,nFix,nAtom,nInter,Coor)
      use Slapaf_Info, only: Grd, Atom, nSup
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom)
      Logical CofM
      Real*8, Allocatable:: TR(:)
*
      iRout=134
      iPrint=nPrint(iRout)
*
      CofM = Iter.eq.1 .and. lNmHss
      Call mma_allocate(TR,18*nAtom,Label='TR')
      TR(:)=Zero
      Call TRPGen(nDimBC,nAtom,Coor,mTR,CofM,TR)
      Call mma_deallocate(TR)
      If (lNmHss) Then
         If (Iter.eq.1) mTROld=mTR
         If (iter.le.2*(nDimBC-mTROld)+1.and.iter.ne.1)
     &    mTR=mTROld
      Else
         mTROld=mTR
      End If
*
*-----Operate according to two modes
*     nLines.gt.0 : user supplied internal coordinates
*     nLines.le.0 : Cartesian or Internal Coordinates
*
      nRowH = 0
      nInter = nDimBC - mTR
      If (nLines.gt.0) Then
*
*--------Find the number of active and frozen internal coordinates.
*
         Call Rd_UDIC(nLines,iInt,nFix,nRowH)
         nQQ=iInt+nFix
         If (nRowH.GT.0) then
            lRowH=.True.
            Call Rd_UDIC_RowH(nQQ,nRowH,mRowH)
         EndIf
         If (nQQ.gt.nInter) Redundant=.True.
*
      Else
*
         nFix=0
      End If
*
*-----Initiate the force constant matrix in internal coordinate
*     basis, excluding rotation and translation.
*     Write to runfile only on the first iteration and that there
*     was not an already defined Hessian.
*
      If (iter.eq.1) Call IntFcm(lOld,lOld_Implicit)
      If (.Not.lOld.and.lOld_Implicit) lOld=.True.
*
*-----Symmetrize forces
*
      If (LSup) Then
         Call SupSym(Grd,nAtom,cMass,Coor,nSupSy,nSup,Atom)
         Call mma_deallocate(Atom)
         Call mma_deallocate(nSup)
      End If
*
      Return
      End
