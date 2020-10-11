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
* Copyright (C) 1992, Martin Schuetz                                   *
*               2017, Roland Lindh                                     *
************************************************************************
      SubRoutine vOO2OV(v1,n1,v2,n2,nD)
************************************************************************
*                                                                      *
*     purpose: converts vector of dim nOO (e.g. gradient) to vector    *
*              of lower dimension nOV (occ-virt block only)            *
*              if n1==nOO && n2==nOV  -> compress                      *
*              else if if n1==nOV && n2==nOO -> decompress             *
*                                                                      *
*     output:                                                          *
*       v2      : compressed or decompressed vector, dep. on input     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*     Roland Lindh                                                     *
*     Uppsala University, Sweden, 2017                                 *
*     Fortran pointer structure                                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
#include "compiler_features.h"
#ifndef POINTER_REMAP
      Use, Intrinsic :: ISO_C_BINDING
#endif
      Implicit Real*8 (a-h,o-z)
*
*     declaration subroutine parameters
      Integer n1,n2
      Real*8, Target:: v1(n1,nD),v2(n2,nD)
      Real*8, Dimension(:,:), Pointer:: pv1, pv2
*
*     declaration local variables
      Integer iSym,ii,ia,ioffs,ivoffs
*
#include "real.fh"

#include "mxdm.fh"
#include "infscf.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Do iD = 1, nD
         iCase = iD - 1
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
         Write (6,*) 'n1,n2,nOO,nOV=',n1,n2,nOO,nOV
         iOff = 1
         Do iSym = 1, nSym
            If (n1.eq.nOO) Then
               nO = nOrb(iSym)
               Call RecPrt('vOO2OV: v1',' ',v1(ioff,iD),nO,nO)
               iOff = iOff + nO**2
            Else If (n1.eq.nOV) Then
               nO = Max(nOcc(iSym,1),nOcc(iSym,2))-nFro(iSym)
               nV = nOrb(iSym)-Min(nOcc(iSym,1),nOcc(iSym,2))-nFro(iSym)
               Call RecPrt('vOO2OV: v1',' ',v1(ioff,iD),nO,nV)
               iOff = iOff + nO*nV
            End If
         End Do
#endif
         ioffs=0
         ivoffs=1
         Do iSym=1,nSym
            nii1=nFro(iSym)+1
            nii2=nOcc(iSym,iD)
            nia1=nOcc(iSym,iD)+1
            nia2=nOrb(iSym)
            nv1=nia2**2
            nv2=(nia2-nia1+1)*(nii2-nii1+1)
*
            If ((n1.eq.nOO).and.(n2.eq.nOV)) Then
*
*           compress
*
#ifdef POINTER_REMAP
               pv1(1:nia2,1:nia2) => v1(ioffs+1:ioffs+nv1,iD)
               pv2(nia1:nia2,nii1:nii2) => v2(ivoffs:ivoffs+nv2-1,iD)
*
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
                     If (pv1(ia,ii).ne.-pv1(ii,ia)) Then
                        Write (6,*) 'inconsistency in gradient'
                        Call Abend()
                     End If
                     pv2(ia,ii)=pv1(ia,ii)
                  End Do
               End Do
#else
               Call C_F_POINTER(C_LOC(v1(ioffs+1,iD)), pv1, [nia2,nia2])
               Call C_F_POINTER(C_LOC(v2(ivoffs,iD)), pv2,
     &                          [nia2-nia1+1,nii2-nii1+1])
*
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
                     If (pv1(ia,ii).ne.-pv1(ii,ia)) Then
                        Write (6,*) 'inconsistency in gradient'
                        Call Abend()
                     End If
                     pv2(ia-nia1+1,ii-nii1+1)=pv1(ia,ii)
                  End Do
               End Do
#endif
*
            Else If ((n1.eq.nOV).and.(n2.eq.nOO)) Then
*
*           decompress
*
#ifdef POINTER_REMAP
               pv1(nia1:nia2,nii1:nii2) => v1(ivoffs:ivoffs+nv2-1,iD)
               pv2(1:nia2,1:nia2) => v2(ioffs+1:ioffs+nv1,iD)
*
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
                     pv2(ia,ii) = pv1(ia,ii)
                     pv2(ii,ia) =-pv1(ia,ii)
                  End Do
               End Do
#else
               Call C_F_POINTER(C_LOC(v1(ivoffs,iD)), pv1,
     &                          [nia2-nia1+1,nii2-nii1+1])
               Call C_F_POINTER(C_LOC(v2(ioffs+1,iD)), pv2, [nia2,nia2])
*
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
                     pv2(ia,ii) = pv1(ia-nia1+1,ii-nii1+1)
                     pv2(ii,ia) =-pv1(ia-nia1+1,ii-nii1+1)
                  End Do
               End Do
#endif
            End If

            Nullify(pv1,pv2)
            ivoffs=ivoffs+nv2
            ioffs=ioffs+(nOrb(iSym)*nOrb(iSym))
         End Do
*
#ifdef _DEBUGPRINT_
         Write (6,*) 'n1,n2,nOO,nOV=',n1,n2,nOO,nOV
         iOff = 1
         Do iSym = 1, nSym
            If (n2.eq.nOO) Then
               nO = nOrb(iSym)
               Call RecPrt('vOO2OV: v2',' ',v2(ioff,iD),nO,nO)
               iOff = iOff + nO**2
            Else If (n2.eq.nOV) Then
               nO = Max(nOcc(iSym,1),nOcc(iSym,2))-nFro(iSym)
               nV = nOrb(iSym)-Min(nOcc(iSym,1),nOcc(iSym,2))-nFro(iSym)
               Call RecPrt('vOO2OV: v2',' ',v2(ioff,iD),nO,nV)
               iOff = iOff + nO*nV
            End If
         End Do
#endif
*
      End Do ! iD
*
      Return
      End
