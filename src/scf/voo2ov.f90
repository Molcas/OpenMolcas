!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine vOO2OV(v1,nOO,v2,mOV,nD,kOV)
      use Constants, only: Zero
      Implicit None
      Integer nOO, mOV, nD, iSt, iEnd, iD
      Integer kOV(nD)
      Real*8 v1(nOO,nD), v2(mOV)

      Interface
        Subroutine vOO2OV_internal(v1,n1,v2,n2,iD)
        Integer n1, n2, iD
        Real*8, Target ::  v1(n1), v2(n2)
        End Subroutine vOO2OV_internal
      End Interface

      iEnd = 0
      v2(:)=Zero
      Do iD = 1, nD
         iSt = iEnd + 1
         iEnd = iEnd + kOV(iD)
         Call vOO2OV_internal(v1(:,iD),nOO,v2(iSt:iEnd),kOV(iD),iD)
      End Do

      End SubRoutine vOO2OV
      SubRoutine vOO2OV_internal(v1,n1,v2,n2,iD)
!***********************************************************************
!                                                                      *
!     purpose: converts vector of dim nOO (e.g. gradient) to vector    *
!              of lower dimension nOV (occ-virt block only)            *
!              if n1==nOO && n2==nOV  -> compress                      *
!              else if if n1==nOV && n2==nOO -> decompress             *
!                                                                      *
!     output:                                                          *
!       v2      : compressed or decompressed vector, dep. on input     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. Schuetz                                                       *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!     Roland Lindh                                                     *
!     Uppsala University, Sweden, 2017                                 *
!     Fortran pointer structure                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
!
#include "compiler_features.h"
#ifndef POINTER_REMAP
      Use, Intrinsic :: ISO_C_BINDING
#endif
      use InfSCF, only: nOO, nSym, nFro, nOcc, nOrb, kOV
      Implicit None
!
!     declaration subroutine parameters
      Integer n1,n2,iD
      Real*8, Target:: v1(n1), v2(n2)
      Real*8, Dimension(:,:), Pointer:: pv1, pv2
!
!     declaration local variables
      Integer iSym,ii,ia,ioffs,ivoffs
      Integer nia1, nia2, nii1, nii2, nv1, nv2
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
#ifdef _DEBUGPRINT_
      Write (6,*) 'n1,n2,nOO,kOV(:)=',n1,n2,nOO,kOV(:)
      iOff = 1
      Do iSym = 1, nSym
         If (n1.eq.nOO) Then
            nO = nOrb(iSym)
            Call RecPrt('vOO2OV: v1',' ',v1(ioff),nO,nO)
            iOff = iOff + nO**2
         Else If (n1.eq.kOV(iD)) Then
            nO = nOcc(iSym,iD)-nFro(iSym)
            nV = nOrb(iSym)-nOcc(iSym,iD)-nFro(iSym)
            Call RecPrt('vOO2OV: v1',' ',v1(ioff),nO,nV)
            iOff = iOff + nO*nV
         End If
      End Do
#endif
         ioffs=0
         ivoffs=1
         Do iSym=1,nSym
! range for occupied non-frozen orbitals
            nii1=nFro(iSym)+1
            nii2=nOcc(iSym,iD)
! range for virtual orbitals
            nia1=nOcc(iSym,iD)+1
            nia2=nOrb(iSym)
! size of the full block
            nv1=nia2**2
! size of the virtual-occupied block
            nv2=(nia2-nia1+1)*(nii2-nii1+1)
!
            If ((n1.eq.nOO).and.(n2.eq.kOV(iD))) Then
!
!           compress
!
#ifdef POINTER_REMAP
               pv1(1:nia2,1:nia2) => v1(ioffs+1:ioffs+nv1)
               pv2(nia1:nia2,nii1:nii2) => v2(ivoffs:ivoffs+nv2-1)
!
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
!                    If (pv1(ia,ii).ne.-pv1(ii,ia)) Then
!                       Write (6,*) 'inconsistency in gradient'
!                       Call Abend()
!                    End If
                     pv2(ia,ii)=pv1(ia,ii)
                  End Do
               End Do
#else
               Call C_F_POINTER(C_LOC(v1(ioffs+1)), pv1, [nia2,nia2])
               Call C_F_POINTER(C_LOC(v2(ivoffs)), pv2,[nia2-nia1+1,nii2-nii1+1])
!
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
!                    If (pv1(ia,ii).ne.-pv1(ii,ia)) Then
!                       Write (6,*) 'inconsistency in gradient'
!                       Call Abend()
                     End If
                     pv2(ia-nia1+1,ii-nii1+1)=pv1(ia,ii)
                  End Do
               End Do
#endif
!
            Else If ((n1.eq.kOV(iD)).and.(n2.eq.nOO)) Then
!
!           decompress
!
#ifdef POINTER_REMAP
               pv1(nia1:nia2,nii1:nii2) => v1(ivoffs:ivoffs+nv2-1)
               pv2(1:nia2,1:nia2) => v2(ioffs+1:ioffs+nv1)
!
               Do ii=nii1,nii2
                  Do ia=nia1,nia2
                     pv2(ia,ii) = pv1(ia,ii)
                     pv2(ii,ia) =-pv1(ia,ii)
                  End Do
               End Do
#else
               Call C_F_POINTER(C_LOC(v1(ivoffs)), pv1,
     &                          [nia2-nia1+1,nii2-nii1+1])
               Call C_F_POINTER(C_LOC(v2(ioffs+1)), pv2, [nia2,nia2])
!
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
!
#ifdef _DEBUGPRINT_
         Write (6,*) 'n1,n2,nOO,kOV=',n1,n2,nOO,kOV(:)
         iOff = 1
         Do iSym = 1, nSym
            If (n2.eq.nOO) Then
               nO = nOrb(iSym)
               Call RecPrt('vOO2OV: v2',' ',v2(ioff),nO,nO)
               iOff = iOff + nO**2
            Else If (n2.eq.kOV(iD)) Then
               nO = nOcc(iSym,iD)-nFro(iSym)
               nV = nOrb(iSym)-nOcc(iSym,iD)-nFro(iSym)
               Call RecPrt('vOO2OV: v2',' ',v2(ioff),nO,nV)
               iOff = iOff + nO*nV
            End If
         End Do
#endif
!
      Return
      End
