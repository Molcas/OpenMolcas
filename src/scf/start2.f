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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine Start2(FName,LuOrb,CMO,mBB,nD,Ovrlp,mBT,
     &                  EOrb,OccNo,mmB)
************************************************************************
*                                                                      *
*     purpose: Get starting orbitals INPORB                            *
*                                                                      *
*     called from: SOrb                                                *
*                                                                      *
*     calls to: Ortho                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "file.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#endif
      Real*8 CMO(mBB,nD), Ovrlp(mBT), EOrb(mmB,nD), OccNo(mmB,nD)
      Character FName*(*)
      Integer nTmp(8)
      Character*6 OrbName
* Pam 2012 Changed VECSORT arg list, need dummy array:
      Integer NewOrd(2)
      Integer, Dimension(:,:), Allocatable:: IndT
      Dimension Dummy(1),iDum(7,8)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Read vectors (vectors MUST !!!!! be written in symmetry blocks
*     of dimension nBas(i)*nBas(i))
*
      Call mma_allocate(IndT,nnB,nD,Label='IndT')
*
      Lu_=LuOrb
      nD = iUHF + 1
      If(iUHF.eq.0) Then
         If (isHDF5) Then
            Call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,
     &                      CMO,OccNo,EOrb,IndT)
         Else
            Call RdVec_(FName,Lu_,'COEI',iUHF,nSym,nBas,nOrb,
     &                  CMO,Dummy,
     &                  OccNo,Dummy,
     &                  EOrb(1,1),Dummy,
     &                  IndT(1,1),VTitle,1,iErr,iWFtype)
         End If
         Call VecSort(nSym,nBas,nBas,
     &               CMO,OccNo,IndT(1,1),0,NewOrd,iErr)
         indx=1
         Do iSym=1,nSym
            nTmp(iSym)=0
            Do iBas=1,nBas(iSym)
               If(IndT(indx,1).eq.7) nTmp(iSym)=nTmp(iSym)+1
               indx=indx+1
            End Do
            If(nOrb(iSym).gt.nBas(iSym)-nTmp(iSym)) Then
               nOrb(iSym)=nBas(iSym)-nTmp(iSym)
               nDel(iSym)=nTmp(iSym)
            End If
         End Do
         Call TrimCMO(CMO,CMO,nSym,nBas,nOrb)
         Call TrimEor(EOrb(1,1),EOrb(1,1),nSym,nBas,nOrb)

         Call Setup
         If(.not.Aufb .and. .not.OnlyProp) Then
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,1)
                  OccNo(iOrb+iOff,1)=2.0d0
               End Do
               Do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,1)=0.0d0
               End Do
               iOff=iOff+nOrb(iSym)
            End Do
         End If
      Else
         If (isHDF5) Then
            isUHF=0
#ifdef _HDF5_
            If (mh5_exists_dset(fileorb_id,'MO_ALPHA_VECTORS')) isUHF=1
#endif
         Else
            Call Chk_Vec_UHF(FNAME,Lu_,isUHF)
         End If
         If(isUHF.eq.1) Then
            If (isHDF5) Then
              Call RdVec_HDF5(fileorb_id,'COEIA',nSym,nBas,
     &                        CMO(1,1),OccNo(1,1),EOrb(1,1),IndT(1,1))
              Call RdVec_HDF5(fileorb_id,'COEIB',nSym,nBas,
     &                        CMO(1,2),OccNo(1,2),EOrb(1,2),IndT(1,2))
            Else
               Call RdVec_(FName,Lu_,'COEI',iUHF,nSym,nBas,nOrb,
     &                     CMO(1,1),CMO(1,2),
     &                     OccNo(1,1),OccNo(1,2),
     &                     EOrb(1,1),EOrb(1,2),
     &                     IndT(1,1),VTitle,1,iErr,iWFtype)
               Call iCopy(nnB,IndT(1,1),1,IndT(1,2),1)
            End If
            Call VecSort(nSym,nBas,nBas,CMO(1,1),OccNo(1,1),
     &                   IndT(1,1),0,NewOrd,iErr)
            Call VecSort(nSym,nBas,nBas,CMO(1,2),OccNo(1,2),
     &                   IndT(1,2),0,NewOrd,iErr)
            indx=1
            Do iSym=1,nSym
               nTmp(iSym)=0
               Do iBas=1,nBas(iSym)
                  If(IndT(indx,1).eq.7) nTmp(iSym)=nTmp(iSym)+1
                  indx=indx+1
               End Do
               If(nOrb(iSym).gt.nBas(iSym)-nTmp(iSym)) Then
                  nOrb(iSym)=nBas(iSym)-nTmp(iSym)
                  nDel(iSym)=nTmp(iSym)
               End If
            End Do
            Call TrimCMO(CMO(1,1),CMO(1,1),nSym,nBas,nOrb)
            Call TrimEor(EOrb(1,1),EOrb(1,2),nSym,nBas,nOrb)
            Call TrimCMO(CMO(1,2),CMO(1,2),nSym,nBas,nOrb)
            Call TrimEor(EOrb(1,2),EOrb(1,2),nSym,nBas,nOrb)
            Call Setup
         Else
            If (isHDF5) Then
              Call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,
     &                        CMO,OccNo,EOrb,IndT)
            Else
               Call RdVec_(FName,Lu_,'COEI',0,nSym,nBas,nOrb,
     &                     CMO,Dummy,
     &                     OccNo,Dummy,
     &                     EOrb(1,1),Dummy,
     &                     IndT(1,1),VTitle,1,iErr,iWFtype)
            End If
            Call VecSort(nSym,nBas,nBas,CMO,OccNo,
     &                   IndT(1,1),0,NewOrd,iErr)
            indx=1
            Do iSym=1,nSym
               nTmp(iSym)=0
               Do iBas=1,nBas(iSym)
                  If(IndT(indx,1).eq.7) nTmp(iSym)=nTmp(iSym)+1
                  indx=indx+1
               End Do
               If(nOrb(iSym).gt.nBas(iSym)-nTmp(iSym)) Then
                  nOrb(iSym)=nBas(iSym)-nTmp(iSym)
                  nDel(iSym)=nTmp(iSym)
               End If
            End Do
            Call TrimCMO(CMO,CMO,nSym,nBas,nOrb)
            Call TrimEor(EOrb(1,1),EOrb(1,1),nSym,nBas,nOrb)
            Call Setup
            Call dCopy_(nBO,CMO(1,1),1,CMO(1,2),1)
            Call dCopy_(nnB,OccNo(1,1),1,OccNo(1,2),1)
            Call dCopy_(nnB,EOrb(1,1),1,EOrb(1,2),1)
            Call dScal_(nnB,0.5d0,OccNo(1,1),1)
            Call dScal_(nnB,0.5d0,OccNo(1,2),1)
         End If
         If(.not.Aufb) Then
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,1)
                  OccNo(iOrb+iOff,1)=1.0d0
               End Do
               Do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,1)=0.0d0
               End Do
               iOff=iOff+nOrb(iSym)
            End Do
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,2)
                  OccNo(iOrb+iOff,2)=1.0d0
               End Do
               Do iOrb=nOcc(iSym,2)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,2)=0.0d0
               End Do
               iOff=iOff+nOrb(iSym)
            End Do
         End If
      End If
      If(Allocated(IndT)) Call mma_deallocate(IndT)
*
      If (MSYMON) Then
#ifdef _MSYM_
         Write(6,*) 'Symmetrizing start orbitals'
         Call fmsym_create_context(msym_ctx)
         Call fmsym_set_elements(msym_ctx)
         Call fmsym_find_symmetry(msym_ctx)
         Do iD = 1, nD
            Call fmsym_symmetrize_orbitals(msym_ctx,CMO(1,iD))
         End Do
#else
         Write(6,*) 'No msym support, skipping symmetrization '
     $           //'of start orbitals...'
#endif
      End If
*
      Do iD = 1, nD
         Call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
      End Do
*
#ifdef _MSYM_
      If (MSYMON) Then
         Call fmsym_release_context(msym_ctx)
      End If
#endif
*
* Dump orbitals
*
      If(iUHF.eq.0) then
         OrbName='SCFORB'
         Call WrVec_(OrbName,LuOut,'COE',iUHF,nSym,nBas,nBas,
     &               CMO,Dummy,OccNo,Dummy,
     &               EOrb(1,1),Dummy,iDum,VTitle,iWFtype)
      Else
         OrbName='UHFORB'
         Call WrVec_(OrbName,LuOut,'COE',iUHF,nSym,nBas,nBas,
     &               CMO(1,1),CMO(1,2),OccNo(1,1),OccNo(1,2),
     &               EOrb(1,1),EOrb(1,2),iDum,VTitle,iWFtype)
      End If
*
      Return
      End
