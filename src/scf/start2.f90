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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine Start2(FName,LuOrb,CMO,mBB,nD,Ovrlp,mBT,EOrb,OccNo,mmB)
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals INPORB                            *
!                                                                      *
!     called from: SOrb                                                *
!                                                                      *
!     calls to: Ortho                                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
#ifdef _MSYM_
      Use, Intrinsic :: iso_c_binding, only: c_ptr
#endif
#ifdef _HDF5_
      Use mh5, Only: mh5_exists_dset
#endif
      use InfSCF, only: Aufb, FileOrb_id, isHDF5, nBO, nBT, nSym, OnlyProp, VTitle, nOcc, nOrb, nBas, nnB, nDel
      use InfSCF, only: mSymON
      use Files, only: LuOut
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, Half, One, Two
      Implicit None
      Character(LEN=*) FName
      Integer LuOrb, mBB,nD,mBT,mmB
      Real*8 CMO(mBB,nD), Ovrlp(mBT), EOrb(mmB,nD), OccNo(mmB,nD)

      Integer iBas, iD, iErr, indx, iOff, iOrb, isUHF, iSym, Lu_
      Integer nTmp(8), iWFtype
      Character(LEN=6) OrbName
! Pam 2012 Changed VECSORT arg list, need dummy array:
      Integer iDummy(1)
      Integer, Dimension(:,:), Allocatable:: IndT
#ifdef _MSYM_
      Type(c_ptr) msym_ctx
#endif
      Integer iDum(7,8)
      Real*8 Dummy(1)
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!---- Read vectors (vectors MUST !!!!! be written in symmetry blocks
!     of dimension nBas(i)*nBas(i))
!
      Call mma_allocate(IndT,nnB,nD,Label='IndT')
!
      Lu_=LuOrb
      If(nD==1) Then
         If (isHDF5) Then
            Call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
         Else
            Call RdVec_(FName,Lu_,'COEI',nD-1,nSym,nBas,nOrb,CMO,Dummy,OccNo,Dummy,      &
                        EOrb(1,1),Dummy,IndT(1,1),VTitle,1,iErr,iWFtype)
         End If
         Call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT(1,1),0,iDummy,iErr)
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

         Call Setup_SCF()
         If(.not.Aufb .and. .not.OnlyProp) Then
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,1)
                  OccNo(iOrb+iOff,1)=Two
               End Do
               Do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,1)=Zero
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
              Call RdVec_HDF5(fileorb_id,'COEIA',nSym,nBas,CMO(1,1),OccNo(1,1),EOrb(1,1),IndT(1,1))
              Call RdVec_HDF5(fileorb_id,'COEIB',nSym,nBas,CMO(1,2),OccNo(1,2),EOrb(1,2),IndT(1,2))
            Else
               Call RdVec_(FName,Lu_,'COEI',nD-1,nSym,nBas,nOrb,CMO(1,1),CMO(1,2),OccNo(1,1),OccNo(1,2),   &
                           EOrb(1,1),EOrb(1,2),IndT(1,1),VTitle,1,iErr,iWFtype)
               Call iCopy(nnB,IndT(1,1),1,IndT(1,2),1)
            End If
            Call VecSort(nSym,nBas,nBas,CMO(1,1),OccNo(1,1),IndT(1,1),0,iDummy,iErr)
            Call VecSort(nSym,nBas,nBas,CMO(1,2),OccNo(1,2),IndT(1,2),0,iDummy,iErr)
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
            Call Setup_SCF()
         Else
            If (isHDF5) Then
              Call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
            Else
               Call RdVec_(FName,Lu_,'COEI',0,nSym,nBas,nOrb,CMO,Dummy,OccNo,Dummy,       &
                           EOrb(1,1),Dummy,IndT(1,1),VTitle,1,iErr,iWFtype)
            End If
            Call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT(1,1),0,iDummy,iErr)
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
            Call Setup_SCF()
            Call dCopy_(nBO,CMO(1,1),1,CMO(1,2),1)
            Call dCopy_(nnB,OccNo(1,1),1,OccNo(1,2),1)
            Call dCopy_(nnB,EOrb(1,1),1,EOrb(1,2),1)
            Call dScal_(nnB,Half,OccNo(1,1),1)
            Call dScal_(nnB,Half,OccNo(1,2),1)
         End If
         If(.not.Aufb) Then
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,1)
                  OccNo(iOrb+iOff,1)=One
               End Do
               Do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,1)=Zero
               End Do
               iOff=iOff+nOrb(iSym)
            End Do
            iOff=0
            Do iSym=1,nSym
               Do iOrb=1,nOcc(iSym,2)
                  OccNo(iOrb+iOff,2)=One
               End Do
               Do iOrb=nOcc(iSym,2)+1,nOrb(iSym)
                  OccNo(iOrb+iOff,2)=Zero
               End Do
               iOff=iOff+nOrb(iSym)
            End Do
         End If
      End If
      If(Allocated(IndT)) Call mma_deallocate(IndT)
!
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
         Write(6,*) 'No msym support, skipping symmetrization of start orbitals...'
#endif
      End If
!
      Do iD = 1, nD
         Call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
      End Do
!
#ifdef _MSYM_
      If (MSYMON) Then
         Call fmsym_release_context(msym_ctx)
      End If
#endif
!
! Dump orbitals
!
      If(nD==1) then
         OrbName='SCFORB'
         Call WrVec_(OrbName,LuOut,'COE',nD-1,nSym,nBas,nBas,CMO,Dummy,OccNo,Dummy,     &
                     EOrb(1,1),Dummy,iDum,VTitle,iWFtype)
      Else
         OrbName='UHFORB'
         Call WrVec_(OrbName,LuOut,'COE',nD-1,nSym,nBas,nBas,CMO(1,1),CMO(1,2),OccNo(1,1),OccNo(1,2),  &
                     EOrb(1,1),EOrb(1,2),iDum,VTitle,iWFtype)
      End If
!
      Return
      End SubRoutine Start2
