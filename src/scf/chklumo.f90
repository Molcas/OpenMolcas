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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!#define _DEBUGPRINT_
Subroutine ChkLumo(OccSet,FermSet,SpinSet)
!***********************************************************************
!                                                                      *
! This routine figure out the status of the lumo file, i.e. should it  *
! trigger keywords OCCUpied or FERMi?                                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************
#ifdef _HDF5_
      Use mh5, Only: mh5_exists_dset
#endif
      use InfSCF, only: nSym, nD, SCF_FileOrb, isHDF5, Tot_EL_Charge, iAU_AB, nOcc, nBas, nOrb, vTitle, FileOrb_ID, nSym
#ifdef _DEBUGPRINT_
      use InfSCF, only: Tot_Charge, Tot_Nuc_Charge
#endif
      use Constants, only: Zero, Half, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
      Logical OccSet
      Logical FermSet
      Logical SpinSet
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Character(LEN=512) FNAME
      Logical      Idem
      Real*8, Dimension(:,:), Allocatable:: OccVec, EpsVec
      Real*8 Dummy(1), qA, qB, Tmp, Tmp1
      Integer nVec, iSym, LU, isUHF, LU_, I, iDiff, iOff, N, iBas, iDummy(1), iErr, iWFType
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
      Call Peek_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      nVec=0
      Do iSym=1,nSym
         nVec=nVec+nBas(iSym)
      End Do
!----------------------------------------------------------------------*
! Allocate fields                                                      *
!----------------------------------------------------------------------*
      Call mma_Allocate(OccVec,nVec,nD,Label='OccVec')
      Call mma_Allocate(EpsVec,nVec,nD,Label='EpsVec')
!----------------------------------------------------------------------*
! Read occupation numbers and orbital energies                         *
!----------------------------------------------------------------------*
      Lu=17
      FName=SCF_FileOrb
      If(nD==1) Then
         If (isHDF5) Then
            Call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
         Else
            Call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),Dummy,EpsVec(1,1), &
                        Dummy,iDummy,VTitle,1,iErr,iWFtype)
         End If
      Else
         isUHF=0
         If (isHDF5) Then
#ifdef _HDF5_
            If (mh5_exists_dset(fileorb_id,'MO_ALPHA_VECTORS')) isUHF=1
#endif
         Else
            Lu_=18
            isUHF=-1
            Call Chk_Vec_UHF(FNAME,Lu_,isUHF)
         End If
         If(isUHF.eq.1) Then
            If (isHDF5) Then
               Call RdVec_HDF5(fileorb_id,'OEA',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
               Call RdVec_HDF5(fileorb_id,'OEB',nSym,nBas,Dummy,OccVec(1,2),EpsVec(1,2),iDummy)
            Else
               Call RdVec_(FNAME,Lu,'OE',nD-1,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),OccVec(1,2),EpsVec(1,1), &
                           EpsVec(1,2),iDummy,VTitle,1,iErr,iWFtype)
            End If
         Else
            If (isHDF5) Then
               Call RdVec_HDF5(fileorb_id,'OE',nSym,nBas,Dummy,OccVec(1,1),EpsVec(1,1),iDummy)
            Else
               Call RdVec_(FNAME,Lu,'OE',0,nSym,nBas,nOrb,Dummy,Dummy,OccVec(1,1),Dummy,EpsVec(1,1),Dummy, &
                           iDummy,VTitle,1,iErr,iWFtype)
            End If
            Call dCopy_(nVec,OccVec(1,1),1,OccVec(1,2),1)
            Call dCopy_(nVec,EpsVec(1,1),1,EpsVec(1,2),1)
            Call dScal_(nVec*nD,half,OccVec,1)
         End If
      End If
#ifdef _DEBUGPRINT_
      If(nD==1) Then
         Write(6,*) 'Orbital energies'
         Write(6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
         Write(6,*) 'Occupation numbers'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
      Else
         Write(6,*) 'Alpha orbital energies'
         Write(6,'(10f12.6)') (EpsVec(i,1),i=1,nVec)
         Write(6,*) 'Alpha occupation numbers'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
         Write(6,*) 'Beta orbital energies'
         Write(6,'(10f12.6)') (EpsVec(i,2),i=1,nVec)
         Write(6,*) 'Beta occupation numbers'
         Write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
      End If
#endif
!----------------------------------------------------------------------*
! What are the charges                                                 *
!----------------------------------------------------------------------*
      qa=Zero
      qb=Zero
      If(nD==1) Then
         If (nSym/=1) Then
            tmp1 = Zero
            Do i=1,nVec
               tmp1 = tmp1 + OccVec(i,1)
               If (tmp1>=Two) Then
                  qa=qa+Two
                  tmp1=tmp1-Two
                  OccVec(i,1)=Two
               Else
                  OccVec(i,1)=Zero
               End If
            End Do
         Else
            tmp1 = Zero
            Do i=1,nVec
               tmp1 = tmp1 + OccVec(i,1)
            End Do
            OccVec(:,:)=Zero
            tmp1 = DBLE(NINT(tmp1))
            Do i=1,(NINT(tmp1)+1)/2
               If (tmp1>=Two) Then
                  qa=qa+Two
                  OccVec(i,1)=Two
                  tmp1=tmp1-Two
               End If
            End Do
         End If
         qa=Half*qa
         qb=qa
      Else
         If (nSym/=1) Then
            tmp1 = Zero
            Do i=1,nVec
               tmp1 = tmp1 + OccVec(i,1) + OccVec(i,2)
               If (tmp1>=Two) Then
                  qa=qa+One
                  qb=qb+One
                  tmp1=tmp1-Two
                  OccVec(i,1)=One
                  OccVec(i,2)=One
               Else If (tmp1>=One) Then
                  qa=qa+One
                  qb=qb+Zero
                  OccVec(i,1)=One
                  OccVec(i,2)=Zero
                  tmp1=tmp1-One
               Else
                  OccVec(i,1)=Zero
                  OccVec(i,2)=Zero
               End If
            End Do
         Else
           tmp1 = Zero
           Do i=1,nVec
              tmp1 = tmp1 + OccVec(i,1) + OccVec(i,2)
           End Do
           OccVec(:,:)=Zero
           tmp1 = DBLE(NINT(tmp1))
           Do i=1,(NINT(tmp1)+1)/2
              If (tmp1>=Two) Then
                 qa=qa+One
                 qb=qb+One
                 OccVec(i,1)=One
                 OccVec(i,2)=One
                 tmp1=tmp1-Two
               Else If (tmp1>=One) Then
                 qa=qa+One
                 qb=qb+Zero
                 OccVec(i,1)=One
                 OccVec(i,2)=Zero
                 tmp1=tmp1-One
               Else
                  OccVec(i,1)=Zero
                  OccVec(i,2)=Zero
               End If
           End Do
         End If
      End If
#ifdef _DEBUGPRINT_
      If(nD==1) Then
         Write(6,*) 'chklumo: Idempotency'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
         Write(6,'(a,f12.6)') 'Tot charge         ',Tot_charge
         Write(6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
         Write(6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
         Write(6,'(a,f12.6)') 'Electron count     ',Two*qa
      Else
         Write(6,*) 'chklumo: Alpha idempotency'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
         Write(6,*) 'chklumo: Beta idempotency'
         Write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
         Write(6,'(a,f12.6)') 'Tot charge         ',Tot_charge
         Write(6,'(a,f12.6)') 'Tot nuc. charge    ',Tot_nuc_charge
         Write(6,'(a,f12.6)') 'Tot el. charge     ',Tot_el_charge
         Write(6,'(a,f12.6)') 'Alpha count        ',qa
         Write(6,'(a,f12.6)') 'Beta count         ',qb
      End If
#endif
!----------------------------------------------------------------------*
! Is it the same charge.                                               *
!----------------------------------------------------------------------*
      If(Abs(qa+qb+Tot_el_charge).gt.half) Then
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: System has changed charge!'
#endif
         Occset=.false.
         FermSet=.true.
         Goto 999
      End If
!----------------------------------------------------------------------*
! Is it the same spin.                                                 *
!----------------------------------------------------------------------*
      If(SpinSet) Then
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: System might have changed spin!'
         Write(6,*) '   iAu_ab=',iAu_ab
         Write(6,*) '   qa-qb=',qa-qb
#endif
         idiff=iAu_ab-Int(qa-qb)
         If(idiff.ne.0) Then
#ifdef _DEBUGPRINT_
            Write(6,*) '   yes indeed, spin has changed!'
#endif
            Occset=.false.
            FermSet=.true.
            Goto 999
         End If
      End If
!----------------------------------------------------------------------*
! Is it idempotent     D^2 = 2 D                                       *
!----------------------------------------------------------------------*
      If(nD==1) Then
         Idem=.true.
         Do i=1,nVec
            tmp=half*OccVec(i,1)*(One-half*OccVec(i,1))
            If(Abs(tmp).gt.0.25d0) Idem=.false.
         End Do
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: Idempotency'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
#endif
      Else
         Idem=.true.
         Do i=1,nVec
            tmp=OccVec(i,1)*(One-OccVec(i,1))
            If(Abs(tmp).gt.0.25d0) Idem=.false.
         End Do
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: Alpha idempotency'
         Write(6,'(10f12.6)') (OccVec(i,1),i=1,nVec)
#endif
         Do i=1,nVec
            tmp=OccVec(i,2)*(One-OccVec(i,2))
            If(Abs(tmp).gt.0.25d0) Idem=.false.
         End Do
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: Beta idempotency'
         Write(6,'(10f12.6)') (OccVec(i,2),i=1,nVec)
#endif
      End If
!----------------------------------------------------------------------*
! Was it idempotent?                                                   *
!----------------------------------------------------------------------*
      If(Idem) Then
#ifdef _DEBUGPRINT_
         Write(6,*) 'chklumo: Idempotent'
#endif
         If(nD==1) Then
            iOff=0
            Do iSym=1,nSym
               n=0
               Do iBas=1,nBas(iSym)
                  n = n + INT(OccVec(iOff+iBas,1))
               End Do
               nOcc(iSym,1)=n/2
               iOff=iOff+nBas(iSym)
            End Do
#ifdef _DEBUGPRINT_
            Write(6,'(a,8i5)') 'Occupation       ',(nOcc(i,1),i=1,nSym)
#endif
         Else
            iOff=0
            Do iSym=1,nSym
               n=0
               Do iBas=1,nBas(iSym)
                  n = n + INT(OccVec(iOff+iBas,1))
               End Do
               nOcc(iSym,1)=n
               iOff=iOff+nBas(iSym)
            End Do
            iOff=0
            Do iSym=1,nSym
               n=0
               Do iBas=1,nBas(iSym)
                  n = n + INT(OccVec(iOff+iBas,2))
               End Do
               nOcc(iSym,2)=n
               iOff=iOff+nBas(iSym)
            End Do
#ifdef _DEBUGPRINT_
            Write(6,'(a,8i5)') 'Alpha occupation ',(nOcc(i,1),i=1,nSym)
            Write(6,'(a,8i5)') 'Beta occupation  ',(nOcc(i,2),i=1,nSym)
#endif
         End If
         Occset=.true.
         FermSet=.false.
      Else
#ifdef _DEBUGPRINT_
         Write(6,*) 'Not idempotent'
#endif
         Occset=.false.
         FermSet=.true.
      End If
!----------------------------------------------------------------------*
! Deallocate fields                                                    *
!----------------------------------------------------------------------*
999   Continue
      Call mma_deallocate(EpsVec)
      Call mma_deallocate(OccVec)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
      Return
      End Subroutine ChkLumo
