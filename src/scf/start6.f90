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
! Copyright (C) 2017, Roland Lindh                                     *
!***********************************************************************
      SubRoutine Start6(FName,LuOrb,CMO,mBB,nD,EOrb,OccNo,mmB)
!***********************************************************************
!                                                                      *
!     purpose: Generate constrained orbitals from INPORB               *
!                                                                      *
!***********************************************************************
      use OneDat, only: sNoNuc, sNoOri
      use SpinAV, only: DSC, Do_SpinAV
      use InfSCF, only: nBas, nOrb, nOcc, nFro, nDel, nConstr, IndxC, DoCholesky, E_nondyn, FileOrb_id, isHDF5, MaxBas, &
                        MxConstr, nBB, nBT, nnB, nSym, VTitle
      use Cholesky, only: ChFracMem
      use DCSCF, only: Erest_xc, s2CNO
      use Constants, only: Zero, Half, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Character(LEN=*) FName
      Integer mBB, nD, mmB, LuOrb
      Real*8 CMO(mBB,nD), EOrb(mmB,nD), OccNo(mmB,nD)
!
      Integer i, ibas, ic1, ic2, iCMO, iComp, iDaa, iDbb, iDSc, iErr, Indx, iOcc, iOff, iOpt, iOrb, ipDaa, ipDbb, ipDScc, ipMK,   &
              ipML, iRC, iSym, iSymLbl, iWFType, j, jc, ji, jOcc, jOff, k, kc, kc1, kc2, kDSc, kk, kkc, kks, l, lc, lc1, lc2, llc,&
              lls, lOcc, lOff, lsq, ltri, Lu_, mAdCMOO, mOff, nDiff_ab
      Real*8 ThrD, xNorm, xOkk, yOkk
      Character(LEN=62) Line
      Integer nTmp(8), nIF(8), nRASO(8), nBD(8), nZero(8), nHoles(8)
      Integer nSsh(8), nSsh_ab(8)
! Pam 2012 Changed VECSORT arg list, need dummy array:
      Integer iDummy(1)
      Integer, Dimension(:), Allocatable:: IndT, ID_vir
      Real*8, Allocatable:: Da(:,:)
      Integer, Allocatable:: Match(:,:)
      Real*8, Dimension(:), Allocatable:: Corb, SAV, SLT, SQ
      Real*8 Dummy(1)
      character(len=8) :: Label
!***********************************************************************
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
      Erest_xc=Zero
!
      If(.not.DoCholesky) then
       write(6,*)
       write(6,*) ' ERROR in Constrained SCF: problem in start6.'
       write(6,*) '*** Constrained NOs implemented only with CD or RI.'
       write(6,*) '*** Use Cholesky or RICD in Seward and rerun! *****'
       Call Abend
      Endif
!
      Do iSym=1,nSym
         nDiff_ab=nOcc(iSym,1)-nOcc(iSym,2)
         If (nDiff_ab .lt. 0) Then
            write(6,*)
            write(6,*) ' ERROR in Constrained SCF: problem in start6.'
            write(6,*) '*** #alpha < #beta not permitted in CNOs    ***'
            write(6,*) '*** Change SCF input accordingly and rerun! ***'
            Call Abend
         EndIf
         nHoles(iSym)=nDiff_ab
      End Do
!
      write(6,*) ' -------------------------------------------------'
      If (Do_SpinAV) Then
       write(6,*)' Spin-averaged wavelets (+/-) '
      Else
       write(6,*)' Configuration of the constrained spins (up/down) '
      EndIf
      write(6,*) ' -------------------------------------------------'
      Do iSym=1,nSym
         write(6,'(1X,A,I1)') ' sym: ',iSym
         Line(1:14)='         (+) '
         k=15
         Do j=1,nConstr(iSym)
            If (indxC(j,1,iSym).eq.1) Then
               If (Do_SpinAV) Then
                  Line(k:k+2)=' + '
               Else
                  Line(k:k+2)=' u '
               EndIf
            ElseIf (indxC(j,1,iSym).eq.2) Then
               If (Do_SpinAV) Then
                  Line(k:k+2)=' - '
               Else
                  Line(k:k+2)=' d '
               EndIf
            Else
               Line(k:k+2)='   '
            EndIf
            k=k+3
         End Do
         write(6,*) Line(1:k-1)
         Line(1:14)='         (-) '
         k=15
         Do j=1,nConstr(iSym)
            If (indxC(j,2,iSym).eq.1) Then
               If (Do_SpinAV) Then
                  Line(k:k+2)=' + '
               Else
                  Line(k:k+2)=' u '
               EndIf
            ElseIf (indxC(j,2,iSym).eq.2) Then
               If (Do_SpinAV) Then
                  Line(k:k+2)=' - '
               Else
                  Line(k:k+2)=' d '
               EndIf
            Else
               Line(k:k+2)='   '
            EndIf
            k=k+3
         End Do
         write(6,*) Line(1:k-1)
      End Do
      write(6,*) ' -------------------------------------------------'
!
      Lu_=LuOrb
      Call mma_Allocate(IndT,nnB,Label='IndT')
      If (isHDF5) Then
         Call RdVec_HDF5(fileorb_id,'COEI',nSym,nBas,CMO,OccNo,EOrb,IndT)
      Else
         Call RdVec_(FName,Lu_,'COEI',0,nSym,nBas,nOrb,    &
                     CMO,Dummy,                            &
                     OccNo,Dummy,                          &
                     EOrb,Dummy,                           &
                     IndT,VTitle,1,iErr,iWFtype)
      End If
      Call RdTwoEnrg(Lu_,E_nondyn)
      Call VecSort(nSym,nBas,nBas,CMO,OccNo,IndT,0,iDummy,iErr)
      indx=1
      Do iSym=1,nSym
         nZero(iSym)=0
         nTmp(iSym)=0
         nIF(iSym)=0
         nRASO(iSym)=0
         nDiff_ab=0
         Do iBas=1,nBas(iSym)
            If(IndT(indx).eq.1 .or. IndT(indx).eq.2) Then
               nIF(iSym)=nIF(iSym)+1 ! froz + inac orbitals
            End If
            If(IndT(indx).eq.3) Then ! a trick for odd number of
               nIF(iSym)=nIF(iSym)+1  ! electrons (place them in RAS1)
               nDiff_ab=nDiff_ab+1
            End If
            If(IndT(indx).gt.3 .and. IndT(indx).lt.6) Then
               nRASO(iSym)=nRASO(iSym)+1
            End If
            If(IndT(indx).eq.7) Then
               nTmp(iSym)=nTmp(iSym)+1
            End If
            indx=indx+1
         End Do
         If (nRASO(iSym).ne.2*nConstr(iSym)) Then
            write(6,*) ' ERROR in Constrained SCF: problem in start6.'
            write(6,*) ' Detected inconsistency between # of partially occupied orbitals and # of constraints. Sym: ',iSym
            Call Abend
         EndIf
         If (nHoles(iSym).ne.nDiff_ab) Then
            write(6,*) ' ERROR in Constrained SCF: problem in start6.'
            write(6,*) ' Detected inconsistency between # of excess alpha orbitals and # of RAS1 orbitals. Sym: ',iSym
            Call Abend
         EndIf
         If(nOrb(iSym).gt.nBas(iSym)-nTmp(iSym)) Then
            nOrb(iSym)=nBas(iSym)-nTmp(iSym)
            nDel(iSym)=nTmp(iSym)
         End If
      End Do
!
      Call TrimCMO(CMO,CMO,nSym,nBas,nOrb)
      Call TrimEor(EOrb,EOrb,nSym,nBas,nOrb)
      Call mma_deallocate(IndT)
!
      Call Setup_SCF()
!
      Call Izero(nBD,nSym)
      Do iSym=2,nSym
         nBD(iSym) = nBD(iSym-1) + nBas(iSym-1)*(nBas(iSym-1)+1)/2
      End Do
      Call mma_allocate(Da,nBT,2,Label='Da')
      Da(:,:)=Zero
      Call mma_allocate(Match,2,MxConstr,Label='Match')
      Call mma_allocate(Corb,MaxBas,Label='Corb')
!
      If (Do_SpinAV) call mma_allocate(SAV,2*MaxBas**2,Label='SAV')
!
      iOff=1
      jOff=0
      Do iSym=1,nSym
         call dcopy_(nBas(iSym)*nOrb(iSym),CMO(iOff,1),1,CMO(iOff,2),1)
         call dcopy_(nOrb(iSym),EOrb(1+jOff,1),1,EOrb(1+jOff,2),1)
         lOcc=1+jOff+nIF(iSym)
         call dcopy_(nRASO(iSym),OccNo(lOcc,1),1,OccNo(1,2),1)
         Call BestMatch(nConstr(iSym),nRASO(iSym),OccNo(1,2),Match,MxConstr)
         Do i=1,nConstr(iSym)
            k=Match(1,i) ! (+) wavelet
            jOcc=jOff+nIF(iSym)+k
            xOkk=OccNo(jOcc,1)/Two
            kc=iOff+nBas(iSym)*(nIF(iSym)+k-1)
            l=Match(2,i)  ! (-) wavelet
            iOcc=jOff+nIF(iSym)+l
            yOkk=OccNo(iOcc,1)/Two
            xnorm=sqrt(abs(xOkk)+abs(yOkk)) ! ensures correct normaliz
            lc=iOff+nBas(iSym)*(nIF(iSym)+l-1)
            xOkk=sqrt(abs(xOkk))/xnorm
            yOkk=sqrt(abs(yOkk))/xnorm
            If (Do_SpinAV) Then
               kkc=1+nBas(iSym)*(k-1)
               llc=1+nBas(iSym)*(nConstr(iSym)+l-1)
               call dcopy_(nBas(iSym),CMO(kc,1),1,SAV(kkc),1)
               call dcopy_(nBas(iSym),CMO(lc,1),1,SAV(llc),1)
               call dscal_(nBas(iSym),yOkk,SAV(kkc),1)
               call dscal_(nBas(iSym),xOkk,SAV(llc),1)
               call dcopy_(nBas(iSym),SAV(kkc),1,Corb,1)
               call dcopy_(nBas(iSym),CMO(kc,1),1,SAV(kkc),1)
               Call daxpy_(nBas(iSym),-One,SAV(llc),1,SAV(kkc),1)
               Call daxpy_(nBas(iSym),One,Corb,1,SAV(llc),1)
            EndIf
            call dscal_(nBas(iSym),xOkk,CMO(kc,1),1)
            call dscal_(nBas(iSym),yOkk,CMO(lc,1),1)
            call dcopy_(nBas(iSym),CMO(lc,1),1,Corb,1)
            Call daxpy_(nBas(iSym), One,CMO(kc,1),1,Corb,1)
            Call daxpy_(nBas(iSym),-One,CMO(kc,1),1,CMO(lc,1),1)
            call dscal_(nBas(iSym),-One,CMO(lc,1),1)
            call dcopy_(nBas(iSym),Corb,1,CMO(kc,1),1)
         End Do
         jc=1
         kc=nConstr(iSym)+1
         Do i=1,nConstr(iSym)
            l=Match(indxC(i,2,iSym),i)
            lc1=iOff+nBas(iSym)*(nIF(iSym)+l-1)
            lc2=iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+jc-1)
            call dcopy_(nBas(iSym),CMO(lc1,1),1,CMO(lc2,2),1)
            k=Match(indxC(i,1,iSym),i)
            kc1=iOff+nBas(iSym)*(nIF(iSym)+k-1)
            kc2=iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+kc-1)
            call dcopy_(nBas(iSym),CMO(kc1,1),1,CMO(kc2,2),1)
            jc=jc+1
            kc=kc+1
         End Do
         kc=nConstr(iSym)+1
         Do i=1,nConstr(iSym)
            ic1=iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+i-1)
            ic2=iOff+nBas(iSym)*(nIF(iSym)+kc-1)
            call dcopy_(nBas(iSym),CMO(ic1,2),1,CMO(ic2,1),1)
            kc1=iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym)+kc-1)
            kc2=iOff+nBas(iSym)*(nIF(iSym)+i-1)
            call dcopy_(nBas(iSym),CMO(kc1,2),1,CMO(kc2,1),1)
            kc=kc+1
         End Do
         kc=nConstr(iSym)+1  ! wavelets in virt space
         If (Do_SpinAV) Then
            Do i=1,nConstr(iSym)
               k=Match(1,i)
               kks=1+nBas(iSym)*(k-1)
               l=Match(2,i)
               lls=1+nBas(iSym)*(nConstr(iSym)+l-1)
               mOff=iOff+nBas(iSym)*(nIF(iSym)+kc-1)
               kk=indxC(i,1,iSym)
               If (kk.eq.1) Then ! => (+) wavelet is in alpha
                  ipMK=mOff
                  ipML=mOff-nBas(iSym)*nHoles(iSym)
                  call dcopy_(nBas(iSym),SAV(kks),1,CMO(ipMK,1),1)
                  call dcopy_(nBas(iSym),SAV(lls),1,CMO(ipML,2),1)
               ElseIf (kk.eq.2) Then
                  ipMK=mOff-nBas(iSym)*nHoles(iSym)
                  ipML=mOff
                  call dcopy_(nBas(iSym),SAV(kks),1,CMO(ipMK,2),1)
                  call dcopy_(nBas(iSym),SAV(lls),1,CMO(ipML,1),1)
               Else
                  ipMK=666666  ! avoid compiler wrngs
                  ipML=666666
                  write(6,*) ' Start6: wrong indxC value: ',kk
                  Call Abend()
               EndIf
               kc=kc+1
            End Do
         EndIf
         iOff=iOff+nBas(iSym)*nOrb(iSym)
         jOff=jOff+nOrb(iSym)
      End Do
!
      If (Do_SpinAV) Then
         Call mma_deallocate(SAV)
         Call mma_Allocate(DSc,nBB,Label='DSc')
         DSC(:)=Zero
      EndIf
!
      iOff=1
      lOff=0
      Do iSym=1,nSym
         ipDaa=1+nBD(iSym)
         mAdCMOO=iOff+nBas(iSym)*nIF(iSym)
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),            &
                          One,CMO(mAdCMOO,1),nBas(iSym),                      &
                                CMO(mAdCMOO,1),nBas(iSym),                      &
                          Zero,Da(ipDaa,1),nBas(iSym))
         ipDbb=1+nBD(iSym)
         mAdCMOO=iOff+nBas(iSym)*(nIF(iSym)-nHoles(iSym))
         Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),            &
                          One,CMO(mAdCMOO,2),nBas(iSym),                      &
                                CMO(mAdCMOO,2),nBas(iSym),                      &
                          Zero,Da(ipDbb,2),nBas(iSym))
!
         If (Do_SpinAV) Then
            ipDScc=lOff
            Do j=1,nBas(iSym)
               Do i=1,j
                  ji=j*(j-1)/2+i
                  iDaa=ipDaa-1+ji
                  iDbb=ipDbb-1+ji
                  iDSc=ipDScc+nBas(iSym)*(j-1)+i
                  DSc(iDSc)=Half*(Da(iDaa,1)-Da(iDbb,2))
                  kDSc=ipDScc+nBas(iSym)*(i-1)+j
                  DSc(kDSc)=DSc(iDSc)
               End Do
            End Do
            lOff=lOff+nBas(iSym)**2
         EndIf
!
         Do j=1,nBas(iSym)
            Do i=1,j-1
               ji=j*(j-1)/2+i
               iDaa=ipDaa-1+ji
               Da(iDaa,1)=Two*Da(iDaa,1)
               iDbb=ipDbb-1+ji
               Da(iDbb,2)=Two*Da(iDbb,2)
            End Do
         End Do
         iOff=iOff+nBas(iSym)*nOrb(iSym)
      End Do
!
      Call Cho_X_init(irc,ChFracMem)
      if (irc.ne.0) then
         Call WarningMessage(2,'Start6. Non-zero rc in Cho_X_init.')
         Call Abend
      endif
!----------------------------------------------------------------------*
      Call Get_Fmat_nondyn(Da(:,1),Da(:,2),nBT,.false.)
!----------------------------------------------------------------------*

      Call Cho_X_Final(irc)
      if (irc.ne.0) then
         Call WarningMessage(2,'Start6. Non-zero rc in Cho_X_Final.')
         CALL Abend
      endif
!
      Call mma_deallocate(Da)
      Call mma_deallocate(Corb)
      Call mma_deallocate(Match)
!
      iOff=0
      jOff=0
      Do iSym=1,nSym
         Do iOrb=1,nOcc(iSym,1)
            OccNo(iOrb+iOff,1)=One
         End Do
         Do iOrb=nOcc(iSym,1)+1,nOrb(iSym)
            OccNo(iOrb+iOff,1)=Zero
         End Do
!
         Do iOrb=1,nOcc(iSym,2)
            OccNo(iOrb+iOff,2)=One
         End Do
         Do iOrb=nOcc(iSym,2)+1,nOrb(iSym)
            OccNo(iOrb+iOff,2)=Zero
         End Do
         iOff=iOff+nOrb(iSym)
      End Do
!
      Call mma_allocate(SLT,nBT,Label='SLT')
      isymlbl=1
      iOpt=ibset(ibset(0,sNoOri),sNoNuc)
      Label='Mltpl  0'
      iComp=1
      Call RdOne(irc,iOpt,Label,iComp,SLT,isymlbl)
      If(irc.ne.0) Then
       write(6,*) ' Start6 : error in getting overlap matrix '
       Call Abend
      Endif
      Call s2calc(CMO(1,1),CMO(1,2),SLT,nOcc(1,1),nOcc(1,2),nBas,nOrb,nSym,s2CNO)

      If (.not.Do_SpinAV) Then
      write(6,'(A,f9.6)')'  Initial value of Total Spin, S(S+1): ',s2CNO
      write(6,*)' -------------------------------------------------'
      EndIf
      write(6,*)
!
!----------------------------------------------------------------------*
!  Virtual space must be orthogonal to the occupied space              *
!----------------------------------------------------------------------*
      If (Do_SpinAV) Then
         Do i=1,nSym
            nOcc(i,1)=nOcc(i,1)+nConstr(i)
            nOcc(i,2)=nOcc(i,2)+nConstr(i)
         End Do
      EndIf
!
      Thrd=1.0d-6
      Do i=1,nSym
         nSsh(i)=nOrb(i)-nOcc(i,1)-nFro(i)
         nSsh_ab(i)=nOrb(i)-nOcc(i,2)-nFro(i)
      End Do
!
      Call mma_allocate(SQ,nBB,Label='SQ')
      ltri=1
      lsq=1
      Do iSym=1,nSym
         Call Square(SLT(ltri),SQ(lsq),1,nBas(iSym),nBas(iSym))
         ltri=ltri+nBas(iSym)*(nBas(iSym)+1)/2
         lsq=lsq+nBas(iSym)**2
      End Do
      Call mma_allocate(ID_vir,nnB,Label='ID_vir')
      Call Cho_ov_Loc(irc,Thrd,nSym,nBas,nOcc(1,1),nZero,nZero,nSsh,CMO(1,1),SQ,ID_vir)

      If(irc.ne.0) then
       write(6,*) ' Start6 : error in getting alpha virt MOs '
       Call Abend
      Endif
!
      iOff=1
      Do iSym=1,nSym
         iCMO=iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym,1))
         Call Ortho_Orb(CMO(iCMO,1),SQ,nBas(iSym),nSsh(iSym),2,.true.)
         iOff=iOff+nBas(iSym)*nOrb(iSym)
      End Do
!      Call ChkOrt(1,Whatever) ! silent
!

      Call Cho_ov_Loc(irc,Thrd,nSym,nBas,nOcc(1,2),nZero,nZero,nSsh_ab,CMO(1,2),SQ,iD_vir)

      If(irc.ne.0) then
       write(6,*) ' Start6 : error in getting beta virt MOs '
       Call Abend
      Endif
!
      iOff=1
      Do iSym=1,nSym
         iCMO=iOff+nBas(iSym)*(nFro(iSym)+nOcc(iSym,2))
         Call Ortho_Orb(CMO(iCMO,2),SQ,nBas(iSym),nSsh_ab(iSym),2,.true.)
         iOff=iOff+nBas(iSym)*nOrb(iSym)
      End Do
!      Call ChkOrt(2,Whatever) ! silent
!
      If (Do_SpinAV) Then ! reset # of occupied
         Do i=1,nSym
            nOcc(i,1)=nOcc(i,1)-nConstr(i)
            nOcc(i,2)=nOcc(i,2)-nConstr(i)
         End Do
      EndIf
!
      Call mma_deallocate(ID_vir)
      Call mma_deallocate(SQ)
      Call mma_deallocate(SLT)
!
      Return
      End Subroutine Start6
!***********************************************************************
!                                                                      *
!***********************************************************************
      Subroutine Get_Fmat_nondyn(Dma,Dmb,nBDT,DFTX)
      Use Fock_util_global, only: Deco
      use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
      use SpinAV, only: Do_SpinAV, DSC
      use InfSCF, only: E_nondyn, KSDFT, nBB, nSym, nBas
      use ChoSCF, only: dmpk, nScreen
      use DCSCF, only: Erest_xc
      use Constants, only: Zero, Half, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nBDT
      Real*8  Dma(nBDT), Dmb(nBDT)
      Logical DFTX

      Real*8 dFMat, FactXI
      Integer i, iOff, ipDai, ipDbi, iRC, nDMat
      Real*8, External:: DDot_
      Integer nForb(8,2), nIorb(8,2)
      Real*8, Dimension(:,:), Allocatable:: Dm
      Real*8 E2act(1)
!
      Real*8   Get_ExFac
      External Get_ExFac
      Type (DSBA_Type) FLT(2), KLT(2), POrb(2), PLT(2)
!
      nDMat=2
      Do i=1,nSym
         nForb(i,1)=0
         nForb(i,2)=0
      End Do
      If (DFTX) Then
         FactXI=Get_ExFac(KSDFT)-One ! note this trick
      Else
         FactXI=One
      EndIf
!
      Call Allocate_DT(PLT(1),nBas,nBas,nSym,aCase='TRI')
      Call Allocate_DT(PLT(2),nBas,nBas,nSym,aCase='TRI',Ref=PLT(1)%A0)
      If (DFTX) Then
         PLT(1)%A0(:)=Zero
      Else
         PLT(1)%A0(:)= Dma(:) + Dmb(:)
      EndIf
!
      Call Allocate_DT(POrb(1),nBas,nBas,nSym)
      Call Allocate_DT(POrb(2),nBas,nBas,nSym)

      Call mma_allocate(Dm,nBB,2,Label='Dm')
      Call UnFold(Dma,nBDT,Dm(:,1),nBB,nSym,nBas)
      Call UnFold(Dmb,nBDT,Dm(:,2),nBB,nSym,nBas)
!
      If (Do_SpinAV) Then
         If (.not.DECO) Then
           write(6,*) ' Keywords NODE and SAVE are incompatible. '
           write(6,*) ' NODE will be reset to default. '
           DECO=.true.
         EndIf
         Call daxpy_(NBB,-One,DSc,1,Dm(:,1),1)
         Call daxpy_(NBB, One,DSc,1,Dm(:,2),1)
      EndIf
!
      iOff=0
      Do i=1,nSym
         ipDai=1+iOff
         Call CD_InCore(Dm(ipDai,1),nBas(i),Porb(1)%SB(i)%A2,nBas(i),nIorb(i,1),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Alpha density. Sym= ',i,'   rc= ',irc
            Call RecPrt('Dm',' ',Dm(ipDai,1),nBas(i),nBas(i))
            Call Abend()
         EndIf
         ipDbi=1+iOff
         Call CD_InCore(Dm(ipDbi,2),nBas(i),Porb(2)%SB(i)%A2,nBas(i),nIorb(i,2),1.0d-12,irc)
         If (irc.ne.0) Then
            write(6,*) ' Beta density. Sym= ',i,'   rc= ',irc
            Call RecPrt('Dm',' ',Dm(ipDbi,1),nBas(i),nBas(i))
            Call Abend()
         EndIf
         iOff=iOff+nBas(i)**2
      End Do
!
      Call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI')
      Call Allocate_DT(FLT(2),nBas,nBas,nSym,aCase='TRI')
      FLT(1)%A0(:)=Zero
      FLT(2)%A0(:)=Zero

      Call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI')
      Call Allocate_DT(KLT(2),nBas,nBas,nSym,aCase='TRI')
      KLT(1)%A0(:)=Zero
      KLT(2)%A0(:)=Zero
!
      dFmat=Zero
      Call CHO_LK_SCF(irc,nDMat,FLT,KLT,nForb,nIorb,Porb,PLT,FactXI,nSCReen,dmpk,dFmat)
      if (irc.ne.0) then
         Call WarningMessage(2,'Start6. Non-zero rc in Cho_LK_scf.')
         CALL Abend
      endif
!
      If (Do_SpinAV) Then
         Call UnFold(Dma,nBDT,Dm(1,1),nBB,nSym,nBas)
         Call UnFold(Dmb,nBDT,Dm(1,2),nBB,nSym,nBas)
         Call daxpy_(NBB,-One,DSc,1,Dm(1,1),1)
         Call daxpy_(NBB, One,DSc,1,Dm(1,2),1)
         Call Fold(nSym,nBas,Dm(1,1),Dma)
         Call Fold(nSym,nBas,Dm(1,2),Dmb)
      EndIf
!
      E2act(1) = Half*(ddot_(nBDT,Dma,1,FLT(1)%A0,1)+ddot_(nBDT,Dmb,1,FLT(2)%A0,1))
      Call GADSum(E2act(1),1)
!
      If (DFTX) Then
         Erest_xc=Erest_xc-E2act(1)
      Else
         E_nondyn=E_nondyn-E2act(1)
      EndIf
!
      Call Deallocate_DT(KLT(2))
      Call Deallocate_DT(KLT(1))
      Call Deallocate_DT(FLT(2))
      Call Deallocate_DT(FLT(1))
      Call mma_deallocate(Dm)
      Call Deallocate_DT(POrb(2))
      Call Deallocate_DT(POrb(1))
      Call Deallocate_DT(PLT(2))
      Call Deallocate_DT(PLT(1))
!
      Return
      End Subroutine Get_Fmat_nondyn
!***********************************************************************
!                                                                      *
!***********************************************************************
      Subroutine RdTwoEnrg(LU,E2act)

      Implicit None
      Integer LU
      Real*8  E2act

      Logical Exist
      CHARACTER(LEN=80) LINE

      Call OpnFl('INPORB',LU,Exist)
      If (.Not.Exist) Then
        Write (6,*) 'RdTwoEnrg: INPORB not found!'
        Call Abend()
      End If
      Rewind(LU)
 55   READ(LU,'(A80)',END=888,ERR=888) Line
      If(Line(1:22).ne.'* ACTIVE TWO-EL ENERGY') goto 55
      READ(LU,'(ES19.12)',err=888,end=888) E2act

      Close(LU)
      Return
 888  Call SysWarnFileMsg('RdTwoEnrg','INPORB','Error during reading INPORB\n','Field not there')
      Call Abend()
      End Subroutine RdTwoEnrg
