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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!               2023, Roland Lindh                                     *
!***********************************************************************
      Subroutine ChoSCF_Drv(nBSQT,nD,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,FLT,FLT_ab,nFLT,ExFac,FSQ,nOcc,nOcc_ab)
!
!     Thomas Bondo Pedersen, September 2010.
!
!     This routine calls the original ChoSCF_Drv routine (now
!     ChoSCF_Drv_) in case of Cholesky or full DF.
!
      Implicit None

      Integer nBSQT, nD, nSym, nFLT
      Integer nBas(nSym), nOcc(nSym), nOcc_ab(nSym)
      Real*8  DSQ(*), DLT(*)
      Real*8  DSQ_ab(*), DLT_ab(*)
      Real*8  FLT(*), FLT_ab(*)
      Real*8  FSQ(nBSQT,nD)
      Real*8  ExFac
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
      SUBROUTINE CHOSCF_DRV_Internal(nD,nSym,nBas,W_DSQ,W_DLT,         &
                                     W_DSQ_ab,W_DLT_ab,W_FLT,          &
                                     W_FLT_ab,nFLT,ExFac,              &
                                     W_FSQ,W_FSQ_ab,                   &
                                     nOcc,nOcc_ab)
      Integer nD, nSym
      Integer nBas(nSym)
      Real*8 W_FLT(*),W_FLT_ab(*)
      Real*8 W_FSQ(*),W_FSQ_ab(*)
      Real*8 W_DSQ(*),W_DSQ_ab(*)
      Real*8 W_DLT(*),W_DLT_ab(*)
      Real*8 ExFac
      Integer, Target :: nOcc(nSym),nOcc_ab(nSym)
      End SUBROUTINE CHOSCF_DRV_Internal

      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *

#define ChoCorr

      Call ChoSCF_Drv_Internal(nD,nSym,nBas,DSQ,DLT,                 &
                               DSQ_ab,DLT_ab,FLT,                    &
                               FLT_ab,nFLT,ExFac,                    &
                               FSQ(:,1),FSQ(:,2),                    &
                               nOcc,nOcc_ab)

      End Subroutine ChoSCF_Drv



      SUBROUTINE CHOSCF_DRV_Internal(nD,nSym,nBas,W_DSQ,W_DLT,          &
                                     W_DSQ_ab,W_DLT_ab,W_FLT,           &
                                     W_FLT_ab,nFLT,ExFac,               &
                                     W_FSQ,W_FSQ_ab,                    &
                                     nOcc,nOcc_ab)

      use Scf_Arrays, only: CMO
      Use Fock_util_global, only: Deco, Lunit
      use Data_Structures, only: Allocate_DT, Deallocate_DT
      use Data_Structures, only: DSBA_Type, Integer_Pointer
      use stdalloc, only: mma_allocate, mma_deallocate
      use SpinAV, only: Do_SpinAV
      use ChoSCF, only: Algo, dFKmat, dmpk, nScreen, ReOrd
      use Constants, only: Zero, Half, One, Two
      use ChoAuf, only: Cho_Aufb

      Implicit None

      Integer nD, nSym
      Integer nBas(nSym)
      Real*8 W_DSQ(*),W_DSQ_ab(*)
      Real*8 W_DLT(*),W_DLT_ab(*)
      Real*8 W_FLT(*),W_FLT_ab(*)
      Integer nFLT
      Real*8 ExFac
      Real*8 W_FSQ(*),W_FSQ_ab(*)
      Integer, Target :: nOcc(nSym),nOcc_ab(nSym)

      Integer MinMem(nSym),rc
      Integer, Parameter :: MaxDs = 3
      Logical DoCoulomb(MaxDs),DoExchange(MaxDs)
      Real*8 FactC(MaxDs),FactX(MaxDs)
      Integer nnBSF(8,8),n2BSF(8,8)
      Integer nForb(8,2),nIorb(8,2)
      character(LEN=512)  ww
      Integer :: iTri, i, j
      Integer :: ikk, iSym, ja, k, kj, loff1, nmat, numV, numV1, numV2
      Integer :: nDen
      Real*8 :: Thr, xFac, YMax
      Logical :: ReOrd_Set=.False.

      Type (Integer_Pointer) :: pNocc(3)

      Integer, Allocatable, Target:: nVec(:,:)

      Type (DSBA_Type) Cka(2), FLT(2), KLT(2), MSQ(3), DLT(1), FSQ(3), DSQ(3), DDec(2), Vec(2)
!
!  **************************************************
        iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
!  **************************************************

      rc=0
      Lunit(:) = -1
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      IF (REORD.and..NOT.ReOrd_Set) THEN
         Call Cho_X_ReOVec(rc)
         ReOrd_Set=.TRUE.
      END If

      IF(nD==1) THEN
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
         Call mma_allocate(nVec,nSym,1,Label='nVec')
         nDen = 1
         DoCoulomb(1)  = .true.
         DoExchange(1) = ExFac.ne.Zero ! no SCF-exchange in pure DFT
         FactC(1)      = One
         FactX(1)      = One*ExFac ! ExFac used for hybrid functionals

         xFac = ExFac

         Call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
         Call Allocate_DT(DSQ(1),nBas,nBas,nSym,Ref=W_DSQ)
         Call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
         ! trick to use already allocated memory
         Call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ)
         Call Allocate_DT(FSQ(1),nBas,nBas,nSym,            Ref=W_FSQ)

         If (ExFac.eq.Zero) Then
            CALL CHO_FOCK_DFT_RED(rc,DLT,FLT)
            If (rc.ne.0) Go To 999
            goto 997
          EndIf

      IF (DECO) THEN !use decomposed density

       xFac = ExFac*half

       CALL set_nnBSF(nSym,nBas,nnBSF,n2BSF)

       Call Allocate_DT(Vec(1),nBas,nBas,nSym)
       Call Allocate_DT(DDec(1),nBas,nBas,nSym)
       DDec(1)%A0(:) = DSQ(1)%A0(:)

       Do i=1,nSym
          if(nBas(i).gt.0)then
            Ymax=Zero
            do ja=1,nBas(i)
               Ymax=Max(Ymax,DDec(1)%SB(i)%A2(ja,ja))
            end do
            Thr = 1.0d-8*Ymax
            CALL CD_InCore(DDec(1)%SB(i)%A2,nBas(i),Vec(1)%SB(i)%A2,nBas(i),NumV,Thr,rc)
            If (rc.ne.0) GOTO 999
            nVec(i,1) = NumV
            if ( NumV .ne. nOcc(i) .and. .not.Do_SpinAV .and. .not.Cho_Aufb ) then
            write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')                   &
              'Warning! The number of occupied from the '//            &
              'decomposition of the density matrix is ',numV,          &
              ' in symm. ',i, '; Expected value = ',nOcc(i),           &
              '; Max diagonal of the density in symm. ',i,             &
              ' is equal to ',Ymax
            call WarningMessage(1,ww)
            endif
          else
            nVec(i,1) = 0
          endif
       End Do
       Call Deallocate_DT(DDec(1))

       pNocc(1)%I1(1:) => nVec(1:,1) ! occup. numbers

       Call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=Vec(1)%A0)

       FactX(1) = half*ExFac ! ExFac used for hybrid functionals

      ELSE

         pNocc(1)%I1(1:) => nOcc(1:) ! occup. numbers

         Call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=CMO(:,1))

      ENDIF

      Call CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)

      Select Case(ALGO)
!                                                                      *
!***********************************************************************
!                                                                      *
        Case (1)
!                                                                      *
!***********************************************************************
!                                                                      *
        FactX(1)=Half*ExFac

        if (REORD) then
           Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
        else
           CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
        end if
!                                                                      *
!***********************************************************************
!                                                                      *
        Case (2)
!                                                                      *
!***********************************************************************
!                                                                      *
        if (DECO) then
           FactX(1) = Half*ExFac ! vectors are scaled by construction
        else
           FactX(1) = One*ExFac ! MOs coeff. are not scaled
        end if

        if (REORD)then
           Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
        else
           CALL CHO_FMO_RED(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
        endif
!                                                                      *
!***********************************************************************
!                                                                      *
        Case (3)

          Do iSym=1,nSym
             nIorb(iSym,1) = pNocc(1)%I1(iSym)
          End Do
          Call Allocate_DT(Cka(1),nIorb(:,1),nBas,nSym)

          Do iSym=1,nSym
           If (nBas(iSym)*nIorb(iSym,1).ne.0) Then
             do ikk=1,nIorb(iSym,1)
                Cka(1)%SB(iSym)%A2(ikk,:) = MSQ(1)%SB(iSym)%A2(:,ikk)
             end do
           EndIf
           nForb(iSym,1) = 0
          End Do

          CALL CHO_FSCF(rc,nDen,FLT,nForb,nIorb,Cka(1),DLT,xFac)

          Call Deallocate_DT(Cka(1))
!                                                                      *
!***********************************************************************
!                                                                      *
        Case (4)

             Do iSym=1,nSym
                nForb(iSym,1) = 0
                nIorb(iSym,1) = pNocc(1)%I1(iSym)
             End Do

             CALL CHO_LK_SCF(rc,nDen,FLT,KLT,nForb,nIorb,MSQ,DLT,FactX(1),nSCReen,dmpk,dFKmat)

!                                                                      *
!***********************************************************************
!                                                                      *
        Case default
!                                                                      *
!***********************************************************************
!                                                                      *
          rc=99
          write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          CALL QUIT(rc)
!                                                                      *
!***********************************************************************
!                                                                      *
      End Select
!                                                                      *
!***********************************************************************
!                                                                      *
      If (rc.ne.0) GOTO 999

      IF (DECO) CALL Deallocate_DT(Vec(1))

      If (ALGO.lt.3.and.ExFac.ne.Zero) Then
         CALL CHO_SUM(rc,nSym,nBas,nD,DoExchange,FLT,FSQ)

      EndIf
!----------------------------------------------------
 997  Continue
      Call GADSum(FLT(1)%A0,nFLT)

      pNocc(1)%I1=>Null()
      Call Deallocate_DT(MSQ(1))
      Call Deallocate_DT(FSQ(1))
      Call Deallocate_DT(DSQ(1))
      Call Deallocate_DT(KLT(1))
      Call Deallocate_DT(FLT(1))
      Call Deallocate_DT(DLT(1))
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      ELSE   !  UHF calculation
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
         Call mma_allocate(nVec,nSym,2,Label='nVec')
         nDen = 3
! ========== Assign a truth table ==================

! --- Density(1) is Dalpha + Dbeta in a LT storage
         DoCoulomb(1)  = .true.
         DoExchange(1) = .false.
         FactC(1)      = One

! --- Density(2) is Dalpha in a SQ storage
         DoCoulomb(2)  = .false.
         DoExchange(2) = ExFac.ne.Zero

! --- Density(3) is Dbeta in a SQ storage
         DoCoulomb(3)  = .false.
         DoExchange(3) = ExFac.ne.Zero

! --- Occupation numbers
!      call get_iarray('nIsh',nOcc,nSym)
!      call get_iarray('nIsh beta',nOcc_ab,nSym)


! Compute the total density Dalpha + Dbeta
      CALL DAXPY_(nFLT,One,W_DLT(1),1,W_DLT_ab(1),1)

      Call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_DLT_ab)
      ! alpha density SQ
      Call Allocate_DT(DSQ(2),nBas,nBas,nSym,Ref=W_DSQ)
      ! beta  density SQ
      Call Allocate_DT(DSQ(3),nBas,nBas,nSym,Ref=W_DSQ_ab)

      ! Coulomb (... Falpha LT)
      Call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
      ! (... Fbeta LT)
      Call Allocate_DT(FLT(2),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT_ab)

      ! alpha exchange (... Falpha SQ)
      Call Allocate_DT(FSQ(2),nBas,nBas,nSym,Ref=W_FSQ)
      ! beta exchange (... Fbeta SQ)
      Call Allocate_DT(FSQ(3),nBas,nBas,nSym,Ref=W_FSQ_ab)

      ! trick to use already allocated work
      Call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ)
      Call Allocate_DT(KLT(2),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ_ab)

      FactX(2) = One*ExFac ! UHF SQ-density is not scaled
      FactX(3) = One*ExFac

      If (ExFac.eq.Zero) Then
         CALL CHO_FOCK_DFT_RED(rc,DLT,FLT)
         If (rc.ne.0) Go To 999
         goto 998
      EndIf

      IF (DECO) THEN !use decomposed density

       CALL set_nnBSF(nSym,nBas,nnBSF,n2BSF)

       Call Allocate_DT(Vec(1),nBas,nBas,nSym)
       Call Allocate_DT(Vec(2),nBas,nBas,nSym)

       Call Allocate_DT(DDec(1),nBas,nBas,nSym)
       Call Allocate_DT(DDec(2),nBas,nBas,nSym)
       DDec(1)%A0(:)=DSQ(2)%A0(:)
       DDec(2)%A0(:)=DSQ(3)%A0(:)

       Do i=1,nSym
          if (nBas(i).gt.0)then
             Ymax=Zero
             do ja=1,nBas(i)
                Ymax=Max(Ymax,DDec(1)%SB(i)%A2(ja,ja))
             end do
             Thr = 1.0d-8*Ymax
             CALL CD_InCore(DDec(1)%SB(i)%A2,nBas(i),Vec(1)%SB(i)%A2,nBas(i),NumV1,Thr,rc)
             If (rc.ne.0) GOTO 999
             nVec(i,1) = NumV1
             if ( NumV1 .ne. nOcc(i) .and. .not.Do_SpinAV .and. .not.Cho_Aufb) then
                write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')                 &
                 'Warning! The number of occupied from the '//           &
                 'decomposition of the ALPHA dens. matrix is ',numV1,    &
                 ' in symm. ',i,                                         &
                 ';Expected value = ',nOcc(i),                           &
                 ';Max diagonal of the alpha density in symmetry ',      &
                 i,' is equal to ',Ymax
                 call WarningMessage(1,ww)
             endif

             Ymax=Zero
             do ja=1,nBas(i)
                Ymax=Max(Ymax,DDec(2)%SB(i)%A2(ja,ja))
             end do
             Thr = 1.0d-8*Ymax
             CALL CD_InCore(DDec(2)%SB(i)%A2,nBas(i),Vec(2)%SB(i)%A2,nBas(i),NumV2,Thr,rc)
             If (rc.ne.0) GOTO 999
             nVec(i,2) = NumV2
             if ( NumV2 .ne. nOcc_ab(i) .and. .not.Do_SpinAV .and. .not.Cho_Aufb) then
                write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)')                 &
                 'Warning! The number of occupied from the '//           &
                 'decomposition of the BETA dens. matrix is ',numV2,     &
                 ' in symm. ',i,                                         &
                 ';Expected value = ',nOcc_ab(i),                        &
                 ';Max diagonal of the beta density in symmetry ',       &
                 i,' is equal to ',Ymax
                call WarningMessage(1,ww)
             endif
          else
            nVec(i,1) = 0
            nVec(i,2) = 0
          endif
       End Do

       Call deallocate_DT(DDec(2))
       Call deallocate_DT(DDec(1))

       pNocc(1)%I1(1:) => nVec(1:,1) ! dummy
       pNocc(2)%I1(1:) => nVec(1:,1) ! alpha occup. numbers
       pNocc(3)%I1(1:) => nVec(1:,2) ! beta occup. numbers

       Call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=Vec(1)%A0)
       Call Allocate_DT(MSQ(2),nBas,nBas,nSym,Ref=Vec(1)%A0)
       Call Allocate_DT(MSQ(3),nBas,nBas,nSym,Ref=Vec(2)%A0)
      ELSE
       pNocc(1)%I1(1:) => nOcc(1:) ! dummy assignement
       pNocc(2)%I1(1:) => nOcc(1:) ! occup. numbers alpha MOs
       pNocc(3)%I1(1:) => nOcc_ab(1:) ! occup. numbers beta MOs

       Call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=CMO(:,1))
       Call Allocate_DT(MSQ(2),nBas,nBas,nSym,Ref=CMO(:,1))
       Call Allocate_DT(MSQ(3),nBas,nBas,nSym,Ref=CMO(:,2))

      ENDIF

      Call CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)


      Select Case(ALGO)

         Case(1)

         If (REORD) Then
            Call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
         Else
            CALL CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
         End If

         Case(2)

         if (REORD)then

          Call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)

         else

          CALL CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)

         endif

         Case(3)

          Do iSym=1,nSym
             nIorb(iSym,1) = pNocc(2)%I1(iSym)
             nIorb(iSym,2) = pNocc(3)%I1(iSym)
          End Do
          Call Allocate_DT(Cka(1),nIorb(:,1),nBas,nSym)
          Call Allocate_DT(Cka(2),nIorb(:,2),nBas,nSym)

          Do iSym=1,nSym
           If (nBas(iSym)*nIorb(iSym,1).ne.0) Then
             do ikk=1,nIorb(iSym,1)
                Cka(1)%SB(iSym)%A2(ikk,:) = MSQ(2)%SB(iSym)%A2(:,ikk)
             end do
           EndIf
           If (nBas(iSym)*nIorb(iSym,2).ne.0) Then
             do ikk=1,nIorb(iSym,2)
                Cka(2)%SB(iSym)%A2(ikk,:) = MSQ(3)%SB(iSym)%A2(:,ikk)
             end do
           EndIf
           nForb(iSym,1) = 0
           nForb(iSym,2) = 0
          End Do

          nMat=2  ! alpha and beta Fock matrices

          CALL CHO_FSCF(rc,nMat,FLT,nForb,nIorb,Cka,DLT,ExFac)


          Call Deallocate_DT(Cka(2))
          Call Deallocate_DT(Cka(1))

         Case(4)

             nMat=2  ! alpha and beta Fock matrices

             Do iSym=1,nSym
                nForb(iSym,1) = 0
                nForb(iSym,2) = 0
                nIorb(iSym,1) = pNocc(2)%I1(iSym)
                nIorb(iSym,2) = pNocc(3)%I1(iSym)
             End Do

             CALL CHO_LK_SCF(rc,nMat,FLT,KLT,nForb,nIorb,MSQ(2:3),DLT,FactX(2),nSCReen,dmpk,dFKmat)

          Case Default

          rc=99
          write(6,*)'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          CALL QUIT(rc)
      End  Select

      If (rc.ne.0) GOTO 999

      IF(DECO) CALL deallocate_DT(Vec(2))
      IF(DECO) CALL deallocate_DT(Vec(1))

998   Continue

! --- To get the Fbeta in LT storage ----

      If (ALGO.lt.3 .or. ExFac.eq.Zero) then

        FLT(2)%A0(:) = FLT(1)%A0(:)

      EndIf


! --- Accumulates Coulomb and Exchange contributions
      If (ALGO.lt.3.and.ExFac.ne.Zero) then
         CALL CHO_SUM(rc,nSym,nBas,nD,DoExchange,FLT,FSQ)
      Endif

!----------------------------------------------------
      Call GADSum(FLT(1)%A0,nFLT)
      Call GADSum(FLT(2)%A0,nFLT)


! --- Restore the Beta-density matrix ----
! --- Copy the lower triangular of DSQ(3))
! --- and pack the off-diagonal elements

        do isym=1,nsym
           do j=1,nBas(isym)
              do k=j,nBas(isym)
               kj =  iTri(k,j)
               if (j.eq.k) then
                DLT(1)%SB(iSym)%A1(kj) =     DSQ(3)%SB(iSym)%A2(k,j)
               else ! packing of the matrix
                DLT(1)%SB(iSym)%A1(kj) = Two*DSQ(3)%SB(iSym)%A2(k,j)
               endif
              end do
           end do
        end do

         pNocc(1)%I1 => Null()
         pNocc(2)%I1 => Null()
         pNocc(3)%I1 => Null()
         Call Deallocate_DT(MSQ(3))
         Call Deallocate_DT(MSQ(2))
         Call Deallocate_DT(MSQ(1))
         Call Deallocate_DT(DSQ(3))
         Call Deallocate_DT(DSQ(2))
         Call Deallocate_DT(FSQ(3))
         Call Deallocate_DT(FSQ(2))
         Call Deallocate_DT(KLT(2))
         Call Deallocate_DT(KLT(1))
         Call Deallocate_DT(FLT(2))
         Call Deallocate_DT(FLT(1))
         Call Deallocate_DT(DLT(1))
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
      ENDIF
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *

999   If (rc.ne.0) then
         write(6,*)'CHOSCF_DRV. Non-zero return code.'
         CALL QUIT(rc)
      endif

      Call mma_deallocate(nVec)

      Return
      End Subroutine ChoSCF_Drv_Internal
