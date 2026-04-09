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

subroutine ChoSCF_Drv(nBSQT,nD,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,FLT,FLT_ab,nFLT,ExFac,FSQ,nOcc,nOcc_ab)
!
! Thomas Bondo Pedersen, September 2010.
!
! This routine calls the original ChoSCF_Drv routine (now
! ChoSCF_Drv_) in case of Cholesky or full DF.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBSQT, nD, nSym, nBas(nSym), nFLT, nOcc(nSym), nOcc_ab(nSym)
real(kind=wp), intent(in) :: DSQ(*), DLT(*), DSQ_ab(*), ExFac
real(kind=wp), intent(inout) :: DLT_ab(*), FLT(*), FLT_ab(*), FSQ(nBSQT,nD)

!                                                                      *
!***********************************************************************
!                                                                      *

call ChoSCF_Drv_Inner(nD,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,FLT,FLT_ab,nFLT,ExFac,FSQ(:,1),FSQ(:,nD),nOcc,nOcc_ab)

! This is to allow type punning without an explicit interface
contains

subroutine CHOSCF_DRV_Inner(nD,nSym,nBas,W_DSQ,W_DLT,W_DSQ_ab,W_DLT_ab,W_FLT,W_FLT_ab,nFLT,ExFac,W_FSQ,W_FSQ_ab,nOcc,nOcc_ab)

  use Index_Functions, only: iTri
  use InfSCF, only: Algo, Cho_Aufb, CMO, dFKmat, dmpk, nScreen, ReOrd
  use Fock_util_global, only: Deco, Lunit
  use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type, Integer_Pointer
  use SpinAV, only: Do_SpinAV
  use Constants, only: Zero, One, Two, Half
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: nD, nSym, nBas(nSym), nFLT
  real(kind=wp), intent(in) :: W_DSQ(*), W_DLT(*), W_DSQ_ab(*), ExFac
  real(kind=wp), intent(inout) :: W_DLT_ab(*), W_FLT(*), W_FLT_ab(*), W_FSQ(*), W_FSQ_ab(*)
  integer(kind=iwp), target, intent(in) :: nOcc(nSym), nOcc_ab(nSym)
  integer(kind=iwp), parameter :: MaxDs = 3
  integer(kind=iwp) :: i, ikk, iSym, j, ja, k, kj, loff1, MinMem(nSym), n2BSF(8,8), nDen, nForb(8,2), nIorb(8,2), nmat, &
                       nnBSF(8,8), numV, numV1, numV2, rc
  real(kind=wp) :: FactC(MaxDs), FactX(MaxDs), Thr, xFac, YMax
  logical(kind=iwp) :: DoCoulomb(MaxDs), DoExchange(MaxDs), ReOrd_Set = .false.
  character(len=512) :: ww
  type(DSBA_Type) :: Cka(2), DDec(2), DLT(1), DSQ(3), FLT(2), FSQ(3), KLT(2), MSQ(3), Vec(2)
  type(Integer_Pointer) :: pNocc(3)
  integer(kind=iwp), allocatable, target :: nVec(:,:)

  rc = 0
  Lunit(:) = -1
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  if (REORD .and. (.not. ReOrd_Set)) then
    call Cho_X_ReOVec(rc)
    ReOrd_Set = .true.
  end if

  if (nD == 1) then
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    call mma_allocate(nVec,nSym,1,Label='nVec')
    nDen = 1
    DoCoulomb(1) = .true.
    DoExchange(1) = (ExFac /= Zero) ! no SCF-exchange in pure DFT
    FactC(1) = One
    FactX(1) = ExFac ! ExFac used for hybrid functionals

    xFac = ExFac

    call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
    call Allocate_DT(DSQ(1),nBas,nBas,nSym,Ref=W_DSQ)
    call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
    ! trick to use already allocated memory
    call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ)
    call Allocate_DT(FSQ(1),nBas,nBas,nSym,Ref=W_FSQ)

    if (ExFac == Zero) then
      call CHO_FOCK_DFT_RED(rc,DLT,FLT)
      call Error_check(rc)
    else

      if (DECO) then !use decomposed density

        xFac = ExFac*Half

        call set_nnBSF(nSym,nBas,nnBSF,n2BSF)

        call Allocate_DT(Vec(1),nBas,nBas,nSym)
        call Allocate_DT(DDec(1),nBas,nBas,nSym)
        DDec(1)%A0(:) = DSQ(1)%A0(:)

        do i=1,nSym
          if (nBas(i) > 0) then
            Ymax = Zero
            do ja=1,nBas(i)
              Ymax = max(Ymax,DDec(1)%SB(i)%A2(ja,ja))
            end do
            Thr = 1.0e-8_wp*Ymax
            call CD_InCore(DDec(1)%SB(i)%A2,nBas(i),Vec(1)%SB(i)%A2,nBas(i),NumV,Thr,rc)
            call Error_check(rc)
            nVec(i,1) = NumV
            if ((NumV /= nOcc(i)) .and. (.not. Do_SpinAV) .and. (.not. Cho_Aufb)) then
              write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)') &
                'Warning! The number of occupied from the decomposition of the density matrix is ',numV,' in symm. ',i, &
                ';Expected value = ',nOcc(i),'; Max diagonal of the density in symm. ',i,' is equal to ',Ymax
              call WarningMessage(1,trim(ww))
            end if
          else
            nVec(i,1) = 0
          end if
        end do
        call Deallocate_DT(DDec(1))

        pNocc(1)%I1(1:) => nVec(:,1) ! occup. numbers

        call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=Vec(1)%A0)

        FactX(1) = Half*ExFac ! ExFac used for hybrid functionals

      else

        pNocc(1)%I1(1:) => nOcc(:) ! occup. numbers

        call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=CMO(:,1))

      end if

      call CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      select case (ALGO)

        case (1)
          !                                                            *
          !*************************************************************
          !                                                            *
          FactX(1) = Half*ExFac

          if (REORD) then
            call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
          else
            call CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
          end if
          !                                                            *
          !*************************************************************
          !                                                            *
        case (2)
          !                                                            *
          !*************************************************************
          !                                                            *
          if (DECO) then
            FactX(1) = Half*ExFac ! vectors are scaled by construction
          else
            FactX(1) = ExFac ! MOs coeff. are not scaled
          end if

          if (REORD) then
            call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
          else
            call CHO_FMO_RED(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
          end if
          !                                                            *
          !*************************************************************
          !                                                            *
        case (3)

          nIorb(1:nSym,1) = pNocc(1)%I1(1:nSym)
          call Allocate_DT(Cka(1),nIorb(:,1),nBas,nSym)

          do iSym=1,nSym
            if (nBas(iSym)*nIorb(iSym,1) /= 0) then
              do ikk=1,nIorb(iSym,1)
                Cka(1)%SB(iSym)%A2(ikk,:) = MSQ(1)%SB(iSym)%A2(:,ikk)
              end do
            end if
            nForb(iSym,1) = 0
          end do

          call CHO_FSCF(rc,nDen,FLT,nForb,nIorb,Cka(1),DLT,xFac)

          call Deallocate_DT(Cka(1))
          !                                                            *
          !*************************************************************
          !                                                            *
        case (4)

          nForb(1:nSym,1) = 0
          nIorb(1:nSym,1) = pNocc(1)%I1(1:nSym)

          call CHO_LK_SCF(rc,nDen,FLT,KLT,nForb,nIorb,MSQ,DLT,FactX(1),nSCReen,dmpk,dFKmat)

          !                                                            *
          !*************************************************************
          !                                                            *
        case default
          !                                                            *
          !*************************************************************
          !                                                            *
          rc = 99
          write(u6,*) 'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          call QUIT(rc)
          !                                                            *
          !*************************************************************
          !                                                            *
      end select
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Error_check(rc)

      if (DECO) call Deallocate_DT(Vec(1))

      if ((ALGO < 3) .and. (ExFac /= Zero)) call CHO_SUM(rc,nSym,nBas,nD,DoExchange,FLT,FSQ)

    end if

    !-------------------------------------------------------------------
    call GADSum(FLT(1)%A0,nFLT)

    nullify(pNocc(1)%I1)
    call Deallocate_DT(MSQ(1))
    call Deallocate_DT(FSQ(1))
    call Deallocate_DT(DSQ(1))
    call Deallocate_DT(KLT(1))
    call Deallocate_DT(FLT(1))
    call Deallocate_DT(DLT(1))
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  else   !  UHF calculation
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
    call mma_allocate(nVec,nSym,2,Label='nVec')
    nDen = 3
    ! ========== Assign a truth table ==================

    ! Density(1) is Dalpha + Dbeta in a LT storage
    DoCoulomb(1) = .true.
    DoExchange(1) = .false.
    FactC(1) = One

    ! Density(2) is Dalpha in a SQ storage
    DoCoulomb(2) = .false.
    DoExchange(2) = (ExFac /= Zero)

    ! Density(3) is Dbeta in a SQ storage
    DoCoulomb(3) = .false.
    DoExchange(3) = (ExFac /= Zero)

    ! Occupation numbers
    !call get_iarray('nIsh',nOcc,nSym)
    !call get_iarray('nIsh beta',nOcc_ab,nSym)

    ! Compute the total density Dalpha + Dbeta
    W_DLT_ab(1:nFLT) = W_DLT_ab(1:nFLT)+W_DLT(1:nFLT)

    call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_DLT_ab)
    ! alpha density SQ
    call Allocate_DT(DSQ(2),nBas,nBas,nSym,Ref=W_DSQ)
    ! beta  density SQ
    call Allocate_DT(DSQ(3),nBas,nBas,nSym,Ref=W_DSQ_ab)

    ! Coulomb (... Falpha LT)
    call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
    ! (... Fbeta LT)
    call Allocate_DT(FLT(2),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT_ab)

    ! alpha exchange (... Falpha SQ)
    call Allocate_DT(FSQ(2),nBas,nBas,nSym,Ref=W_FSQ)
    ! beta exchange (... Fbeta SQ)
    call Allocate_DT(FSQ(3),nBas,nBas,nSym,Ref=W_FSQ_ab)

    ! trick to use already allocated work
    call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ)
    call Allocate_DT(KLT(2),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ_ab)

    FactX(2) = ExFac ! UHF SQ-density is not scaled
    FactX(3) = ExFac

    if (ExFac == Zero) then
      call CHO_FOCK_DFT_RED(rc,DLT,FLT)
      call Error_check(rc)
    else

      if (DECO) then !use decomposed density

        call set_nnBSF(nSym,nBas,nnBSF,n2BSF)

        call Allocate_DT(Vec(1),nBas,nBas,nSym)
        call Allocate_DT(Vec(2),nBas,nBas,nSym)

        call Allocate_DT(DDec(1),nBas,nBas,nSym)
        call Allocate_DT(DDec(2),nBas,nBas,nSym)
        DDec(1)%A0(:) = DSQ(2)%A0(:)
        DDec(2)%A0(:) = DSQ(3)%A0(:)

        do i=1,nSym
          if (nBas(i) > 0) then
            Ymax = Zero
            do ja=1,nBas(i)
              Ymax = max(Ymax,DDec(1)%SB(i)%A2(ja,ja))
            end do
            Thr = 1.0e-8_wp*Ymax
            call CD_InCore(DDec(1)%SB(i)%A2,nBas(i),Vec(1)%SB(i)%A2,nBas(i),NumV1,Thr,rc)
            call Error_check(rc)
            nVec(i,1) = NumV1
            if ((NumV1 /= nOcc(i)) .and. (.not. Do_SpinAV) .and. (.not. Cho_Aufb)) then
              write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)') &
                'Warning! The number of occupied from the decomposition of the ALPHA dens. matrix is ',numV1,' in symm. ',i, &
                ';Expected value = ',nOcc(i),';Max diagonal of the alpha density in symmetry ',i,' is equal to ',Ymax
              call WarningMessage(1,trim(ww))
            end if

            Ymax = Zero
            do ja=1,nBas(i)
              Ymax = max(Ymax,DDec(2)%SB(i)%A2(ja,ja))
            end do
            Thr = 1.0e-8_wp*Ymax
            call CD_InCore(DDec(2)%SB(i)%A2,nBas(i),Vec(2)%SB(i)%A2,nBas(i),NumV2,Thr,rc)
            call Error_check(rc)
            nVec(i,2) = NumV2
            if ((NumV2 /= nOcc_ab(i)) .and. (.not. Do_SpinAV) .and. (.not. Cho_Aufb)) then
              write(ww,'(a,i6,a,i6,a,i6,a,i6,a,f6.4)') &
                'Warning! The number of occupied from the decomposition of the BETA dens. matrix is ',numV2,' in symm. ',i, &
                ';Expected value = ',nOcc_ab(i),';Max diagonal of the beta density in symmetry ',i,' is equal to ',Ymax
              call WarningMessage(1,trim(ww))
            end if
          else
            nVec(i,1) = 0
            nVec(i,2) = 0
          end if
        end do

        call deallocate_DT(DDec(2))
        call deallocate_DT(DDec(1))

        pNocc(1)%I1(1:) => nVec(:,1) ! dummy
        pNocc(2)%I1(1:) => nVec(:,1) ! alpha occup. numbers
        pNocc(3)%I1(1:) => nVec(:,2) ! beta occup. numbers

        call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=Vec(1)%A0)
        call Allocate_DT(MSQ(2),nBas,nBas,nSym,Ref=Vec(1)%A0)
        call Allocate_DT(MSQ(3),nBas,nBas,nSym,Ref=Vec(2)%A0)
      else
        pNocc(1)%I1(1:) => nOcc(:) ! dummy assignement
        pNocc(2)%I1(1:) => nOcc(:) ! occup. numbers alpha MOs
        pNocc(3)%I1(1:) => nOcc_ab(:) ! occup. numbers beta MOs

        call Allocate_DT(MSQ(1),nBas,nBas,nSym,Ref=CMO(:,1))
        call Allocate_DT(MSQ(2),nBas,nBas,nSym,Ref=CMO(:,1))
        call Allocate_DT(MSQ(3),nBas,nBas,nSym,Ref=CMO(:,2))

      end if

      call CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,loff1)

      select case (ALGO)

        case (1)

          if (REORD) then
            call CHO_FOCKTWO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
          else
            call CHO_FOCKTWO_RED(rc,nBas,nDen,DoCoulomb,DoExchange,FactC,FactX,DLT,DSQ,FLT,FSQ,pNocc,MinMem)
          end if

        case (2)

          if (REORD) then
            call CHO_FTWO_MO(rc,nSym,nBas,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
          else
            call CHO_FMO_red(rc,nDen,DoCoulomb,DoExchange,lOff1,FactC,FactX,DLT,DSQ,FLT,FSQ,MinMem,MSQ,pNocc)
          end if

        case (3)

          nIorb(1:nSym,1) = pNocc(2)%I1(1:nSym)
          nIorb(1:nSym,2) = pNocc(3)%I1(1:nSym)
          call Allocate_DT(Cka(1),nIorb(:,1),nBas,nSym)
          call Allocate_DT(Cka(2),nIorb(:,2),nBas,nSym)

          do iSym=1,nSym
            if (nBas(iSym)*nIorb(iSym,1) /= 0) then
              do ikk=1,nIorb(iSym,1)
                Cka(1)%SB(iSym)%A2(ikk,:) = MSQ(2)%SB(iSym)%A2(:,ikk)
              end do
            end if
            if (nBas(iSym)*nIorb(iSym,2) /= 0) then
              do ikk=1,nIorb(iSym,2)
                Cka(2)%SB(iSym)%A2(ikk,:) = MSQ(3)%SB(iSym)%A2(:,ikk)
              end do
            end if
            nForb(iSym,1) = 0
            nForb(iSym,2) = 0
          end do

          nMat = 2  ! alpha and beta Fock matrices

          call CHO_FSCF(rc,nMat,FLT,nForb,nIorb,Cka,DLT,ExFac)

          call Deallocate_DT(Cka(2))
          call Deallocate_DT(Cka(1))

        case (4)

          nMat = 2  ! alpha and beta Fock matrices

          nForb(1:nSym,:) = 0
          nIorb(1:nSym,1) = pNocc(2)%I1(1:nSym)
          nIorb(1:nSym,2) = pNocc(3)%I1(1:nSym)

          call CHO_LK_SCF(rc,nMat,FLT,KLT,nForb,nIorb,MSQ(2:3),DLT,FactX(2),nSCReen,dmpk,dFKmat)

        case default

          rc = 99
          write(u6,*) 'Illegal Input. Specified Cholesky Algorithm= ',ALGO
          call QUIT(rc)
      end select

      call Error_check(rc)

      if (DECO) call deallocate_DT(Vec(2))
      if (DECO) call deallocate_DT(Vec(1))

    end if

    ! To get the Fbeta in LT storage ----

    if ((ALGO < 3) .or. (ExFac == Zero)) FLT(2)%A0(:) = FLT(1)%A0(:)

    ! Accumulates Coulomb and Exchange contributions
    if ((ALGO < 3) .and. (ExFac /= Zero)) call CHO_SUM(rc,nSym,nBas,nD,DoExchange,FLT,FSQ)

    !-------------------------------------------------------------------
    call GADSum(FLT(1)%A0,nFLT)
    call GADSum(FLT(2)%A0,nFLT)

    ! --- Restore the Beta-density matrix ---
    ! Copy the lower triangular of DSQ(3))
    ! and pack the off-diagonal elements

    do isym=1,nsym
      do j=1,nBas(isym)
        do k=j,nBas(isym)
          kj = iTri(k,j)
          if (j == k) then
            DLT(1)%SB(iSym)%A1(kj) = DSQ(3)%SB(iSym)%A2(k,j)
          else ! packing of the matrix
            DLT(1)%SB(iSym)%A1(kj) = Two*DSQ(3)%SB(iSym)%A2(k,j)
          end if
        end do
      end do
    end do

    nullify(pNocc(1)%I1,pNocc(2)%I1,pNocc(3)%I1)
    call Deallocate_DT(MSQ(3))
    call Deallocate_DT(MSQ(2))
    call Deallocate_DT(MSQ(1))
    call Deallocate_DT(DSQ(3))
    call Deallocate_DT(DSQ(2))
    call Deallocate_DT(FSQ(3))
    call Deallocate_DT(FSQ(2))
    call Deallocate_DT(KLT(2))
    call Deallocate_DT(KLT(1))
    call Deallocate_DT(FLT(2))
    call Deallocate_DT(FLT(1))
    call Deallocate_DT(DLT(1))
    !                                                                  *
    !*******************************************************************
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *

  call Error_check(rc)

  call mma_deallocate(nVec)

  return

end subroutine ChoSCF_Drv_Inner

subroutine Error_check(rc)

  integer(kind=iwp), intent(in) :: rc

  if (rc /= 0) then
    write(u6,*) 'CHOSCF_DRV. Non-zero return code.'
    call QUIT(rc)
  end if

end subroutine Error_check

end subroutine ChoSCF_Drv
