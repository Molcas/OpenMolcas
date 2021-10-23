!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CHO_CAS_DRV(rc,W_CMO,DI,FI,DA1,FA,W_PWXY,TraOnly)

use Data_Structures, only: DSBA_Type, Allocate_DSBA, Deallocate_DSBA

implicit real*8(a-h,o-z)
integer rc
real*8 W_PWXY(*)
real*8 DA1(*), DI(*), FI(*), FA(*), W_CMO(*)
integer nForb(8), nIorb(8), nAorb(8), nChM(8), nChI(8)
logical TraOnly
#include "real.fh"
#include "chotodo.fh"
#include "chlcas.fh"
#include "cholk.fh"
character(LEN=11), parameter :: SECNAM = 'CHO_CAS_DRV'
#include "rasdim.fh"
#include "wadr.fh"
#include "general.fh"
#include "rasscf.fh"
#include "stdalloc.fh"
type(DSBA_Type) CVa(2), POrb(3), Ddec, ChoIn, CMO, DLT(2), FLT(2), MSQ, FLT_MO(2)
real*8, allocatable :: Tmp1(:), Tmp2(:)
!                                                                      *
!***********************************************************************
!                                                                      *
interface

  subroutine DGEMM_(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

    character*1 TRANSA, TRANSB
    integer M, N, K, LDA, LDB, LDC
    real*8 ALPHA, BETA
    real*8 A(LDA,*), B(LDB,*), C(LDC,*)
  end subroutine DGEMM_

  subroutine MXMT(A,ICA,IRA,B,ICB,IRB,C,NROW,NSUM)
    integer ICA, IRA, ICB, IRB, NROW, NSUM
    real*8 A(*), B(*), C(*)
  end subroutine MXMT

end interface

!                                                                      *
!***********************************************************************
!                                                                      *
rc = 0

call Allocate_DSBA(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=FI)
call Allocate_DSBA(FLT(2),nBas,nBas,nSym,aCase='TRI',Ref=FA)
call Allocate_DSBA(CMO,nBas,nBas,nSym,Ref=W_CMO)

!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (TraOnly) then
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Let us construct a second pointer structure based on nOrb

  call Allocate_DSBA(FLT_MO(1),nOrb,nOrb,nSym,aCase='TRI',Ref=FI)
  call Allocate_DSBA(FLT_MO(2),nOrb,nOrb,nSym,aCase='TRI',Ref=FA)
  !
  ! It only performs the MO transformation of FI and FA
  ! ---------------------------------------------------
  !
  ! transform FI/FA from AO to MO basis  (LT-storage)
  do i=1,2
    do iSym=1,nSym
      iBas = nBas(iSym)
      iOrb = nOrb(iSym)
      iFro = nFro(iSym)
      call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
      call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp1')
      call Square(FLT(i)%SB(iSym)%A1,Tmp1,1,iBas,iBas)
      call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Tmp1,iBas,CMO%SB(iSym)%A1(1+iFro*iBas:),iBas,0.0d0,Tmp2,iBas)
      call MXMT(Tmp2,iBas,1,CMO%SB(iSym)%A1(1+iFro*iBas:),1,iBas,FLT_MO(i)%SB(iSym)%A1,iOrb,iBas)
      call mma_deallocate(Tmp2)
      call mma_deallocate(Tmp1)
    end do
  end do

  call deallocate_DSBA(FLT_MO(1))
  call deallocate_DSBA(FLT_MO(2))
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! It only computes FI and FA  in AO-basis and returns
  ! the active integrals (tw|xy)
  ! If specified in input, the routine also computes
  ! the auxiliary Q-matrix stored as Q(av), where a is an AO index
  ! and v refers to the active orbitals only

  nForb(:) = nFro(:)
  nIorb(:) = nIsh(:)
  nAorb(:) = nAsh(:)

  ! Build the packed densities from the Squared ones
  call Allocate_DSBA(DLT(1),nBas,nBas,nSym,aCase='TRI')
  call Allocate_DSBA(DLT(2),nBas,nBas,nSym,aCase='TRI')

  call Fold(nSym,nBas,DI,DLT(1)%A0)
  call Fold(nSym,nBas,DA1,DLT(2)%A0)

  FactXI = -1.0d0

  !AMS - should this be set differently for ExFac.ne.1?
  !      FactXI = 0-ExFac

  if (Deco) then

    FactXI = -0.5d0

    !AMS - should this be set differently for ExFac.ne.1?
    !      FactXI = 0-(ExFac*0.5)

    ! --- decompose the Inactive density on request
    call Allocate_DSBA(ChoIn,nBas,nBas,nSym)
    call Allocate_DSBA(DDec,nBas,nBas,nSym)
    DDec%A0(1:NTot2) = DI(1:NTot2)

    call Allocate_DSBA(MSQ,nBas,nBas,nSym,Ref=ChoIn%A0)

    Thr = 1.0d-12
    incs = 0
    do i=1,nSym
      if ((nForb(i)+nIorb(i)) > 0) then
        call CD_InCore(DDec%SB(i)%A2,nBas(i),ChoIn%SB(i)%A2,nBas(i),NumV,Thr,rc)
        if (rc /= 0) then
          write(6,*) SECNAM//': ill-defined dens decomp for Inact'
          write(6,*) 'rc value produced = ',rc
          call abend()
        end if
        nChI(i) = NumV
        if (NumV /= nIsh(i)+nForb(i)) then
          write(6,*) 'Warning! The number of occupied from the decomposition of the Inactive density matrix is ',numV,' in symm. ',i
          write(6,*) 'Expected value = ',nIsh(i)+nForb(i)
          incs = incs+1
          Ymax = 0.0d0
          do ja=1,nBas(i)
            Ymax = max(Ymax,DDec%SB(i)%A2(ja,ja))
          end do
          write(6,*) 'Max diagonal of the density in symm. ',i,' is equal to ',Ymax
        end if
      else
        nChI(i) = 0
      end if
    end do

    if ((incs > 0) .and. DoLocK) then
      dmpk_old = dmpk
      dmpk = 1.0d-2*dmpk
      write(6,*) 'LK-damping decreased from ',dmpk_old,' to ',dmpk
    end if

    call Deallocate_DSBA(DDEc)

    ! --- to get the right input arguments for CHO_FCAS_AO and CHO_FMCSCF
    if (.not. DoLocK) then
      nForb(:) = 0
      nIorb(:) = nChI(:)
    end if

  else

    call Allocate_DSBA(MSQ,nBas,nBas,nSym,Ref=W_CMO)

    nChI(:) = nForb(:)+nIorb(:)

  end if

  ! Reordering of the MOs coefficients to fit cholesky needs

  if (.not. DoLocK) then

    call Allocate_DSBA(POrb(1),nChI,nBas,nSym)
    call Allocate_DSBA(POrb(3),nAOrb,nBas,nSym)

    do iSym=1,nSym

      do ikk=1,nChI(iSym)
        POrb(1)%SB(iSym)%A2(ikk,:) = MSQ%SB(iSym)%A2(:,ikk)
      end do

      do ikk=1,nAorb(iSym)
        jkk = nForb(iSym)+nIorb(iSym)+ikk
        POrb(3)%SB(iSym)%A2(ikk,:) = CMO%SB(iSym)%A2(:,jkk)
      end do

    end do

  else

    ! *** Only the active orbitals MO coeff need reordering
    call Allocate_DSBA(CVa(1),nAorb,nBas,nSym)

    do iSym=1,nSym
      do ikk=1,nAorb(iSym)
        jkk = nForb(iSym)+nIorb(iSym)+ikk
        CVa(1)%SB(iSym)%A2(ikk,:) = CMO%SB(iSym)%A2(:,jkk)
      end do
    end do

  end if

  if (DoActive) then
    ! Decompose the active density  -----------------------------

#   ifdef _DEBUGPRINT_
    do i=1,nSym
      call CD_TESTER(rc,DLT(2)%SB(i)%A1,nBas(i),.true.)
      write(6,*) 'DALT for sym=',i
      call TRIPRT('DALT',' ',DLT(2)%SB(i)%A1,nBas(i))
    end do
#   endif

    call Allocate_DSBA(CVa(2),nBas,nBas,nSym)
    call Allocate_DSBA(DDec,nBas,nBas,nSym)
    DDec%A0(1:NTot2) = DA1(1:NTot2)

    Thr = 1.0d-12
    do i=1,nSym
      if (nAorb(i) > 0) then
        ! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
        !                 This will avoid warnings for negative-definit
        call CD_InCore(DDec%SB(i)%A2,nBas(i),CVa(2)%SB(i)%A2,nBas(i),NumV,Thr,rc)
        if (rc /= 0) then
          write(6,*) SECNAM//': ill-defined dens decomp for active'
          write(6,*) 'rc value produced = ',rc
          call abend()
        end if
        nChM(i) = NumV
      else
        nChM(i) = 0
      end if
    end do

    call Deallocate_DSBA(DDec)

  else

    ! Dummy allocation
    call Allocate_DSBA(CVa(2),[1],[1],1)
    nChM(:) = 0

  end if

  if ((.not. DoLocK) .and. DoActive) then

    ! reorder "Cholesky MOs" to Cva storage

    call Allocate_DSBA(POrb(2),nChM,nBas,nSym)
    do iSym=1,nSym
      if (nBas(iSym)*nChM(iSym) /= 0) then
        do ikk=1,nChM(iSym)
          POrb(2)%SB(iSym)%A2(ikk,:) = CVa(2)%SB(iSym)%A2(:,ikk)
        end do
      end if
    end do

  else

    call Allocate_DSBA(POrb(2),[1],[1],1)

  end if
  ! --------------------------------------------------------------------
  FLT(1)%A0(:) = Zero
  FLT(2)%A0(:) = Zero

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  if ((ALGO == 1) .and. (.not. DoLocK)) then

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

    call CHO_FMCSCF(rc,FLT,nForb,nIorb,nAorb,FactXI,DLT,DoActive,POrb,nChM,W_PWXY,CMO,ExFac)

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  else if ((ALGO == 1) .and. DoLocK) then

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

    call CHO_LK_CASSCF(DLT,FLT,MSQ,W_PWXY,FactXI,nChI,nAorb,nChM,CVa,DoActive,nScreen,dmpK,abs(CBLBM),CMO,ExFac)

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  else

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

    write(6,*) SECNAM//': wrong input parameter. ALGO= ',ALGO
    rc = 55
    return

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  end if

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  call Deallocate_DSBA(POrb(3))
  call Deallocate_DSBA(POrb(2))
  call Deallocate_DSBA(POrb(1))
  call Deallocate_DSBA(CVa(1))
  call Deallocate_DSBA(CVa(2))

  if (Deco) call Deallocate_DSBA(ChoIn)

  call Deallocate_DSBA(DLT(2))
  call Deallocate_DSBA(DLT(1))

  call deallocate_DSBA(MSQ)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
call deallocate_DSBA(CMO)
call deallocate_DSBA(FLT(2))
call deallocate_DSBA(FLT(1))

return

end subroutine CHO_CAS_DRV
