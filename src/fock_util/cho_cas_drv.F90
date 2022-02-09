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

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Fock_util_global, only: ALGO, Deco, dmpk, DoActive, DoLocK, Nscreen
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
real(kind=wp), intent(in) :: W_CMO(*), DI(*), DA1(*)
real(kind=wp), intent(inout) :: FI(*), FA(*)
real(kind=wp), intent(_OUT_) :: W_PWXY(*)
logical(kind=iwp), intent(in) :: TraOnly
#include "rasdim.fh"
#include "wadr.fh"
#include "general.fh"
#include "rasscf.fh"
integer(kind=iwp) :: i, iBas, iFro, ikk, incs, iOrb, iSym, ja, jkk, nAorb(8), nChI(8), nChM(8), nForb(8), nIorb(8), NumV
real(kind=wp) :: dmpk_old, FactXI, Ymax
type(DSBA_Type) :: ChoIn, CMO, CVa(2), Ddec, DLT(2), FLT(2), FLT_MO(2), MSQ, POrb(3)
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)
real(kind=wp), parameter :: Thr = 1.0e-12_wp
character(len=*), parameter :: SECNAM = 'CHO_CAS_DRV'

!                                                                      *
!***********************************************************************
!                                                                      *
rc = 0

call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=FI)
call Allocate_DT(FLT(2),nBas,nBas,nSym,aCase='TRI',Ref=FA)
call Allocate_DT(CMO,nBas,nBas,nSym,Ref=W_CMO)

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

  call Allocate_DT(FLT_MO(1),nOrb,nOrb,nSym,aCase='TRI',Ref=FI)
  call Allocate_DT(FLT_MO(2),nOrb,nOrb,nSym,aCase='TRI',Ref=FA)
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
      call DGEMM_('N','N',iBas,iOrb,iBas,One,Tmp1,iBas,CMO%SB(iSym)%A1(1+iFro*iBas:),iBas,Zero,Tmp2,iBas)
      call DGEMM_Tri('T','N',iOrb,iOrb,iBas,One,Tmp2,iBas,CMO%SB(iSym)%A1(1+iFro*iBas:),iBas,Zero,FLT_MO(i)%SB(iSym)%A1,iOrb)
      call mma_deallocate(Tmp2)
      call mma_deallocate(Tmp1)
    end do
  end do

  call Deallocate_DT(FLT_MO(1))
  call Deallocate_DT(FLT_MO(2))
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
  call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI')
  call Allocate_DT(DLT(2),nBas,nBas,nSym,aCase='TRI')

  call Fold(nSym,nBas,DI,DLT(1)%A0)
  call Fold(nSym,nBas,DA1,DLT(2)%A0)

  FactXI = -One

  !AMS - should this be set differently for ExFac /= 1?
  !      FactXI = 0-ExFac

  if (Deco) then

    FactXI = -Half

    !AMS - should this be set differently for ExFac /= 1?
    !      FactXI = 0-(ExFac*Half)

    ! --- decompose the Inactive density on request
    call Allocate_DT(ChoIn,nBas,nBas,nSym)
    call Allocate_DT(DDec,nBas,nBas,nSym)
    DDec%A0(1:NTot2) = DI(1:NTot2)

    call Allocate_DT(MSQ,nBas,nBas,nSym,Ref=ChoIn%A0)

    incs = 0
    do i=1,nSym
      if ((nForb(i)+nIorb(i)) > 0) then
        call CD_InCore(DDec%SB(i)%A2,nBas(i),ChoIn%SB(i)%A2,nBas(i),NumV,Thr,rc)
        if (rc /= 0) then
          write(u6,*) SECNAM//': ill-defined dens decomp for Inact'
          write(u6,*) 'rc value produced = ',rc
          call abend()
        end if
        nChI(i) = NumV
        if (NumV /= nIsh(i)+nForb(i)) then
          write(u6,*) 'Warning! The number of occupied from the decomposition of the Inactive density matrix is ',numV, &
                      ' in symm. ',i
          write(u6,*) 'Expected value = ',nIsh(i)+nForb(i)
          incs = incs+1
          Ymax = Zero
          do ja=1,nBas(i)
            Ymax = max(Ymax,DDec%SB(i)%A2(ja,ja))
          end do
          write(u6,*) 'Max diagonal of the density in symm. ',i,' is equal to ',Ymax
        end if
      else
        nChI(i) = 0
      end if
    end do

    if ((incs > 0) .and. DoLocK) then
      dmpk_old = dmpk
      dmpk = 1.0e-2_wp*dmpk
      write(u6,*) 'LK-damping decreased from ',dmpk_old,' to ',dmpk
    end if

    call Deallocate_DT(DDEc)

    ! --- to get the right input arguments for CHO_FCAS_AO and CHO_FMCSCF
    if (.not. DoLocK) then
      nForb(:) = 0
      nIorb(:) = nChI(:)
    end if

  else

    call Allocate_DT(MSQ,nBas,nBas,nSym,Ref=W_CMO)

    nChI(:) = nForb(:)+nIorb(:)

  end if

  ! Reordering of the MOs coefficients to fit cholesky needs

  if (.not. DoLocK) then

    call Allocate_DT(POrb(1),nChI,nBas,nSym)
    call Allocate_DT(POrb(3),nAOrb,nBas,nSym)

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
    call Allocate_DT(CVa(1),nAorb,nBas,nSym)

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
      write(u6,*) 'DALT for sym=',i
      call TRIPRT('DALT',' ',DLT(2)%SB(i)%A1,nBas(i))
    end do
#   endif

    call Allocate_DT(CVa(2),nBas,nBas,nSym)
    call Allocate_DT(DDec,nBas,nBas,nSym)
    DDec%A0(1:NTot2) = DA1(1:NTot2)

    do i=1,nSym
      if (nAorb(i) > 0) then
        ! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
        !                 This will avoid warnings for negative-definite
        call CD_InCore(DDec%SB(i)%A2,nBas(i),CVa(2)%SB(i)%A2,nBas(i),NumV,Thr,rc)
        if (rc /= 0) then
          write(u6,*) SECNAM//': ill-defined dens decomp for active'
          write(u6,*) 'rc value produced = ',rc
          call abend()
        end if
        nChM(i) = NumV
      else
        nChM(i) = 0
      end if
    end do

    call Deallocate_DT(DDec)

  else

    ! Dummy allocation
    call Allocate_DT(CVa(2),[1],[1],1)
    nChM(:) = 0

  end if

  if ((.not. DoLocK) .and. DoActive) then

    ! reorder "Cholesky MOs" to Cva storage

    call Allocate_DT(POrb(2),nChM,nBas,nSym)
    do iSym=1,nSym
      if (nBas(iSym)*nChM(iSym) /= 0) then
        do ikk=1,nChM(iSym)
          POrb(2)%SB(iSym)%A2(ikk,:) = CVa(2)%SB(iSym)%A2(:,ikk)
        end do
      end if
    end do

  else

    call Allocate_DT(POrb(2),[1],[1],1)

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

    write(u6,*) SECNAM//': wrong input parameter. ALGO= ',ALGO
    rc = 55
    return

    !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  end if

  !)()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()

  call Deallocate_DT(POrb(3))
  call Deallocate_DT(POrb(2))
  call Deallocate_DT(POrb(1))
  call Deallocate_DT(CVa(1))
  call Deallocate_DT(CVa(2))

  if (Deco) call Deallocate_DT(ChoIn)

  call Deallocate_DT(DLT(2))
  call Deallocate_DT(DLT(1))

  call Deallocate_DT(MSQ)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
call Deallocate_DT(CMO)
call Deallocate_DT(FLT(2))
call Deallocate_DT(FLT(1))

return

end subroutine CHO_CAS_DRV
