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
! Copyright (C) 1990, Roland Lindh                                     *
!               1995,1996, Anders Bernhardsson                         *
!***********************************************************************

subroutine Drvg2(Hess,nHess,lGrad,lHess)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will controll the type of the two-electron integral, eg.    *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
! Input:                                                               *
!              nHess         : Size of gradient and hessian            *
!              lGrad,lHess   : Boolean on/off for gradient/hessian     *
!                              generation                              *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March 1990                                               *
!             Anders Bernhardsson 1995-1996                            *
!***********************************************************************

use setup, only: MxPrm
use McKinley_global, only: ipDisp, ipDisp2, ipDisp3, ipMO, nFck, nMethod, RASSCF
use Index_Functions, only: nTri_Elem, nTri_Elem1
use iSD_data, only: iSD, nSD
use k2_arrays, only: DeDe, DeDe2, DoHess_, ipDijS, ipDijS2, ipOffD, ipOffDA, MxDij, nDeDe, Sew_Scr
use Disp, only: lDisp
use Etwas, only: nAsh
use pso_stuff, only: nDens
use Basis_Info, only: dbsc, nBas, nCnttp, Shells
use Symmetry_Info, only: iOper, nIrrep
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nHess
real(kind=wp), intent(out) :: Hess(nHess)
logical(kind=iwp), intent(in) :: lGrad, lHess
integer(kind=iwp) :: i, iAng, iBas, iCmp, iCnttp, id, id_Tsk, idd, iDisk, iDisp, iIrr, iIrrep, ij, ijSh, ip, iPrim, iS, iShll, &
                     jBas, jCmp, jDisp, jIrr, js, kIrr, klSh, ks, ls, mDeDe, mIndij, mmdede, moip(0:7), MxBsC, n_Int, nAco, &
                     nBuffer, ndisp, nij, nIndij, nMO, nPairs, nQuad, nSkal, nTemp
real(kind=wp) :: A_int, ThrAO, TMax_all
logical(kind=iwp) :: DoFock, DoGrad, Indexation, lpick, Post_Process
integer(kind=iwp), allocatable :: Pair_Index(:,:)
real(kind=wp), allocatable :: DInAc(:), DTemp(:), iInt(:), TMax(:,:), Buffer(:)
real(kind=wp), pointer :: Temp(:)
integer(kind=iwp), parameter :: Nr_of_Densities = 1
logical(kind=iwp), parameter :: Int_Direct = .true.
integer(kind=iwp), external :: NrOpr
logical(kind=iwp), external :: Rsv_Tsk

!                                                                      *
!***********************************************************************
!                                                                      *
! PROLOGUE

call StatusLine('McKinley: ','Computing 2-electron 2nd order derivatives')

ipDijS = 0
ipDijS2 = 0
lpick = lgrad .and. (nIrrep /= 1)

ndisp = 0
nACO = 0
do iS=0,nIrrep-1
  moip(iS) = nACO
  nDisp = nDisp+ldisp(is)
  nACO = naco+nAsh(is)
end do

Hess(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()

Indexation = .false.
ThrAO = Zero
DoFock = .false.
DoGrad = .false.
DoHess_ = .true.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate working area

MxPrm = 0
MxDij = 0
MxBsC = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
  do iCnttp=1,nCnttp
    iShll = dbsc(iCnttp)%iVal+iAng
    iPrim = Shells(iShll)%nExp
    if (iPrim == 0) cycle
    if (Shells(iShll)%nBasis == 0) cycle
    iBas = Shells(iShll)%nBasis
    iCmp = nTri_Elem1(iAng)
    MxBsC = max(MxBsC,iBas*iCmp)
    MxDij = max(MxDij,(iBas**2+1)*iCmp**2+iPrim**2+1)
  end do
end do
MxDij = 6*nIrrep*MxDij

!                                                                      *
!***********************************************************************
!                                                                      *
if (lGrad) then

  ! Calculate the size of memory needed for storing fock matrices and
  ! MO integrals and allocate it.

  nMO = nTri_Elem(nTri_Elem(naco))

  call mma_allocate(ipDisp,nDisp,label='ipDisp')
  if (nMethod == RASSCF) then
    call mma_allocate(ipMO,nDisp,label='ipMO')
    call mma_allocate(ipDisp2,nDisp,label='ipDisp2')
    call mma_allocate(ipDisp3,nDisp,label='ipDisp3')
  end if

  nIndij = nTri_Elem(S%nShlls)
  n_Int = 0
  jDisp = 0
  do iIrrep=0,nIrrep-1
    do iDisp=1,lDisp(iIrrep)

      jDisp = jDisp+1
      ipDisp(jDisp) = n_Int+1

      do jIrr=0,nIrrep-1
        kIrr = nrOpr(ieor(iOper(iIrrep),iOper(jIrr)))
        if (jIrr == kIrr) then
          n_Int = n_Int+nTri_Elem(nBas(jIrr))
        else if (kIrr < jIrr) then
          n_Int = n_Int+nBas(jIrr)*nBas(kIrr)
        end if
      end do

      if (nMethod == RASSCF) then
        ipMO(jDisp) = n_Int+1
        n_Int = n_Int+nMO
        ipDisp2(jDisp) = n_Int+1

        do jIrr=0,nIrrep-1
          kIrr = nrOpr(ieor(iOper(iIrrep),iOper(jIrr)))
          if (jIrr == jIrr) then
            n_Int = n_Int+nTri_Elem(nBas(jIrr))
          else if (kIrr < jIrr) then
            n_Int = n_Int+nBas(jIrr)*nBas(kIrr)
          end if
        end do

      end if

    end do
  end do
  if (nMethod == RASSCF) then
    jDisp = 0
    do iIrrep=0,nIrrep-1
      do iDisp=1,lDisp(iIrrep)
        jDisp = jDisp+1
        ipDisp3(jDisp) = n_Int+1
        do iS=0,nirrep-1
          js = nrOpr(ieor(iOper(is),iOper(iIrrep)))
          n_Int = n_Int+nBas(iS)*nAsh(jS)
        end do
      end do
    end do
  end if
  call mma_allocate(iInt,n_Int,Label='iInt')
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Desymmetrize  densities.
  ! Observe that the desymmetrized 1st order density matrices are canonical,
  ! i.e. the relative order of the indices are canonically ordered.

  iInt(:) = Zero
  call mma_allocate(DTemp,nDens,Label='DTemp')
  DTemp(:) = Zero
  call mma_allocate(DInAc,nDens,Label='DInAc')
  DInAc(:) = Zero
  if (nIrrep == 1) then
    if (nmethod /= RASSCF) then
      call Get_D1ao_Var(DTemp,nDens)
      DTemp(:) = Half*DTemp
      ij = 0
      do i=1,nBas(0)
        ij = ij+i
        DTemp(ij) = Two*DTemp(ij)
      end do
    else
      call Din(DInAc)
      DInAc(:) = Half*DInAc
      ij = 0
      do i=1,nBas(0)
        ij = ij+i
        DInAc(ij) = Two*DInAc(ij)
      end do
      call Dan(DTemp)
      DTemp(:) = Half*DTemp
      ij = 0
      do i=1,nBas(0)
        ij = ij+i
        DTemp(ij) = Two*DTemp(ij)
      end do
    end if
  else
    mmDeDe = nDeDe
    call mma_allocate(ipOffD,3,nIndij,label='ipOffD')
    call mma_allocate(DeDe,mmDeDe+MxDij,label='DeDe')
    ipDijS = 1+mmDeDe
    if (nMethod /= RASSCF) then
      call Get_D1ao_Var(DTemp,nDens)
      call DeDe_mck(DTemp,nFck(0),ipOffD,nIndij,Dede,mmDeDe,mDeDe,mIndij)
    else
      call mma_allocate(ipOffDA,3,nIndij,Label='ipOffDA')
      call mma_allocate(DeDe2,mmDeDe+MxDij,label='DeDe2')
      ipDijS2 = 1+mmDeDe

      call Dan(DTemp)
      call DeDe_mck(DTemp,nFck(0),ipOffD,nIndij,DeDe,mmDeDe,mDeDe,mIndij)

      call Din(DInAc)
      call DeDe_mck(DInAc,nFck(0),ipOffDA,nIndij,DeDe2,mmDeDe,mDeDe,mIndij)

      if (mDeDe /= nDeDe) then
        write(u6,*) 'DrvG2: mDeDe /= nDeDe'
        write(u6,*) 'mDeDe,nDeDe=',mDeDe,nDeDe
        call Abend()
      end if
    end if
  end if

  call mma_allocate(DeDe,[-1,-1],label='DeDe',safe='*') ! Dummy allocation
  call mma_allocate(DeDe2,[-1,-1],label='DeDe2',safe='*') ! Dummy allocation

end if ! lGrad
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal)
call Setup_iSD()

nPairs = nTri_Elem(nSkal)
nQuad = nTri_Elem(nPairs)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
call Shell_MxSchwz(nSkal,TMax)

TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do

!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Pair_Index,2,nPairs,Label='Ind_ij')
nij = 0
nBuffer = 1  ! Dummy length
do iS=1,nSkal
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  do jS=1,iS
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      Pair_Index(1,nij) = iS
      Pair_Index(2,nij) = jS
      if ((nMethod == RASSCF) .and. lGrad) nBuffer = max(nBuffer,nTri_Elem(nACO)*iCmp*iBas*jCmp*jBas*nDisp*nIrrep)
    end if
  end do
end do
call mma_allocate(Buffer,nBuffer,Label='Buffer')
!                                                                      *
!***********************************************************************
!                                                                      *
call Init_Tsk(id_Tsk,nij)
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes

! make reservation of a task on global task list and get task range
! in return. Function will be false if no more tasks to execute.
do while (Rsv_Tsk(id_Tsk,ijSh))
  iS = Pair_Index(1,ijSh)
  jS = Pair_Index(2,ijSh)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Outer loops (ij) over angular momenta and centers

  Buffer(:) = Zero

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Post_Process = .false.
  do klSh=1,nij
    ks = Pair_Index(1,klSh)
    ls = Pair_Index(2,klSh)

    A_int = TMax(iS,jS)*TMax(kS,lS)
    if (A_Int < CutInt) cycle

    call Eval_g2_ijkl(iS,jS,kS,lS,Hess,nHess,Post_Process,iInt,n_Int,nACO,lGrad,lHess,lPick,nBuffer,Buffer,nDens,DTemp,DInAc,moip, &
                      n_Int,nQuad)

  end do ! klS

  if ((nMethod == RASSCF) .and. Post_Process) then
    nTemp = size(Sew_Scr)
    Temp(1:nTemp) => Sew_Scr(1:nTemp)
    call CLR2(Buffer,iInt,nACO,nSD,iSD(:,iS),iSD(:,jS),nDisp,nTemp,Temp)
    nullify(Temp)
  end if

end do
call mma_deallocate(Buffer)
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
! EPILOGUE
!                                                                      *
!***********************************************************************
!                                                                      *
if (nIrrep == 1) then
  idd = 0
  do iS=0,nirrep-1
    do iD=1,ldisp(is)
      idd = idd+1
      ip = ipDisp(idd)
      iInt(ip:ip+nDens-1) = Half*iInt(ip:ip+nDens-1)
      ij = ip-1
      do i=1,nBas(0)
        ij = ij+i
        iInt(ij) = Two*iInt(ij)
      end do
    end do
  end do
  if (nmethod == RASSCF) then
    idd = 0
    do iS=0,nirrep-1
      do iD=1,ldisp(is)
        idd = idd+1
        ip = ipDisp2(idd)
        iInt(ip:ip+nDens-1) = Half*iInt(ip:ip+nDens-1)
        ij = ip-1
        do i=1,nBas(0)
          ij = ij+i
          iInt(ij) = Two*iInt(ij)
        end do
      end do
    end do

  end if
end if
call mma_deallocate(Sew_Scr,safe='*')
call Free_Tsk(id_Tsk)

! YIPPIEEEE Finished OK fill it UP!!

call GADSum(iInt,n_Int)
jDisp = 0
do iIrr=0,nIrrep-1
  do iDisk=1,lDisp(iIrr)
    jDisp = jDisp+1
    call WrDisk(iInt,n_Int,jDisp,iIrr)
  end do
end do

call mma_deallocate(Pair_Index)
call mma_deallocate(TMax)
call Free_iSD()

call mma_deallocate(DeDe)
call mma_deallocate(DeDe2)
if (nIrrep /= 1) then
  call mma_deallocate(ipOffD)
  if (nMethod == RASSCF) then
    call mma_deallocate(ipOffDA)
  end if
end if

call mma_deallocate(DInAc)
call mma_deallocate(DTemp)
call mma_deallocate(iInt)

call Term_Ints()

call mma_deallocate(ipDisp,safe='*')
call mma_deallocate(ipDisp2,safe='*')
call mma_deallocate(ipDisp3,safe='*')
call mma_deallocate(ipMO,safe='*')

end subroutine Drvg2
