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

use setup, only: MxPrm, nAux
use McKinley_global, only: ipDisp, ipDisp2, ipDisp3, ipMO, nFck, nMethod, nTwoDens, RASSCF
use Index_Functions, only: iTri, nTri_Elem, nTri_Elem1
use iSD_data, only: iSD, nSD
use k2_arrays, only: Aux, Create_BraKet_Base, DeDe, DeDe2, Destroy_BraKet_Base, ipDijS, ipDijS2, &
                     ipOffD, ipOffDA, MxDij, nDeDe, Sew_Scr
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
integer(kind=iwp) :: i, iBas, iCmp, iCnttp, &
                     id, id_Tsk, idd, ider, iDisk, iDisp, iIrr, iIrrep, ij, ijSh,  &
                     ip, iPrim, &
                     iS, iShll, jBas, jCmp, jDisp, jIrr, js, kIrr, klSh, iAng, ks, ls, &
                     mDeDe, MemBuffer, &
                     mIndij, mmdede, moip(0:7), MxBsC, n_Int, nAco, nb, ndisp, &
                     nijS, nIndij, nMO, nPairs, nQuad, nSkal, nTwo, nTwo2, iSD4(0:nSD,4), nTemp, &
                     ipDum
real(kind=wp) :: A_int, Prem, Pren, TMax_all
logical(kind=iwp) :: lpick, new_fock, Post_Process
#ifdef _DEBUGPRINT_
character(len=40) :: frmt
#endif
integer(kind=iwp), allocatable :: Ind_ij(:,:)
real(kind=wp), allocatable :: DInAc(:), DTemp(:), iInt(:), TMax(:,:), Buffer(:)
real(kind=wp), pointer :: Temp(:)
integer(kind=iwp), parameter :: Nr_of_Densities = 1
logical(kind=iwp), parameter :: Int_Direct = .true.
integer(kind=iwp), external :: MemSO2_P, NrOpr
logical(kind=iwp), external :: Rsv_Tsk

!                                                                      *
!***********************************************************************
!                                                                      *
! PROLOGUE

call StatusLine('McKinley: ','Computing 2-electron 2nd order derivatives')

ipDijS = 0
ipDijS2 = 0

ndisp = 0
nACO = 0
New_Fock = nirrep == 1
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
!                                                                      *
!***********************************************************************
!                                                                      *
! Precompute k2 entities.

lpick = lgrad .and. (.not. New_Fock)
Pren = Zero
Prem = Zero

call Drvk2_mck(new_Fock)

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate auxiliary array for symmetry transformation

nAux = nIrrep**3
if (nIrrep == 1) nAux = 1
call mma_allocate(Aux,nAux,Label='Aux')
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

call Create_BraKet_Base(MxPrm**2)
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
  nTwo = 0
  do iIrrep=0,nIrrep-1
    nTwo = max(nTwo,nFck(iIrrep))
  end do
  if (Int_Direct) then
    nTwo2 = n_Int
  else
    nTwo2 = nTwo
  end if
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
  if (New_Fock) then
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
    mmdede = ndede
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

  nb = 0
  do is=0,nIrrep-1
    nb = nb+nBas(iS)
  end do

  call mma_allocate(DeDe,0,label='DeDe',safe='*') ! Dummy allocation
  call mma_allocate(DeDe2,0,label='DeDe2',safe='*') ! Dummy allocation

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

call mma_allocate(Ind_ij,2,nPairs,Label='Ind_ij')
nijS = 0
MemBuffer = 1  ! Dummy length
do iS=1,nSkal
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  do jS=1,iS
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nijS = nijS+1
      Ind_ij(1,nijS) = iS
      Ind_ij(2,nijS) = jS
      if ((nMethod == RASSCF) .and. lGrad)   &
         MemBuffer = Max(MemBuffer,nTri_Elem(nACO)*iCmp*iBas*jCmp*jBas*nDisp*nIrrep)
    end if
  end do
end do
call Init_Tsk(id_Tsk,nijS)
!                                                                    *
!*********************************************************************
!                                                                    *
! Cltrls for MO transformation
!                                                                    *
!*********************************************************************
!                                                                    *
Call mma_allocate(Buffer,MemBuffer,Label='Buffer')
!                                                                      *
!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes

! make reservation of a task on global task list and get task range
! in return. Function will be false if no more tasks to execute.
do while (Rsv_Tsk(id_Tsk,ijSh))
  iS = Ind_ij(1,ijSh)
  jS = Ind_ij(2,ijSh)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Outer loops (ij) over angular momenta and centers

  Buffer(:) = Zero

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Post_Process = .false.
  do klSh=1,nijS
    ks = Ind_ij(1,klSh)
    ls = Ind_ij(2,klSh)

    A_int = TMax(iS,jS)*TMax(kS,lS)
    if (A_Int < CutInt) cycle

    Call Eval_g2_ijkl(iS,jS,kS,lS,Hess,nHess,Post_Process,iInt,n_Int,nACO,lHess,lPick,MemBuffer, &
                      Buffer,nDens, DTemp, DInAc)

  end do ! klS

  if ((nMethod == RASSCF) .and. Post_Process) then
    nTemp = Size(Sew_Scr)
    Temp(1:nTemp) => Sew_Scr(1:nTemp)
    call CLR2(Buffer,iInt,nACO,nSD,iSD(:,iS),iSD(:,jS),nDisp,nTemp,Temp)
    nullify(Temp)
  end if

end do
Call mma_deallocate(Buffer)
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
! EPILOGUE
!                                                                      *
!***********************************************************************
!                                                                      *
if (New_Fock) then
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
#ifdef _DEBUGPRINT_
call GADSum_SCAL(Pren)
call GADSum_SCAL(Prem)
write(frmt,'(A,I2,A,I2,A)') '(A,F',3+int(log10(Pren)),'.0,A,F',3+int(log10(Prem)),'.0,A)'
write(u6,frmt) ' A total of',Pren,' entities were prescreened and',Prem,' were kept.'
#endif
call mma_deallocate(Sew_Scr)
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

call mma_deallocate(Ind_ij)
call mma_deallocate(TMax)
call Free_iSD()

call mma_deallocate(DeDe)
call mma_deallocate(DeDe2)
if (.not. New_Fock) then
  call mma_deallocate(ipOffD)
  if (nMethod == RASSCF) then
    call mma_deallocate(ipOffDA)
  end if
end if

call Destroy_BraKet_Base()

call mma_deallocate(DInAc)
call mma_deallocate(DTemp)
call mma_deallocate(iInt)

call mma_deallocate(Aux)

call Term_Ints()

call mma_deallocate(ipDisp,safe='*')
call mma_deallocate(ipDisp2,safe='*')
call mma_deallocate(ipDisp3,safe='*')
call mma_deallocate(ipMO,safe='*')

contains

subroutine Dens_Infos(nMethod)

  use Dens_stuff, only: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, &
                        ipDij, ipDij2, ipDik, ipDik2, ipDil, ipDil2, ipDjk, ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, mDCRij, mDCRik, &
                        mDCRil, mDCRjk, mDCRjl, mDCRkl
  use k2_arrays, only: ipDijS, ipDijS2
  use Index_Functions, only: iTri

  integer(kind=iwp), intent(in) :: nMethod
  integer(kind=iwp) :: ijS, ikS, ilS, ipTmp, ipTmp2, iS, jkS, jlS, jS, klS, kS, lS
  integer(kind=iwp), parameter :: Nr_of_D = 1

  iS = iSD4(11,1)
  jS = iSD4(11,2)
  kS = iSD4(11,3)
  lS = iSD4(11,4)

  ijS = iTri(iS,jS)
  klS = iTri(kS,lS)
  ikS = iTri(iS,kS)
  ilS = iTri(iS,lS)
  jkS = iTri(jS,kS)
  jlS = iTri(jS,lS)
  ijS = iTri(iS,jS)
  klS = iTri(kS,lS)
  ikS = iTri(iS,kS)
  ilS = iTri(iS,lS)
  jkS = iTri(jS,kS)
  jlS = iTri(jS,lS)
  ipTmp = ipDijs
  if (nMethod == RASSCF) ipTmp2 = ipDijs2
  call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDij2,ipDDij2)
  call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDkl2,ipDDkl2)
  call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDik2,ipDDik2)
  call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDil2,ipDDil2)
  call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDjk2,ipDDjk2)
  call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,nr_of_Densities,nMethod,ipTmp2,ipDjl2,ipDDjl2)

end subroutine Dens_Infos

subroutine Eval_g2_ijkl(iS,jS,kS,lS,Hess,nHess,Post_Process,iInt,n_Int,nACO,lHess,lPick,MemBuffer, &
                        Buffer,nDens, DTemp, DInAc)
use setup, only: nAux
use McKinley_global, only: nMethod, RASSCF
use Index_Functions, only: iTri
use Definitions, only: wp, iwp, u6
use iSD_data, only: iSD, nSD
use k2_arrays, only: Create_Braket, Destroy_Braket, Sew_Scr, nFT, Aux
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero

Implicit None
integer(kind=iwp), intent(in):: iS, jS, kS, lS, nHess, n_Int, nACO, MemBuffer, nDens
real(kind=wp), intent(inout) :: Hess(nHess), iInt(n_Int), Buffer(MemBuffer), DTemp(nDens), DInAc(nDens)
logical(kind=iwp), intent(inout):: Post_Process
logical(kind=iwp), intent(in):: lHess, lPick

real(kind=wp) :: Coor(3,4), PMax
logical(kind=iwp) :: lTri, lDot, lDot2
logical(kind=iwp), parameter :: n8=.true.
integer(kind=iwp) :: nSO, nRys, iFnc(4)
integer(kind=iwp) :: iBasAO, jBasAO, kBasAO, lBasAO
integer(kind=iwp) :: iBasi, jBasj, kBask, lBasl
integer(kind=iwp) :: iBsInc, jBsInc, kBsInc, lBsInc
integer(kind=iwp) :: iBasn, jBasn, kBasn, lBasn
integer(kind=iwp) :: JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7)
logical(kind=iwp) :: JfG(4), JfGrd(3,4), JfHss(4,3,4,3)
integer(kind=iwp) :: MemMax, ipMOC, MemCMO
integer(kind=iwp) :: Mem1, Mem2, Mem3, Mem4
integer(kind=iwp) :: ipPSO, ipFin, ipMem2, ipMem3, ipMem4, ipMemX
integer(kind=iwp) :: MemFck, MemFin, MemPrm, MemPSO, MemX
integer(kind=iwp) :: kCmp, lCmp, nijkl
real(kind=wp), pointer :: Fin(:), MOC(:), PSO(:,:), Work2(:), Work3(:), Work4(:), WorkX(:)

iFnc(:)=-99
PMax=Zero
if (.not. allocated(Sew_Scr)) Then
   call mma_MaxDBLE(MemMax)
   if (MemMax > 8000) MemMax = MemMax-8000
   call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
else
   MemMax=Size(Sew_Scr)
endif
ipMOC = 1

call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

Call Coor_setup(iSD4,nSD,Coor)
!                                                                  *
!*******************************************************************
!                                                                  *
! The code is working in such away that the MO needs upper and lower
! triangular parts of ij kl but hessian needs only lower, check if the
! integralbatch is lower or upper!!

lTri = iTri(iS,jS) >= iTri(kS,lS)
if ((.not. lTri) .and. (nMethod /= RASSCF)) Return
lDot = (lTri .and. lHess)

!                                                                  *
!*******************************************************************
!                                                                  *
! Allocate memory for zeta, eta, kappa, P and Q.
! Allocate also for Alpha, Beta , Gamma and Delta in expanded form.

call Create_BraKet(iSD4(5,1)*iSD4(5,2),iSD4(5,3)*iSD4(5,4))
!                                                                  *
!*******************************************************************
!                                                                  *
! Fix the 1st order density matrix

! Pick up pointers to desymmetrized 1st order density matrices.
! Observe that the desymmetrized 1st order density matrices
! follow the contraction index.

if (lTri .and. lPick) call Dens_Infos(nMethod)

!                                                                  *
!*******************************************************************
!                                                                  *
! Compute total size of the second order density matrix in SO basis.
!------------------------------------------------------------------*
nSO = MemSO2_P(nSD,iSD4)
ldot2 = ldot
if (nSO == 0) ldot2 = .false.

! Compute memory request for the primitives.

iDer = 2
if (.not. ldot2) iDer = 1
call MemRg2(iSD4(1,:),nRys,MemPrm,iDer)

!------------------------------------------------------------------*
!
! Calculate which derivatives should be made.
!
!------------------------------------------------------------------*

call DerCtr(ldot2,JfGrd,JndGrd,JfHss,JndHss,JfG,nSD,iSD4)

!------------------------------------------------------------------*
!
! Decide on the partioning of the shells based on the
! available memory and the requested memory.
!
!------------------------------------------------------------------*

call PSOAO2(nSO,MemPrm,MemMax,iFnc,nAco,Mem1,Mem2,Mem3,Mem4,MemX,MemPSO,MemFck,nFT,MemFin,MemBuffer,nSD,iSD4)

iBasi = iSD4(3,1)
jBasj = iSD4(3,2)
kBask = iSD4(3,3)
lBasl = iSD4(3,4)

iBsInc = iSD4(4,1)
jBsInc = iSD4(4,2)
kBsInc = iSD4(4,3)
lBsInc = iSD4(4,4)

kCmp = iSD(2,kS)
lCmp = iSD(2,lS)

!------------------------------------------------------------------*
!
! Loop over basis function if we do not have enough of memory to
! calculate them in one step.
!
!------------------------------------------------------------------*
do iBasAO=1,iBasi,iBsInc
  iBasn = min(iBsInc,iBasi-iBasAO+1)
  iSD4(8,1) = iBasAO-1
  iSD4(19,1) = iBasn

  !----------------------------------------------------------------*
  !
  ! Move appropriate portions of the desymmetrized 1st order density matrix.
  !
  !----------------------------------------------------------------*
  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iSD4(8,2) = jBasAO-1
    iSD4(19,2) = jBasn

    if (lpick) call Picky_Mck(nSD,iSD4,1,2,nMethod)

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iSD4(8,3) = kBasAO-1
      iSD4(19,3) = kBasn

      if (lpick) then
        call Picky_Mck(nSD,iSD4,1,3,nMethod)
        call Picky_Mck(nSD,iSD4,2,3,nMethod)
      end if

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iSD4(8,4) = lBasAO-1
        iSD4(19,4) = lBasn

        if (lpick) then
          call Picky_Mck(nSD,iSD4,3,4,nMethod)
          call Picky_Mck(nSD,iSD4,1,4,nMethod)
          call Picky_Mck(nSD,iSD4,2,4,nMethod)
        end if

        !----------------------------------------------------------*

        nijkl = iBasn*jBasn*kBasn*lBasn

        ! Mark out the memory allocations explicitly with pointers
        ! MO tranformation buffer
        MemCMO = nACO*(kCmp*kBasn+lCmp*lBasn)
        MOC(1:MemCMO) => Sew_Scr(ipMOC:ipMOC+MemCMO-1)
        ! Area for the AO integrals
        ipFin = ipMOC+MemCMO
        Fin(1:MemFin) => Sew_Scr(ipFin:ipFin+MemFin-1)
        ! Area for 2el density
        ipPSO = ipFin+MemFin
        PSO(1:nijkl,1:nSO) => Sew_Scr(ipPSO:ipPSO+nijkl*nSO-1)
        if (nijkl*nSO > Mem1) then
          write(u6,'(A)') 'nijkl*nSO>Mem1'
          write(u6,*) 'njikl,nSO=',nijkl,nSO
          write(u6,*) 'Mem1=',Mem1
          call Abend()
        end if
        ipMem2 = ipPSO+Mem1  ! Work
        Work2(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem2-1)
        ipMem3 = ipMem2+Mem2 ! Work
        Work3(1:Mem3) => Sew_Scr(ipMem3:ipMem3+Mem3-1)
        ipMemX = ipMem3+Mem3 ! Work
        WorkX(1:MemX) => Sew_Scr(ipMemX:ipMemX+MemX-1)

        ! If MO transformation is performed in the standard way
        ! reserve memory for partial transformed integrals

        ! Multilayer

        ipMem4 = ipMem2+Mem2-Mem4
        Work4(1:Mem4) => Sew_Scr(ipMem4:ipMem4+Mem4-1)
        nTemp = Mem2+Mem3+MemX
        Temp(1:nTemp) => Sew_Scr(ipMem2:ipMem2+nTemp-1)

        !----------------------------------------------------------*
        !
        ! Get the 2nd order density matrix in SO basis.
        !
        !----------------------------------------------------------*

        if (n8) call PickMO(MOC,MemCMO,nSD,iSD4)
        if (ldot2) call PGet0(nijkl,PSO,nSO,iFnc,MemPSO,Work2,Mem2,nQuad,PMax,iSD4)

        ! Compute gradients of shell quadruplet

        call TwoEl_mck(Coor,nRys,Pren,Prem,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,JfG,PSO,nijkl,nSO,Work2,Mem2,Work3,Mem3,Work4, &
                       Mem4,Aux,nAux,WorkX,MemX,Fin,MemFin,Temp,nTemp,nTwo2,nFT,iInt,Buffer,MemBuffer,lgrad,ldot2,n8,ltri, &
                       DTemp,DInAc,moip,nAco,MOC,MemCMO,new_fock,iSD4)
        Post_Process = .true.

        nullify(MOC)
        nullify(Fin)
        nullify(PSO)
        nullify(Work2)
        nullify(Work3)
        nullify(WorkX)
        nullify(Work4)
        nullify(Temp)

        !----------------------------------------------------------*

      end do
    end do
  end do
end do
call Destroy_Braket()

end subroutine Eval_g2_ijkl

end subroutine Drvg2
