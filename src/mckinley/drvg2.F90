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

subroutine Drvg2(Hess,nHess,l_Grd,l_Hss)
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
!              l_Grd,l_Hss   : Boolean on/off for gradient/hessian     *
!                              generation                              *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March 1990                                               *
!             Anders Bernhardsson 1995-1996                            *
!***********************************************************************

use setup, only: MxPrm, nAux
use McKinley_global, only: CPUStat, ipDisp, ipDisp2, ipDisp3, ipMO, nFck, nMethod, nTwoDens, RASSCF
use Index_Functions, only: iTri, nTri_Elem, nTri_Elem1
use iSD_data, only: iSD, nSD
use k2_arrays, only: Aux, Create_Braket, Create_BraKet_Base, DeDe, Destroy_Braket, Destroy_BraKet_Base, ipDijS, ipOffD, MxDij, &
                     ndede, nFT, Sew_Scr, ipOffDA
use Disp, only: lDisp
use Etwas, only: nAsh
use pso_stuff, only: nDens
use Basis_Info, only: dbsc, nBas, nCnttp, Shells
use Symmetry_Info, only: iOper, nIrrep
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6
use Dens_stuff, only: mDCRij,mDCRkl,mDCRik,mDCRil,mDCRjk,mDCRjl,&
                      ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl,&
                       ipDij, ipDkl, ipDik, ipDil, ipDjk, ipDjl,&
                        mDij,  mDkl,  mDik,  mDil,  mDjk,  mDjl



implicit none
integer(kind=iwp), intent(in) :: nHess
real(kind=wp), intent(out) :: Hess(nHess)
logical(kind=iwp), intent(in) :: l_Grd, l_Hss
integer(kind=iwp) :: i, iBas, iBasAO, ibasI, iBasn, iBsInc, iCmp, iCmpV(4), iCnt, iCnttp, &
                     id, id_Tsk, idd, ider, iDisk, iDisp, iFnc(4), iii, iIrr, iIrrep, ij, ijS, ijSh,  ikS, ilS, &
                     ip, ipPSO, ipDDij2, ipDDik2, ipDDil2, &
                     ipDDjk2, ipDDjl2, ipDDkl2, ipDij2, ipDijS2, ipDik2, ipDil2,  &
                     ipDjk2, ipDjl2, ipDkl2, ipFin, ipMem, ipMem2, ipMem3, ipMem4, ipMemX, ipMOC, iPrim, iPrimi, &
                     ipTmp, ipTmp2, iS, iShell, iShll, jBas, jBasAO, jBasj, jBasn, &
                     jBsInc, jCmp, jCnt, jCnttp, jDisp, jIrr, jkS, jlS, JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), jPrimj, &
                     js, jShell, kBasAO, kBask, kBasn, kBsInc, kCmp, kCnt, kCnttp, kIrr, klS, klSh, kPrimk, iAng, &
                     ks, kShell, lBasAO, lBasl, lBasn, lBsInc, lCmp, lCnt, lCnttp, lPriml, ls, &
                     lShell, mDeDe, Mem1, Mem2, Mem3, Mem4, MemBuffer, MEMCMO, memCMO2, MemFck, MemFin, MemMax, MemPrm, &
                     MemPSO, MemX, mIndij, mmdede, moip(0:7), MxBsC, n_Int, nAco, nb, nDik, nDil, ndisp, nDjk, &
                     nDjl, nDkl, nijkl, nijS, nIndij, nMO, nPairs, nQuad, nRys, nSkal, nSO, nTwo, nTwo2, iSD4(0:nSD,4), nTemp, &
                     ipDum
real(kind=wp) :: A_int, dum1, dum2, dum3, Coor(3,4), PMax, Prem, Pren, TCpu1, TCpu2, Time, TMax_all, TWall1, TWall2
logical(kind=iwp) :: JfG(4), JfGrd(3,4), JfHss(4,3,4,3), ldot, ldot2, lGrad, lpick, ltri, n8, new_fock, Post_Process, Shijij, &
                     Shik, Shjl
#ifdef _DEBUGPRINT_
character(len=40) :: frmt
#endif
logical(kind=iwp), parameter :: Int_Direct = .true.
integer(kind=iwp), allocatable :: Ind_ij(:,:)
real(kind=wp), allocatable :: DeDe2(:), DInAc(:), DTemp(:), iInt(:), TMax(:,:)
integer(kind=iwp), external :: MemSO2_P, NrOpr
logical(kind=iwp), external :: Rsv_Tsk
real(kind=wp), pointer :: Buffer(:)=>Null(), MOC(:)=>Null(), Fin(:)=>Null(), PSO(:,:)=>Null(), Temp(:)=>Null()
real(kind=wp), pointer :: Work2(:)=>Null(), Work3(:)=>Null(), WorkX(:)=>Null(), Work4(:)=>Null()
integer(kind=iwp), parameter :: Nr_of_Densities=1

!                                                                      *
!***********************************************************************
!                                                                      *
! PROLOGUE

call StatusLine('McKinley: ','Computing 2-electron 2nd order derivatives')

ipDij2 = 0
ipDDij = 0
ipDDij2 = 0
ipDkl2 = 0
ipDDkl = 0
ipDDkl2 = 0
ipDik2 = 0
ipDDik = 0
ipDDik2 = 0
ipDil2 = 0
ipDDil = 0
ipDDil2 = 0
ipDjk2 = 0
ipDDjk = 0
ipDDjk2 = 0
ipDjl2 = 0
ipDDjl = 0
ipDDjl2 = 0

ipMOC = 0

iFnc(:) = -99
nDkl = 0
nDik = 0
nDjl = 0
nDil = 0
nDjk = 0
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
n8 = .true.

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

lgrad = l_Grd
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
iii = nDens*10+10

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
do iS=1,nSkal
  do jS=1,iS
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nijS = nijS+1
      Ind_ij(1,nijS) = iS
      Ind_ij(2,nijS) = jS
    end if
  end do
end do
call Init_Tsk(id_Tsk,nijS)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
if (MemMax > 8000) MemMax = MemMax-8000
call mma_allocate(Sew_Scr,MemMax-iii,Label='Sew_Scr')
ipMem = 1
memmax = memmax-iii
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes

! make reservation of a task on global task list and get task range
! in return. Function will be false if no more tasks to execute.
do while (Rsv_Tsk(id_Tsk,ijSh))
  iS = Ind_ij(1,ijSh)
  jS = Ind_ij(2,ijSh)
  call CWTime(TCpu1,TWall1)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Outer loops (ij) over angular momenta and centers
  !

  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  Coor(1:3,1) = dbsc(iCnttp)%Coor(1:3,iCnt)

  jCmp = iSD(2,jS)
  jBas = iSD(3,jS)
  jShell = iSD(11,jS)
  jCnttp = iSD(13,jS)
  jCnt = iSD(14,jS)
  Coor(1:3,2) = dbsc(jCnttp)%Coor(1:3,jCnt)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Cltrls for MO transformation
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if ((nMethod == RASSCF) .and. l_Grd) then
    MemBuffer = nTri_Elem(nACO)*iCmp*iBas*jCmp*jBas*nDisp*nirrep
  else
    MemBuffer = 1  ! Dummy length
  end if
  if (MemBuffer > MemMax) then
     write(u6,*) 'DrvG2: MemBuffer > MemMax'
     write(u6,*) 'MemBuffer=',MemBuffer
     write(u6,*) 'MemMax=',MemMax
     write(u6,*) 'Increase MOLCAS_MEM!'
     call Abend()
  end if

  Buffer(1:MemBuffer)=>Sew_Scr(ipMem:ipMem+MemBuffer-1)
  Buffer(:)=Zero
  ipMOC = ipMem+MEMBUFFER

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  Post_Process = .false.
  do klSh=1,nijS
    ks = Ind_ij(1,klSh)
    ls = Ind_ij(2,klSh)

    A_int = TMax(iS,jS)*TMax(kS,lS)
    !write(u6,*) 'is,js,ks,ls=',is,js,ks,ls
    if (A_Int < CutInt) cycle

    !do kS=1,nSkal
    kCmp = iSD(2,kS)
    kShell = iSD(11,kS)
    kCnttp = iSD(13,kS)
    kCnt = iSD(14,kS)
    Coor(1:3,3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    Shik = iShell == kShell

    !  do lS=1,kS
    lCmp = iSD(2,lS)
    lShell = iSD(11,lS)
    lCnttp = iSD(13,lS)
    lCnt = iSD(14,lS)
    Coor(1:3,4) = dbsc(lCnttp)%Coor(1:3,lCnt)

    call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

    iCmpV(:) = iSD4( 2,:)

    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! The code is working in such away that the MO needs upper and lower
    ! triangular parts of ij kl but hessian needs only lower, check if the
    ! integralbatch is lower or upper!!

    lTri = iTri(iS,jS) >= iTri(kS,lS)
    if ((.not. lTri) .and. (nMethod /= RASSCF)) cycle
    lDot = (lTri .and. l_Hss)

    Shjl = jShell == lShell
    Shijij = Shik .and. Shjl

    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Allocate memory for zeta, eta, kappa, P and Q.
    ! Allocate also for Alpha, Beta , Gamma and Delta in expanded form.

    call Create_BraKet(iSD4(5,1)*iSD4(5,2),iSD4(5,3)*iSD4(5,4))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ijS = iTri(iShell,jShell)
    klS = iTri(kShell,lShell)
    ikS = iTri(iShell,kShell)
    ilS = iTri(iShell,lShell)
    jkS = iTri(jShell,kShell)
    jlS = iTri(jShell,lShell)

    if (ltri) then

      !----------------------------------------------------------------*

      ! Fix the 1st order density matrix

      ! Pick up pointers to desymmetrized 1st order density matrices.
      ! Observe that the desymmetrized 1st order density matrices
      ! follow the contraction index.

      ipTmp = 0
      ipTmp2 = 0
      if (lpick) then

        ipTmp = ipDijs
        if (nMethod == RASSCF) ipTmp2 = ipDijs2

#ifdef _SKIP_
        ipDij = ipOffD(1,ijS)
        mDCRij = ipOffD(2,ijS)
        nDij = ipOffD(3,ijS)
        if (nMethod == RASSCF) ipDij2 = ipOffDA(1,ijS)

        if (mDCRij*nDij /= 0) then
          ipDDij = ipTmp
          ipTmp = ipTmp+nDij*mDCRij
          if (nMethod == RASSCF) then
            ipDDij2 = ipTmp2
            ipTmp2 = ipTmp2+nDij*mDCRij
          end if
        else
          ipDDij = 1
        end if
#endif

        Call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,nr_of_Densities,nMethod, &
                       ipTmp2, ipDij2, ipDDij2)

        ipDkl = ipOffD(1,klS)
        mDCRkl = ipOffD(2,klS)
        nDkl = ipOffD(3,klS)
        if (nMethod == RASSCF) ipDkl2 = ipOffDA(1,klS)

        if (mDCRkl*nDkl /= 0) then
          ipDDkl = ipTmp
          ipTmp = ipTmp+nDkl*mDCRkl
          if (nMethod == RASSCF) then
            ipDDkl2 = ipTmp2
            ipTmp2 = ipTmp2+nDkl*mDCRkl
          end if
        else
          ipDDkl = 1
        end if

        ipDik = ipOffD(1,ikS)
        mDCRik = ipOffD(2,ikS)
        nDik = ipOffD(3,ikS)
        if (nMethod == RASSCF) ipDik2 = ipOffDA(1,ikS)

        if (mDCRik*nDik /= 0) then
          ipDDik = ipTmp
          ipTmp = ipTmp+nDik*mDCRik
          if (nMethod == RASSCF) then
            ipDDik2 = ipTmp2
            ipTmp2 = ipTmp2+nDik*mDCRik
          end if
        else
          ipDDik = 1
        end if

        ipDil = ipOffD(1,ilS)
        mDCRil = ipOffD(2,ilS)
        nDil = ipOffD(3,ilS)
        if (nMethod == RASSCF) ipDil2 = ipOffDA(1,ilS)

        if (mDCRil*nDil /= 0) then
          ipDDil = ipTmp
          ipTmp = ipTmp+nDil*mDCRil
          if (nMethod == RASSCF) then
            ipDDil2 = ipTmp2
            ipTmp2 = ipTmp2+nDil*mDCRil
          end if
        else
          ipDDil = 1
        end if

        ipDjk = ipOffD(1,jkS)
        mDCRjk = ipOffD(2,jkS)
        nDjk = ipOffD(3,jkS)
        if (nMethod == RASSCF) ipDjk2 = ipOffDA(1,jkS)

        if (mDCRjk*nDjk /= 0) then
          ipDDjk = ipTmp
          ipTmp = ipTmp+nDjk*mDCRjk
          if (nMethod == RASSCF) then
            ipDDjk2 = ipTmp2
            ipTmp2 = ipTmp2+nDjk*mDCRjk
          end if
        else
          ipDDjk = 1
        end if

        ipDjl = ipOffD(1,jlS)
        mDCRjl = ipOffD(2,jlS)
        nDjl = ipOffD(3,jlS)
        if (nMethod == RASSCF) ipDjl2 = ipOffDA(1,jlS)

        if (mDCRjl*nDjl /= 0) then
          ipDDjl = ipTmp
          ipTmp = ipTmp+nDjl*mDCRjl
          if (nMethod == RASSCF) then
            ipDDjl2 = ipTmp2
            ipTmp2 = ipTmp2+nDjl*mDCRjl
          end if
        else
          ipDDjl = 1
        end if

      end if  ! if (lpick) then
    end if  ! if (ltri) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute total size of the second order density matrix in SO basis.
    !
    !------------------------------------------------------------------*
    nSO = MemSO2_P(nSD,iSD4)
    ldot2 = ldot
    if (nSO == 0) ldot2 = .false.

    ! Compute memory request for the primitives.

    ider = 2
    if (.not. ldot2) iDer = 1
    call MemRg2(iSD4( 1,:),nRys,MemPrm,ider)

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

    call PSOAO2(nSO,MemPrm,MemMax,iFnc,nAco,Mem1,Mem2,Mem3,Mem4,MemX,MemPSO, &
                MemFck,nFT,memCMO2,MemFin,MemBuffer,nSD,iSD4)


    iBasi = iSD4(3,1)
    jBasj = iSD4(3,2)
    kBask = iSD4(3,3)
    lBasl = iSD4(3,4)

    iBsInc= iSD4(4,1)
    jBsInc= iSD4(4,2)
    kBsInc= iSD4(4,3)
    lBsInc= iSD4(4,4)


!todo
    iPrimi=iSD4( 5,1)
    jPrimj=iSD4( 5,2)
    kPrimk=iSD4( 5,3)
    lPriml=iSD4( 5,4)
    !------------------------------------------------------------------*
    !
    ! Loop over basis function if we do not have enough of memory to
    ! calculate them in one step.
    !
    !------------------------------------------------------------------*
    do iBasAO=1,iBasi,iBsInc
      iBasn = min(iBsInc,iBasi-iBasAO+1)
      iSD4( 8,1) = iBasAO-1
      iSD4(19,1) = iBasn


      !----------------------------------------------------------------*
      !
      ! Move appropriate portions of the desymmetrized 1st order density matrix.
      !
      !----------------------------------------------------------------*
      do jBasAO=1,jBasj,jBsInc
        jBasn = min(jBsInc,jBasj-jBasAO+1)
        iSD4( 8,2) = jBasAO-1
        iSD4(19,2) = jBasn


        if (lpick .and. (mDCRij /= 0)) then
          call Picky_inner(DeDe(ipDij),iBasi,jBasj,iPrimi*jPrimj,iCmpV(1)*iCmpV(2),mDCRij,iBasAO,iBasAO+iBasn-1,jBasAO, &
                           jBasAO+jBasn-1,DeDe(ipDDij))
          if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDij2),iBasi,jBasj,iPrimi*jPrimj,iCmpV(1)*iCmpV(2),mDCRij,iBasAO, &
                                                  iBasAO+iBasn-1,jBasAO,jBasAO+jBasn-1,DeDe2(ipDDij2))
        end if
        mDij = (iBasn*jBasn+1)*iCmpV(1)*iCmpV(2)+iPrimi*jPrimj+1

        do kBasAO=1,kBask,kBsInc
          kBasn = min(kBsInc,kBask-kBasAO+1)
          iSD4( 8,3) = kBasAO-1
          iSD4(19,3) = kBasn

          if (lpick .and. (mDCRik /= 0)) then
            call Picky_inner(DeDe(ipDik),iBasi,kBask,iPrimi*kPrimk,iCmpV(1)*iCmpV(3),mDCRik,iBasAO,iBasAO+iBasn-1,kBasAO, &
                             kBasAO+kBasn-1,DeDe(ipDDik))
            if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDik2),iBasi,kBask,iPrimi*kPrimk,iCmpV(1)*iCmpV(3),mDCRik,iBasAO, &
                                                    iBasAO+iBasn-1,kBasAO,kBasAO+kBasn-1,DeDe2(ipDDik2))
          end if
          mDik = (iBasn*kBasn+1)*iCmpV(1)*iCmpV(3)+iPrimi*kPrimk+1
          if (lpick .and. (mDCRjk /= 0)) then
            call Picky_inner(DeDe(ipDjk),jBasj,kBask,jPrimj*kPrimk,iCmpV(2)*iCmpV(3),mDCRjk,jBasAO,jBasAO+jBasn-1,kBasAO, &
                             kBasAO+kBasn-1,DeDe(ipDDjk))
            if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDjk2),jBasj,kBask,jPrimj*kPrimk,iCmpV(2)*iCmpV(3),mDCRjk,jBasAO, &
                                                    jBasAO+jBasn-1,kBasAO,kBasAO+kBasn-1,DeDe2(ipDDjk2))
          end if
          mDjk = (jBasn*kBasn+1)*iCmpV(2)*iCmpV(3)+jPrimj*kPrimk+1

          do lBasAO=1,lBasl,lBsInc
            lBasn = min(lBsInc,lBasl-lBasAO+1)
            iSD4( 8,4) = lBasAO-1
            iSD4(19,4) = lBasn

            if (lpick .and. (mDCRkl /= 0)) then
              call Picky_inner(DeDe(ipDkl),kBask,lBasl,kPrimk*lPriml,iCmpV(3)*iCmpV(4),mDCRkl,kBasAO,kBasAO+kBasn-1,lBasAO, &
                               lBasAO+lBasn-1,DeDe(ipDDkl))
              if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDkl2),kBask,lBasl,kPrimk*lPriml,iCmpV(3)*iCmpV(4),mDCRkl,kBasAO, &
                                                      kBasAO+kBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDkl2))
            end if
            mDkl = (kBasn*lBasn+1)*iCmpV(3)*iCmpV(4)+kPrimk*lPriml+1
            if (lpick .and. (mDCRil /= 0)) then
              call Picky_inner(DeDe(ipDil),iBasi,lBasl,iPrimi*lPriml,iCmpV(1)*iCmpV(4),mDCRil,iBasAO,iBasAO+iBasn-1,lBasAO, &
                               lBasAO+lBasn-1,DeDe(ipDDil))
              if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDil2),iBasi,lBasl,iPrimi*lPriml,iCmpV(1)*iCmpV(4),mDCRil,iBasAO, &
                                                      iBasAO+iBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDil2))
            end if
            mDil = (iBasn*lBasn+1)*iCmpV(1)*iCmpV(4)+iPrimi*lPriml+1
            if (lpick .and. (mDCRjl /= 0)) then
              call Picky_inner(DeDe(ipDjl),jBasj,lBasl,jPrimj*lPriml,iCmpV(2)*iCmpV(4),mDCRjl,jBasAO,jBasAO+jBasn-1,lBasAO, &
                               lBasAO+lBasn-1,DeDe(ipDDjl))
              if (nMethod == RASSCF) call Picky_inner(DeDe2(ipDjl2),jBasj,lBasl,jPrimj*lPriml,iCmpV(2)*iCmpV(4),mDCRjl,jBasAO, &
                                                      jBasAO+jBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDjl2))
            end if
            mDjl = (jBasn*lBasn+1)*iCmpV(2)*iCmpV(4)+jPrimj*lPriml+1

            !----------------------------------------------------------*

            nijkl = iBasn*jBasn*kBasn*lBasn

            ! Mark out the memory allocations explicitly with pointers
            ! MO tranformation buffer
            MEMCMO = nACO*(kCmp*kBasn+lCmp*lBasn)
            MOC(1:MemCMO)=>Sew_Scr(ipMOC:ipMOC+MemCMO-1)
            ! Area for the AO integrals
            ipFin = ipMOC+MemCMO
            Fin(1:MemFin)=>Sew_Scr(ipFin:ipFin+MemFin-1)
            ! Area for 2el density
            ipPSO = ipFin+MemFin
            PSO(1:nijkl,1:nSO)=>Sew_Scr(ipPSO:ipPSO+nijkl*nSO-1)
            If (nijkl*nSO>Mem1) Then
               Write(u6,'(A)') 'nijkl*nSO>Mem1'
               Write(u6,*) 'njikl,nSO=',nijkl,nSO
               Write(u6,*) 'Mem1=',Mem1
               Call Abend()
            End If
            ipMem2 = ipPSO+Mem1  ! Work
            Work2(1:Mem2)=>Sew_Scr(ipMem2:ipMem2+Mem2-1)
            ipMem3 = ipMem2+Mem2 ! Work
            Work3(1:Mem3)=>Sew_Scr(ipMem3:ipMem3+Mem3-1)
            ipMemX = ipMem3+Mem3 ! Work
            WorkX(1:MemX)=>Sew_Scr(ipMemX:ipMemX+MemX-1)

            ! If MO transformation is performed in the standard way
            ! reserve memory for partial transformed integrals

            ! Multilayer

            ipMem4 = ipMem2+Mem2-Mem4
            Work4(1:Mem4)=>Sew_Scr(ipMem4:ipMem4+Mem4-1)
            nTemp=Mem2+Mem3+MemX
            Temp(1:nTemp)=>Sew_Scr(ipMem2:ipMem2+nTemp-1)

            !----------------------------------------------------------*
            !
            ! Get the 2nd order density matrix in SO basis.
            !
            !----------------------------------------------------------*

            call Timing(dum1,Time,dum2,dum3)
            if (n8) call PickMO(MOC,MemCMO,nSD,iSD4)
            if (ldot2) call PGet0(nijkl,PSO,nSO,iFnc,MemPSO,Work2,Mem2,nQuad,PMax,iSD4)
            call Timing(dum1,Time,dum2,dum3)
            CPUStat(nTwoDens) = CPUStat(nTwoDens)+Time

            ! Compute gradients of shell quadruplet

            call TwoEl_mck(Coor,nRys,Pren,Prem, &
                           Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,JfG,PSO,nijkl,nSO, &
                           Work2,Mem2,Work3,Mem3,Work4,Mem4,Aux,nAux,WorkX,MemX, &
                           Shijij,DeDe(ipDDij),DeDe2(ipDDij2),mDij,mDCRij,DeDe(ipDDkl),DeDe2(ipDDkl2),mDkl,mDCRkl,DeDe(ipDDik), &
                           DeDe2(ipDDik2),mDik,mDCRik,DeDe(ipDDil),DeDe2(ipDDil2),mDil,mDCRil,DeDe(ipDDjk),DeDe2(ipDDjk2),mDjk, &
                           mDCRjk,DeDe(ipDDjl),DeDe2(ipDDjl2),mDjl,mDCRjl,iCmpV,Fin,MemFin,Temp, &
                           nTemp,nTwo2,nFT,iInt,Buffer,MemBuffer,lgrad,ldot2,n8,ltri,DTemp,DInAc,moip,nAco, &
                           MOC,MemCMO,new_fock,iSD4)
            Post_Process = .true.
            MOC=>Null()
            Fin=>Null()
            PSO=>Null()
            Work2=>Null()
            Work3=>Null()
            WorkX=>Null()
            Work4=>Null()
            Temp=>Null()

            !----------------------------------------------------------*

          end do
        end do
      end do
    end do
    call Destroy_Braket()

    !  end do ! lS
    !end do ! kS
  end do ! klS

  if ((nMethod == RASSCF) .and. Post_Process) then
    nTemp=MemMax-MemBuffer
    Temp(1:nTemp)=>Sew_Scr(ipMOC:ipMOC+nTemp-1)
    call CLR2(Buffer,iInt,nACO,nSD,iSD4,nDisp,nTemp,Temp)
    Temp=>Null()
  end if
  Buffer=>Null()

  !  end do ! jS
  !end do ! iS

  call CWTime(TCpu2,TWall2)
end do
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

return

end subroutine Drvg2
