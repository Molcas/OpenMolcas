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
use iSD_data, only: iSD
use k2_arrays, only: Aux, Create_Braket, Create_BraKet_Base, DeDe, Destroy_Braket, Destroy_BraKet_Base, ipDijS, ipOffD, MxDij, &
                     ndede, nFT, Sew_Scr
use k2_structure, only: Indk2, k2Data
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
logical(kind=iwp), intent(in) :: l_Grd, l_Hss
integer(kind=iwp) :: i, iAng, iAngV(4), iAO, iAOst(4), iAOV(4), iBas, iBasAO, ibasI, iBasn, iBsInc, iCmp, iCmpV(4), iCnt, iCnttp, &
                     id, id_Tsk, idd, ider, iDisk, iDisp, iFnc(4), iii, iIrr, iIrrep, ij, ijMax, ijS, ijSh, ik2, ikS, ilS, iMemB, &
                     ip, ip1, ip2, ip3, ip4, ip5, ip6, ip_PP, ipBuffer, ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, &
                     ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, ipDij, ipDij2, ipDijS2, ipDik, ipDik2, ipDil, ipDil2, ipDjk, &
                     ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, ipFin, ipMem, ipMem2, ipMem3, ipMem4, ipMemX, ipMOC, iPrim, iPrimi, &
                     iPrInc, ipTmp, ipTmp2, iS, iShell, iShelV(4), iShll, iShllV(4), j, jAng, jAO, jBas, jBasAO, jBasj, jBasn, &
                     jBsInc, jCmp, jCnt, jCnttp, jDisp, jIrr, jk2, jkS, jlS, JndGrd(3,4,0:7), JndHss(4,3,4,3,0:7), jPrimj, jPrInc, &
                     js, jShell, jShll, kAng, kAO, kBasAO, kBask, kBasn, kBsInc, kCmp, kCnt, kCnttp, kIrr, klS, klSh, kPrimk, &
                     kPrInc, ks, kShell, kShll, lAng, lAO, lBasAO, lBasl, lBasn, lBsInc, lCmp, lCnt, lCnttp, lPriml, lPrInc, ls, &
                     lShell, lShll, mdci, mdcj, mdck, mdcl, mDCRij, mDCRik, mDCRil, mDCRjk, mDCRjl, mDCRkl, mDeDe, mDij, mDik, &
                     mDil, mDjk, mDjl, mDkl, Mem1, Mem2, Mem3, Mem4, MemBuffer, MEMCMO, memCMO2, MemFck, MemFin, MemMax, MemPrm, &
                     MemPSO, MemX, mIndij, mmdede, moip(0:7), MxBsC, n_Int, nAco, nb, nDCRR, nDCRS, nDij, nDik, nDil, ndisp, nDjk, &
                     nDjl, nDkl, nEta, nHrrab, nHrrcd, nijkl, nijS, nIndij, nMO, nPairs, nQuad, nRys, nSkal, nSO, nTwo, nTwo2, nZeta
real(kind=wp) :: A_int, dum1, dum2, dum3, Coor(3,4), PMax, Prem, Pren, TCpu1, TCpu2, Time, TMax_all, TWall1, TWall2
logical(kind=iwp) :: JfG(4), JfGrd(3,4), JfHss(4,3,4,3), ldot, ldot2, lGrad, lpick, ltri, n8, new_fock, Post_Process, Shijij, &
                     Shik, Shjl
#ifdef _DEBUGPRINT_
character(len=40) :: frmt
#endif
logical(kind=iwp), parameter :: Int_Direct = .true.
integer(kind=iwp), allocatable :: Ind_ij(:,:), ipOffDA(:,:)
real(kind=wp), allocatable :: DeDe2(:), DInAc(:), DTemp(:), iInt(:), TMax(:,:)
integer(kind=iwp), external :: MemSO2_P, NrOpr
logical(kind=iwp), external :: Rsv_Tsk

!                                                                      *
!***********************************************************************
!                                                                      *
! PROLOGUE

call StatusLine(' McKinley:',' Computing 2-electron 2nd order derivatives')

ipDij = 0
ipDij2 = 0
ipDDij = 0
ipDDij2 = 0
ipDkl = 0
ipDkl2 = 0
ipDDkl = 0
ipDDkl2 = 0
ipDik = 0
ipDik2 = 0
ipDDik = 0
ipDDik2 = 0
ipDil = 0
ipDil2 = 0
ipDDil = 0
ipDDil2 = 0
ipDjk = 0
ipDjk2 = 0
ipDDjk = 0
ipDDjk2 = 0
ipDjl = 0
ipDjl2 = 0
ipDDjl = 0
ipDDjl2 = 0
ipBuffer = 0
ipMOC = 0
iFnc(1) = -99
iFnc(2) = -99
iFnc(3) = -99
iFnc(4) = -99
nDij = 0
nDkl = 0
nDik = 0
nDjl = 0
nDil = 0
nDjk = 0
mDCRij = 0
mDCRkl = 0
mDCRik = 0
mDCRjl = 0
mDCRil = 0
mDCRjk = 0
ipDijS = 0
ipDijS2 = 0

call CtrlMO(moip,nAco)

ndisp = 0
naco = 0
New_Fock = nirrep == 1
do iS=0,nIrrep-1
  nDisp = nDisp+ldisp(is)
  naco = naco+nAsh(is)
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
nZeta = MxPrm*MxPrm
nEta = MxPrm*MxPrm
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

  if (.not. allocated(DeDe)) call mma_allocate(DeDe,[-1,-1],label='DeDe') ! Dummy allocation
  if (.not. allocated(DeDe2)) call mma_allocate(DeDe2,[-1,-1],label='DeDe2') ! Dummy allocation

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
if (MemMax > 1000) MemMax = MemMax-1000
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
  !do iS=1,nSkal
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  mdci = iSD(10,iS)
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  Coor(1:3,1) = dbsc(iCnttp)%Coor(1:3,iCnt)

  iAngV(1) = iAng
  iShllV(1) = iShll
  iCmpV(1) = iCmp
  iShelV(1) = iShell
  iAOV(1) = iAO

  !  do jS=1,iS
  jShll = iSD(0,jS)
  jAng = iSD(1,jS)
  jCmp = iSD(2,jS)
  jBas = iSD(3,jS)
  jAO = iSD(7,jS)
  mdcj = iSD(10,jS)
  jShell = iSD(11,jS)
  jCnttp = iSD(13,jS)
  jCnt = iSD(14,jS)
  Coor(1:3,2) = dbsc(jCnttp)%Coor(1:3,jCnt)

  iAngV(2) = jAng
  iShllV(2) = jShll
  iCmpV(2) = jCmp
  iShelV(2) = jShell
  iAOV(2) = jAO

  nHrrab = 0
  do i=0,iAng+1
    do j=0,jAng+1
      if (i+j <= iAng+jAng+1) then
        ijMax = min(iAng,jAng)+1
        nHrrab = nHrrab+ijMax*2+1
      end if
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Cltrls for MO transformation
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if ((nMethod == RASSCF) .and. l_Grd) then
    iMemB = nACO**2*iCmp*iBas*jCmp*jBas*nDisp*nirrep
    if (iMemB > MemMax) then
      write(u6,*) 'DrvG2: iMemB > MemMax'
      write(u6,*) 'iMemB=',iMemB
      write(u6,*) 'MemMax=',MemMax
      write(u6,*) 'Increase MOLCAS_MEM!'
      call Abend()
    end if
    Sew_Scr(1:iMemb) = Zero
  else
    iMemb = 0
  end if
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
    kShll = iSD(0,kS)
    kAng = iSD(1,kS)
    kCmp = iSD(2,kS)
    kAO = iSD(7,kS)
    mdck = iSD(10,kS)
    kShell = iSD(11,kS)
    kCnttp = iSD(13,kS)
    kCnt = iSD(14,kS)
    Coor(1:3,3) = dbsc(kCnttp)%Coor(1:3,kCnt)

    iAngV(3) = kAng
    iShllV(3) = kShll
    iCmpV(3) = kCmp
    iShelV(3) = kShell
    iAOV(3) = kAO

    Shik = iShell == kShell

    !  do lS=1,kS
    lShll = iSD(0,lS)
    lAng = iSD(1,lS)
    lCmp = iSD(2,lS)
    lAO = iSD(7,lS)
    mdcl = iSD(10,lS)
    lShell = iSD(11,lS)
    lCnttp = iSD(13,lS)
    lCnt = iSD(14,lS)
    Coor(1:3,4) = dbsc(lCnttp)%Coor(1:3,lCnt)

    iAngV(4) = lAng
    iShllV(4) = lShll
    iCmpV(4) = lCmp
    iShelV(4) = lShell
    iAOV(4) = lAO

    nHrrcd = 0
    do i=0,kAng+1
      do j=0,lAng+1
        if (i+j <= kAng+lAng+1) then
          ijMax = min(kAng,lAng)+1
          nHrrcd = nHrrcd+ijMax*2+1
        end if
      end do
    end do
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
    iCmpV(1) = icmp
    iCmpV(2) = jcmp
    iCmpV(3) = kcmp
    iCmpV(4) = lcmp
    iPrimi = Shells(iShllV(1))%nExp
    jPrimj = Shells(iShllV(2))%nExp
    kPrimk = Shells(iShllV(3))%nExp
    lPriml = Shells(iShllV(4))%nExp
    iBasi = Shells(iShllV(1))%nBasis
    jBasj = Shells(iShllV(2))%nBasis
    kBask = Shells(iShllV(3))%nBasis
    lBasl = Shells(iShllV(4))%nBasis
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Allocate memory for zeta, eta, kappa, P and Q.
    ! Allocate also for Alpha, Beta , Gamma and Delta in expanded form.

    nZeta = iPrimi*jPrimj
    nEta = kPrimk*lPriml

    call Create_BraKet(nZeta,nEta)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ijS = iTri(iShell,jShell)
    klS = iTri(kShell,lShell)
    ikS = iTri(iShell,kShell)
    ilS = iTri(iShell,lShell)
    jkS = iTri(jShell,kShell)
    jlS = iTri(jShell,lShell)
    nDCRR = Indk2(2,ijS)
    ik2 = Indk2(3,ijS)
    nDCRS = Indk2(2,klS)
    jk2 = Indk2(3,klS)

    if (ltri) then

      !----------------------------------------------------------------*

      ! Fix the 1st order density matrix

      ! Pick up pointers to desymmetrized 1st order density matrices.
      ! Observe that the desymmetrized 1st order density matrices
      ! follow the contraction index.

      ipTmp = 0
      ipTmp2 = 0
      if (lpick) then

        ipDij = ipOffD(1,ijS)
        mDCRij = ipOffD(2,ijS)
        nDij = ipOffD(3,ijS)

        ipTmp = ipDijs
        if (nMethod == RASSCF) then
          ipDij2 = ipOffDA(1,ijS)
          ipTmp2 = ipDijs2
        end if

        if (mDCRij /= 0) then
          ipDDij = ipTmp
          ipTmp = ipTmp+nDij*mDCRij
          if (nMethod == RASSCF) then
            ipDDij2 = ipTmp2
            ipTmp2 = ipTmp2+nDij*mDCRij
          end if
        else
          ipDDij = 0
        end if

        ipDkl = ipOffD(1,klS)
        if (nMethod == RASSCF) ipDkl2 = ipOffDA(1,klS)
        mDCRkl = ipOffD(2,klS)
        nDkl = ipOffD(3,klS)
        if (mDCRkl /= 0) then
          ipDDkl = ipTmp
          ipTmp = ipTmp+nDkl*mDCRkl
          if (nMethod == RASSCF) then
            ipDDkl2 = ipTmp2
            ipTmp2 = ipTmp2+nDkl*mDCRkl
          end if
        else
          ipDDkl = 0
        end if

        ipDik = ipOffD(1,ikS)
        if (nMethod == RASSCF) ipDik2 = ipOffDA(1,ikS)
        mDCRik = ipOffD(2,ikS)
        nDik = ipOffD(3,ikS)
        if (mDCRik /= 0) then
          ipDDik = ipTmp
          ipTmp = ipTmp+nDik*mDCRik
          if (nMethod == RASSCF) then
            ipDDik2 = ipTmp2
            ipTmp2 = ipTmp2+nDik*mDCRik
          end if
        else
          ipDDik = 0
        end if

        ipDil = ipOffD(1,ilS)
        if (nMethod == RASSCF) ipDil2 = ipOffDA(1,ilS)
        mDCRil = ipOffD(2,ilS)
        nDil = ipOffD(3,ilS)
        if (mDCRil /= 0) then
          ipDDil = ipTmp
          ipTmp = ipTmp+nDil*mDCRil
          if (nMethod == RASSCF) then
            ipDDil2 = ipTmp2
            ipTmp2 = ipTmp2+nDil*mDCRil
          end if
        else
          ipDDil = 0
        end if

        ipDjk = ipOffD(1,jkS)
        if (nMethod == RASSCF) ipDjk2 = ipOffDA(1,jkS)
        mDCRjk = ipOffD(2,jkS)
        nDjk = ipOffD(3,jkS)
        if (mDCRjk /= 0) then
          ipDDjk = ipTmp
          ipTmp = ipTmp+nDjk*mDCRjk
          if (nMethod == RASSCF) then
            ipDDjk2 = ipTmp2
            ipTmp2 = ipTmp2+nDjk*mDCRjk
          end if
        else
          ipDDjk = 0
        end if

        ipDjl = ipOffD(1,jlS)
        if (nMethod == RASSCF) ipDjl2 = ipOffDA(1,jlS)
        mDCRjl = ipOffD(2,jlS)
        nDjl = ipOffD(3,jlS)
        if (mDCRjl /= 0) then
          ipDDjl = ipTmp
          ipTmp = ipTmp+nDjl*mDCRjl
          if (nMethod == RASSCF) then
            ipDDjl2 = ipTmp2
            ipTmp2 = ipTmp2+nDjl*mDCRjl
          end if
        else
          ipDDjl = 0
        end if

      end if  ! if (lpick) then
    end if  ! if (ltri) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute total size of the second order density matrix in SO basis.
    !
    !------------------------------------------------------------------*
    nSO = MemSO2_P(iCmp,jCmp,kCmp,lCmp,iAOV(1),iAOV(2),iAOV(3),iAOV(4))
    ldot2 = ldot
    if (nSO == 0) ldot2 = .false.

    ! Compute memory request for the primitives.

    ider = 2
    if (.not. ldot2) iDer = 1
    call MemRg2(iAngV,nRys,MemPrm,ider)

    !------------------------------------------------------------------*
    !
    ! Calculate which derivatives should be made.
    !
    !------------------------------------------------------------------*

    call DerCtr(mdci,mdcj,mdck,mdcl,ldot2,JfGrd,JndGrd,JfHss,JndHss,JfG)

    !------------------------------------------------------------------*
    !
    ! Decide on the partioning of the shells based on the
    ! available memory and the requested memory.
    !
    !------------------------------------------------------------------*

    call PSOAO2(nSO,MemPrm,MemMax,iAngV,iCmpV,iAOV,iFnc,iBasi,iBsInc,jBasj,jBsInc,kBask,kBsInc,lBasl,lBsInc,iPrimi,iPrInc,jPrimj, &
                jPrInc,kPrimk,kPrInc,lPriml,lPrInc,nAco,Mem1,Mem2,Mem3,Mem4,MemX,MemPSO,MemFck,nFT,memCMO2,MemFin,MemBuffer,iMemB)

    !------------------------------------------------------------------*
    !
    ! Loop over basis function if we do not have enough of memory to
    ! calculate them in one step.
    !
    !------------------------------------------------------------------*
    do iBasAO=1,iBasi,iBsInc
      iBasn = min(iBsInc,iBasi-iBasAO+1)
      iAOst(1) = iBasAO-1
      !----------------------------------------------------------------*
      !
      ! Move appropriate portions of the desymmetrized 1st order density matrix.
      !
      !----------------------------------------------------------------*
      do jBasAO=1,jBasj,jBsInc
        jBasn = min(jBsInc,jBasj-jBasAO+1)
        iAOst(2) = jBasAO-1
        if (lpick .and. (nDij*mDCRij /= 0)) then
          call Picky(DeDe(ipDij),iBasi,jBasj,iPrimi*jPrimj,iCmpV(1)*iCmpV(2),mDCRij,iBasAO,iBasAO+iBasn-1,jBasAO,jBasAO+jBasn-1, &
                     DeDe(ipDDij))
          if (nMethod == RASSCF) call Picky(DeDe2(ipDij2),iBasi,jBasj,iPrimi*jPrimj,iCmpV(1)*iCmpV(2),mDCRij,iBasAO, &
                                            iBasAO+iBasn-1,jBasAO,jBasAO+jBasn-1,DeDe2(ipDDij2))
        end if
        mDij = (iBasn*jBasn+1)*iCmpV(1)*iCmpV(2)+iPrimi*jPrimj+1
        mDij = min(nDij,mDij)

        do kBasAO=1,kBask,kBsInc
          kBasn = min(kBsInc,kBask-kBasAO+1)
          iAOst(3) = kBasAO-1
          if (lpick .and. (nDik*mDCRik /= 0)) then
            call Picky(DeDe(ipDik),iBasi,kBask,iPrimi*kPrimk,iCmpV(1)*iCmpV(3),mDCRik,iBasAO,iBasAO+iBasn-1,kBasAO,kBasAO+kBasn-1, &
                       DeDe(ipDDik))
            if (nMethod == RASSCF) call Picky(DeDe2(ipDik2),iBasi,kBask,iPrimi*kPrimk,iCmpV(1)*iCmpV(3),mDCRik,iBasAO, &
                                              iBasAO+iBasn-1,kBasAO,kBasAO+kBasn-1,DeDe2(ipDDik2))
          end if
          mDik = (iBasn*kBasn+1)*iCmpV(1)*iCmpV(3)+iPrimi*kPrimk+1
          mDik = min(nDik,mDik)
          if (lpick .and. (nDjk*mDCRjk /= 0)) then
            call Picky(DeDe(ipDjk),jBasj,kBask,jPrimj*kPrimk,iCmpV(2)*iCmpV(3),mDCRjk,jBasAO,jBasAO+jBasn-1,kBasAO,kBasAO+kBasn-1, &
                       DeDe(ipDDjk))
            if (nMethod == RASSCF) call Picky(DeDe2(ipDjk2),jBasj,kBask,jPrimj*kPrimk,iCmpV(2)*iCmpV(3),mDCRjk,jBasAO, &
                                              jBasAO+jBasn-1,kBasAO,kBasAO+kBasn-1,DeDe2(ipDDjk2))
          end if
          mDjk = (jBasn*kBasn+1)*iCmpV(2)*iCmpV(3)+jPrimj*kPrimk+1
          mDjk = min(nDjk,mDjk)

          do lBasAO=1,lBasl,lBsInc
            lBasn = min(lBsInc,lBasl-lBasAO+1)
            iAOst(4) = lBasAO-1
            if (lpick .and. (nDkl*mDCRkl /= 0)) then
              call Picky(DeDe(ipDkl),kBask,lBasl,kPrimk*lPriml,iCmpV(3)*iCmpV(4),mDCRkl,kBasAO,kBasAO+kBasn-1,lBasAO, &
                         lBasAO+lBasn-1,DeDe(ipDDkl))
              if (nMethod == RASSCF) call Picky(DeDe2(ipDkl2),kBask,lBasl,kPrimk*lPriml,iCmpV(3)*iCmpV(4),mDCRkl,kBasAO, &
                                                kBasAO+kBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDkl2))
            end if
            mDkl = (kBasn*lBasn+1)*iCmpV(3)*iCmpV(4)+kPrimk*lPriml+1
            mDkl = min(nDkl,mDkl)
            if (lpick .and. (nDil*mDCRil /= 0)) then
              call Picky(DeDe(ipDil),iBasi,lBasl,iPrimi*lPriml,iCmpV(1)*iCmpV(4),mDCRil,iBasAO,iBasAO+iBasn-1,lBasAO, &
                         lBasAO+lBasn-1,DeDe(ipDDil))
              if (nMethod == RASSCF) call Picky(DeDe2(ipDil2),iBasi,lBasl,iPrimi*lPriml,iCmpV(1)*iCmpV(4),mDCRil,iBasAO, &
                                                iBasAO+iBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDil2))
            end if
            mDil = (iBasn*lBasn+1)*iCmpV(1)*iCmpV(4)+iPrimi*lPriml+1
            mDil = min(nDil,mDil)
            if (lpick .and. (nDjl*mDCRjl /= 0)) then
              call Picky(DeDe(ipDjl),jBasj,lBasl,jPrimj*lPriml,iCmpV(2)*iCmpV(4),mDCRjl,jBasAO,jBasAO+jBasn-1,lBasAO, &
                         lBasAO+lBasn-1,DeDe(ipDDjl))
              if (nMethod == RASSCF) call Picky(DeDe2(ipDjl2),jBasj,lBasl,jPrimj*lPriml,iCmpV(2)*iCmpV(4),mDCRjl,jBasAO, &
                                                jBasAO+jBasn-1,lBasAO,lBasAO+lBasn-1,DeDe2(ipDDjl2))
            end if
            mDjl = (jBasn*lBasn+1)*iCmpV(2)*iCmpV(4)+jPrimj*lPriml+1
            mDjl = min(nDjl,mDjl)
            if (.not. lpick) then
              ipddjl2 = 0
              ipddil2 = 0
              ipddkl2 = 0
              ipddij2 = 0
              ipddik2 = 0
              ipddjk2 = 0
            end if

            !----------------------------------------------------------*

            MEMCMO = nACO*(kCmp*kBasn+lCmp*lBasn)
            ! MO tranformation buffer
            ipBuffer = ipMem
            ipMOC = ipBuffer+MEMBUFFER
            ! Area for the AO integrals
            ipFin = ipMOC+MemCMO
            ! Area for 2el density
            ip_PP = ipFin+MemFin
            ipMem2 = ip_PP+Mem1  ! Work
            ipMem3 = ipMem2+Mem2 ! Work
            ipMemX = ipMem3+Mem3 ! Work

            ! If MO transformation is performed in the standard way
            ! reserve memory for partial transfromed integrals

            ! Multilayer

            ipMem4 = ipMem2+Mem2-Mem4

            !----------------------------------------------------------*
            !
            ! Get the 2nd order density matrix in SO basis.
            !
            !----------------------------------------------------------*

            nijkl = iBasn*jBasn*kBasn*lBasn
            call Timing(dum1,Time,dum2,dum3)
            if (n8) call PickMO(Sew_Scr(ipMOC),MemCMO,iCmpV,iBasAO,iBasn,jBasAO,jBasn,kBasAO,kBasn,lBasAO,lBasn,iAOV)
            if (ldot2) call PGet0(iCmpV,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,Sew_Scr(ip_PP),nSO,iFnc(1)*iBasn, &
                                  iFnc(2)*jBasn,iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
            call Timing(dum1,Time,dum2,dum3)
            CPUStat(nTwoDens) = CPUStat(nTwoDens)+Time

            ! Compute gradients of shell quadruplet

            call TwoEl_mck(Coor,iAngV,iCmpV,iShelV,iShllV,iAOV,iAOst,mdci,mdcj,mdck,mdcl,nRys,nDCRR,nDCRS, &
                           k2data(:,ik2),k2data(:,jk2), &
                           Pren,Prem,iPrimi,jPrimj,jPrInc,kPrimk,lPriml,lPrInc, &
                           Shells(iShllV(1))%pCff(1,iBasAO),iBasn,Shells(iShllV(2))%pCff(1,jBasAO),jBasn, &
                           Shells(iShllV(3))%pCff(1,kBasAO),kBasn,Shells(iShllV(4))%pCff(1,lBasAO),lBasn, &
                           nZeta,nEta,Hess,nHess,JfGrd,JndGrd,JfHss,JndHss,JfG,Sew_Scr(ip_PP),nSO, &
                           Sew_Scr(ipMem2),Mem2,Sew_Scr(ipMem3),Mem3,Sew_Scr(ipMem4),Mem4,Aux,nAux,Sew_Scr(ipMemX),MemX, &
                           Shijij,DeDe(ipDDij),DeDe2(ipDDij2),mDij,mDCRij,DeDe(ipDDkl),DeDe2(ipDDkl2),mDkl,mDCRkl,DeDe(ipDDik), &
                           DeDe2(ipDDik2),mDik,mDCRik,DeDe(ipDDil),DeDe2(ipDDil2),mDil,mDCRil,DeDe(ipDDjk),DeDe2(ipDDjk2),mDjk, &
                           mDCRjk,DeDe(ipDDjl),DeDe2(ipDDjl2),mDjl,mDCRjl,iCmpV,Sew_Scr(ipFin),MemFin,Sew_Scr(ipMem2), &
                           Mem2+Mem3+MemX,nTwo2,nFT,iInt,Sew_Scr(ipBuffer),MemBuffer,lgrad,ldot2,n8,ltri,DTemp,DInAc,moip,nAco, &
                           Sew_Scr(ipMOC),MemCMO,new_fock)
            Post_Process = .true.

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
    ip1 = ipMOC
    ip2 = ip1+iCmp*iBas*naco
    ip3 = ip2+nAco**2
    ip4 = ip3+jcmp*jBas*naco
    ip5 = ip4+iCmp*naco*iBas
    ip6 = ip5+jcmp*jbas*naco
    call CLR2(Sew_Scr(ipBuffer),iInt,ibas,icmp,jbas,jcmp,iAOV(1),iAOV(2),naco,ishelV,Sew_Scr(ip1),Sew_Scr(ip2),Sew_Scr(ip3), &
              Sew_Scr(ip4),Sew_Scr(ip5),Sew_Scr(ip6))
  end if

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

if (allocated(ipDisp)) call mma_deallocate(ipDisp)
if (allocated(ipDisp2)) call mma_deallocate(ipDisp2)
if (allocated(ipDisp3)) call mma_deallocate(ipDisp3)
if (allocated(ipMO)) call mma_deallocate(ipMO)

return

end subroutine Drvg2
