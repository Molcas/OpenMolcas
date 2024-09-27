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
! Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
!               1990, IBM                                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine Drv2El_dscf(Dens,TwoHam,nDens,nDisc,FstItr)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will control the type of the two-electron integral, e.g.    *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
!          Dens is the folded lower triangular of the 1st order        *
!               density matrix.                                        *
!          Twoham is the lower triangular of the two-electron contri-  *
!               bution to the Fock matrix.                             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for direct SCF. January '93                     *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified by M.Schuetz @teokem.lu.se :                    *
!             parallel region split off in drvtwo.f, April '95         *
!             Modified by R. Lindh  @teokem.lu.se :                    *
!             total repacking of code September '96                    *
!***********************************************************************

use IOBUF, only: lBuf
use Gateway_Info, only: CutInt, ThrInt
use RICD_Info, only: Do_DCCD
use iSD_data, only: iSD
use Integral_Interfaces, only: DeDe_SCF, Int_PostProcess, int_wrout
use Int_Options, only: Disc, Disc_Mx, DoFock, DoIntegrals, Exfac, FckNoClmb, FckNoExch, Quad_ijkl, Thize, W2Disc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Eight, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDens, nDisc
real(kind=wp), target, intent(inout) :: Dens(nDens), TwoHam(nDens)
logical(kind=iwp), intent(inout) :: FstItr
integer(kind=iwp), parameter :: nTInt = 1
integer(kind=iwp) :: ijS, iOpt, iS, jS, klS, kS, lS, mDens, nIJ, nSkal
real(kind=wp) :: A_Int, Dtst, P_Eff, PP_Count, PP_Eff, PP_Eff_Delta, S_Eff, ST_Eff, T_Eff, TCPU1, TCPU2, ThrAO, TInt(nTInt), &
                 TMax_All, TskHi, TskLw, TWALL1, TWALL2
logical(kind=iwp) :: DoGrad, Indexation, Semi_Direct, Skip, Triangular
character(len=72) :: SLine
integer(kind=iwp), allocatable :: ip_ij(:,:)
real(kind=wp), allocatable :: DMax(:,:), TMax(:,:)
procedure(int_wrout) :: No_Routine
logical(kind=iwp), external :: Rsv_GTList

!                                                                      *
!***********************************************************************
!                                                                      *
SLine = 'Computing 2-electron integrals'
call StatusLine(' SCF:',SLine)
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
Int_PostProcess => No_Routine
!                                                                      *
!***********************************************************************
!                                                                      *
! Set variables in module Int_Options
DoIntegrals = .false.
DoFock = .true.
FckNoExch = ExFac == Zero
W2Disc = .false.     ! Default value
! Disc_Mx = file size in Real*8 128=1024/8
Disc_Mx = real(nDisc,kind=wp)*128.0_wp
! Subtract for the last buffer
Disc_Mx = Disc_Mx-lBuf
Disc = Zero        ! Default value
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up for partial SO/AO integral storage.

! nDisc = file size in kbyte from input
Semi_Direct = (nDisc /= 0)
if (Semi_Direct) call Init_SemiDSCF(FstItr,Thize,Cutint)
!                                                                      *
!***********************************************************************
!                                                                      *
! Desymmetrize differential densities.
! Observe that the desymmetrized 1st order density matrices are
! canonical, i.e. the relative order of the indices are canonically
! ordered.

call DeDe_SCF(Dens,TwoHam,nDens,mDens)
!                                                                      *
!***********************************************************************
!                                                                      *
Indexation = .false.
ThrAO = Zero           ! Do not modify CutInt
DoGrad = .false.

call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
!                                                                      *
!***********************************************************************
!                                                                      *
TskHi = Zero
TskLw = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    if (Do_DCCD .and. (iSD(10,iS) /= iSD(10,jS))) cycle
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
call mma_allocate(DMax,nSkal,nSkal,Label='DMax')
call Shell_MxDens(Dens,DMax,nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(ip_ij,2,nSkal*(nSkal+1),Label='ip_ij')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if (Do_DCCD .and. (iSD(10,iS) /= iSD(10,jS))) cycle
    if (TMax_All*TMax(iS,jS) >= CutInt) then
      nij = nij+1
      ip_ij(1,nij) = iS
      ip_ij(2,nij) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)

PP_Eff = P_Eff**2
PP_Eff_delta = 0.1_wp*PP_Eff
PP_Count = Zero

!                                                                      *
!***********************************************************************
!                                                                      *
! For distributed parallel SCF initiate (sequential code is special
! case when the number of nodes in the mpp is 1).
!
! 1: Task list (tlist)
! 2: Private priority list (pplist)
! 3: Global task list (gtlist)

if (FstItr) then
  Triangular = .true.
  call Init_TList(Triangular,P_Eff)
  call Init_PPList()
  call Init_GTList()
else
  call ReInit_PPList(Semi_Direct)
  call ReInit_GTList()
end if
iOpt = 0
if ((.not. FstItr) .and. Semi_direct) iOpt = 2

call CWTime(TCpu1,TWall1)

! big loop over individual tasks, distributed over individual nodes

do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.

  if (.not. Rsv_GTList(TskLw,TskHi,iOpt,W2Disc)) exit

  call Mode_SemiDSCF(W2Disc)
  !write(u6,*) 'TskLw,TskHi,W2Disc=',TskLw,TskHi,W2Disc

  ! Now do a quadruple loop over shells

  ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
  iS = ip_ij(1,ijS)
  jS = ip_ij(2,ijS)
  klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
  kS = ip_ij(1,klS)
  lS = ip_ij(2,klS)
  Quad_ijkl = TskLw

  Skip = (Quad_ijkl-TskHi > 1.0e-10_wp) ! Cut off check

  do
    if (Skip) exit
    ! What are these variables
    S_Eff = real(ijS,kind=wp)
    T_Eff = real(klS,kind=wp)
    ST_Eff = S_Eff*(S_Eff-One)*Half+T_Eff

    if (ST_Eff >= PP_Count) then
      write(SLine,'(A,F5.2,A)') 'Computing 2-electron integrals,',ST_Eff/PP_Eff,'% done so far.'
      call StatusLine(' Seward:',SLine)
      PP_Count = PP_Count+PP_Eff_delta
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    A_int = TMax(iS,jS)*TMax(kS,lS)
    if (Semi_Direct) then

      ! No density screening in semi-direct case!
      ! Cutint: Threshold for Integrals. In semi-direct case, this
      !         must be the final threshold used in the last scf
      !         iteration
      ! Thrint: Threshold for Density screening. This the actual
      !         threshold
      !         for the current iteration

      if (A_Int < CutInt) Skip = .true.
    else

      if (FckNoClmb) then
        Dtst = max(DMax(is,ls)/Four,DMax(is,ks)/Four,DMax(js,ls)/Four,DMax(js,ks)/Four)
      else if (FckNoExch) then
        Dtst = max(DMax(is,js),DMax(ks,ls))
      else
        Dtst = max(DMax(is,ls)/Four,DMax(is,ks)/Four,DMax(js,ls)/Four,DMax(js,ks)/Four,DMax(is,js),DMax(ks,ls))
      end if

      if (A_int*Dtst < ThrInt) Skip = .true.

    end if
    if (Do_DCCD .and. (iSD(10,iS) /= iSD(10,kS))) Skip = .true.
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (.not. Skip) call Eval_IJKL(iS,jS,kS,lS,TInt,nTInt)
    Skip = .false.

    Quad_ijkl = Quad_ijkl+One
    if (Quad_ijkl-TskHi > 1.0e-10_wp) exit
    klS = klS+1
    if (klS > ijS) then
      ijS = ijS+1
      klS = 1
    end if
    iS = ip_ij(1,ijS)
    jS = ip_ij(2,ijS)
    kS = ip_ij(1,klS)
    lS = ip_ij(2,klS)
  end do

  ! Task endpoint

  if (Semi_Direct) then
    if (W2Disc) then
      call Put_QLast()
    else
      call Pos_QLast(Disc)
    end if
  end if

end do
! End of big task loop
call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
if (Semi_Direct) call Close_SemiDSCF()
FstItr = .false.

call mma_deallocate(ip_ij)
call mma_deallocate(DMax)
call mma_deallocate(TMax)

call Term_Ints()

call Free_DeDe(Dens,TwoHam,nDens)
Int_PostProcess => null()

! Broadcast contributions to the Fock matrix

call Sync_TH(TwoHam,nDens)
!                                                                      *
!***********************************************************************
!                                                                      *
!MAW start
!call fmm_call_get_J_matrix(nDens,1,dens,TwoHam)
!MAW end
call Free_iSD()
!call Init_Int_Options()    ?

end subroutine Drv2El_dscf
