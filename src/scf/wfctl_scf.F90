!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it ana/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               1995,1996, Martin Schuetz                              *
!               2003, Valera Veryazov                                  *
!               2016,2017,2022, Roland Lindh                           *
!               2024, Daniel Wessling                                  *
!               2024, Ignacio Fdez. Galvan                             *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine WfCtl_SCF(iTerm,Meth,FstItr,SIntTh)
!***********************************************************************
!                                                                      *
!     purpose: Optimize SCF wavefunction.                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
!     University of Lund, Sweden, 1992                                 *
!     QNR/DIIS, M. Schuetz, 1995                                       *
!     NDDO start orbitals, M. Schuetz, 1996                            *
!     University of Lund, Sweden, 1992,95,96                           *
!     UHF, V.Veryazov, 2003                                            *
!     Cleanup, R. Lindh, 2016                                          *
!                                                                      *
!***********************************************************************

#ifdef _MSYM_
use, intrinsic :: iso_c_binding, only: c_ptr
use InfSCF, only: nBO
#endif
use Interfaces_SCF, only: OptClc_X, TraClc_i
use LnkLst, only: GetVec, LLDelt, LLGrad, LLx, LstPtr, PutVec, SCF_V
use InfSCF, only: AccCon, Aufb, CMO, CMO_Ref, CPUItr, Damping, DIIS, DIISTh, DltNrm, DltnTh, DMOMax, DoCholesky, DSCF, DThr, E1V, &
                  E2V, EDiff, Energy, EneV, EOrb, EThr, Expand, FckAuf, FMOMax, FThr, idKeep, iDMin, Iter, Iter_Ref, Iter_Start, &
                  iterGEK, iterSO, iterSO_Max, jPrint, kOptim, kOptim_Max, kOV, KSDFT, Loosen, MaxFlip, MiniDn, mOV, MSYMON, &
                  MxIter, MxOptm, nAufb, nBas, nBB, nBB, nBT, nD, Neg2_Action, nIter, nIterP, nnB, nnB, nOcc, nOrb, nSym, OccNo, &
                  One_Grid, Ovrlp, QNRTh, RGEK, RSRFO, rTemp, S2Uhf, Teee, TemFac, TimFld, TrDD, TrDh, TrDP, TrM, TStop, &
                  Two_Thresholds, WarnCfg, WarnPocc
use Cholesky, only: ChFracMem
use SCFFiles, only: LuOut
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Ten, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iTerm
character(len=*), intent(in) :: Meth
logical(kind=iwp), intent(inout) :: FstItr
real(kind=wp), intent(inout) :: SIntTh
integer(kind=iwp) :: i, iAufOK, iBas, iCMO, iDummy(7,8), Ind(MxOptm), iNode, iOffOcc, iOpt, iOpt_DIIS, iRC, iSym, iter_, &
                     Iter_DIIS, Iter_DIIS_min, Iter_no_DIIS, iTrM, itSave, jpxn, lth, MinDMx, nBs, nCI, nOr, nTr
real(kind=wp) :: ang, DD, DiisTH_Save, dqdq, dqHdq, Dummy(1), EnVOld, EThr_new, GradNorm, LastStatus, LastStep = 0.1_wp, TCP1, &
                 TCP2, TCPU1, TCPU2, TWall1, TWall2
logical(kind=iwp) :: AllowFlip, AufBau_Done, Converged, Diis_Save, FckAuf_save, FrstDs, Loosen_Active, QNR1st, Reset, Reset_GEK, &
                     Reset_NQ, Reset_Thresh, SORange
character(len=128) :: OrbName
character(len=72) :: Note
character(len=32) :: IterText
character(len=10) :: Meth_
character(len=9) :: StrSave
#ifdef _MSYM_
integer(kind=iwp) :: iD
real(kind=wp) :: Whatever
type(c_ptr) :: msym_ctx
#endif
real(kind=wp), allocatable :: CInter(:,:), D1Sao(:), Disp(:), Grd1(:), Prev(:), Xn(:), Xnp1(:)
real(kind=wp), parameter :: E2VTolerance = -1.0e-8_wp, StepMax = Ten
real(kind=wp), external :: DDot_, Seconds

#include "warnings.h"

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
! Allocate memory for some arrays

nTr = MxIter
call mma_allocate(TrDh,nTR,nTR,nD,Label='TrDh')
call mma_allocate(TrDP,nTR,nTR,nD,Label='TrDP')
call mma_allocate(TrDD,nTR,nTR,nD,Label='TrDD')
nCI = MxOptm+1
call mma_allocate(CInter,nCI,nD,Label='CInter')

!----------------------------------------------------------------------*

call CWTime(TCpu1,TWall1)
! Put empty spin density on the runfile.
call mma_allocate(D1Sao,nBT,Label='D1Sao')
D1Sao(:) = Zero
call Put_dArray('D1sao',D1Sao,nBT)
call mma_deallocate(D1Sao)

if (len_trim(Meth) > len(Meth_)) then
  Meth_ = '[...]'
else
  Meth_ = trim(Meth)
end if
iTerm = 0
iDMin = -1
! Choose between normal and minimized differences
MinDMx = 0
if (MiniDn) MinDMx = max(0,nIter(nIterP)-1)

! Optimization options
!
! iOpt=0: DIIS interpolation on density or/and Fock matrices
! iOpt=1: DIIS extrapolation on gradients w.r.t orbital rot.
! iOpt=2: DIIS extrapolation on the anti-symmetrix X matrix.
! iOpt=3: RS-RFO in the space of the anti-symmetric X matrix.
! iOpt=4: s-GEK/RVO in the space of the anti-symmetric X matrix.

iOpt = 0
QNR1st = .true.
FrstDs = .true.

! START INITIALIZATION

! Initialize Cholesky information if requested

if (DoCholesky) then
  call Cho_X_init(irc,ChFracMem)
  if (irc /= 0) then
    call WarningMessage(2,'WfCtl. Non-zero rc in Cho_X_init.')
    call Abend()
  end if
end if
!                                                                      *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!                                                                      *
IterSO = 0        ! number of second order steps.
IterGEK = 0       ! number of data points in S-GEK
kOptim = 1
Iter_Diis_min = 2
Iter_no_Diis = 2
Converged = .false.
SORange = Damping ! within 2nd-order range for S-GEK
Reset_GEK = .false.
Reset_NQ = .false.
IterText = 'Iterations'
!                                                                      *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!                                                                      *
! END INITIALIZATION
!                                                                      *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!                                                                      *
DiisTh = max(DiisTh,QNRTh)

! If DIIS is turned off make threshold for activation impossible

if (.not. DIIS) then
  DiisTh = Zero
  Iter_no_Diis = 1000000
end if

! If no damping to precede the DIIS make threshold being fulfilled at once.

if (.not. Damping) then
  DiisTh = DiisTh*1.0e99_wp
  Iter_no_Diis = 1
  Iter_Diis_min = 0
end if

! turn temporarily off DIIS & QNR/DIIS, if Aufbau is active...

DiisTh_save = DiisTh
Diis_Save = Diis
if (Aufb) then
  DiisTh = Zero
  Diis = .false.
end if

! Print header to iterations

if ((KSDFT == 'SCF') .or. One_Grid) call PrBeg(Meth_)
AufBau_Done = .false.
!                                                                      *
!======================================================================*
!======================================================================*
!                                                                      *
! If direct SCF modify thresholds temporarily

Reset = .false.
Reset_Thresh = .false.

! pow: temporary disabling of threshold switching

if ((DSCF .or. (KSDFT /= 'SCF')) .and. (nIter(nIterP) > 10)) then

  EThr_New = EThr*Ten**2

  if (DSCF .and. (KSDFT == 'SCF') .and. Two_Thresholds) then
    Reset = .true.
    Reset_Thresh = .true.
    call Reduce_Thresholds(EThr_New,SIntTh)
  end if

  if ((KSDFT /= 'SCF') .and. (.not. One_Grid)) then
    Reset = .true.
    call Modify_NQ_Grid()
    call PrBeg(Meth_)
  end if

end if
!                                                                      *
!======================================================================*
!======================================================================*
!                                                                      *
!                     I  t  e  r  a  t  i  o  n  s                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! nIterP=0: NDDO                                                       *
! nIterP=1: SCF                                                        *
!                                                                      *
!======================================================================*
!                                                                      *
! Set some parameters to starting defaults.

AllowFlip = .true.
Loosen_Active = .false.
iAufOK = 0
Iter_DIIS = 0
EDiff = Zero
DMOMax = 1.0e99_wp
FMOMax = Zero
DltNrm = Zero
LastStatus = -One

if (MSYMON) then
# ifdef _MSYM_
  write(u6,*) 'Symmetrizing orbitals during SCF calculation'
  call fmsym_create_context(msym_ctx)
  call fmsym_set_elements(msym_ctx)
  call fmsym_find_symmetry(msym_ctx)
# else
  write(u6,*) 'No msym support, skipping symmetrization of SCF orbitals...'
# endif
end if
!                                                                      *
!======================================================================*
!                                                                      *
! Start of iteration loop
!                                                                      *
!======================================================================*
!                                                                      *
do iter_=1,nIter(nIterP)
  iter = iter_
  WarnCfg = .false.

  if ((.not. Aufb) .and. (iter > MaxFlip)) AllowFlip = .false.

  TCP1 = seconds()
  if (LastStatus < Zero) LastStatus = TCP1

  iDMin = iDMin+1
  if (iDMin > MinDMx) iDMin = MinDMx

# ifdef _MSYM_
  if (MSYMON) then
    do iD=1,nD
      call fmsym_symmetrize_orbitals(msym_ctx,CMO(1,iD))
      call ChkOrt(iD,Whatever)
      call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
      call ChkOrt(iD,Whatever)
    end do
  end if
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Do Aufbau procedure, if active...

  if (Aufb) call Aufbau(nAufb,OccNo,nnB,iAufOK,nD)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Save the variational energy of the previous step

  EnVold = EneV
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  call SCF_Energy(FstItr,E1V,E2V,EneV)
  Energy(iter) = EneV
  if (iter == 1) then
    EDiff = Zero
  else
    EDiff = EneV-EnVold
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  !      Select on the fly the current optimization method
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Test if this is a DIIS extrapolation iteration, alternatively
  ! the iteration is a DIIS interpolation iteration. The former
  ! is activated if the DMOMax is lower than the threshold
  ! after a specific number of iteration, or if the condition
  ! has already been achieved.
  !
  ! 2017-02-03: add energy criterion to make sure that the DIIS
  !             gets some decent to work with.

  if ((DMOMax < DiisTh) .and. (Iter > Iter_no_Diis) .and. (EDiff < 0.1_wp)) then

    if (iOpt == 0) kOptim = 2
    iOpt = max(1,iOpt)
    Iter_DIIS = Iter_DIIS+1
  end if

  ! Test if the DIIS scheme will be operating in an orbital
  ! rotation mode or linear combination of density matrices. This
  ! option is available only in the extrapolation mode of DIIS.
  !
  ! 2017-02-03: Make sure that the density based DIIS is in
  !             action for at least 2 iterations such that
  !             when the orbital rotation DIIS is turned on
  !             we are firmly in the NR region.

  if (RGEK .and. (.not. (Damping .or. Aufb))) then
    iOpt = 4
    if ((.not. SORange) .and. (Iter > 1)) then
      call mma_allocate(Grd1,mOV,Label='Grd1')
      call GetVec(iter-1,LLGrad,inode,Grd1,mOV)
      GradNorm = sqrt(ddot_(mov,Grd1,1,Grd1,1))
      call mma_deallocate(Grd1)
      if (GradNorm < QNRTh) Reset_GEK = .true.
    end if
  else if ((iOpt >= 2) .or. ((iOpt == 1) .and. (DMOMax < QNRTh) .and. (Iter_DIIS >= Iter_DIIS_min))) then
    if (RSRFO .or. RGEK) then
      if (RSRFO) then
        iOpt = 3
      else
        iOpt = 4
      end if
    else
      iOpt = 2
    end if

  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  select case (iOpt)

    case (0)
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! Interpolation DIIS
      !
      ! Compute traces: TrDh, TrDP and TrDD.
      !
      !FckAuf_Save = FckAuf
      !FckAuf = .false.
      call TraClc_i(iter,nD)

      ! DIIS interpolation optimization: EDIIS, ADIIS, LDIIS

      iOpt_DIIS = 1 ! EDIIS option

      call DIIS_i(CInter,nCI,TrDh,TrDP,TrDD,MxIter,nD,iOpt_DIIS,Ind)

      ! Compute optimal density, dft potentials, and TwoHam

      call OptClc(CInter,nCI,nD,Ind,MxOptm)

      ! Update Fock Matrix from OneHam and extrapolated TwoHam & Vxc

      call mk_FockAO(nIter(nIterP))
      ! Diagonalize Fock matrix and obtain new orbitals

      call NewOrb_SCF(AllowFlip)

      ! Transform density matrix to MO basis

      call MODens()
      !FckAuf = FckAuf_Save
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
    case (1)
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! Extrapolation DIIS
      !
      ! The section of the code operates on the derivatives of
      ! the energy w.r.t the elements of the anti-symmetric X matrix.
      !
      ! The optimal density matrix is found with the DIIS and
      ! canonical CMOs are generated by the diagonalization of the
      ! Fock matrix.

      FckAuf_Save = FckAuf
      FckAuf = .false.
      if (FrstDs) then
        Iter_Ref = 1
        Iter_Start = 1
        CMO_Ref(:,:) = CMO(:,:)
      end if
      call GrdClc(FrstDs)

      call DIIS_x(nD,CInter,nCI,.false.,Ind)

      ! Compute optimal density, dft potentials, and TwoHam

      call OptClc(CInter,nCI,nD,Ind,MxOptm)

      ! Update Fock Matrix from OneHam and extrapolated TwoHam & Vxc

      call mk_FockAO(nIter(nIterP))
      ! Diagonalize Fock matrix and obtain new orbitals

      call NewOrb_SCF(AllowFlip)

      ! Transform density matrix to MO basis

      call MODens()
      FckAuf = FckAuf_Save
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
    case (2,3,4)
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! Extrapolation DIIS & QNR
      !
      ! or
      !
      ! Quasi-2nd order scheme (rational function optimization)
      ! with restricted.step.
      !
      ! In this section we operate directly on the anti-symmetric X
      ! matrix.
      !
      ! Initially the energy is determined through a line seach,
      ! followed by the Fock matrix etc. In this respect, this
      ! approach is doing in this iteration what was already done
      ! by the two other approaches in the end of the previous
      ! iteration.
      !
      ! Note also, that since we work directly with the X matrix we
      ! also rotate the orbitals with this matrix (note the call to
      ! RotMOs right before the call to SCF_Energy in the line
      ! search routine, linser). Thus, the orbitals here are not
      ! canonical and would require at the termination of the
      ! optimization that these are computed.
      !
      ! Initiate if the first QNR step

      if (QNR1st) then

        ! 1st QNR step, reset kOptim to 1

        kOptim = 1
        CInter(1,1) = One
        CInter(1,nD) = One

        Iter_Start = Iter
        QNR1st = .false.
      end if

      ! Set the reference set of parameters and the corresponding
      ! CMOs to be the current iteration.

      if (Iter == Iter_Start) then
        Iter_Ref = Iter
        CMO_Ref(:,:) = CMO(:,:)
        ! init 1st orb rot parameter X1 (set it to zero)
        call mma_allocate(Xn,mOV,Label='Xn')
        Xn(:) = Zero
        ! and store it on appropriate LList
        call PutVec(Xn,mOV,iter,'OVWR',LLx)
        call mma_deallocate(Xn)
      end if

      ! Compute the current gradient
      call SCF_Gradient()

      ! Set the reference set of parameters and the corresponding
      ! CMOs to be the current iteration.
      call Move_Ref(Iter)

      ! Update the Fock Matrix from actual OneHam, Vxc & TwoHam AO basis

      call mk_FockAO(nIter(nIterP))

      ! and transform the Fock matrix into the new MO space,

      call TraFck(.false.,FMOMax)

      ! compute initial diagonal Hessian, Hdiag

      call SOIniH()
      AccCon(8:8) = 'H'

      ! update the QNR iteration counter
      IterSO = min(IterSO+1,IterSO_Max)

      ! Allocate memory for the current gradient and displacement vector.

      call mma_allocate(Grd1,mOV,Label='Grd1')
      call mma_allocate(Disp,mOV,Label='Disp')
      call mma_allocate(Xnp1,mOV,Label='Xnp1')

      select case (iOpt)

        case (2) ! qNRC2DIIS
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
          ! Compute extrapolated g_x(n) and X_x(n)

          do
            call DIIS_x(nD,CInter,nCI,.true.,Ind)

            call OptClc_X(CInter,nCI,nD,Grd1,mOV,Ind,MxOptm,kOptim,kOV,LLGrad)
            call OptClc_X(CInter,nCI,nD,Xnp1,mOV,Ind,MxOptm,kOptim,kOV,LLx)

            ! compute new displacement vector delta
            ! dX_x(n) = -H(-1)*g_x(n) ! Temporary storage in Disp

            call SOrUpV(Grd1,mOV,Disp,'DISP','BFGS')

            DD = sqrt(DDot_(mOV,Disp,1,Disp,1))

            if (DD > Pi) then
#             ifdef _DEBUGPRINT_
              write(u6,*) 'WfCtl_SCF: Additional displacement is too large.'
              write(u6,*) 'DD=',DD
#             endif
              if (kOptim /= 1) then
#               ifdef _DEBUGPRINT_
                write(u6,*) 'Reset update depth in BFGS, redo the DIIS'
#               endif
                kOptim = 1
                Iter_Start = Iter
                IterSO = 1
                cycle
              else
#               ifdef _DEBUGPRINT_
                write(u6,*) 'Scale the step to be within the threshold.'
                write(u6,*) 'LastStep=',LastStep
#               endif
                Disp(:) = Disp(:)*(LastStep/DD)
              end if
            end if

            ! from this, compute new orb rot parameter X(n+1)
            !
            ! X(n+1) = X_x(n) - H(-1)g_x(n)
            ! X(n+1) = X_x(n) + dX_x(n)

            Xnp1(:) = Xnp1(:)-Disp(:)

            ! get address of actual X(n) in corresponding LList

            jpXn = LstPtr(iter,LLx)

            ! and compute actual displacement dX(n)=X(n+1)-X(n)

            Disp(:) = Xnp1(:)-SCF_V(jpXn)%A(:)

            DD = sqrt(DDot_(mOV,Disp,1,Disp,1))

            if (DD <= Pi) exit
#           ifdef _DEBUGPRINT_
            write(u6,*) 'WfCtl_SCF: Total displacement is too large.'
            write(u6,*) 'DD=',DD
#           endif
            if (kOptim /= 1) then
#             ifdef _DEBUGPRINT_
              write(u6,*) 'Reset update depth in BFGS, redo the DIIS'
#             endif
              kOptim = 1
              Iter_Start = Iter
              IterSO = 1
            else
#             ifdef _DEBUGPRINT_
              write(u6,*) 'Scale the step to be within the threshold.'
#             endif
              Disp(:) = Disp(:)*(LastStep/DD)
              DD = sqrt(DDot_(mOV,Disp,1,Disp,1))
              exit
            end if
          end do
          LastStep = min(DD,1.0e-2_wp)
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
        case (3) ! RS-RFO
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
          ! Get g(n)

          call GetVec(iter,LLGrad,inode,Grd1,mOV)
          !                                                            *
          !*************************************************************
          !                                                            *
          dqHdq = Zero
          do
            call rs_rfo_scf(Grd1,mOV,Disp,AccCon(1:6),dqdq,dqHdq,StepMax,AccCon(9:9),3)
            DD = sqrt(DDot_(mOV,Disp,1,Disp,1))
            if (DD <= Pi) exit
#           ifdef _DEBUGPRINT_
            write(u6,*) 'WfCtl_SCF: Total displacement is too large.'
            write(u6,*) 'DD=',DD
#           endif
            if (kOptim /= 1) then
#             ifdef _DEBUGPRINT_
              write(u6,*) 'Reset update depth in BFGS, redo the RS-RFO.'
#             endif
              kOptim = 1
              Iter_Start = Iter
              IterSO = 1
            else
              write(u6,*) 'Probably a bug.'
              call Abend()
            end if
          end do
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Pick up X(n) and compute X(n+1)=X(n)+dX(n)

          call GetVec(iter,LLx,inode,Xnp1,mOV)

          Xnp1(:) = Xnp1(:)+Disp(:)
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
        case (4) ! S-GEK
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
          ! update the S-GEK iteration counter
          if (Reset_GEK) then
            ! Reset the GEK surrogate model if the NQ grid has been reset
            if (Reset_NQ) IterGEK = 1
            IterSO = 1
            kOptim = 1
            SORange = .true.
            Reset_GEK = .false.
            Reset_NQ = .false.
#           ifdef _DEBUGPRINT_
            write(u6,*) 'Reset GEK'
#           endif
          else
            IterGEK = min(IterGEK+1,IterSO_Max)
          end if
          Iter_Start = max(Iter_Start,Iter-IterGEK+1)
          IterText = 'Macro Iterations'

          ! Get g(n)

          call GetVec(iter,LLGrad,inode,Grd1,mOV)
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Expand the subspace for GEK
          StrSave = AccCon
          itSave = Iter_Start
          select case (Expand)
            case (1) ! Use DIIS
              ! Compute extrapolated g_x(n) and X_x(n)

              call DIIS_x(nD,CInter,nCI,.true.,Ind)

              call OptClc_X(CInter,nCI,nD,Grd1,mOV,Ind,MxOptm,kOptim,kOV,LLGrad)
              call OptClc_X(CInter,nCI,nD,Xnp1,mOV,Ind,MxOptm,kOptim,kOV,LLx)

              ! compute new displacement vector delta
              ! dX_x(n) = -H(-1)*g_x(n) ! Temporary storage in Disp

              ! Reduce the BFGS update depth until the step is reasonable
              do i=1,IterSO
                call SOrUpV(Grd1,mOV,Disp,'DISP','BFGS')
                DD = sqrt(DDot_(mOV,Disp,1,Disp,1))
                if (DD <= Pi) exit
#               ifdef _DEBUGPRINT_
                if (i == 1) then
                  write(u6,*) 'WfCtl_SCF: Total displacement is large.'
                  write(u6,*) 'DD=',DD
                end if
#               endif
                if (IterSO > 1) then
                  if (i == 1) write(u6,*) 'Reset update depth in BFGS'
                  IterSO = IterSO-1
                  Iter_Start = Iter-IterSO+1
                  kOptim = min(kOptim,IterSO)
                end if
              end do
#             ifdef _DEBUGPRINT_
              if (i > 1) write(u6,*) 'IterSO=',IterSO
#             endif

              ! from this, compute new orb rot parameter X(n+1)
              !
              ! X(n+1) = X_x(n) - H(-1)g_x(n)
              ! X(n+1) = X_x(n) + dX_x(n)

              Xnp1(:) = Xnp1(:)-Disp(:)

              ! get address of actual X(n) in corresponding LList

              jpXn = LstPtr(iter,LLx)

              ! and compute actual displacement dX(n)=X(n+1)-X(n)

              Disp(:) = Xnp1(:)-SCF_V(jpXn)%A(:)
            case (2) ! Use BFGS
              ! Reduce the BFGS update depth until the step is reasonable
              do i=1,IterSO
                call SOrUpV(Grd1,mOV,Disp,'DISP','BFGS')
                DD = sqrt(DDot_(mOV,Disp,1,Disp,1))
                if (DD <= Pi) exit
#               ifdef _DEBUGPRINT_
                if (i == 1) then
                  write(u6,*) 'WfCtl_SCF: Total displacement is large.'
                  write(u6,*) 'DD=',DD
                end if
#               endif
                if (IterSO > 1) then
#                 ifdef _DEBUGPRINT_
                  if (i == 1) write(u6,*) 'Reset update depth in BFGS'
#                 endif
                  IterSO = IterSO-1
                  Iter_Start = Iter-IterSO+1
                  kOptim = min(kOptim,IterSO)
                end if
              end do
#             ifdef _DEBUGPRINT_
              if (i > 1) write(u6,*) 'IterSO=',IterSO
#             endif
            case (3) ! Use RS-RFO
              dqHdq = Zero
              call rs_rfo_scf(Grd1,mOV,Disp,AccCon(1:6),dqdq,dqHdq,StepMax,AccCon(9:9),1)
          end select
          AccCon = StrSave
          Iter_Start = itSave

          Loosen_Active = (Loosen%Factor > One)
          call S_GEK_Optimizer(Disp,mOV,dqdq,AccCon(1:6),AccCon(9:9),SORange)
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Pick up X(n) and compute X(n+1)=X(n)+dX(n)

          call GetVec(iter,LLx,inode,Xnp1,mOV)

          !                                                            *
          !*************************************************************
          !                                                            *
          ! Undershoot avoidance,
          ! when consecutive steps have a large overlap

          if ((Loosen%Step > One) .and. (iterGEK > 1)) then
            call mma_allocate(Prev,mOV,Label='Prev')
            call GetVec(iter-1,LLDelt,inode,Prev,mOV)
            dqdq = DDot_(mOV,Disp,1,Disp,1)*DDot_(mOV,Prev,1,Prev,1)
            ang = DDot_(mOV,Prev,1,Disp,1)/sqrt(dqdq)
            if ((ang < Loosen%Thrs2) .or. (AccCon(9:9) /= ' ')) then
              Loosen%Factor = One
            else if (ang > Loosen%Thrs) then
              Loosen%Factor = Loosen%Factor*Loosen%Step
            end if
            call mma_deallocate(Prev)
          end if

          Xnp1(:) = Xnp1(:)+Disp(:)
          !                                                            *
          !*************************************************************
          !*************************************************************
          !                                                            *
      end select
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! Store X(n+1) and dX(n)

      dqdq = sqrt(DDot_(mOV,Disp,1,Disp,1))
      call PutVec(Xnp1,mOV,iter+1,'OVWR',LLx)
      call PutVec(Disp,mOV,iter,'OVWR',LLDelt)

      ! compute Norm of dX(n)

      DltNrm = real(nD,kind=wp)*dqdq

      ! Generate the CMOs, rotate MOs accordingly to new point

      call RotMOs(Disp,mOV)

      ! and release memory...
      call mma_deallocate(Xnp1)
      call mma_deallocate(Disp)
      call mma_deallocate(Grd1)
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
    case default
      write(u6,*) 'WfCtl_SCF: Illegal option'
      call Abend()
  end select
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Update DIIS interpolation and BFGS depth kOptim

  if (idKeep == 0) then
    kOptim = 1
  else if (idKeep == 1) then
    kOptim = 2
  else
    kOptim = min(kOptim_Max,kOptim+1)
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! End optimization section
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Save the new orbitals in case the SCF program aborts

  iTrM = 1
  iCMO = 1
  do iSym=1,nSym
    nBs = nBas(iSym)
    nOr = nOrb(iSym)
    lth = nBs*nOr
    TrM(iTrM:iTrM+lth-1,1:nD) = CMO(iCMO:iCMO+lth-1,1:nD)
    TrM(iTrm+nBs*nOr:iTrm+nBs*nBs-1,1:nD) = Zero
    iTrM = iTrM+nBs*nBs
    iCMO = iCMO+nBs*nOr
  end do

  if (nD == 1) then
    OrbName = 'SCFORB'
    Note = '*  intermediate SCF orbitals'

    call WrVec_(OrbName,LuOut,'CO',nD-1,nSym,nBas,nBas,TrM(:,1),Dummy,OccNo(:,1),Dummy,Dummy,Dummy,iDummy,Note,2)
    call Put_darray('SCF orbitals',TrM(:,1),nBB)
    call Put_darray('OrbE',Eorb(:,1),nnB)
    if (.not. Aufb) call Put_iarray('SCF nOcc',nOcc(:,1),nSym)
  else
    OrbName = 'UHFORB'
    Note = '*  intermediate UHF orbitals'
    call WrVec_(OrbName,LuOut,'CO',nD-1,nSym,nBas,nBas,TrM(:,1),TrM(:,2),OccNo(:,1),OccNo(:,2),Dummy,Dummy,iDummy,Note,3)
    call Put_darray('SCF orbitals',TrM(:,1),nBB)
    call Put_darray('SCF orbitals_ab',TrM(:,2),nBB)
    call Put_darray('OrbE',Eorb(:,1),nnB)
    call Put_darray('OrbE_ab',Eorb(:,2),nnB)
    if (.not. Aufb) then
      call Put_iarray('SCF nOcc',nOcc(:,1),nSym)
      call Put_iarray('SCF nOcc_ab',nOcc(:,2),nSym)
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  !====================================================================*
  !                                                                    *
  !      A U F B A U  section                                          *
  !                                                                    *
  !====================================================================*
  !                                                                    *
  if (Aufb .and. (.not. Teee)) then

    if ((iter /= 1) .and. (abs(EDiff) < 0.002_wp) .and. (.not. AufBau_Done)) then
      Aufbau_Done = .true.
      Diis = Diis_Save
      DiisTh = DiisTh_save
    end if

  else if (Aufb .and. Teee) then

    ! If one iteration has been performed with the end temperature and we
    ! have a stable aufbau Slater determinant, terminate.
    ! The temperature is scaled down each iteration until it reaches the end
    ! temperature

    ! !!! continue here
    if ((abs(RTemp-TStop) <= 1.0e-3_wp) .and. (iAufOK == 1) .and. (abs(Ediff) < 1.0e-2_wp)) then
      Aufbau_Done = .true.
      Diis = Diis_Save
      DiisTh = DiisTh_save
    end if
    RTemp = TemFac*RTemp
    !RTemp = max(TStop,RTemp)
    TSTop = min(TStop,RTemp)
  end if
  !                                                                    *
  !====================================================================*
  !====================================================================*
  !                                                                    *
  ! Print out necessary things

  TCP2 = seconds()
  CpuItr = TCP2-TCP1
  ! Update status every 10 seconds at most
  if (TCP2-LastStatus > Ten) then
    write(Note,'(A,I0)') 'Iteration ',Iter
    call StatusLine('SCF: ',Note)
    LastStatus = TCP2
  end if

  call PrIte(iOpt >= 2,CMO,nBB,nD,Ovrlp,nBT,OccNo,nnB)

  !--------------------------------------------------------------------*
  call Scf_Mcontrol(iter)
  !--------------------------------------------------------------------*

  ! Check that two-electron energy is positive.
  ! The integrals may not be positive definite, leading to
  ! negative eigenvalues of the ERI matrix - i.e. attractive forces
  ! between the electrons! When that occurs the SCF iterations
  ! still might find a valid SCF solution, but it will obviously be
  ! unphysical and we might as well just stop here.

  if (Neg2_Action /= 'CONT') then
    if (E2V < E2VTolerance) then
      call WarningMessage(2,'WfCtl_SCF: negative two-electron energy')
      write(u6,'(A,ES25.10)') 'Two-electron energy E2V=',E2V
      call xFlush(u6)
      if (Neg2_Action == 'STOP') call Abend()
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Perform another iteration if necessary
  !
  ! Convergence criteria are
  !
  ! Either,
  !
  ! 1) it is not the first iteration
  ! 2) the absolute energy change is smaller than EThr
  ! 3) the absolute Fock matrix change is smaller than FThr
  ! 4) the absolute density matrix change is smaller than DThr,
  !    and
  ! 5) step is NOT a Quasi NR step
  !
  ! or
  !
  ! 1) step is a Quasi NR step, and
  ! 2) DltNrm <= DltNth
  !
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if ((EDiff > Zero) .and. (.not. Reset)) EDiff = Ten*EThr
  if ((iter /= 1) .and. (abs(EDiff) <= EThr) .and. (abs(FMOMax) <= FThr) .and. &
      (((abs(DMOMax) <= DThr) .and. (iOpt < 2)) .or. ((DltNrm <= DltNTh) .and. iOpt >= 2))) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (Aufb) WarnPocc = .true.

    ! If calculation with two different sets of parameters reset
    ! the convergence parameters to default values. This is
    ! possibly used for direct SCF and DFT.

    if (Reset) then
      ! Reset thresholds for direct SCF procedure

      Reset = .false.

      if (iOpt == 2) then
        iOpt = 1        ! True if step is QNR
        QNR1st = .true.
      end if
      IterSO = 0
      if (Reset_Thresh) call Reset_Thresholds()
      if (KSDFT /= 'SCF') then
        if (.not. One_Grid) then
          call Reset_NQ_grid()
          Reset_GEK = .true.
          Reset_NQ = .true.
        end if
        if (iOpt == 0) kOptim = 1
      end if
      cycle
    end if

    ! Here if we converged!

    if (Loosen_Active) then
      ! If in loosen mode, disable it
      Loosen%Factor = One
      Loosen%Step = One
    else
      ! Branch out of the iterative loop! Done!!!
      Converged = .true.
      exit
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else if (Aufb .and. AufBau_Done) then
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Here if Aufbau stage is finished and the calculation will
    ! continue normally.

    Aufb = .false.
    if (jPrint >= 2) then
      write(u6,*)
      write(u6,'(6X,A)') ' Fermi aufbau procedure completed!'
    end if
    if (nD == 1) then
      iOffOcc = 0
      do iSym=1,nSym
        do iBas=1,nOrb(iSym)
          iOffOcc = iOffOcc+1
          if (iBas <= nOcc(iSym,1)) then
            OccNo(iOffOcc,1) = Two
          else
            OccNo(iOffOcc,1) = Zero
          end if
        end do
      end do
      if (jPrint >= 2) then
        write(u6,'(6X,A,8I5)') 'nOcc=',(nOcc(iSym,1),iSym=1,nSym)
        write(u6,*)
      end if
    else
      iOffOcc = 0
      do iSym=1,nSym
        do iBas=1,nOrb(iSym)
          iOffOcc = iOffOcc+1
          if (iBas <= nOcc(iSym,1)) then
            OccNo(iOffOcc,1) = One
          else
            OccNo(iOffOcc,1) = Zero
          end if
          if (iBas <= nOcc(iSym,2)) then
            OccNo(iOffOcc,2) = One
          else
            OccNo(iOffOcc,2) = Zero
          end if
        end do
      end do
      if (jPrint >= 2) then
        write(u6,'(6X,A,8I5)') 'nOcc(alpha)=',(nOcc(iSym,1),iSym=1,nSym)
        write(u6,'(6X,A,8I5)') 'nOcc(beta) =',(nOcc(iSym,2),iSym=1,nSym)
        write(u6,*)
      end if

    end if

    ! Now when nOcc is known compute standard sizes of arrays.

    call Setup_SCF()

  end if
  !                                                                    *
  !====================================================================*
  !                                                                    *
  ! End of iteration loop

end do ! iter_
!                                                                      *
!======================================================================*
!                                                                      *

if (Converged) then

  if (jPrint >= 2) then
    write(u6,*)
    write(u6,'(6X,A,I3,1X,A)') ' Convergence after ',iter,trim(IterText)
  end if

else

  ! Here if we didn't converge or if this was a forced one
  ! iteration calculation.

  iter = iter-1
  if (nIter(nIterP) > 1) then
    if (jPrint >= 1) then
      write(u6,*)
      write(u6,'(6X,A,I3,1X,A)') ' No convergence after ',iter,trim(IterText)
    end if
    iTerm = _RC_NOT_CONVERGED_
  else
    if (jPrint >= 2) then
      write(u6,*)
      write(u6,'(6X,A)') ' Single iteration finished!'
    end if
  end if

  if (jPrint >= 2) write(u6,*)

  if (Reset) then
    if (DSCF .and. (KSDFT == 'SCF')) call Reset_Thresholds()
    if ((KSDFT /= 'SCF') .and. (.not. One_Grid)) then
      call Reset_NQ_grid()
      Reset_GEK = .true.
      Reset_NQ = .true.
    end if
  end if

end if

!**********************************************************
!                      S   T   O   P                      *
!**********************************************************

if (jPrint >= 2) then
  call CollapseOutput(0,'Convergence information')
  write(u6,*)
end if
! Compute total spin (if UHF)
if (nD == 1) then
  s2uhf = Zero
else
  call s2calc(CMO(:,1),CMO(:,2),Ovrlp,nOcc(:,1),nOcc(:,2),nBas,nOrb,nSym,s2uhf)
end if

call Add_Info('SCF_ITER',[real(Iter,kind=wp)],1,8)
!call Scf_XML(0)

call KiLLs()

! If the orbitals are generated with orbital rotations in
! RotMOs we need to generate the canonical orbitals.

if (iOpt >= 2) then

  ! Generate canonical orbitals

  call TraFck(.true.,FMOMax)

  ! Transform density matrix to MO basis

  call MODens()

end if

! Compute correct orbital energies and put on the runfile.

call Mk_EOrb()

! Put orbital coefficients on the runfile.

call Put_darray('SCF orbitals',CMO(1,1),nBB)
if (nD == 2) call Put_darray('SCF orbitals_ab',CMO(1,2),nBB)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Finalize Cholesky information if initialized
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
if (DoCholesky) then
  call Cho_X_Final(irc)
  if (irc /= 0) then
    call WarningMessage(2,'WfCtl. Non-zero rc in Cho_X_Final.')
    call Abend()
  end if
end if
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
#ifdef _MSYM_
if (MSYMON) call fmsym_release_context(msym_ctx)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Exit                                                             *
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
TimFld(2) = TimFld(2)+(TCpu2-TCpu1)

call mma_deallocate(CInter)
call mma_deallocate(TrDD)
call mma_deallocate(TrDP)
call mma_deallocate(TrDh)

!                                                                      *
!***********************************************************************
!                                                                      *
end subroutine WfCtl_SCF
