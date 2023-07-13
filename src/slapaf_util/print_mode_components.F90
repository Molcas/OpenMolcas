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
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
! Print_Mode_Components
!
!> @brief
!> Print the contributions from the primitive internal coordinates to the
!> vibrational modes.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Compute and print the contributions from the primitive internal coordinates
!> (stretches, angles, dihedrals, out-of-planes) to the vibrational modes
!>
!> @param[in] Modes Vibrational modes, as computed by e.g. ::freqanal
!> @param[in] Freq Vibrational frequencies
!> @param[in] nModes Number of modes
!> @param[in] lModes Size of \p Modes
!> @param[in] lDisp Number of displacements per irrep
!***********************************************************************

subroutine Print_Mode_Components(Modes,Freq,nModes,lModes,lDisp)

use Symmetry_Info, only: nIrrep, VarR, VarT
use Slapaf_Info, only: Analytic_Hessian, ANr, ApproxNADC, AtomLbl, Baker, Beta, Beta_Disp, BM, BMx, BSet, CallLast, CnstWght, &
                       Coor, Cubic, Curvilinear, Cx, dBM, ddV_Schlegel, Degen, Delta, DipM, dMass, dMEPStep, dqInt, EDiffZero, &
                       eMEPTest, Energy, Energy0, FindTS, GNrm, GNrm_Threshold, Grd, Gx, Gx0, Header, HrmFrq_Show, HSet, HUpMet, &
                       HWRS, iBM, iCoSet, idBM, iInt, iOptC, iOptH, IRC, iRef, iRow, iRow_c, isFalcon, iState, iter, jStab, &
                       Lambda, Lbl, lCtoF, lDoubleIso, Line_Search, lNmHss, lOld, lOld_Implicit, lSoft, lTherm, Max_Center, &
                       MaxItr, mB_Tot, mdB_Tot, MEP, MEP_Algo, MEP_Type, MEPCons, MEPNum, MF, Mode, mq, mRowH, mTROld, mTtAtm, &
                       MxItr, NAC, NADC, nBVec, nDimBC, nFix, nLambda, nMEP, NmIter, nqBM, nsRot, nStab, Numerical, nUserPT, &
                       nWndw, PrQ, Q_nuclear, qInt, Redundant, RefGeo, Request_Alaska, Request_RASSI, rFuzz, rHidden, rMEP, &
                       RootMap, RtRnc, Shift, SlStop, Smmtrc, ThrCons, ThrEne, ThrGrd, ThrMEP, Track, TSConstraints, TwoRunFiles, &
                       UpMeth, User_Def, UserP, UserT, WeightedConstraints, Weights
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Modes(*), Freq(*)
integer(kind=iwp), intent(in) :: nModes, lModes, lDisp(nIrrep)
#include "Molcas.fh"
integer(kind=iwp) :: Bk_iInt, Bk_iOptC, Bk_iOptH, Bk_IRC, Bk_iRef, Bk_iRow, Bk_iRow_c, Bk_iState(2), Bk_iter, Bk_Max_Center, &
                     Bk_mB_Tot, Bk_mdB_Tot, Bk_MEPnum, Bk_mode, Bk_mq, Bk_mTROld, Bk_mTtAtm, Bk_MxItr, Bk_nBVec, Bk_nDimBC, &
                     Bk_nFix, Bk_nLambda, Bk_nMEP, Bk_NmIter, Bk_nsRot, Bk_nUserPT, Bk_nWndw, i, iB, iErr, ii, im, iq, j, LuIC, &
                     LuInput, nAll_Atoms, nB, nK, nQQ, nsAtom, nUnique_Atoms, nX, PLback
real(kind=wp) :: Bk_Beta, Bk_Beta_Disp, Bk_CnstWght, Bk_Delta, Bk_dMEPStep, Bk_GNrm_Threshold, Bk_rFuzz, Bk_rHidden, Bk_Rtrnc, &
                 Bk_ThrCons, Bk_ThrEne, Bk_ThrGrd, Bk_ThrMEP, Bk_UserP, Bk_UserT(64), Mx, MinComp, rDum(1)
logical(kind=iwp) :: Bk_Analytic_Hessian, Bk_ApproxNADC, Bk_Baker, Bk_BSet, Bk_CallLast, Bk_Cubic, Bk_CurviLinear, &
                     Bk_DDV_Schlegel, Bk_EDiffZero, Bk_eMEPtest, Bk_FindTS, Bk_HrmFrq_Show, Bk_HSet, Bk_HWRS, Bk_isFalcon, &
                     Bk_lCtoF, Bk_lDoubleIso, Bk_Line_Search, Bk_lNmHss, Bk_lOld, Bk_lOld_Implicit, Bk_lSoft, Bk_lTherm, Bk_MEP, &
                     Bk_MEPCons, Bk_NADC, Bk_Redundant, Bk_Request_Alaska, Bk_Request_RASSI, Bk_rMEP, Bk_SlStop, Bk_Track, &
                     Bk_TSConstraints, Bk_TwoRunFiles, Bk_User_Def, Bk_VarR, Bk_VarT, Bk_WeightedConstraints, Found
character(len=180) :: Line
character(len=16) :: StdIn
character(len=10) :: Bk_MEP_Type
character(len=8) :: Filename
character(len=6) :: Bk_HUpMet, Bk_UpMeth
character(len=2) :: Bk_MEP_Algo
character :: Bk_Header(144)
integer(kind=iwp), allocatable :: Bk_ANr(:), Bk_iBM(:), Bk_iCoSet(:,:), Bk_idBM(:), Bk_jStab(:,:), Bk_mRowH(:), Bk_nqBM(:), &
                                  Bk_nStab(:), Bk_RootMap(:), Sort(:)
real(kind=wp), allocatable :: Bk_BM(:), Bk_BMx(:,:), Bk_Coor(:,:), Bk_Cx(:,:,:), Bk_dBM(:), Bk_Degen(:,:), Bk_DipM(:,:), &
                              Bk_dMass(:), Bk_dqInt(:,:), Bk_Energy(:), Bk_Energy0(:), Bk_GNrm(:), Bk_Grd(:,:), Bk_Gx(:,:,:), &
                              Bk_Gx0(:,:,:), Bk_Lambda(:,:), Bk_MF(:,:), Bk_NAC(:,:,:), Bk_Q_nuclear(:), Bk_qInt(:,:), &
                              Bk_RefGeo(:,:), Bk_Shift(:,:), Bk_Weights(:), IntMod(:,:), KKtB(:,:), KMtrx(:,:), KTrsp(:,:), &
                              NMod(:,:)
logical(kind=iwp), allocatable :: Bk_Smmtrc(:,:)
character(len=LenIn), allocatable :: Bk_AtomLbl(:)
character(len=24), allocatable :: Label(:)
character(len=8), allocatable :: Bk_Lbl(:)
integer(kind=iwp), parameter :: nLbl = 10*MxAtom
integer(kind=iwp), external :: AixRm, iPrintLevel, IsFreeUnit
real(kind=wp), external :: DDot_
character(len=180), external :: Get_Ln_EOF

! Ugly hack: backup all "global" slapaf variables in case this is
!            called from inside slapaf

! Note, this routine might be called from outside the Slapaf
! environment. In which case there is no backup to be made.

if (allocated(Cx)) then
  nsAtom = size(Cx,2)
  call mma_allocate(Bk_Cx,3,nsAtom,MaxItr+1,Label='Bk_Cx')
  Bk_Cx(:,:,:) = Cx(:,:,:)
  call mma_deallocate(Cx)
  if (allocated(Gx)) then
    call mma_allocate(Bk_Gx,3,nsAtom,MaxItr+1,Label='Bk_Gx')
    Bk_Gx(:,:,:) = Gx(:,:,:)
    call mma_deallocate(Gx)
  end if
  if (allocated(Gx0)) then
    call mma_allocate(Bk_Gx0,3,nsAtom,MaxItr+1,Label='Bk_Gx0')
    Bk_Gx0(:,:,:) = Gx0(:,:,:)
    call mma_deallocate(Gx0)
  end if
  if (allocated(NAC)) then
    call mma_allocate(Bk_NAC,3,nsAtom,MaxItr+1,Label='Bk_NAC')
    Bk_NAC(:,:,:) = NAC(:,:,:)
    call mma_deallocate(NAC)
  end if
  if (allocated(Q_nuclear)) then
    call mma_allocate(Bk_Q_nuclear,nsAtom,Label='Bk_Q_nuclear')
    Bk_Q_nuclear(:) = Q_nuclear(:)
    call mma_deallocate(Q_nuclear)
  end if
  if (allocated(dMass)) then
    call mma_allocate(Bk_dMass,nsAtom,Label='Bk_dMass')
    Bk_dMass(:) = dMass(:)
    call mma_deallocate(dMass)
  end if
  if (allocated(Coor)) then
    call mma_allocate(Bk_Coor,3,nsAtom,Label='Bk_Coor')
    Bk_Coor(:,:) = Coor(:,:)
    call mma_deallocate(Coor)
  end if
  if (allocated(Grd)) then
    call mma_allocate(Bk_Grd,3,nsAtom,Label='Bk_Grd')
    Bk_Grd(:,:) = Grd(:,:)
    call mma_deallocate(Grd)
  end if
  if (allocated(ANr)) then
    call mma_allocate(Bk_ANr,nsAtom,Label='Bk_ANr')
    Bk_ANr(:) = ANr(:)
    call mma_deallocate(ANr)
  end if
  if (allocated(Weights)) then
    call mma_allocate(Bk_Weights,size(Weights),Label='Bk_Weights')
    Bk_Weights(:) = Weights(:)
    call mma_deallocate(Weights)
  end if
  if (allocated(Shift)) then
    call mma_allocate(Bk_Shift,size(Shift,1),MaxItr,Label='Bk_Shift')
    Bk_Shift(:,:) = Shift(:,:)
    call mma_deallocate(Shift)
  end if
  if (allocated(GNrm)) then
    call mma_allocate(Bk_GNrm,MaxItr+1,Label='Bk_GNrm')
    Bk_GNrm(:) = GNrm(:)
    call mma_deallocate(GNrm)
  end if
  if (allocated(Lambda)) then
    call mma_allocate(Bk_Lambda,size(Lambda,1),MaxItr+1,Label='Bk_Lambda')
    Bk_Lambda(:,:) = Lambda(:,:)
    call mma_deallocate(Lambda)
  end if
  if (allocated(Energy)) then
    call mma_allocate(Bk_Energy,MaxItr+1,Label='Bk_Energy')
    Bk_Energy(:) = Energy(:)
    call mma_deallocate(Energy)
  end if
  if (allocated(Energy0)) then
    call mma_allocate(Bk_Energy0,MaxItr+1,Label='Bk_Energy0')
    Bk_Energy0(:) = Energy0(:)
    call mma_deallocate(Energy0)
  end if
  if (allocated(MF)) then
    call mma_allocate(Bk_MF,3,nsAtom,Label='Bk_MF')
    Bk_MF(:,:) = MF(:,:)
    call mma_deallocate(MF)
  end if
  if (allocated(DipM)) then
    call mma_allocate(Bk_DipM,3,MaxItr+1,Label='Bk_DipM')
    Bk_DipM(:,:) = DipM(:,:)
    call mma_deallocate(DipM)
  end if
  if (allocated(qInt)) then
    call mma_allocate(Bk_qInt,size(qInt,1),MaxItr,Label='Bk_qInt')
    Bk_qInt(:,:) = qInt(:,:)
    call mma_deallocate(qInt)
  end if
  if (allocated(dqInt)) then
    call mma_allocate(Bk_dqInt,size(dqInt,1),MaxItr,Label='Bk_dqInt')
    Bk_dqInt(:,:) = dqInt(:,:)
    call mma_deallocate(dqInt)
  end if
  if (allocated(RefGeo)) then
    call mma_allocate(Bk_RefGeo,3,nsAtom,Label='Bk_RefGeo')
    Bk_RefGeo(:,:) = RefGeo(:,:)
    call mma_deallocate(RefGeo)
  end if
  if (allocated(BM)) then
    call mma_allocate(Bk_BM,size(BM),Label='Bk_BM')
    Bk_BM(:) = BM(:)
    call mma_deallocate(BM)
  end if
  if (allocated(dBM)) then
    call mma_allocate(Bk_dBM,size(dBM),Label='Bk_dBM')
    Bk_dBM(:) = dBM(:)
    call mma_deallocate(dBM)
  end if
  if (allocated(iBM)) then
    call mma_allocate(Bk_iBM,size(iBM),Label='Bk_iBM')
    Bk_iBM(:) = iBM(:)
    call mma_deallocate(iBM)
  end if
  if (allocated(idBM)) then
    call mma_allocate(Bk_idBM,size(idBM),Label='Bk_idBM')
    Bk_idBM(:) = idBM(:)
    call mma_deallocate(idBM)
  end if
  if (allocated(nqBM)) then
    call mma_allocate(Bk_nqBM,size(nqBM),Label='Bk_nqBM')
    Bk_nqBM(:) = nqBM(:)
    call mma_deallocate(nqBM)
  end if
  if (allocated(BMx)) then
    call mma_allocate(Bk_BMx,size(BMx,1),size(BMx,2),Label='Bk_BMx')
    Bk_BMx(:,:) = BMx(:,:)
    call mma_deallocate(BMx)
  end if
  if (allocated(Degen)) then
    call mma_allocate(Bk_Degen,size(Degen,1),size(Degen,2),Label='Bk_Degen')
    Bk_Degen(:,:) = Degen(:,:)
    call mma_deallocate(Degen)
  end if
  if (allocated(jStab)) then
    call mma_allocate(Bk_jStab,[0,7],[1,size(jStab,2)],Label='Bk_jStab')
    Bk_jStab(:,:) = jStab(:,:)
    call mma_deallocate(jStab)
  end if
  if (allocated(iCoSet)) then
    call mma_allocate(Bk_iCoSet,[0,7],[1,size(iCoSet,2)],Label='Bk_iCoSet')
    Bk_iCoSet(:,:) = iCoSet(:,:)
    call mma_deallocate(iCoSet)
  end if
  if (allocated(nStab)) then
    call mma_allocate(Bk_nStab,size(nStab),Label='Bk_nStab')
    Bk_nStab(:) = nStab(:)
    call mma_deallocate(nStab)
  end if
  if (allocated(AtomLbl)) then
    call mma_allocate(Bk_AtomLbl,size(AtomLbl),Label='Bk_AtomLbl')
    Bk_AtomLbl(:) = AtomLbl(:)
    call mma_deallocate(AtomLbl)
  end if
  if (allocated(Smmtrc)) then
    call mma_allocate(Bk_Smmtrc,3,size(Smmtrc,2),Label='Bk_Smmtrc')
    Bk_Smmtrc(:,:) = Smmtrc(:,:)
    call mma_deallocate(Smmtrc)
  end if
  if (allocated(Lbl)) then
    call mma_allocate(Bk_Lbl,size(Lbl),Label='Bk_Lbl')
    Bk_Lbl(:) = Lbl(:)
    call mma_deallocate(Lbl)
  end if
  if (allocated(mRowH)) then
    call mma_allocate(Bk_mRowH,size(mRowH),Label='Bk_mRowH')
    Bk_mRowH(:) = mRowH(:)
    call mma_deallocate(mRowH)
  end if
  if (allocated(RootMap)) then
    call mma_allocate(Bk_RootMap,size(RootMap),Label='Bk_RootMap')
    Bk_RootMap(:) = RootMap(:)
    call mma_deallocate(RootMap)
  end if
end if

Bk_Header(:) = Header(:)
Bk_iRef = iRef
Bk_NmIter = NmIter
Bk_iter = iter
Bk_nDimBC = nDimBC
Bk_MxItr = MxItr
Bk_Max_Center = Max_Center
Bk_iOptC = iOptC
Bk_mode = mode
Bk_mTROld = mTROld
Bk_nWndw = nWndw
Bk_iOptH = iOptH
Bk_nLambda = nLambda
Bk_nMEP = nMEP
Bk_nBVec = nBVec
Bk_IRC = IRC
Bk_mTtAtm = mTtAtm
Bk_MEPnum = MEPnum
Bk_SlStop = SlStop
Bk_lOld = lOld
Bk_CurviLinear = CurviLinear
Bk_HSet = HSet
Bk_BSet = BSet
Bk_lNmHss = lNmHss
Bk_Cubic = Cubic
Bk_Baker = Baker
Bk_DDV_Schlegel = DDV_Schlegel
Bk_Line_Search = Line_Search
Bk_HWRS = HWRS
Bk_Analytic_Hessian = Analytic_Hessian
Bk_FindTS = FindTS
Bk_MEP = MEP
Bk_User_Def = User_Def
Bk_rMEP = rMEP
Bk_lOld_Implicit = lOld_Implicit
Bk_HrmFrq_Show = HrmFrq_Show
Bk_eMEPtest = eMEPtest
Bk_Redundant = Redundant
Bk_lCtoF = lCtoF
Bk_lSoft = lSoft
Bk_CallLast = CallLast
Bk_TwoRunFiles = TwoRunFiles
Bk_TSConstraints = TSConstraints
Bk_MEPCons = MEPCons
Bk_Track = Track
Bk_Request_Alaska = Request_Alaska
Bk_Request_RASSI = Request_RASSI
Bk_ThrEne = ThrEne
Bk_ThrGrd = ThrGrd
Bk_Beta = Beta
Bk_Beta_Disp = Beta_Disp
Bk_Delta = Delta
Bk_Rtrnc = Rtrnc
Bk_rHidden = rHidden
Bk_ThrCons = ThrCons
Bk_GNrm_Threshold = GNrm_Threshold
Bk_CnstWght = CnstWght
Bk_dMEPStep = dMEPStep
Bk_rFuzz = rFuzz
Bk_ThrMEP = ThrMEP
Bk_lTherm = lTherm
Bk_lDoubleIso = lDoubleIso
Bk_nUserPT = nUserPT
Bk_nsRot = nsRot
Bk_UserT(:) = UserT(:)
Bk_UserP = UserP
Bk_HUpMet = HUpMet
Bk_UpMeth = UpMeth
Bk_MEP_Type = MEP_Type
Bk_MEP_Algo = MEP_Algo
Bk_isFalcon = isFalcon
Bk_mB_Tot = mB_Tot
Bk_mdB_Tot = mdB_Tot
Bk_mq = mq
Bk_WeightedConstraints = WeightedConstraints
Bk_NADC = NADC
Bk_EDiffZero = EDiffZero
Bk_ApproxNADC = ApproxNADC
Bk_iState(:) = iState(:)
Bk_iRow = iRow
Bk_iRow_c = iRow_c
Bk_nFix = nFix
Bk_iInt = iInt
Bk_VarR = VarR
Bk_VarT = VarT

! Make a backup of the runfile, since we are going to change the
! internal coordinates definition.

call fCopy('RUNFILE','RUNBCK2',iErr)
if (iErr /= 0) call Abend()

! Remove the Hessian and disable translational and rotational invariance

call Put_dArray('Hess',rDum,0)
VarR = .true.
VarT = .true.
!                                                                      *
!***********************************************************************
! Call Slapaf to build the B matrix and get the displacement vectors   *
! corresponding to the primitive internal coordinates (bonds and       *
! angles)                                                              *
!***********************************************************************
!                                                                      *
! Process the input

PLback = iPrintLevel(-1)
i = iPrintLevel(0)
LuInput = 21
LuInput = IsFreeUnit(LuInput)
call StdIn_Name(StdIn)
call Molcas_open(LuInput,StdIn)

iRow = 0
iRow_c = 0
nFix = 0
iInt = 0
call RdCtl_Slapaf(LuInput,.true.)
Curvilinear = .true.
Numerical = .false.

close(LuInput)
i = iPrintLevel(PLback)
!                                                                      *
!***********************************************************************
!                                                                      *
BSet = .true.
HSet = .false.
PrQ = .false.
nWndw = iter
iRef = 0
call BMtrx(size(Coor,2),Coor,iter,mTtAtm,nWndw)
nQQ = size(Shift,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! First get the (transposed) K matrix, (dQ/dq)^T

call mma_allocate(KMtrx,mq,nQQ,label='KMtrx')
call mma_allocate(KTrsp,nQQ,mq,label='KTrsp')
call Qpg_dArray('K',Found,nK)
if ((.not. Found) .or. (nK /= mq*nQQ)) call Abend()
call Get_dArray('K',KMtrx,nK)
KTrsp(:,:) = transpose(KMtrx)
call mma_deallocate(KMtrx)

! Form the full B matrix in the redundant internal coordinates
! (primitives) and multiply by KK(t) on the fly

nX = 3*mTtAtm
call mma_allocate(KKtB,mq,nX,label='KKtB')
KKtB(:,:) = Zero
i = 1
do iq=1,mq
  nB = nqBM(iq)
  do iB=i,i+nB-1
    j = iBM(iB)
    do ii=1,mq
      KKtB(ii,j) = KKtB(ii,j)+BM(iB)*DDot_(nQQ,KTrsp(1,ii),1,KTrsp(1,iq),1)
    end do
  end do
  i = i+nB
end do
call mma_deallocate(KTrsp)

! Get the Cartesian normal modes

call Get_iScalar('Unique atoms',nUnique_Atoms)
call Get_nAtoms_All(nAll_Atoms)
if (3*nAll_Atoms /= nX) call Abend()
call mma_allocate(NMod,nX,nModes,label='NMod')
NMod(:,:) = Zero
call Get_NMode_All(Modes,lModes,nModes,nUnique_Atoms,NMod,nAll_Atoms,lDisp)

! Compute the overlaps with the primitive displacements

call mma_allocate(IntMod,mq,nModes,label='IntMod')
call DGeMM_('N','N',mq,nModes,nX,One,KKtB,mq,NMod,nX,Zero,IntMod,mq)
call mma_deallocate(KKtB)
call mma_deallocate(NMod)

! "Normalize" the maximum value to 1

do i=1,nModes
  Mx = Zero
  do j=1,mq
    Mx = max(Mx,abs(IntMod(j,i)))
  end do
  if (Mx > 1.0e-10_wp) IntMod(:,i) = IntMod(:,i)/Mx
end do

! Print the overlaps

call mma_allocate(Label,nLbl,label='Label')
Filename = 'INTCOR'
LuIC = 21
LuIC = IsFreeUnit(LuIC)
call Molcas_open(LuIC,Filename)
i = 1
Line = Get_Ln_EOF(LuIC)
do while (Line /= 'EOF')
  j = index(Line,'=')+1
  Label(i) = adjustl(Line(j:))
  i = i+1
  Line = Get_Ln_EOF(LuIC)
end do
close(LuIC)

MinComp = Half
call CollapseOutput(1,'Principal components of the normal modes')
write(u6,'(3X,A)') '----------------------------------------'
write(u6,*)
write(u6,'(3X,A,F4.2,A)') '(Only contributions larger than ',MinComp,' times the maximum are printed)'
write(u6,*)
call mma_allocate(Sort,mq,label='Sort')
do i=1,nModes
  write(u6,*)
  write(u6,'(6X,A,1X,I6)') 'Mode',i
  write(Line,'(5X,F10.2)') Freq(i)
  j = index(Line,'-')
  if (j > 0) Line(j:j) = 'i'
  write(u6,'(8X,A)') 'Frequency: '//trim(Line)//' cm-1'
  write(u6,'(6X,A)') '---------------------------------'
  do j=1,mq
    Sort(j) = j
  end do
  do j=1,mq
    do ii=j+1,mq
      if (abs(IntMod(Sort(ii),i)) > abs(IntMod(Sort(j),i))) then
        im = Sort(j)
        Sort(j) = Sort(ii)
        Sort(ii) = im
      end if
    end do
    if (abs(IntMod(Sort(j),i)) < MinComp) exit
    write(u6,'(8X,A,F7.4)') Label(Sort(j)),IntMod(Sort(j),i)
  end do
  write(u6,'(6X,A)') '---------------------------------'
end do
call CollapseOutput(0,'Principal components of the normal modes')

! Clean up

call mma_deallocate(Label)
call mma_deallocate(IntMod)
call mma_deallocate(Sort)
!                                                                      *
!***********************************************************************
!                                                                      *
! Restore the runfile and the "global" variables

call fCopy('RUNBCK2','RUNFILE',iErr)
if (iErr /= 0) call Abend()
if (AixRm('RUNBCK2') /= 0) call Abend()

Header(:) = Bk_Header(:)
iRef = Bk_iRef
NmIter = Bk_NmIter
iter = Bk_iter
nDimBC = Bk_nDimBC
MxItr = Bk_MxItr
Max_Center = Bk_Max_Center
iOptC = Bk_iOptC
mode = Bk_mode
mTROld = Bk_mTROld
nWndw = Bk_nWndw
iOptH = Bk_iOptH
nLambda = Bk_nLambda
nMEP = Bk_nMEP
nBVec = Bk_nBVec
IRC = Bk_IRC
mTtAtm = Bk_mTtAtm
MEPnum = Bk_MEPnum
SlStop = Bk_SlStop
lOld = Bk_lOld
CurviLinear = Bk_CurviLinear
HSet = Bk_HSet
BSet = Bk_BSet
lNmHss = Bk_lNmHss
Cubic = Bk_Cubic
Baker = Bk_Baker
DDV_Schlegel = Bk_DDV_Schlegel
Line_Search = Bk_Line_Search
HWRS = Bk_HWRS
Analytic_Hessian = Bk_Analytic_Hessian
FindTS = Bk_FindTS
MEP = Bk_MEP
User_Def = Bk_User_Def
rMEP = Bk_rMEP
lOld_Implicit = Bk_lOld_Implicit
HrmFrq_Show = Bk_HrmFrq_Show
eMEPtest = Bk_eMEPtest
Redundant = Bk_Redundant
lCtoF = Bk_lCtoF
lSoft = Bk_lSoft
CallLast = Bk_CallLast
TwoRunFiles = Bk_TwoRunFiles
TSConstraints = Bk_TSConstraints
MEPCons = Bk_MEPCons
Track = Bk_Track
Request_Alaska = Bk_Request_Alaska
Request_RASSI = Bk_Request_RASSI
ThrEne = Bk_ThrEne
ThrGrd = Bk_ThrGrd
Beta = Bk_Beta
Beta_Disp = Bk_Beta_Disp
Delta = Bk_Delta
Rtrnc = Bk_Rtrnc
rHidden = Bk_rHidden
ThrCons = Bk_ThrCons
GNrm_Threshold = Bk_GNrm_Threshold
CnstWght = Bk_CnstWght
dMEPStep = Bk_dMEPStep
rFuzz = Bk_rFuzz
ThrMEP = Bk_ThrMEP
lTherm = Bk_lTherm
lDoubleIso = Bk_lDoubleIso
nUserPT = Bk_nUserPT
nsRot = Bk_nsRot
UserT(:) = Bk_UserT(:)
UserP = Bk_UserP
HUpMet = Bk_HUpMet
UpMeth = Bk_UpMeth
MEP_Type = Bk_MEP_Type
MEP_Algo = Bk_MEP_Algo
isFalcon = Bk_isFalcon
mB_Tot = Bk_mB_Tot
mdB_Tot = Bk_mdB_Tot
mq = Bk_mq
WeightedConstraints = Bk_WeightedConstraints
NADC = Bk_NADC
EDiffZero = Bk_EDiffZero
ApproxNADC = Bk_ApproxNADC
iState(:) = Bk_iState(:)
iRow = Bk_iRow
iRow_c = Bk_iRow_c
nFix = Bk_nFix
iInt = Bk_iInt
VarR = Bk_VarR
VarT = Bk_VarT

! Process arrays that are always allocated.

if (allocated(Bk_Cx)) then
  Cx(:,:,:) = Bk_Cx(:,:,:)
  call mma_deallocate(Bk_Cx)
else
  call mma_deallocate(Cx)
end if
if (allocated(Bk_Gx)) then
  Gx(:,:,:) = Bk_Gx(:,:,:)
  call mma_deallocate(Bk_Gx)
else
  call mma_deallocate(Gx)
end if
if (allocated(Bk_Gx0)) then
  Gx0(:,:,:) = Bk_Gx0(:,:,:)
  call mma_deallocate(Bk_Gx0)
else
  call mma_deallocate(Gx0)
end if
if (allocated(Bk_NAC)) then
  NAC(:,:,:) = Bk_NAC(:,:,:)
  call mma_deallocate(Bk_NAC)
else
  call mma_deallocate(NAC)
end if
if (allocated(Bk_Q_nuclear)) then
  Q_nuclear(:) = Bk_Q_nuclear(:)
  call mma_deallocate(Bk_Q_nuclear)
else
  call mma_deallocate(Q_nuclear)
end if
if (allocated(Bk_dMass)) then
  dMass(:) = Bk_dMass(:)
  call mma_deallocate(Bk_dMass)
else
  call mma_deallocate(dMass)
end if
if (allocated(Bk_Coor)) then
  Coor(:,:) = Bk_Coor(:,:)
  call mma_deallocate(Bk_Coor)
else
  if (allocated(Coor)) call mma_deallocate(Coor)
end if
if (allocated(Bk_Grd)) then
  Grd(:,:) = Bk_Grd(:,:)
  call mma_deallocate(Bk_Grd)
else
  if (allocated(Grd)) call mma_deallocate(Grd)
end if
if (allocated(Bk_ANr)) then
  ANr(:) = Bk_ANr(:)
  call mma_deallocate(Bk_ANr)
else
  call mma_deallocate(ANr)
end if
if (allocated(Bk_Weights)) then
  Weights(:) = Bk_Weights(:)
  call mma_deallocate(Bk_Weights)
else
  call mma_deallocate(Weights)
end if
if (allocated(Bk_Shift)) then

  if (size(Shift,1) /= size(Bk_Shift,1)) then
    call mma_deallocate(Shift)
    call mma_allocate(Shift,size(Bk_Shift,1),MaxItr,Label='Shift')
  end if
  Shift(:,:) = Bk_Shift(:,:)
  call mma_deallocate(Bk_Shift)
else
  call mma_deallocate(Shift)
end if
if (allocated(Bk_GNrm)) then
  GNrm(:) = Bk_GNrm(:)
  call mma_deallocate(Bk_GNrm)
else
  call mma_deallocate(GNrm)
end if
if (allocated(Bk_Energy)) then
  Energy(:) = Bk_Energy(:)
  call mma_deallocate(Bk_Energy)
else
  call mma_deallocate(Energy)
end if
if (allocated(Bk_Energy0)) then
  Energy0(:) = Bk_Energy0(:)
  call mma_deallocate(Bk_Energy0)
else
  call mma_deallocate(Energy0)
end if
if (allocated(Bk_MF)) then
  MF(:,:) = Bk_MF(:,:)
  call mma_deallocate(Bk_MF)
else
  call mma_deallocate(MF)
end if
if (allocated(Bk_DipM)) then
  DipM(:,:) = Bk_DipM(:,:)
  call mma_deallocate(Bk_DipM)
else
  if (allocated(DipM)) call mma_deallocate(DipM)
end if
if (allocated(Bk_RefGeo)) then
  RefGeo(:,:) = Bk_RefGeo(:,:)
  call mma_deallocate(Bk_RefGeo)
else
  if (allocated(RefGeo)) call mma_deallocate(RefGeo)
end if

! Process arrays that are allocated optionally.

if (allocated(Bk_Lambda)) then
  if (.not. allocated(Lambda)) call mma_allocate(Lambda,size(Bk_Lambda,1),MaxItr+1,Label='Lambda')
  Lambda(:,:) = Bk_Lambda(:,:)
  call mma_deallocate(Bk_Lambda)
else
  if (allocated(Lambda)) call mma_deallocate(Lambda)
end if

if (allocated(Bk_qInt)) then
  if (allocated(qInt)) call mma_deallocate(qInt)
  call mma_allocate(qInt,size(Bk_qInt,1),MaxItr,Label='qInt')
  qInt(:,:) = Bk_qInt(:,:)
  call mma_deallocate(Bk_qInt)
else
  if (allocated(qInt)) call mma_deallocate(qInt)
end if
if (allocated(Bk_dqInt)) then
  if (allocated(dqInt)) call mma_deallocate(dqInt)
  call mma_allocate(dqInt,size(Bk_dqInt,1),MaxItr,Label='dqInt')
  dqInt(:,:) = Bk_dqInt(:,:)
  call mma_deallocate(Bk_dqInt)
else
  if (allocated(dqInt)) call mma_deallocate(dqInt)
end if

if (allocated(Bk_BM)) then
  if (allocated(BM)) call mma_deallocate(BM)
  call mma_allocate(BM,size(Bk_BM),Label='BM')
  BM(:) = Bk_BM(:)
  call mma_deallocate(Bk_BM)
else
  if (allocated(BM)) call mma_deallocate(BM)
end if
if (allocated(Bk_dBM)) then
  if (allocated(dBM)) call mma_deallocate(dBM)
  call mma_allocate(dBM,size(Bk_dBM),Label='dBM')
  dBM(:) = Bk_dBM(:)
  call mma_deallocate(Bk_dBM)
else
  if (allocated(dBM)) call mma_deallocate(dBM)
end if
if (allocated(Bk_iBM)) then
  if (allocated(iBM)) call mma_deallocate(iBM)
  call mma_allocate(iBM,size(Bk_iBM),Label='iBM')
  iBM(:) = Bk_iBM(:)
  call mma_deallocate(Bk_iBM)
else
  if (allocated(iBM)) call mma_deallocate(iBM)
end if
if (allocated(Bk_idBM)) then
  if (allocated(idBM)) call mma_deallocate(idBM)
  call mma_allocate(idBM,size(Bk_idBM),Label='idBM')
  idBM(:) = Bk_idBM(:)
  call mma_deallocate(Bk_idBM)
else
  if (allocated(idBM)) call mma_deallocate(idBM)
end if
if (allocated(Bk_nqBM)) then
  if (allocated(nqBM)) call mma_deallocate(nqBM)
  call mma_allocate(nqBM,size(Bk_nqBM),Label='nqBM')
  nqBM(:) = Bk_nqBM(:)
  call mma_deallocate(Bk_nqBM)
else
  if (allocated(nqBM)) call mma_deallocate(nqBM)
end if
if (allocated(Bk_BMx)) then
  if (allocated(BMx)) call mma_deallocate(BMx)
  call mma_allocate(BMx,size(Bk_BMx,1),size(Bk_BMx,2),Label='BMx')
  BMx(:,:) = Bk_BMx(:,:)
  call mma_deallocate(Bk_BMx)
else
  if (allocated(BMx)) call mma_deallocate(BMx)
end if
if (allocated(Bk_Degen)) then
  if (allocated(Degen)) call mma_deallocate(Degen)
  call mma_allocate(Degen,size(Bk_Degen,1),size(Bk_Degen,2),Label='Degen')
  Degen(:,:) = Bk_Degen(:,:)
  call mma_deallocate(Bk_Degen)
else
  if (allocated(Degen)) call mma_deallocate(Degen)
end if
if (allocated(Bk_jStab)) then
  if (allocated(jStab)) call mma_deallocate(jStab)
  call mma_allocate(jStab,[0,7],[1,size(Bk_jStab,2)],Label='jStab')
  jStab(:,:) = Bk_jStab(:,:)
  call mma_deallocate(Bk_jStab)
else
  if (allocated(jStab)) call mma_deallocate(jStab)
end if
if (allocated(Bk_iCoSet)) then
  if (allocated(iCoSet)) call mma_deallocate(iCoSet)
  call mma_allocate(iCoSet,[0,7],[1,size(Bk_iCoSet,2)],Label='iCoSet')
  iCoSet(:,:) = Bk_iCoSet(:,:)
  call mma_deallocate(Bk_iCoSet)
else
  if (allocated(iCoSet)) call mma_deallocate(iCoSet)
end if
if (allocated(Bk_nStab)) then
  if (allocated(nStab)) call mma_deallocate(nStab)
  call mma_allocate(nStab,size(Bk_nStab,1),Label='nStab')
  nStab(:) = Bk_nStab(:)
  call mma_deallocate(Bk_nStab)
else
  if (allocated(nStab)) call mma_deallocate(nStab)
end if
if (allocated(Bk_AtomLbl)) then
  if (allocated(AtomLbl)) call mma_deallocate(AtomLbl)
  call mma_allocate(AtomLbl,size(Bk_AtomLbl,1),Label='AtomLbl')
  AtomLbl(:) = Bk_AtomLbl(:)
  call mma_deallocate(Bk_AtomLbl)
else
  if (allocated(AtomLbl)) call mma_deallocate(AtomLbl)
end if
if (allocated(Bk_Smmtrc)) then
  if (allocated(Smmtrc)) call mma_deallocate(Smmtrc)
  call mma_allocate(Smmtrc,3,size(Bk_Smmtrc,2),Label='Smmtrc')
  Smmtrc(:,:) = Bk_Smmtrc(:,:)
  call mma_deallocate(Bk_Smmtrc)
else
  if (allocated(Smmtrc)) call mma_deallocate(Smmtrc)
end if
if (allocated(Bk_Lbl)) then
  if (allocated(Lbl)) call mma_deallocate(Lbl)
  call mma_allocate(Lbl,size(Bk_Lbl),Label='Lbl')
  Lbl(:) = Bk_Lbl(:)
  call mma_deallocate(Bk_Lbl)
else
  if (allocated(Lbl)) call mma_deallocate(Lbl)
end if
if (allocated(Bk_mRowH)) then
  if (allocated(mRowH)) call mma_deallocate(mRowH)
  call mma_allocate(mRowH,size(Bk_mRowH),Label='mRowH')
  mRowH(:) = Bk_mRowH(:)
  call mma_deallocate(Bk_mRowH)
else
  if (allocated(mRowH)) call mma_deallocate(mRowH)
end if
if (allocated(Bk_RootMap)) then
  if (allocated(RootMap)) call mma_deallocate(RootMap)
  call mma_allocate(RootMap,size(Bk_RootMap),Label='RootMap')
  RootMap(:) = Bk_RootMap(:)
  call mma_deallocate(Bk_RootMap)
else
  if (allocated(RootMap)) call mma_deallocate(RootMap)
end if

end subroutine Print_Mode_Components
