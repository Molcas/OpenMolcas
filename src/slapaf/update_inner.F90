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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine Update_inner(kIter,Beta,Beta_Disp,Step_Trunc,nWndw,mIter,Kriging_Hessian,qBeta,iOpt_RS,First_MicroIteration,Iter, &
                        qBeta_Disp,Hide)
!***********************************************************************
!     Object: to update coordinates                                    *
!                                                                      *
!    Input:                                                            *
!      kIter          : iteration counter                              *
!      Beta           : damping factor step length                     *
!      Beta_Disp      : damping factor variance                        *
!                                                                      *
!    OutPut:                                                           *
!      Step_Trunc     : character label to denote truncation type      *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh                                             *
!             2000                                                     *
!***********************************************************************

use Slapaf_Info, only: BMx, Curvilinear, Degen, dqInt, Energy, FindTS, GNrm, GNrm_Threshold, GrdMax, HrmFrq_Show, iInt, iNeg, &
                       iOptC, iOptH, iRow_c, Lambda, Lbl, Mode, nBVec, nDimBC, nFix, nLambda, nStab, qInt, Shift, Smmtrc, StpMax, &
                       TSConstraints
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Five, Ten, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: kIter, nWndw, mIter, iOpt_RS, Iter
real(kind=wp), intent(inout) :: Beta
real(kind=wp), intent(in) :: Beta_Disp, qBeta_Disp
character, intent(out) :: Step_Trunc
logical(kind=iwp), intent(in) :: Kriging_Hessian, First_MicroIteration, Hide
real(kind=wp), intent(out) :: qBeta
integer(kind=iwp) :: i, iAd, iAtom, iDo_DipM, iDum(1), iEnd, iIter, iLambda, iOptH_, ix, ixyz, j, jAtom, jPrint, jx, jxyz, k, &
                     kStart, lIter, LudRdx, M, Mode_Save, N, n1, nA, nLoop, nQQ, NRHS, nRP, nsAtom
real(kind=wp) :: Disp(1), dxdx, dxdx_last, fact, fCart, gBeta, gg, gg_last, rCart, rInter, Sf, tBeta, Thr, Thr_RS, xBeta
logical(kind=iwp) :: Found, lWrite
character(len=8) :: File1, File2
character :: Step_Trunc_
integer(kind=iwp), allocatable :: iFlip(:), Indx(:)
real(kind=wp), allocatable :: AMat(:), BM(:,:), BVec(:,:), CInt(:), CInt0(:), d2L(:,:,:), dBM(:,:,:), dBVec(:), dEdq(:,:), &
                              dEdx(:,:), dRdq(:,:,:), du(:), dx(:,:), dy(:), EMtrx(:,:), Energy_L(:), ErrVec(:,:), Hessian(:,:), &
                              Mult(:), R(:,:), RHS(:), Scr1(:), Scr2(:), T(:,:), Tmp(:), Val(:), Value0(:), Wess(:,:), x(:,:)
character(len=8), allocatable :: Lbl_Tmp(:)
real(kind=wp), external :: DDot_
!                                                                      *
!***********************************************************************
!                                                                      *
! This option should be more carefully debugged.
!                                              RL December 2020
#define _NEW_CODE_
#ifdef _NEW_CODE_
real(kind=wp), allocatable :: QC(:,:,:)
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
nQQ = size(qInt,1)
nsAtom = size(Degen,2)
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Lu = u6
write(Lu,*) 'Update_inner:iOpt_RS,Beta,Beta_Disp=',iOpt_RS,Beta,Beta_Disp
call RecPrt('Update_inner: qInt',' ',qInt,nQQ,kIter)
call RecPrt('Update_inner: Shift',' ',Shift,nQQ,kIter-1)
call RecPrt('Update_inner: GNrm',' ',GNrm,kIter,1)
call RecPrt('Update_inner: Energy',' ',Energy,kIter,1)
n1 = 3*nsAtom
call RecPrt('Update_inner: dQ/dx(BMx)',' ',BMx,n1,nQQ)
#endif

GrdMax = Zero
StpMax = Zero

Step_Trunc = 'N'
nA = (max(nQQ,kIter)+1)**2
call mma_Allocate(AMat,nA,Label='AMat')
call mma_Allocate(Indx,kIter,Label='Indx')
call mma_Allocate(ErrVec,nQQ,(kIter+1),Label='ErrVec')
call mma_Allocate(EMtrx,kIter+1,kIter+1,Label='EMtrx')
call mma_Allocate(RHS,kIter+1)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Pick up the Hessian in internal coordinates and then update it
! according to some Hessian update method (BFGS, MSP, etc.)

if (Kriging_Hessian) then
  iOptH_ = ibset(0,3)
else
  iOptH_ = iOptH
end if
call mma_Allocate(Hessian,nQQ,nQQ,Label='Hessian')
call Get_dArray('Hss_Q',Hessian,nQQ**2)

! Perform the Hessian update, in case of GEK it essentially will
! modify the Hessian if it is needed to guide 2nd order
! optimization towards a minimum or a TS.

#ifdef _DEBUGPRINT_
write(Lu,*)
write(Lu,*)
write(Lu,*) ' *** Updating the molecular Hessian ***'
write(Lu,*)
jPrint = 99
#else
jPrint = 5
#endif
call Update_H(nWndw,Hessian,nQQ,mIter,iOptC,Shift(:,kIter-mIter+1),dqInt(:,kIter-mIter+1),iOptH_,jPrint,GNrm(kIter),nsAtom,.true., &
              First_MicroIteration)

!call RecPrt('Update_inner: Hessian',' ',Hessian,nQQ,nQQ)
!write(u6,*) 'After corrections'
!call DiagMtrx(Hessian,nQQ,jNeg)

! Save the number of internal coordinates on the runfile.

call Put_iScalar('No of Internal coordinates',nQQ)

! Save the force constant matrix on the runfile.

call Put_dArray('Hess',Hessian,nQQ**2)

! Optional frequency analysis

if (HrmFrq_Show) call GF_on_the_fly(iDo_DipM)

!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! If this is a constrained optimization with the goal of finding a
! transition state and we have one negative eigenvalue of the
! Hessian then shift over to Mode Following RF to find the TS.

! If TSConstraints were given, merge them with global constraints if
! they exist, unless we are in TS regime.
! If TS constraints were not given, remove global constraints when in
! TS regime.

call qpg_darray('TanVec',Found,nRP)
if (FindTS .and. First_MicroIteration) then
  File1 = 'UDC'
  File2 = 'TSC'
  if (.not. TSConstraints) File2 = ''
  if (iNeg(1) >= 1) then
    if ((GNrm(kIter) <= GNrm_Threshold) .or. Found) then
      ! Change to MFRF optimization.
      iOptC = ibclr(iOptC,7)
      iOptC = ibclr(iOptC,12)
      call Put_lScalar('TS Search',.true.)
      if (TSConstraints) then
        File2 = ''
      else
        File1 = ''
      end if
    end if
  end if
  call Merge_Constraints(File1,File2,'UDC',nLambda,iRow_c)
  call Fix_UDC(iRow_c,nLambda,nsAtom,nStab,.false.)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (nLambda == 0) then
  M = 3*nsAtom
  N = nQQ
  NRHS = 1
  call mma_Allocate(Tmp,M,Label='Tmp:M')

  ! Iterate to avoid too large displacement in cartesians.
  ! If the maximum displacement is more than 2*Beta, reduce the step.

  ! Initial setup to ensure fCart=1.0 at first iteration.
  rInter = Beta
  fCart = Ten
  rCart = fCart*rInter
  nLoop = 0
  do while (rCart >= Two*Beta)
    nLoop = nLoop+1
    if (nLoop > 10) exit
    if (rCart > rInter) then
      fCart = fCart*rInter/rCart
    else
      fCart = fCart*0.9_wp
    end if
    nA = (max(nQQ,kIter)+1)**2
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    gBeta = One
    xBeta = One
    gg_last = Beta
    dxdx_last = Beta
    Sf = sqrt(Two)
    if (iOpt_RS == 0) then
      kStart = max(1,kIter-4)
      iEnd = kIter
    else
      kStart = max(1,Iter-4)
      iEnd = Iter
    end if
    Thr = 1.0e-6_wp
    do iIter=kStart,iEnd

      if (iIter /= kIter) then
        dxdx = sqrt(DDot_(nQQ,Shift(:,iIter),1,Shift(:,iIter),1))

        if ((dxdx < 0.75_wp*dxdx_last) .and. (dxdx < Beta-Thr)) then
          ! Increase trust radius
          !xBeta = xBeta*Sf
          xBeta = min(Two,xBeta*Sf)
        else if ((dxdx > 1.25_wp*dxdx_last) .or. (dxdx >= Beta+Thr)) then
          ! Reduce trust radius
          xBeta = max(One/Five,xBeta/Sf)
        end if
        dxdx_last = dxdx
      end if

      gg = sqrt(DDot_(nQQ-nLambda,dqInt(:,iIter),1,dqInt(:,iIter),1))
      if ((gg < 0.75_wp*gg_last) .and. (gg < Beta-Thr)) then
        ! Increase trust radius
        !gBeta = gBeta*Sf
        gBeta = min(Two,gBeta*Sf)
      else if ((gg > 1.25_wp*gg_last) .or. (gg >= Beta+Thr)) then
        ! Reduce trust radius
        gBeta = max(One/Five,gBeta/Sf)
      end if
      gg_last = gg

    end do
    tBeta = max(Beta*min(xBeta,gBeta),Beta/Ten)
    !write(u6,*) 'tBeta=',tBeta
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute updated geometry in Internal coordinates

    fact = One
    qBeta = fCart*tBeta
    Thr_RS = 1.0e-7_wp
    do
      Step_Trunc_ = Step_Trunc
      call Newq(qInt,nQQ,kIter,Shift,Hessian,dqInt,ErrVec,EMtrx,RHS,AMat,nA,qBeta,nFix,Indx,Energy,Step_Trunc_,Thr_RS)
      if (Step_Trunc == 'N') Step_Trunc = ' '
      if (iOpt_RS == 0) then
        if (Step_Trunc_ == 'N') Step_Trunc_ = ' '
        Step_Trunc = Step_Trunc_
        exit
      end if
      if (Step_Trunc//Step_Trunc_ == ' *') Step_Trunc = '.'

      qInt(:,kIter+1) = qInt(:,kIter)+Shift(:,kIter)
      call Dispersion_Kriging_Layer(qInt(:,kIter+1),Disp,nQQ)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Disp,Beta_Disp=',Disp(1),Beta_Disp
#     endif
      fact = Half*fact
      qBeta = Half*qBeta
      if (One-disp(1)/Beta_Disp > 1.0e-3_wp) exit
      if ((fact < 1.0e-5_wp) .or. (disp(1) < Beta_Disp)) exit
      Step_Trunc = '*'
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Step_Trunc=',Step_Trunc
#   endif

    call MxLbls(nQQ,dqInt(:,kIter),Shift(:,kIter),Lbl)

    ! Set shift vector to zero for frozen internal coordinates.

    if (nFix > 0) Shift(iInt+1:iInt+nFix,kIter) = Zero

    ! Rough conversion to Cartesians

    call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,Shift(:,kIter),Tmp)
    rInter = sqrt(dDot_(N,Shift(:,kIter),1,Shift(:,kIter),1))
    rCart = Zero
    do i=1,nsAtom
      rCart = max(rCart,sqrt(dDot_(3,Tmp(i),nsAtom,Tmp(i),nsAtom)))
    end do
  end do
  call mma_deallocate(Tmp)
  call Put_dScalar('Max error',Zero)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Optimization with constraints using Lagrangian technique.

  ! Allocate memory for the new arrays

  call mma_allocate(R,nLambda,kIter,Label='R')
  call mma_allocate(dRdq,nQQ,nLambda,kIter,Label='dRdq')
  dRdq(:,:,:) = Zero
  call mma_allocate(T,nQQ,nQQ,Label='T')
  T(:,:) = Zero
  call mma_allocate(dy,nLambda,Label='dy')
  call mma_allocate(dx,(nQQ-nLambda),kIter,Label='dx')
  dx(:,:) = Zero
  call mma_allocate(dEdq,nQQ,kIter,Label='dEdQ')
  call mma_allocate(du,nQQ,Label='du')
  call mma_allocate(X,(nQQ-nLambda),(kIter+1),Label='X')
  X(:,:) = Zero
  call mma_allocate(dEdx,nQQ-nLambda,kIter,Label='dEdx')
  call mma_allocate(Wess,nQQ-nLambda,nQQ-nLambda,Label='Wess')
  Wess(:,:) = Zero
  call mma_allocate(Energy_L,kIter,Label='Energy_L')

  Energy_L(:) = Energy(1:kIter)

  if (nQQ+nLambda > size(Lbl)) then
    call mma_allocate(Lbl_tmp,nQQ+nLambda,Label='Lbl_tmp')
    Lbl_tmp(:) = ''
    Lbl_tmp(1:size(Lbl)) = Lbl(:)
    call mma_deallocate(Lbl)
    call mma_allocate(Lbl,nQQ+nLambda,Label='Lbl')
    Lbl(:) = Lbl_tmp(:)
    call mma_deallocate(Lbl_tmp)
  end if

  nBVec = iRow_c-nLambda-1

  n1 = 3*nsAtom
  call mma_allocate(BM,n1,nLambda,Label='BM')
  call mma_allocate(dBM,n1,n1,nLambda,Label='dBM')

  call mma_allocate(BVec,3*nsAtom,nBVec,Label='BVec')
  call mma_allocate(Val,nBVec,Label='Val')
  call mma_allocate(Value0,nBVec,Label='Value0')
  call mma_allocate(CInt,nLambda,Label='CInt')
  call mma_allocate(CInt0,nLambda,Label='CInt0')
  call mma_allocate(Mult,nBVec**2,Label='Mult')
  call mma_allocate(iFlip,nBVec,Label='iFlip')
  call mma_allocate(dBVec,nBVec*(3*nsAtom)**2,Label='dBVec')
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the constraints
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do lIter=1,kIter
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    lWrite = (lIter == kIter) .and. First_MicroIteration
    lWrite = lWrite .and. (.not. Hide)
    call DefInt2(BVec,dBVec,nBVec,BM,nLambda,nsAtom,iRow_c,Val,cInt,cInt0,Lbl(nQQ+1),lWrite,Mult,dBM,Value0,lIter,iFlip)

    ! Assemble r

    R(:,lIter) = cInt(:)-cInt0(:)

    ! Assemble dC/dQ: Solve  B dC/dQ = dC/dx

    ! where B = dQ/dx

    dRdq(:,:,lIter) = Zero
#   ifdef _DEBUGPRINT_
    write(Lu,*) 'Update_inner: lIter=',lIter
    call RecPrt('Update_inner: dQ/dx(BMx)',' ',BMx,n1,nQQ)
    call RecPrt('Update_inner: dC/dx(BM)',' ',BM,n1,nLambda)
#   endif

    M = n1
    N = nQQ
    NRHS = nLambda

    ! Temporary fix of the dC/dx vector which always
    ! is propted up with the full degeneracy factor.

    if (.not. Curvilinear) then
      do iLambda=1,nLambda
        do i=1,n1
          iAtom = (i+2)/3
          ixyz = i-(iAtom-1)*3
          BM(i,iLambda) = BM(i,iLambda)/Degen(ixyz,iAtom)
        end do
      end do
    end if
    if (lIter == kIter) then
      LudRdX = 30
      call DaName(LudRdX,'dRdX')
      iAd = 0
      iDum(1) = nLambda
      call iDaFile(LudRdX,1,iDum,1,iAd)
      iDum(1) = n1
      call iDaFile(LudRdX,1,iDum,1,iAd)
      call dDaFile(LudRdX,1,BM,nLambda*n1,iAd)
      call DaClos(LudRdX)
    end if
    call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,BM,dRdq(:,:,lIter))
#   ifdef _DEBUGPRINT_
    call RecPrt('Update_inner: dRdq(:,:,lIter)',' ',dRdq(:,:,lIter),nQQ,nLambda)
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do     ! lIter
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_deallocate(dBVec)
  call mma_deallocate(iFlip)
  call mma_deallocate(Mult)
  call mma_deallocate(cInt0)
  call mma_deallocate(cInt)
  call mma_deallocate(Value0)
  call mma_deallocate(Val)
  call mma_deallocate(BVec)

# ifdef _DEBUGPRINT_
  call RecPrt('Update_inner: R',' ',R,nLambda,kIter)
  call RecPrt('Update_inner: dRdq',' ',dRdq,nQQ*nLambda,kIter)
  do i=1,nLambda
    write(u6,*) ' iLambda=',i
    call RecPrt('Update_inner: d2C/dx2',' ',dBM(:,:,i),n1,n1)
  end do
  if (Curvilinear) call dBPrint(nQQ,nDimBC)
# endif
  call mma_deallocate(BM)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Assemble d^2C/dQ^2
  !
  !  The expression is written as
  !
  !  dQ/dx * d^2C/dQ^2 * (dQ/dx)^T = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
  !
  !  Dimensions
  !  dQ/dx     :  3*nsAtom x nQQ
  !  d^2C/dx^2 :  3*nsatom x 3*nsatom
  !  d^2C/dQ^2 :  nQQ   x nQQ

# ifdef _NEW_CODE_
  ! Compute Sum (d^2Q/dx^2 * dC/dQ)

  call mma_allocate(QC,nDimBC,nDimBC,nLambda,Label='QC')
  QC(:,:,:) = Zero
  if (Curvilinear) call dBMult(dRdq(:,:,kIter),QC,nQQ,nDimBC,nLambda)
# ifdef _DEBUGPRINT_
  write(Lu,*) 'Update_inner: kIter=',kIter
  call RecPrt('dRdq(:,1,kIter)',' ',dRdq(:,1,kIter),nQQ,1)
  do iLambda=1,nLambda
    write(u6,*) 'Update_inner: iLambda=',iLambda
    call RecPrt('Update_inner: QC',' ',QC(:,:,iLambda),nDimBC,nDimBC)
  end do
# endif

  ! Subtract the term from d^2C/dx^2

  do k=1,nLambda
    j = 0
    do jx=1,n1
      jAtom = (jx+2)/3
      jxyz = jx-(jAtom-1)*3
      if (Smmtrc(jxyz,jAtom)) then
        j = j+1

        i = 0
        do ix=1,n1
          iAtom = (ix+2)/3
          ixyz = ix-(iAtom-1)*3
          if (Smmtrc(ixyz,iAtom)) then
            i = i+1
            dBM(ix,jx,k) = dBM(ix,jx,k)-QC(i,j,k)
          end if
        end do
      end if
    end do
  end do
  call mma_deallocate(QC)
# endif

  call mma_allocate(Scr2,nQQ*n1*nLambda,Label='Scr2')
  call mma_allocate(Scr1,nQQ*n1*nLambda,Label='Scr1')
# ifdef _DEBUGPRINT_
  call RecPrt('Update_inner: d^2C/dx^2(dBM)',' ',dBM,n1**2,nLambda)
# endif

  ! Temporary fix of the d^2C/dx^2 vector which always
  ! is propted up with the full degeneracy factor.

  if (.not. Curvilinear) then
    do k=1,nLambda
      do j=1,n1
        jAtom = (j+2)/3
        jxyz = j-(jAtom-1)*3
        do i=1,n1
          iAtom = (i+2)/3
          ixyz = i-(iAtom-1)*3
          dBM(i,j,k) = dBM(i,j,k)/(Degen(ixyz,iAtom)*Degen(jxyz,jAtom))
        end do
      end do
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('Update_inner: d^2C/dx^2(dBM)',' ',dBM,n1**2,nLambda)
# endif

  ! Solve first for y in
  !         dQ/dx y = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
  !
  ! where y   = d^2C/dQ^2 * (dQ/dx)^T
  !
  ! and   y^T = dQ/dx * d^2C/dQ^2

  M = 3*nsAtom
  N = nQQ
  NRHS = n1*nLambda
  call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,dBM,Scr1)

  ! Generate y^T in Scr2

  call TRNSPS(nQQ,n1*nLambda,Scr1,Scr2)
# ifdef _DEBUGPRINT_
  call RecPrt('d^2C/dQ^2 * (dQ/dx)^T',' ',Scr1,nQQ,n1*nLambda)
  call RecPrt('dQ/dx * d^2C/dQ^2',' ',Scr2,n1*nLambda,nQQ)
# endif
  call mma_deallocate(dBM)

  ! Followed by solve for x in B x = y^T for which B = dQ/dx

  NRHS = nLambda*nQQ
  call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,Scr2,Scr1)
  call mma_allocate(d2L,nQQ,nQQ,nLambda,Label='d2L')

  call TRNSPS(nQQ*nLambda,nQQ,Scr1,d2L)
# ifdef _DEBUGPRINT_
  call RecPrt('Scr1',' ',Scr1,nQQ*nLambda,nQQ)
  do i=1,nLambda
    write(u6,*) ' iLambda=',i
    call RecPrt('Update_inner: d2L',' ',d2L(:,:,i),nQQ,nQQ)
  end do
# endif
  call mma_deallocate(Scr2)
  call mma_deallocate(Scr1)
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Compute updated geometry in Internal coordinates

  Mode_Save = Mode
  Mode = 0
  M = 3*nsAtom
  N = nQQ
  NRHS = 1
  call mma_Allocate(Tmp,M,Label='Tmp:M')

  ! Iterate to avoid too large displacement in cartesians.
  ! If the maximum displacement is more than 2*Beta, reduce the step.

  ! Initial setup to ensure fCart=1.0 at first iteration.
  Thr_RS = 1.0e-7_wp
  rInter = Beta
  fCart = Ten
  rCart = fCart*rInter
  nLoop = 0
  do while (rCart >= Two*Beta)
    nLoop = nLoop+1
    if (nLoop > 100) exit
    if (rCart > rInter) then
      fCart = fCart*rInter/rCart
    else
      fCart = fCart*0.9_wp
    end if
    call Con_Opt(R,dRdq,T,dqInt,Lambda,qInt,Shift,dy,dx,dEdq,du,x,dEdx,Wess,GNrm(kIter),nWndw,Hessian,nQQ,kIter,iOptH_,jPrint, &
                 Energy_L,nLambda,ErrVec,EMtrx,RHS,AMat,nA,Beta,fCart,qBeta_Disp,nFix,Indx,Step_Trunc,Lbl,d2L,nsAtom,iOpt_RS, &
                 Thr_RS,iter,First_Microiteration)
    qBeta = fCart*Beta
    if (iOpt_RS == 1) exit

    ! Rough conversion to Cartesians

    call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,Shift(:,kIter),Tmp)
    rInter = sqrt(dDot_(N,Shift(:,kIter),1,Shift(:,kIter),1))
    rCart = Zero
    do i=1,nsAtom
      rCart = max(rCart,sqrt(dDot_(3,Tmp(i),nsAtom,Tmp(i),nsAtom)))
    end do
  end do
  Mode = Mode_Save
  call mma_deallocate(Tmp)

# ifdef _DEBUGPRINT_
  write(Lu,*)
  write(Lu,*) '********************************************'
  write(Lu,*) '* Lagrange multipliers for the constraints *'
  write(Lu,*) '********************************************'
  write(Lu,'(1X,A,2X,ES13.6)') (Lbl(nQQ+iInt),-One*Lambda(iInt,mIter),iInt=1,nLambda)
  write(Lu,*)
# endif

  call mma_deallocate(d2L)
  call mma_deallocate(Energy_L)
  call mma_deallocate(Wess)
  call mma_deallocate(dEdx)
  call mma_deallocate(x)
  call mma_deallocate(du)
  call mma_deallocate(dEdq)
  call mma_deallocate(dx)
  call mma_deallocate(dy)
  call mma_deallocate(T)
  call mma_deallocate(dRdq)
  call mma_deallocate(R)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_Deallocate(Hessian)
call mma_Deallocate(RHS)
call mma_Deallocate(EMtrx)
call mma_Deallocate(ErrVec)
call mma_Deallocate(Indx)
call mma_Deallocate(AMat)
#ifdef _WARNING_WORKAROUND_
if (Smmtrc(1,1)) i = nDimBC
#endif

return

end subroutine Update_inner
