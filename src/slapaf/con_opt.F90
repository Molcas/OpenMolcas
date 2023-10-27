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
! Copyright (C) 2003,2020, Roland Lindh                                *
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Con_Opt(r,drdq,T,dEdq,rLambda,q,dq,dy,dx,hql,du,x,dEdx,W,GNrm,nWndw,Hess,nInter,nIter,iOptH,jPrint,Energy,nLambda,Err, &
                   EMx,RHS,A,nA,cBeta,fCart,Beta_Disp_Seed,nFix,iP,Step_Trunc,Lbl,d2rdq2,nsAtom,iOpt_RS,Thr_RS,iter_, &
                   First_Microiteration)
!***********************************************************************
!                                                                      *
!     Object: to perform an constrained optimization. The constraints  *
!             are optimized directly as a linear problem in the        *
!             subspace defined by the gradients of the constraints.    *
!             The minimization is performed in the complemental        *
!             subspace with the ordinary routines.                     *
!                                                                      *
!             L(q,l) = E(q) - l^T r(q)                                 *
!                                                                      *
!             Ref:                                                     *
!             "A Reduced-Restricted-Quasi-Newton-Raphson Method for    *
!              Locating and Optimizing Energy Crossing Points Between  *
!              Two Potential Energy Surfaces",                         *
!             J. M. Anglada and J. M. Bofill,                          *
!             J. Comput. Chem., 18, 992-1003 (1997)                    *
!                                                                      *
!             For more details consult:                                *
!             "New General Tools for Constrained Geometry Optimization"*
!             L. De Vico, M. Olivucci, and R. Lindh,                   *
!             JCTC, 1:1029-1037, 2005                                  *
!             doi:10.1021/ct500949                                     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!             July 2003                                                *
!***********************************************************************

use kriging_mod, only: Max_MicroIterations, nSet
use Slapaf_Info, only: CnstWght, E_Delta, GrdMax, iOptC, IRC, MF, NADC, StpLbl, StpMax
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Five, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nWndw, nInter, nIter, iOptH, jPrint, nA, nFix, nsAtom, iOpt_RS, iter_
integer(kind=iwp), intent(inout) :: nLambda
real(kind=wp), intent(inout) :: r(nLambda,nIter), drdq(nInter,nLambda,nIter), dEdq(nInter,nIter), rLambda(nLambda,nIter+1), &
                                q(nInter,nIter+1), dq(nInter,nIter), Energy(nIter), cBeta, Thr_RS
real(kind=wp), intent(out) :: T(nInter,nInter), dy(nLambda), dx(nInter-nLambda,nIter), hql(nInter,nIter), du(nInter), &
                              x(nInter-nLambda,nIter+1), dEdx(nInter-nLambda,nIter), W(nInter-nLambda,nInter-nLambda), GNrm, &
                              Err(nInter,nIter+1), EMx((nIter+1)**2), RHS(nIter+1), A(nA)
real(kind=wp), intent(in) :: Hess(nInter,nInter), fCart, Beta_Disp_Seed, d2rdq2(nInter,nInter,nLambda)
integer(kind=iwp), intent(out) :: iP(nInter)
character, intent(out) :: Step_Trunc
character(len=8), intent(inout) :: Lbl(nInter+nLambda)
logical(kind=iwp), intent(in) :: First_MicroIteration
integer(kind=iwp) :: i, iCount, iCount_Max, iInter, iIter, iLambda, iMEP, iOff_Iter, iOptC_Temp, ipTb, ipTti, iSt, j, jLambda, &
                     nTrans
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iNeg
#endif
real(kind=wp) :: Beta, Beta_Disp, Beta_Disp_Save = Zero, Det, disp(nSet), disp_long, Disp_Save = Zero, disp_short, Dummy, dxdx, &
                 dxdx_last, dydy, dydy_last, dydymax, Fact, Fact_long, Fact_short, gBeta, gg, gg_last, RG, RR_, Sf, StpMax_Save, &
                 tBeta, Thr, tmp, xBeta, yBeta
logical(kind=iwp) :: Found, IRC_setup, Recompute_disp, RVO
character(len=8) :: StpLbl_Save
character :: Step_Trunc_
real(kind=wp), allocatable :: dq_xy(:), Hessian(:,:,:), RR(:,:), RRInv(:,:), RRR(:,:), RT(:,:), RTInv(:,:), Tdy(:), Tmp1(:), &
                              Tmp2(:,:), Tr(:), Trans(:), WTr(:)
character(len=8), allocatable :: LblSave(:)
real(kind=wp), parameter :: Beta_Disp_Min = 1.0e-10_wp
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) '****************************************************'
write(u6,*) '**** Con_Opt: Input data ***************************'
write(u6,*) '****************************************************'
write(u6,*)
write(u6,*) 'iOpt_RS=',iOpt_RS
write(u6,*) 'cBeta=',cBeta
write(u6,*) 'Beta_Disp_Seed=',Beta_Disp_Seed
write(u6,*)
write(u6,*) 'First_Microiteration=',First_Microiteration
write(u6,*)
call RecPrt('Con_Opt: Energy',' ',Energy,nIter,1)
call RecPrt('Con_Opt: q',' ',q,nInter,nIter)
call RecPrt('Con_Opt: dEdq',' ',dEdq,nInter,nIter)
if (nIter > 1) call RecPrt('Con_Opt: dq',' ',dq,nInter,nIter-1)
call RecPrt('Con_Opt: Hess(in)',' ',Hess,nInter,nInter)
call RecPrt('Con_Opt: r',' ',r,nLambda,nIter)
do iIter=1,nIter
  write(u6,*) ' iIter=',iIter
  call RecPrt('Con_Opt: drdq(orig)',' ',drdq(:,:,iIter),nInter,nLambda)
end do
do iLambda=1,nLambda
  write(u6,*) 'iLambda=',iLambda
  call RecPrt('Con_Opt: d2rdq2(iLambda)',' ',d2rdq2(:,:,iLambda),nInter,nInter)
end do
write(u6,*) '****************************************************'
write(u6,*) ' End of input arrays'
write(u6,*) '****************************************************'
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
RVO = iOpt_RS /= 0

ipTb = 1
ipTti = ipTb+nLambda

yBeta = One
gBeta = One
xBeta = One
Beta = fCart*cBeta
dydy_last = Beta
gg_last = Beta
dxdx_last = Beta
Sf = sqrt(Two)
dxdx = Zero
dydy = Zero
Thr = 1.0e-6_wp
call Get_iScalar('iOff_Iter',iOff_Iter)
call mma_allocate(dq_xy,nInter,Label='dq_xy')
call mma_allocate(Hessian,nInter,nInter,nSet,Label='Hessian')
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
iSt = max(iOff_Iter+1,nIter-nWndw+1)

#ifdef _DEBUGPRINT_
write(u6,*) '****************************************************'
write(u6,*) ' Compute the lambda values'
write(u6,*) '****************************************************'
#endif
do iIter=iSt,nIter
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) '>>>>>> iIter=',iIter
  write(u6,*)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute the lambda values (Eqn. 17)
  !
  ! dEdq=drdq l ; from  h(q_0,l)=dEdq - drdq l = 0
  !
  ! drdq^T dEdq= (drdq^T drdq) l
  !
  ! l =  (drdq^T drdq)^{-1} drdq^T dEdq

  call mma_allocate(RRR,nLambda,nInter,Label='RRR')
  call mma_allocate(RRInv,nLambda,nLambda,Label='RRInv')
  call mma_allocate(RR,nLambda,nLambda,Label='RR')

  ! drdq^T drdq

  call DGEMM_('T','N',nLambda,nLambda,nInter,One,drdq(:,:,iIter),nInter,drdq(:,:,iIter),nInter,Zero,RR,nLambda)

  ! (drdq^T drdq)^{-1}

  call MInv(RR,RRInv,Det,nLambda)
  call mma_deallocate(RR)

  ! (drdq^T drdq)^{-1} drdq^T

  call DGEMM_('N','T',nLambda,nInter,nLambda,One,RRInv,nLambda,drdq(:,:,iIter),nInter,Zero,RRR,nLambda)
  call mma_deallocate(RRInv)

  ! l = (T_b^T drdq)^{-1} drdq^T dEdq

  ! Note the sign conflict due to dEdq stored as a force.

  call DGEMM_('N','N',nLambda,1,nInter,-One,RRR,nLambda,dEdq(:,iIter),nInter,Zero,rLambda(:,iIter),nLambda) ! Sign conflict
  call mma_deallocate(RRR)
# ifdef _DEBUGPRINT_
  call RecPrt('rLambda(iIter)',' ',rLambda(:,iIter),nLambda,1)
# endif
end do
#ifdef _DEBUGPRINT_
write(u6,*) '****************************************************'
write(u6,*) ' End of computing the lambda values'
write(u6,*) '****************************************************'
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) '****************************************************'
write(u6,*) ' Check that grads of constraint are not null vectors'
write(u6,*) '****************************************************'
#endif
!do iIter = iOff_Iter+1,nIter
do iIter=iSt,nIter
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) '>>>>>> iIter,nIter=',iIter,nIter
  write(u6,*)
# endif
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! If we have that the derivative of the constraint is a null
  ! vector replace it with an arbitrary vector which is orthogonal
  ! to the gradient.
  !
  ! In case of the first macro iteration of an IRC search replace
  ! it with the reaction vector.

# ifdef _DEBUGPRINT_
  write(u6,*) '****************************************************'
  write(u6,*) ' Check that grads of constraint are not null vectors'
  write(u6,*) '****************************************************'
# endif
  do iLambda=1,nLambda
    RR_ = sqrt(DDot_(nInter,drdq(:,iLambda,iIter),1,drdq(:,iLambda,iIter),1))
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iLambda=',iLambda
    call RecPrt('drdq',' ',drdq(:,iLambda,iIter),1,nInter)
#   endif

    ! Make sure that we don't mess up gradients which are zero vectors.

    if (RR_ < 1.0e-12_wp) then
      xBeta = xBeta*Half

      if (abs(IRC) /= 0) then
        iMEP = 0
        call Qpg_iScalar('nMEP',Found)
        if (Found) call Get_iScalar('nMEP',iMEP)
        IRC_Setup = (iIter == 1) .and. (iMEP == 0)
      else
        IRC_Setup = .false.
      end if

      if (iIter == nIter) then

        write(u6,*) 'Warning: constraint ',iLambda,' has a null vector, I''ll fix it!'

        if (IRC_Setup .and. (IRC == 1)) then
          write(u6,*) ' IRC forward direction.'
        else if (IRC_Setup .and. (IRC == -1)) then
          write(u6,*) ' IRC backward direction.'
        end if
      end if
      r(iLambda,iIter) = Zero

      if (IRC_SetUp) call ReacQ(MF,3*nsAtom,dEdq(:,iIter),nInter)

      ! Try to use the transverse direction if the gradient is too small

      RR_ = sqrt(DDot_(nInter,dEdq(:,iIter),1,dEdq(:,iIter),1))
      if (RR_ < 1.0e-12_wp) then
        call qpg_dArray('Transverse',Found,nTrans)
        if (Found .and. (nTrans == 3*nsAtom)) then
          call mma_Allocate(Trans,3*nsAtom,Label='Trans')
          call Get_dArray('Transverse',Trans,nTrans)
          call ReacQ(Trans,3*nsAtom,dEdq(:,iIter),nInter)
          RR_ = sqrt(DDot_(nInter,dEdq(:,iIter),1,dEdq(:,iIter),1))
          dEdq(:,iIter) = dEdq(:,iIter)/RR_
          call mma_deallocate(Trans)
        end if
      else
        dEdq(:,iIter) = dEdq(:,iIter)/RR_
      end if

      do iInter=1,nInter
        if (abs(dEdq(iInter,iIter)) <= 1.0e-6_wp) then
          drdq(iInter,iLambda,iIter) = Zero
        else
          drdq(iInter,iLambda,iIter) = sign(One,dEdq(iInter,iIter))
        end if
      end do
#     ifdef _DEBUGPRINT_
      call RecPrt('Con_Opt: drdq',' ',drdq,nInter,nLambda*nIter)
#     endif

      ! Orthogonalize against the gradient and all other constraints

      RG = DDot_(nInter,drdq(:,iLambda,iIter),1,dEdq(:,iIter),1)
      drdq(:,iLambda,iIter) = drdq(:,iLambda,iIter)-RG*dEdq(:,iIter)
      RR_ = sqrt(DDot_(nInter,drdq(:,iLambda,iIter),1,drdq(:,iLambda,iIter),1))
      do jLambda=1,nLambda
        if (jLambda /= iLambda) then
          RG = DDot_(nInter,drdq(:,iLambda,iIter),1,drdq(:,jLambda,iIter),1)
          drdq(:,iLambda,iIter) = drdq(:,iLambda,iIter)-RG*drdq(:,jLambda,iIter)
          RR_ = sqrt(DDot_(nInter,drdq(:,iLambda,iIter),1,drdq(:,iLambda,iIter),1))
          drdq(:,iLambda,iIter) = drdq(:,iLambda,iIter)/RR_
        end if
      end do

#     ifdef _DEBUGPRINT_
      call RecPrt('Con_Opt; drdq(2)',' ',drdq,nInter,nLambda*nIter)
#     endif
    end if
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) '****************************************************'
  write(u6,*) ' End Check'
  write(u6,*) '****************************************************'
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! NOTE: for historical reasons the code stores the force rather than *
  !       the gradient. Hence, be careful of the sign whenever the     *
  !       terms dEdq, hql show up.                                     *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Set up the T-matrix which separates the space into the
  ! subspace in which the constraint is accomplished and
  ! the complemental subspace.
  !
  ! [T_b,T_{ti}]
  !
  ! With the properties
  !
  ! drdq^T T_b =/= 0 (but not = 1, see note below)
  !
  ! drdq^T T_{ti} = T_b^T T_{ti} = 0
  !
  ! T_{ti}^T T_{ti} = 1

  call GS(drdq(:,:,iIter),nLambda,T,nInter,.false.,.false.)
# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: T-Matrix',' ',T,nInter,nInter)
  call RecPrt('Con_Opt: T_b',' ',T(:,ipTB),nInter,nLambda)
  call RecPrt('Con_Opt: T_ti',' ',T(:,ipTti),nInter,nInter-nLambda)
# endif
  !                                                                    *
  ! Note that the paper by Anglada and Bofill has some errors on the   *
  ! properties of T. Especially with respect to the properties of T_b. *
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute delta y  (Eqn. 11) -- in case of full relaxation --
  !
  ! r(q)+drdq^T dq=0
  !
  ! dq = [T_b, T_{ti}] (dy,dx)^T
  !
  ! r(q) = - drdq^T T_b dy - drdq^T T_{ti} dx
  !
  ! r(q) = - drdq^T T_b dy
  !
  ! dy = - (drdq^T T_b)^{-1} r(q)

  call mma_Allocate(RTInv,nLambda,nLambda,Label='RTInv')
  call mma_Allocate(RT,nLambda,nLambda,Label='RT   ')

  ! drdq^T T_b

  call DGEMM_('T','N',nLambda,nLambda,nInter,One,drdq(:,:,iIter),nInter,T(:,ipTb),nInter,Zero,RT,nLambda)
# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: RT = drdq^T T_b',' ',RT,nLambda,nLambda)
# endif

  call MInv(RT,RTInv,Det,nLambda)
  call mma_deallocate(RT)
# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: RTInv = (drdq^T T_b)^{-1}',' ',RTInv,nLambda,nLambda)
  call RecPrt('Con_Opt: r(:,iIter)',' ',r(:,iIter),nLambda,1)
# endif

  ! dy = - (drdq^T T_b)^{-1} r(q)

  call DGEMM_('N','N',nLambda,1,nLambda,-One,RTInv,nLambda,r(:,iIter),nLambda,Zero,dy,nLambda)
  call mma_deallocate(RTInv)
# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: dy(full) = - (drdq^T T_b)^{-1} r(q)',' ',dy,nLambda,1)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Add contributions from constraints according to the last
  ! iterations. See Eqn. 7,
  !
  ! The Hessian of the Lagrangian expressed in the full space.
  !
  ! W^0 = H^0 - Sum(i) l_{0,i} d2rdq2(q_0)

  if (.not. RVO) then
    Hessian(:,:,1) = Hess(:,:)
    if (iIter == nIter) then
      do iLambda=1,nLambda
        Hessian(:,:,1) = Hessian(:,:,1)-rLambda(iLambda,iIter)*d2rdq2(:,:,iLambda)
      end do
    end if
  else
    call Hessian_Kriging_Layer(q(:,iIter),Hessian,nInter)
    if ((nSet > 2) .and. NADC) Hessian(:,:,1) = Half*(Hessian(:,:,1)+Hessian(:,:,2))
    do iLambda=1,nLambda
      Hessian(:,:,1) = Hessian(:,:,1)-rLambda(iLambda,iIter)*d2rdq2(:,:,iLambda)
    end do
  end if
# ifdef _DEBUGPRINT_
  call RecPrt('W^0 = H^0 - Sum(i) l_{0,i} d2rdq2(q_0)',' ',Hessian(:,:,1),nInter,nInter)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Transform various vectors and matrices to the nInter-nLambda
  ! subspace. See Eqn. 10 and text following.
  !
  ! Compute x, the coordinates of the reduced space.

  if (nInter-nLambda > 0) then

    ! q = T_b y + T_{ti} x = T [y,x]^T
    !
    ! x = T_{ti}^T q, since T_{ti}^T T_b = 0

    call DGEMM_('T','N',nInter-nLambda,1,nInter,One,T(:,ipTti),nInter,q(:,iIter),nInter,Zero,x(:,iIter),nInter-nLambda)

    ! Compute dx
    !
    ! dx = T_{ti}^T dq

    call DGEMM_('T','N',nInter-nLambda,1,nInter,One,T(:,ipTti),nInter,dq(:,iIter),nInter,Zero,dx(:,iIter),nInter-nLambda)

    ! Compute step restriction based on information from the
    ! minimization in the x subspace. This restriction is based
    ! on the step length in the x subspace.

    if (iIter /= nIter) then
#     ifdef _DEBUGPRINT_
      write(u6,*) 'dx = T_{ti}^T dq'
      call RecPrt('dx(:,iIter)',' ',dq(:,iIter),nInter,1)
#     endif

      ! Compute dq step in the x subspace
      !
      ! du = [0,dx]^T
      !
      ! dq = T [0,dx]^T

      du(:) = Zero
      du(nLambda+1:) = dx(1:nInter-nLambda,iIter)
      call Backtrans_T(du,dq_xy)
      dxdx = sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))

#     ifdef _DEBUGPRINT_
      write(u6,*) 'dxdx=',dxdx
      write(u6,*) 'dxdx_last=',dxdx_last
#     endif

      if ((dxdx < 0.75_wp*dxdx_last) .and. (dxdx < Beta-Thr)) then
        ! Increase trust radius
        xBeta = min(One,xBeta*Sf)
      else if ((dxdx > 1.25_wp*dxdx_last) .or. (dxdx >= Beta+Thr)) then
        ! Reduce trust radius
        xBeta = max(One/Five,xBeta/Sf)
      end if
      dxdx_last = dxdx
#     ifdef _DEBUGPRINT_
      write(u6,*) 'xBeta=',xBeta
#     endif
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute the reduced gradient, observe again that we store
    ! the forces rather than the gradients.
    !
    ! See text after Eqn. 13.
    !
    ! dEdx = T_{ti}^T (dEdq + W_{ex} T_b dy)

    call mma_allocate(Tmp1,nInter,Label='Tmp1')
    call mma_allocate(Tmp2,nInter,nLambda,Label='Tmp2')

#   ifdef _DEBUGPRINT_
    call RecPrt('Con_Opt: dEdq',' ',dEdq(:,iIter),1,nInter)
    call RecPrt('Con_Opt: W',' ',Hessian(:,:,1),nInter,nInter)
    call RecPrt('Con_Opt: T',' ',T,nInter,nInter)
    call RecPrt('Con_Opt: dy',' ',dy,1,nLambda)
#   endif

    ! W_{ex} T_b

    call DGEMM_('N','N',nInter,nLambda,nInter,One,Hessian,nInter,T(:,ipTb),nInter,Zero,Tmp2,nInter)
#   ifdef _DEBUGPRINT_
    call RecPrt('W_{ex} T_b',' ',Tmp2,nInter,nLambda)
#   endif

    ! dEdq + W_{ex} T_b dy
    !
    ! Since we are computing the force we have the conflicting sign below.

    Tmp1(:) = dEdq(:,iIter)
    call DGEMM_('N','N',nInter,1,nLambda,-One,Tmp2,nInter,dy,nLambda,One,Tmp1,nInter) ! Sign conflict
    call mma_deallocate(Tmp2)
#   ifdef _DEBUGPRINT_
    call RecPrt('dEdq + W_{ex} T_b dy',' ',Tmp1,1,nInter)
#   endif

    ! dEdx = T_{ti}^T (dEdq + W_{ex} T_b dy)

    call DGEMM_('T','N',nInter-nLambda,1,nInter,One,T(:,ipTti),nInter,Tmp1,nInter,Zero,dEdx(:,iIter),nInter-nLambda)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iIter=',iIter
    call RecPrt('dEdx(:,iIter)',' ',dEdx(:,iIter),1,nInter-nLambda)
#   endif
    call mma_deallocate(Tmp1)

    ! Only now do we have the constrained gradient so we can
    ! reduce the step size (see update_kriging)

    if (RVO .and. First_Microiteration .and. (iIter == nIter)) then
      GNrm = sqrt(DDot_(nInter-nLambda,dEdx(:,nIter),1,dEdx(:,nIter),1))
      cBeta = min(1.0e3_wp*GNrm,cBeta)
      Beta = fCart*cBeta
    end if

    ! Compute step restriction based on information from the
    ! minimization in the x subspace. This restriction is based
    ! on the gradient in the x subspace.

    gg = sqrt(DDot_(nInter-nLambda,dEdx(:,iIter),1,dEdx(:,iIter),1))
    if ((gg < 0.75_wp*gg_last) .and. (gg < Beta-Thr)) then
      ! Increase trust radius
      gBeta = min(One,gBeta*Sf)
      !gBeta = gBeta*Sf
    else if ((gg > 1.25_wp*gg_last) .or. (gg >= Beta+Thr)) then
      ! Reduce trust radius
      gBeta = max(One/Five,gBeta/Sf)
    end if
    gg_last = gg
#   ifdef _DEBUGPRINT_
    write(u6,*) 'gg=',gg
    write(u6,*) 'gBeta=',gBeta
#   endif

  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Add contributions to the Lagrangian energy change

  if (iIter == nIter) then

    ! Term due to constraint

    call mma_allocate(Tmp1,nLambda,Label='Tmp1')
    call mma_allocate(Tdy,nInter,Label='Tdy')

    ! T_b dy

    call DGEMM_('N','N',nInter,1,nLambda,One,T(:,ipTb),nInter,dy,nLambda,Zero,Tdy,nInter)

    ! drdq^T dq = drdq^T T_b dy

    call DGEMM_('T','N',nLambda,1,nInter,One,drdq(:,:,iIter),nInter,Tdy,nInter,Zero,Tmp1,nLambda)
    call mma_deallocate(Tdy)

    ! r + drdq^T dq

    Tmp1(1:nLambda) = Tmp1(1:nLambda)+r(1:nLambda,iIter)

    ! l^T (r + drdq^T dq)

    E_Delta = E_Delta-DDot_(nLambda,rLambda(:,iIter),1,Tmp1,1)
    call mma_deallocate(Tmp1)

    ! Term due to coupling

    call mma_allocate(Tr,nInter,Label='Tr')

    ! T_b dy

    call DGEMM_('N','N',nInter,1,nLambda,One,T(:,ipTb),nInter,dy,nLambda,Zero,Tr,nInter)

    ! dEdq^T T_b dy
    !
    ! Note the sign conflict.

    E_Delta = E_Delta-DDot_(nInter,Tr,1,dEdq,1)

    call mma_allocate(WTr,nInter,Label='WTr')

    ! W T_b dy

    call DGEMM_('N','N',nInter,1,nInter,One,Hessian,nInter,Tr,nInter,Zero,WTr,nInter)

    ! dy^T T_b^T W T_b dy

    E_Delta = E_Delta+Half*DDot_(nInter,Tr,1,WTr,1)
    call mma_deallocate(WTr)
    call mma_deallocate(Tr)
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Restrict actual step, dy
  !
  ! Step in direction which fulfills the restriction is given priority.
  ! This step is only reduced if it is larger than half the overall step
  ! restriction.

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Restricting the dy step.'
  write(u6,*) 'Start: dy(:)=',dy(:)
  write(u6,*)
# endif
  if (.not. RVO) then

    ! Compute dq step in the y subspace
    !
    ! dq_y = T [dy, 0]^T

    du(:) = Zero
    du(1:nLambda) = dy(:)
    call Backtrans_T(du,dq_xy)
    dydy = sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))

    ! Reduce y step size if larger than some maximum size
    ! (the x step or half the total max step times a factor)

    dydymax = CnstWght*max(dxdx,Half*Beta)
    if (dydy > dydymax) then
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Reduce dydy!',dydy,' -> ',dydymax
      write(u6,*) 'Scaling Factor=',(dydymax/dydy)
#     endif
      dy(:) = (dydymax/dydy)*dy(:)
      dydy = dydymax
      Step_Trunc = '*'
    else
#     ifdef _DEBUGPRINT_
      write(u6,*) 'No reduce dydy!',dydy,' < ',dydymax
#     endif
      Step_Trunc = ' '
    end if

    ! The step reduction in the space which we minimize is such
    ! that while the fulfillment of the constraint is not
    ! improving from step to step we reduce the step length in
    ! the subspace in which we do the minimization.

    if ((dydy < 0.75_wp*dydy_last) .or. (dydy < 1.0e-2_wp)) then
      ! Recent step in the space for the restriction is smaller
      ! than the previous, or the recent step is smaller than
      ! some threshold. Then increase step length for x, however
      ! not more than the overall step restriction.
      yBeta = min(Two,yBeta*Sf)
    else if ((dydy > 1.25_wp*dydy_last) .and. (dydy >= 1.0e-5_wp)) then
      ! Otherwise decrease step direction.
      yBeta = max(One/Ten,yBeta/Sf)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  else
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Here in the case of kriging and restricted-variance
    ! optimization.
    !
    ! We only need this for the last point.

    Fact = One
#   ifdef _DEBUGPRINT_
    if (iIter /= nIter) write(u6,*) 'Skip'
#   endif
    if (iIter == nIter) then

      ! Pick up the largest element of the gradient of the constraints.
      tmp = Zero
      do i=1,nLambda
        do j=1,nInter
          tmp = max(tmp,abs(drdq(j,i,iIter)))
        end do
      end do
      tmp = min(tmp,0.3_wp) ! Some constraints can have huge gradients. So be a bit careful.
      ! Allowed dispersion in the y subspace
      Beta_Disp = max(Beta_Disp_Min,tmp*CnstWght/(CnstWght+One)*Beta_Disp_Seed)
      if (.not. First_MicroIteration) Beta_Disp = min(Disp_Save+Beta_Disp,Beta_Disp_Save)

#     ifdef _DEBUGPRINT_
      write(u6,*) 'Step_trunc=',Step_trunc
      write(u6,*) 'CnstWght=',CnstWght
      write(u6,*) 'Beta=',Beta
      write(u6,*) 'tmp=',tmp
      write(u6,*) 'Beta_Disp=',Beta_Disp
#     endif

      dydy = DDot_(nLambda,dy,1,dy,1)
      if (Disp_Save/Beta_Disp > 0.99_wp) then
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Skip'
#       endif
        dy(:) = Zero
      else if (dydy >= 1.0e-12_wp) then
        ! Restrict dy step during micro iterations
        Fact = max(sqrt(dydy)/(CnstWght/(CnstWght+One)*Beta),One)

        iCount = 1
        iCount_Max = 100
        if (Step_Trunc == 'N') Step_Trunc = ' '
        do
          du(1:nLambda) = (One/Fact)*dy(:)
          du(nLambda+1:nInter) = Zero
          call Backtrans_T(du,dq_xy)
          q(:,iIter+1) = q(:,iIter)+dq_xy(:)

          call Dispersion_Kriging_Layer(q(:,iIter+1),disp,nInter)

          if (iCount == 1) then
            Fact_long = Fact
            disp_long = disp(1)
            Fact_short = Zero
            disp_short = disp_long+One
          end if
#         ifdef _DEBUGPRINT_
          write(u6,*) 'disp,Fact,iCount=',disp(1),Fact,iCount
#         endif
          if ((disp(1) <= Beta_Disp) .and. (iCount == 1)) exit
          if (abs(Beta_Disp-disp(1)) < Thr_RS) exit
          iCount = iCount+1
          if (iCount > iCount_Max) then
            write(u6,*) 'iCount > iCount_Max'
            call Abend()
          end if
          call Find_RFO_Root(Fact_long,disp_long,Fact_short,disp_short,Fact,disp(1),Beta_Disp)
          Step_Trunc = '*'
        end do
      end if
    end if
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Final scaling Factor:',(One/Fact)
#   endif
    dy(:) = (One/Fact)*dy(:)
#   ifdef _FindTS_

    ! If a FindTS optimization hold back a bit. The purpose of
    ! the constrained optimization phase is actually not to
    ! find the constrained structure.

    if (btest(iOptC,12)) then
      dydy = sqrt(DDot_(nLambda,dy,1,dy,1))
      Thrdy = 0.075_wp
      if (dydy > Thrdy) then
        dy(:) = (Thrdy/dydy)*dy(:)
        Step_Trunc = '*'
      end if
    end if
#   endif
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Step_trunc=',Step_trunc
    write(u6,*) 'Final: dy(:)=',dy(:)
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !      Twist for MEP optimizations.

  if ((iIter == iOff_iter+1) .and. (dydy < 1.0e-4_wp) .and. (iIter /= 1)) xBeta = xBeta*Half
  dydy_last = dydy

# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: dy(actual)',' ',dy,nLambda,1)
# endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
end do ! Do iIter = iSt, nIter
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: dEdx',' ',dEdx,nInter-nLambda,nIter)
call RecPrt('Con_Opt: Lambda',' ',rLambda,nLambda,nIter)
call RecPrt('Con_Opt: x',' ',x,nInter-nLambda,nIter)
call RecPrt('Con_Opt: dx',' ',dx,nInter-nLambda,nIter-1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the h vectors.   (Eqn. 16)
!
! h(q,l) = dEdq - drdq l_0^T
!
! Note the sign conflict due to storage of the force rather than the
! gradient.
#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: dEdq',' ',dEdq,nInter,nIter)
write(u6,*) 'iOff_iter=',iOff_iter
#endif
hql(:,:) = dEdq(:,:)
do iIter=iOff_iter+1,nIter
  do iLambda=1,nLambda
    hql(:,iIter) = hql(:,iIter)+rLambda(iLambda,nIter)*drdq(:,iLambda,iIter) ! Sign conflict
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: h(q,l) = dEdq - drdq l_0^T',' ',hql,nInter,nIter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the value of the Lagrangian
!
! L = E + l r(q)

do iIter=iOff_iter+1,nIter
  do iLambda=1,nLambda
    Energy(iIter) = Energy(iIter)+rLambda(iLambda,iIter)*r(iLambda,iIter)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: Lagrangian, L = E + l r(q)',' ',Energy,nIter,1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Update the Hessian in the m subspace
! There should be no negative eigenvalue, if so change the sign.

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*)
write(u6,*) ' *** Updating the reduced Hessian ***'
write(u6,*)
#endif

if (btest(iOptC,12)) then

  ! If FINDTS option force minimization option during the Hessian update.

  iOptC_Temp = 128
else
  iOptC_Temp = iOptC
end if

#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: Hessian before the update',' ',Hessian(:,:,1),nInter,nInter)
write(u6,*) 'iOptH=',iOptH
#endif
if (Step_Trunc == 'N') Step_Trunc = ' '
Dummy = Zero
call Update_H(nWndw,Hessian,nInter,nIter,iOptC_Temp,dq,hql,iOptH,jPrint,Dummy,nsAtom,.false.,.false.)

#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: W(updated)',' ',Hessian(:,:,1),nInter,nInter)
write(u6,*) 'Step_Trunc=',Step_trunc
call DiagMtrx(Hessian(:,:,1),nInter,iNeg)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the reduced Hessian
!
! See text after Eqn. 13.
!
! T_{ti}^T W T_{ti}

if (nInter-nLambda > 0) then

  call mma_allocate(Tmp2,nInter-nLambda,nInter,Label='Tmp2')
  call DGEMM_('T','N',nInter-nLambda,nInter,nInter,One,T(:,ipTti),nInter,Hessian,nInter,Zero,Tmp2,nInter-nLambda)
  call DGEMM_('N','N',nInter-nLambda,nInter-nLambda,nInter,One,Tmp2,nInter-nLambda,T(:,ipTti),nInter,Zero,W,nInter-nLambda)
  call mma_deallocate(Tmp2)
# ifdef _DEBUGPRINT_
  call RecPrt('Con_Opt: the reduced Hessian, T_{ti}^T W T_{ti}',' ',W,nInter-nLambda,nInter-nLambda)
  call DiagMtrx(W,nInter-nLambda,iNeg)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Update dx
  !
  ! Set threshold depending on if restriction is w.r.t. step-size
  ! or variance.

  if (.not. RVO) then
    Beta_Disp = One ! Dummy assign
  else

    ! Note that we use the dEdx data for the last point on the real PES.

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Beta_Disp=',Beta_Disp
#   endif
    tmp = Zero
    do i=1,nInter-nLambda
      tmp = max(tmp,abs(dEdx(i,iter_)))
    end do
    ! Add the allowed dispersion in the x subspace.
    if (First_MicroIteration) Beta_Disp_Save = max(Beta_Disp+tmp*Beta_Disp_Seed,Beta_Disp_Min)
    Beta_Disp = Beta_Disp_Save
#   ifdef _DEBUGPRINT_
    write(u6,*) 'tmp,Beta_Disp=',tmp,Beta_Disp
#   endif
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Compute updated geometry in internal coordinates

  fact = One
  tBeta = max(Beta*yBeta*min(xBeta,gBeta),Beta/Ten)
  Thr_RS = 1.0e-7_wp
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Compute dx'
  write(u6,*)
  write(u6,*) 'Step_Trunc(0)=',Step_Trunc
  write(u6,*) 'Beta=',Beta
  write(u6,*) 'yBeta=',yBeta
  write(u6,*) 'xBeta=',xBeta
  write(u6,*) 'gBeta=',gBeta
  write(u6,*) 'tBeta=',tBeta
# endif
  do
    Step_Trunc_ = Step_Trunc
    call Newq(x,nInter-nLambda,nIter,dx,W,dEdx,Err,EMx,RHS,A,nA,tBeta,nFix,iP,Energy,Step_Trunc_,Thr_RS)
    if (Step_Trunc == 'N') Step_Trunc = ' '

    if (.not. RVO) then
      if (Step_Trunc_ == 'N') Step_Trunc_ = ' '
      Step_Trunc = Step_Trunc_
      exit
    end if

    if (Step_Trunc//Step_Trunc_ == ' *') Step_Trunc = '.'

    du(1:nLambda) = dy(:)
    du(nLambda+1:nInter) = dx(:,nIter)
    call Backtrans_T(du,dq_xy)
    q(:,nIter+1) = q(:,nIter)+dq_xy(:)

    call Dispersion_Kriging_Layer(q(:,nIter+1),disp,nInter)
    Disp_Save = disp(1)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'disp=',disp(1)
#   endif
    fact = Half*fact
    tBeta = Half*tBeta
    if (One-disp(1)/Beta_Disp > 1.0e-3_wp) exit
    if ((fact < 1.0e-5_wp) .or. (disp(1) < Beta_Disp)) exit
    Step_Trunc = '*'
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) 'Step_Trunc(n)=',Step_Trunc
# endif
  GNrm = sqrt(DDot_(nInter-nLambda,dEdx(:,nIter),1,dEdx(:,nIter),1))

else
  ! Negative norm should serve as sign that there is no gradient
  GNrm = -One
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'Step_Trunc=',Step_trunc
call RecPrt('Con_Opt: dx(actual)',' ',dx,nInter-nLambda,nIter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Back transform dy and dx to dq.
!
! See Eqn. 10.

#ifdef _DEBUGPRINT_

! dy only, constraint

du(:) = Zero
du(1:nLambda) = dy(:)
call Backtrans_T(du,dq_xy)
dydy = sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))
call RecPrt('dq(dy)',' ',dq_xy,nInter,1)
write(u6,*) '<R(q_0)|dy>=',DDot_(nInter,dRdq(:,1,nIter),1,dq_xy,1)

! dx only, constrained minimization

du(:) = Zero
du(nLambda+1:nInter) = dx(:,nIter)
call Backtrans_T(du,dq_xy)
dxdx = sqrt(DDot_(nInter,dq_xy,1,dq_xy,1))
call RecPrt('dq(dx)',' ',dq_xy,nInter,1)
write(u6,*) '<R(q_0)|dx>=',DDot_(nInter,dRdq(:,1,nIter),1,dq_xy,1)
#endif

! dy+dx, full step
!
! dq = T [dy, dx]

du(:) = Zero
du(1:nLambda) = dy(1:nLambda)
du(nLambda+1:nInter) = dx(:,nIter)
#ifdef _DEBUGPRINT_
call RecPrt('du',' ',du,1,nInter)
#endif

! For kriging, in the last 10 micro iterations, give up trying to
! optimize and just focus on fulfilling the constraints.

Recompute_disp = .false.
if ((Max_MicroIterations >= 50) .and. (nIter-iter_+1 > Max_Microiterations-10)) then
  dydy = sqrt(DDot_(nLambda,dy,1,dy,1))
  if (dydy > 1.0e-12_wp) du(nLambda+1:nInter) = Zero
  Recompute_disp = .true.
end if
call Backtrans_T(du,dq(:,nIter))
#ifdef _DEBUGPRINT_
call RecPrt('dq(dx+dy)',' ',dq(:,nIter),1,nInter)
#endif

! Compute q for the next iteration

q(:,nIter+1) = q(:,nIter)+dq(:,nIter)
if (Recompute_disp) then
  call Dispersion_Kriging_Layer(q(:,nIter+1),Disp,nInter)
  Disp_Save = Disp(1)
end if

#ifdef _DEBUGPRINT_
call RecPrt('Con_Opt: q',' ',q,nInter,nIter+1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! StpMax from q

call MxLbls(nInter,hql(:,nIter),dq(:,nIter),Lbl)

! GrdMax for dEdx

if (nInter-nLambda > 0) then

  StpMax_Save = StpMax
  StpLbl_Save = StpLbl
  call mma_allocate(LblSave,size(Lbl),Label='LblSave')
  LblSave(:) = Lbl
  GrdMax = Zero
  do i=1,nInter-nLambda
    write(Lbl(i),'(A,I3.3)') 'dEdx',i
  end do
  call MxLbls(nInter-nLambda,dEdx(:,nIter),dx(:,nIter),LblSave)
  call mma_deallocate(LblSave)
  StpMax = StpMax_Save
  StpLbl = StpLbl_Save

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(dq_xy)
call mma_deallocate(Hessian)

return

contains

subroutine Backtrans_T(du,dq)
  implicit none
  real(kind=wp), intent(in) :: du(nInter)
  real(kind=wp), intent(out) :: dq(nInter)
  call DGEMM_('N','N',nInter,1,nInter,One,T,nInter,du,nInter,Zero,dq,nInter)
end subroutine Backtrans_T

end subroutine Con_Opt
