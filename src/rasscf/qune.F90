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

subroutine QUNE(NCALLS,ENOW,BK,XSX,VL,VM,XQN,XOLD,V1,V2,NDIM,LUQUNE,TMIN,QNSTEP,QNUPDT,KSDFT)
! RASSCF PROGRAM VERSION IBM-3090: SX SECTION
!
! PURPOSE: THIS SUBROUTINE MAKES A QUASI NEWTON UPDATE
!          OF THE ROTATION MATRIX X, OR A LINE SEARCH.

use RASDim, only: MxIter
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(inout) :: NCALLS
real(kind=wp), intent(in) :: ENOW
integer(kind=iwp), intent(in) :: NDIM, LUQUNE
real(kind=wp), intent(_IN_) :: BK(NDIM)
real(kind=wp), intent(inout) :: XSX(NDIM)
real(kind=wp), intent(out) :: VL(NDIM), VM(NDIM), XQN(NDIM), XOLD(NDIM), V1(NDIM), V2(NDIM), TMIn
character(len=2), intent(out) :: QNSTEP
character(len=3), intent(out) :: QNUPDT
character(len=80), intent(in) :: KSDFT
integer(kind=iwp) :: IAD, IVEC, NLM
real(kind=wp) :: C0, C1, C2, C3, E0, E1, ELM, EMIN, EPRED_LS, EPRED_QN, EPRED_SX, FP, P, P1, P2, Q, SCLFCT, T0, T1, TLM, X, &
                 XSXNRM, Y
integer(kind=iwp), save :: NLS, NVEC
real(kind=wp), save :: ALPHA(mxiter+2), BETA(mxiter+2), ELAST, FPLAST
real(kind=wp), external :: DDot_, DNrm2_

NCALLS = NCALLS+1

if (NCALLS == 1) then
  NVEC = 0
  ! NLS: Nr of consecutive line searches.
  NLS = 0
  ALPHA(:) = Zero
  BETA(:) = Zero
  IAD = 0
  call DDAFILE(LUQUNE,1,BK,NDIM,IAD)
  call DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
  call DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
  ELAST = ENOW
  FPLAST = Two*DDOT_(NDIM,BK,1,XSX,1)
  TMIN = Zero
  QNSTEP = 'SX'
  QNUPDT = ' NO'
  return
end if
! -- 1. CALC L AND M ARRAYS
! Read from the beginning of LUQUNE. There is the BLB gradient, the proposed optimal
! (QN) step at the previous iteration, and the actual step taken in that iteration.
! VM should be computed as the difference in QN steps, where the new one is computed
! without the new update vectors (which are still unknown).
IAD = 0
call DDAFILE(LUQUNE,2,V1,NDIM,IAD)
VL(:) = BK(:)-V1(:)
! VL is the difference between the old and the new BLB gradient array:
! Read old QN step into VM:
call DDAFILE(LUQUNE,2,VM,NDIM,IAD)
call DDAFILE(LUQUNE,2,XOLD,NDIM,IAD)
! -- FP WILL BE USED IN LINE SEARCH ANALYSIS LATER.
FP = Two*DDOT_(NDIM,BK,1,XOLD,1)
! -- QN VECTOR UPDATED WITHOUT NEWEST UPDATE VECTORS: XQN=-HINV*BK
XQN(:) = XSX(:)
do IVEC=1,NVEC
  call DDAFILE(LUQUNE,2,V1,NDIM,IAD)
  call DDAFILE(LUQUNE,2,V2,NDIM,IAD)
  X = -DDOT_(NDIM,V1,1,BK,1)
  Y = -DDOT_(NDIM,V2,1,BK,1)
  P1 = ALPHA(IVEC)*X+BETA(IVEC)*Y
  P2 = BETA(IVEC)*X
  XQN(:) = XQN(:)+P1*V1(:)+P2*V2(:)
end do
! Subtract. VM will now contain the difference in QN steps, where the new one
! is computed without the update (which we still do not know).
VM(:) = VM(:)-XQN(:)
! -- NOTE: M ARRAY = HK(INV)*LK
! -- 2. DECIDE ON WHICH ROTATION TO USE.
! -- A LINE SEARCH ANALYSIS.
C0 = ELAST
C1 = FPLAST
C2 = Three*(ENOW-ELAST)-Two*FPLAST-FP
C3 = -Two*(ENOW-ELAST)+FPLAST+FP
!write(u6,*) repeat('*',60)
!write(u6,*) ' SEARCH FOR MINIMUM:'
!write(u6,*) ' POLYNOMIAL COEFFICIENTS:'
!write(u6,'(A,4F16.8)') ' ENOW,ELAST,FP,FPLAST',ENOW,ELAST,FP,FPLAST
!write(u6,'(A,4F16.8)') ' C0,C1,C2,C3',C0,C1,C2,C3
P = Three*C1*C3
Q = C2**2
! -- NLM=NR OF LOCAL MINIMA
NLM = 0
TLM = Zero
if (abs(P) > 0.001_wp*Q) then
  if (Q > P) then
    !write(u6,*) ' THIS IS A 3RD DEGREE POLY WITH 2 STAT. POINTS.'
    NLM = 1
    TLM = (sqrt(Q-P)-C2)/(Three*C3)
    !write(u6,*) ' THERE IS A LOCAL MINIMUM AT'
    !write(u6,*) ' TLM=',TLM
  else
    !write(u6,*) ' THIS IS A MONOTONOUS 3RD DEGREE POLY.'
  end if
else if (abs(C2) > 0.001_wp*C1) then
  if (C2 > Zero) then
    !write(u6,*) ' THIS IS A 2ND DEGREE POLY WITH A MINIMUM.'
    NLM = 1
    TLM = -C1/(Two*C2)
    !write(u6,*) ' THERE IS A LOCAL MINIMUM AT'
    !write(u6,*) ' TLM=',TLM
  else
    !write(u6,*) ' THIS IS A 2ND DEGREE POLY WITH A MAXIMUM.'
  end if
end if
!if (NLM == 0) write(u6,*) ' NO LOCAL MINIMUM.'
T0 = -Half
T1 = +2.5_wp
if (NLM == 1) then
  ELM = C0+TLM*(C1+TLM*(C2+TLM*C3))
  !write(u6,*) ' LOCAL MINIMUM AT TLM=',TLM
  !write(u6,*) '     ENERGY AT TLM IS=',ELM
  ! --  DISREGARD LOCAL MINIMUM IF TOO FAR AWAY:
  if (TLM > T1) NLM = 0
  if (TLM < T0) NLM = 0
  !if (NLM == 0) write(u6,*) ' TOO FAR AWAY -- REJECTED.'
end if
! -- RULES: 1. IF THERE IS A LOCAL MINIMUM IN THE TRUST REGION
!              -0.4<TLM<1.4, THIS IS USED.
! --        2. ELSE, SELECT THE LOWEST MINIMUM IN -0.5..+2.5.
if ((NLM == 1) .and. (abs(TLM-Half) < 0.9_wp)) then
  TMIN = TLM
  EMIN = ELM
  !write(u6,*) ' THE LOCAL MINIMUM IS USED.'
else
  E0 = C0+T0*(C1+T0*(C2+T0*C3))
  E1 = C0+T1*(C1+T1*(C2+T1*C3))
  !write(u6,*) ' EXTRAP. ENERGIES IN T=-0.5 AND 2.5 ARE'
  !write(u6,*) E0,E1
  if (E0 < E1) then
    TMIN = T0
    EMIN = E0
  else
    TMIN = T1
    EMIN = E1
  end if
  if (NLM == 1) then
    ELM = C0+TLM*(C1+TLM*(C2+TLM*C3))
    if (ELM < EMIN) then
      TMIN = TLM
      EMIN = ELM
    else
      NLM = 0
    end if
  end if
end if
!write(u6,*) '  SELECTED TMIN:',TMIN
!write(u6,*) ' PREDICTED EMIN:',EMIN
!write(u6,*) repeat('*',60)

! Predicted energy lowering, depending on step taken.
! EPRED_LS should be quite exact, if step length is not too long.
! But EPRED_SX and EPRED_QN are much more uncertain.
EPRED_LS = EMIN-ENOW
EPRED_SX = Half*DDOT_(NDIM,BK,1,XSX,1)
EPRED_QN = Half*DDOT_(NDIM,BK,1,XQN,1)
!write(u6,*) ' Predicted QN energy change:',EPRED_QN
!write(u6,*) ' Predicted SX energy change:',EPRED_SX
!write(u6,*) ' Predicted LS energy change:',EPRED_LS

! Here follows decision whether to update inverse Hessian, or not:
if (NLM == 1) then
  X = C3*(One-TMIN)/sqrt(Q-P)
  if ((abs(X) < 0.2_wp) .and. (TMIN > Half)) then
    !PAM THEN THE ERROR ARISING FROM NONLINEARITY OF GRADIENT ALONG THE
    !PAM SEARCH DIRECTION IS SMALLER THAN ABOUT 25 PERCENT, SO IT IS
    !PAM MEANINGFUL INFORMATION ABOUT GRADIENT HESSIAN. ALSO, THE
    !PAM INFORMATION GIVEN BY THE LAST GRADIENT IS BETTER THAN THE ONE
    !PAM FROM THE NEXT-TO-LAST. THEREFORE, UPDATE INVERSE HESSIAN.
    ! -- WE ARE PROBABLY IN A QUADRATIC MINIMUM REGION, SO UPDATE
    !write(u6,*) ' QUAD MINIMUM REGION -- UPDATE INVERSE HESSIAN.'
    QNUPDT = 'YES'
  else
    ! -- WE ARE OUTSIDE QUADRATIC MINIMUM.
    !write(u6,*) ' NOT IN QUADRATIC MIN. INV HESSIAN NOT UPDATED.'
    QNUPDT = ' NO'
  end if
else
  ! -- DEFINITELY NON-QUADRATIC REGION.
  !write(u6,*) ' VERY NON-QUADRATIC REGION.'
  !write(u6,*) ' SHOULD WE SCRAP ALL THE UPDATE VECTORS??'
  !write(u6,*) ' NO, BUT DO NOT UPDATE.'
  QNUPDT = ' NO'
  !if (KSDFT(1:3) /= 'SCF') NVEC = 0
end if

! I strongly mistrust the above analysis. For the moment, decide to update
! always. The analysis should instead be based on if whether a large rotation
! remains (do not update) or has been done (scrap the vectors).
!if (QNUPDT /= 'YES') write(u6,*) ' DECIDE TO UPDATE ANYWAY! (PAM Dec 2006)'
QNUPDT = 'YES'

if (QNUPDT == 'YES') then
  ! Determine new pair of update vectors, and also use it to correct XQN:
  ! Recall: VL is the difference between the old and the new BLB gradient array:
  ! VM is the difference between the old, and the new provisional, XQN arrays.
  NVEC = NVEC+1
  ! METHOD: BFGS
  V1(:) = XOLD(:)
  V2(:) = VM(:)
  X = DDOT_(NDIM,V1,1,VL,1)
  Y = DDOT_(NDIM,V2,1,VL,1)
  if (X == Zero) then
    ALPHA(NVEC) = 1.0e99_wp
    BETA(NVEC) = -1.0e49_wp
  else
    ALPHA(NVEC) = (One+Y/X)/X
    BETA(NVEC) = -One/X
  end if
  ! METHOD: SYMMETRIZED POWELL 1-RANK:
  !V1(:) = VL(:)
  !V2(:) = XOLD(:)-VM(:)
  !X = DDOT_(NDIM,V1,1,V1,1)
  !Y = DDOT_(NDIM,V2,1,V1,1)
  !ALPHA(NVEC) = -Y/(X**2)
  !BETA(NVEC) = One/X
  ! --  ADD THE NEWEST UPDATE CONTRIBUTION INTO XQN:
  X = -DDOT_(NDIM,V1,1,BK,1)
  Y = -DDOT_(NDIM,V2,1,BK,1)
  P1 = ALPHA(NVEC)*X+BETA(NVEC)*Y
  P2 = BETA(NVEC)*X
  XQN(:) = XQN(:)+P1*V1(:)+P2*V2(:)
  call DDAFILE(LUQUNE,1,V1,NDIM,IAD)
  call DDAFILE(LUQUNE,1,V2,NDIM,IAD)
  ! We have added to XQN the two-rank update,
  !  (alpha*|V1><V1| + beta*|v1><V2| + beta |V2><V1|)|BK>
end if

EPRED_QN = Half*DDOT_(NDIM,BK,1,XQN,1)
!write(u6,*) ' Predicted QN energy change (revised):', EPRED_QN

! -- IF LINE SEARCH IS NEARLY CONVERGED, USE THE QN OR SX STEP
! -- TO GET NEW DIRECTION, ELSE CONTINUE LINE SEARCH:
X = TMIN-One
if ((abs(X) < 0.4_wp) .or. (abs(EPRED_LS) < 1.0e-8_wp) .or. (KSDFT(1:3) /= 'SCF')) then
  !write(u6,*) ' THE LINE SEARCH MINIMUM IS PREDICTED TO BE'
  !write(u6,*) ' RATHER CLOSE TO THE CURRENT POINT.'
  !write(u6,*) ' THEREFORE, DO NOT CONTINUE BY LS.'
  if (NVEC == 0) QNSTEP = 'SX'
  if (NVEC > 0) QNSTEP = 'QN'
else if (NLS >= 2) then
  !write(u6,*) ' WE HAVE ALREADY USED LINE SEARCHES CONSECUTIVELY'
  !write(u6,*) '    FOR SEVERAL ITERATIONS. WE MUST TRY TO BREAK'
  !write(u6,*) '    OUT FROM THIS LINE SEARCH: DISABLE IT.'
  if (NVEC == 0) QNSTEP = 'SX'
  if (NVEC > 0) QNSTEP = 'QN'
else
  ! Else we should use a line search. But maybe check predicted effect:
  ! If predicted LS-energy minimum is not an improvement
  ! then do not use line search.
  !write(u6,*) ' THE LINE-SEARCH ANALYSIS PREDICTS A MINIMUM'
  !write(u6,*) '    PRETTY FAR FROM WHERE WE GOT WITH THE LAST'
  !write(u6,*) '    STEP. THEREFORE; WE MAY CONSIDER CONTINUING'
  !write(u6,*) '    BY LS, ALONG THE SAME DIRECTION AS LAST STEP.'
  !write(u6,*) ' BUT WILL WE GAIN FROM THAT? LETS SEE!'
  if (EPRED_LS > EPRED_QN) then
    !write(u6,*) ' LS REJECTED -- GIVES TOO LITTLE.'
    if (NVEC == 0) QNSTEP = 'SX'
    if (NVEC > 0) QNSTEP = 'QN'
  else if (EPRED_LS > EPRED_SX) then
    !write(u6,*) ' USE SX.'
    QNSTEP = 'SX'
  else
    QNSTEP = 'LS'
  end if
end if
if (QNSTEP == 'LS') then
  !write(u6,*) ' USE LINE SEARCH.'
  XSX(:) = (TMIN-One)*XOLD(:)
end if
if (QNSTEP == 'QN') then
  !write(u6,*) ' USE QN STEP.'
  XSX(:) = XQN(:)
end if
!if (QNSTEP == 'SX') write(u6,*) ' USE SX STEP.'

! Often, use of qune results in getting caught in a cyclic or caotic
! attractor towards the end. As an attempt to break such a pattern,
! apply a scaling whenever energy goes up:
if (ENOW > ELAST) XSX(:) = 0.7_wp*XSX(:)

! Before committing the finally suggested step (which is now in XSX),
! also apply a step size limitation. The squared 2-norm of the vector XSX
! sum of squares of rotation angles, so the 2-norm is a strict limit on
! largest rotation angle, which we (arbitrarily) limit to 0.5 (say):
XSXNRM = DNRM2_(NDIM,XSX,1)
SCLFCT = One/(One+Two*XSXNRM)
XSX(:) = SCLFCT*XSX(:)

ELAST = ENOW
FPLAST = Two*DDOT_(NDIM,BK,1,XSX,1)
IAD = 0
call DDAFILE(LUQUNE,1,BK,NDIM,IAD)
call DDAFILE(LUQUNE,1,XQN,NDIM,IAD)
call DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
NLS = NLS+1
if (QNSTEP /= 'LS') NLS = 0

end subroutine QUNE
