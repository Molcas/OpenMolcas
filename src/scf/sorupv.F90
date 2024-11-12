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
! Copyright (C) 1994, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _DIAGONAL_ONLY_
subroutine SOrUpV(V,lvec,W,Mode,UpTp)
!***********************************************************************
!     for Ref., see T.H. Fischer and J. Almloef, JPC 96, 9768 (1992)   *
!               doi:10.1021/j100203a036                                *
!                                                                      *
!     purpose: Second Order Updated vector V using BFGS and            *
!              diagonal Hessian of Orbital Rotations                   *
!     input  : V          -> input vector                              *
!              HDiag      -> initial diagonal Hessian                  *
!              lvec       -> lengths of vectors delta, grad, HDiag & V *
!              Mode       -> update mode, see below                    *
!     The routine expects, that grad(n) and delta(n-1) are already     *
!     stored on the appropriate linked lists LLdGrd & LLDelt.          *
!     output:  W          -> H(n)V with H(n) SOrUp (inverse) Hessian   *
!                            after n-1 updates                         *
!                            mem has to be allocated from caller side  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     Adapted to work both for d=-H(-1)g and g = -Hd                   *
!     depending on the value of MODE                                   *
!     Mode='DISP':  d=-H(-1)g    BFGS update for the inverse Hessian   *
!     Mode='GRAD':  g=-Hd        DFP  update for the Hessian           *
!     2017, Roland Lindh, Harvard, Cambridge                           *
!                                                                      *
!     Adapted to do both BFGS and DFP updates depending on the value   *
!     of UPTP                                                          *
!     UpTp='BFGS':  BFGS update                                        *
!     UpTp='DFP ':  DFP update                                         *
!     Mode='DISP':  d=-H(-1)g    update for the Hessian                *
!     Mode='GRAD':  g=-Hd        update for the inverse Hessian        *
!***********************************************************************

use LnkLst, only: GetNod, iVPtr, LLDelt, LLdGrd, LLLen, LLy, LstPtr, PutVec, SCF_V
use Interfaces_SCF, only: yHx
use InfSCF, only: HDiag, Iter, IterSO, TimFld
use SCFFiles, only: LuDel, LuDGd
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lvec
real(kind=wp), intent(in) :: V(lVec)
real(kind=wp), intent(out) :: W(lvec)
character(len=4), intent(in) :: Mode, UpTp
integer(kind=iwp) :: i, inode, ipdel, ipdgd, ipynm1, it, leny, LL1, LL2, Lu1
real(kind=wp) :: Cpu1, Cpu2, S(6), T(4), Tim1, Tim2, Tim3
logical(kind=iwp) :: diag, Inverse_H, updy
character(len=8) :: Mode_Old = ''
real(kind=wp), allocatable :: SOGrd(:), SODel(:), SOScr(:)
real(kind=wp), parameter :: Thr = 1.0e-9_wp
real(kind=wp), external :: ddot_

if (DDot_(lVec,V,1,V,1) == Zero) then
  W(:) = Zero
  return
end if
call Timing(Cpu1,Tim1,Tim2,Tim3)
S(:) = Zero
T(:) = Zero

! Dummy initializations

Inverse_H = .false.
Lu1 = -1
LL1 = -1
LL2 = -1

! This section will control the mode

if (Mode == 'DISP') then

  ! Mode for computation d = -H(-1)g

  Inverse_H = .true.
  Lu1 = LuDel
  LL1 = LLDelt
  LL2 = LLdGrd

else if (Mode == 'GRAD') then

  ! Mode for computation g = -Hd

  Inverse_H = .false.
  Lu1 = LudGd
  LL1 = LLdGrd
  LL2 = LLDelt

end if

if (Lu1 < 0) then
  write(u6,*) 'SOrUpV: Illegal mode'
  call Abend()
end if

if ((UpTp /= 'BFGS') .and. (UpTp /= 'DFP ')) then
  write(u6,*) 'SOrUpV: Illegal update type',UpTp
  call Abend()
end if

if ((iterso > 1) .and. (Mode//UpTp /= Mode_Old)) then
  write(u6,*) 'IterSO=',IterSO
  write(u6,*) 'Mode_Old:',Mode_Old
  write(u6,*) 'Mode//UpTp:',Mode//UpTp
  write(u6,*) 'SOrUpV: Illegal mode switch'
  call Abend()
end if

! Steps denoted as in the reference paper!

! (1): initialize w=HDiag*v

if (Inverse_H) then
  do i=1,lvec
    if (abs(HDiag(i)) < Thr) then
      W(i) = 1.0e2_wp*V(i)
    else
      W(i) = V(i)/HDiag(i)
    end if
  end do
else
  !do i=1,lvec
  !  W(i) = HDiag(i)*V(i)
  !end do
  call yHx(V,W,lvec)
end if
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*)
call Check_Vec(W,lVec,'H_{n-1}v')
call NrmClc(V,lVec,'SORUPV','V')
call NrmClc(HDiag,lVec,'SORUPV','HDiag')
call NrmClc(W,lVec,'SORUPV','W')
#endif

diag = (iterso == 1)
#ifdef _DIAGONAL_ONLY_
write(u6,*) ' SorUpV: Only diagonal approximation'
diag = .true.
#endif
if (diag) then
  ! (2): no further work required this time
  call Timing(Cpu2,Tim1,Tim2,Tim3)
  TimFld(7) = TimFld(7)+(Cpu2-Cpu1)
  Mode_Old = Mode//UpTp
  return
end if

! need three vectors of length lvec as scratch space

call mma_allocate(SOGrd,lvec,Label='SOGrd')
call mma_allocate(SODel,lvec,Label='SODel')
call mma_allocate(SOScr,lvec,Label='SOScr')

call GetNod(iter-1,LL2,inode)
if (inode == 0) call Error_Handling()
call iVPtr(SOGrd,lvec,inode)

! (3b): initialize y(n-1)=HDiag*dGrd(n-1) ...

if (Inverse_H) then
  do i=1,lvec
    if (abs(HDiag(i)) < Thr) then
      SOScr(i) = 1.0e2_wp*SOGrd(i)
    else
      SOScr(i) = SOGrd(i)/HDiag(i)
    end if
  end do
else
  !do i=1,lvec
  !  SOScr(i) = HDiag(i)*SOGrd(i)
  !end do
  call yHx(SOGrd,SOScr,lvec)
end if
#ifdef _DEBUGPRINT_
!call RecPrt('Init y(n-1)',' ',SOScr,1,lVec)
call NrmClc(SOScr,lVec,'SOrUpV','Init y(n-1)')
call Check_vec(SOScr,lVec,'y_{n-1}=H_{0}Delta_{n-1}')
#endif

! and store it on appropriate linked list

leny = LLLen(LLy)
call PutVec(SOScr,lvec,iter-1,'NOOP',LLy)
if (leny == LLLen(LLy)) then
  ! already there, so we don't have to recalculate later
  updy = .false.
else
  ! new listhead generated in ipy, we have to update
  updy = .true.
end if

! (4): now loop over 1..n-2 iterations.

#ifdef _DEBUGPRINT_
write(u6,*) 'IterSO=',IterSO
#endif
do it=iter-iterso+1,iter-2

  ! fetch delta(i), dGrd(i) and y(i) from corresponding LLists

  call GetNod(it,LL1,inode)
  if (inode == 0) call Error_Handling()
  call iVPtr(SODel,lvec,inode)

  call GetNod(it,LL2,inode)
  if (inode == 0) call Error_Handling()
  call iVPtr(SOGrd,lvec,inode)

  call GetNod(it,LLy,inode)
  if (inode == 0) call Error_Handling()
  call iVPtr(SOScr,lvec,inode)

  ! calculate S_k and T_k dot products.
  ! (note that S(2) is the inverse of the one in the paper

  S(1) = ddot_(lvec,SODel,1,SOGrd,1)
  if (abs(S(1)) < Thr) then
    S(1) = Zero
    !S(1) = One/Thr
  else
    S(1) = One/S(1)
  end if
  S(2) = ddot_(lvec,SOGrd,1,SOScr,1)
  S(3) = ddot_(lvec,SODel,1,V,1)
  S(4) = ddot_(lvec,SOScr,1,V,1)
  if (updy) then

    ! here we have to reload dGrd(n-1) from llist, but this
    ! for sure is a memory hit, since it was put there last

    ipdgd = LstPtr(iter-1,LL2)

    S(5) = ddot_(lvec,SODel,1,SCF_V(ipdgd)%A,1)
    S(6) = ddot_(lvec,SOScr,1,SCF_V(ipdgd)%A,1)
  else
    S(5) = Zero
    S(6) = Zero
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'it=',it
  write(u6,*) '(S(i),i=1,6)=',(S(i),i=1,6)
# endif

  if ((Mode_Old == 'DISPBFGS') .or. (Mode_Old == 'GRADDFP ')) then
    if (abs(S(2)) < Thr) then
      S(2) = One
      !S(2) = One+S(1)/Thr
    else
      S(2) = One+S(1)*S(2)
    end if
    T(2) = S(1)*S(3)
    T(4) = S(1)*S(5)
    T(1) = S(2)*T(2)-S(1)*S(4)
    T(3) = S(2)*T(4)-S(1)*S(6)
  else if ((Mode_Old == 'DISPDFP ') .or. (Mode_Old == 'GRADBFGS')) then
    ! DFP update (or BFGS for the Hessian) use
    ! a slightly different and simpler formula
    if (abs(S(2)) < Thr) then
      S(2) = Zero
      S(2) = One/Thr
    else
      S(2) = One/S(2)
    end if
    T(1) = S(1)*S(3)
    T(2) = S(2)*S(4)
    T(3) = S(1)*S(5)
    T(4) = S(2)*S(6)
  end if

  ! Compute w and y(n-1)

# ifdef _DEBUGPRINT_
  write(u6,*) '(T(i),i=1,4)=',(T(i),i=1,4)
  call Check_Vec(W,lVec,'W(0)')
  W(:) = W(:)+T(1)*SODel(:)
  call Check_Vec(W,lVec,'W(1)')
  W(:) = W(:)-T(2)*SOScr(:)
  call Check_Vec(W,lVec,'W(2)')
# else
  W(:) = W(:)+T(1)*SODel(:)-T(2)*SOScr(:)
# endif
  if (updy) then

    ! here we have to reload y(n-1) from llist, but this
    ! for sure is a memory hit, since it was put there last
    ! -> we operate directly on the memory cells of the LList,
    ! where y(n-1) resides.

    ipynm1 = LstPtr(iter-1,LLy)
    SCF_V(ipynm1)%A(:) = SCF_V(ipynm1)%A(:)+T(3)*SODel(:)-T(4)*SOScr(:)
  end if

end do

! (5): reload y(n-1), delta(n-1) & dGrd(n-1) from linked list.
! these all are memory hits, of course

ipynm1 = LstPtr(iter-1,LLy)
ipdel = LstPtr(iter-1,LL1)
ipdgd = LstPtr(iter-1,LL2)
#ifdef _DEBUGPRINT_
write(u6,*)
!call RecPrt('y(n-1)',' ',SCF_V(ipynm1)%A,1,lVec)
call NrmClc(SCF_V(ipynm1)%A,lVec,'SOrUpV','y(n-1)')
if (Mode == 'DISP') then
  !call RecPrt('dX(n-1)',' ',SCF_V(ipdel)%A,1,lVec)
  !call RecPrt('dg(n-1)',' ',SCF_V(ipdgd)%A,1,lVec)
  call NrmClc(SCF_V(ipdel)%A,lVec,'SOrUpV','dX(n-1)')
  call NrmClc(SCF_V(ipdgd)%A,lVec,'SOrUpV','dg(n-1)')
else
  !call RecPrt('dg(n-1)',' ',SCF_V(ipdel)%A,1,lVec)
  !call RecPrt('dX(n-1)',' ',SCF_V(ipdgd)%A,1,lVec)
  call NrmClc(SCF_V(ipdel)%A,lVec,'SOrUpV','dg(n-1)')
  call NrmClc(SCF_V(ipdgd)%A,lVec,'SOrUpV','dX(n-1)')
end if
write(u6,*)
call Check_Vec(SCF_V(ipdel)%A,lvec,'delta_{n-1}')
call Check_Vec(SCF_V(ipdgd)%A,lvec,'Delta_{n-1}')
#endif

! calculate diverse dot products...

S(1) = ddot_(lvec,SCF_V(ipdel)%A,1,SCF_V(ipdgd)%A,1)
!write(u6,*) 'S(1)=',S(1)
if (abs(S(1)) < Thr) then
  S(1) = Zero
  !S(1) = One/Thr
else
  S(1) = One/S(1)
end if
!write(u6,*) 'S(1)=Alpha=',S(1)
!call Check_vec(SCF_V(ipynm1)%A,lVec,'y_{n-1}=H_{n-1}Delta_{n-1}')
S(2) = ddot_(lvec,SCF_V(ipdgd)%A,1,SCF_V(ipynm1)%A,1)
!write(u6,*) 'S(2)=Delta_{n-1}^T*y_{n-1}=',S(2)
S(3) = ddot_(lvec,SCF_V(ipdel)%A,1,V,1)
!write(u6,*) 'S(3)=delta_{n-1}^T*v=',S(3)
S(4) = ddot_(lvec,SCF_V(ipynm1)%A,1,V,1)
!write(u6,*) 'S(4)=y_{n-1}^T*v=',S(4)

#ifdef _DEBUGPRINT_
write(u6,*) '(S(i),i=1,4)=',(S(i),i=1,4)
#endif

if ((Mode_Old == 'DISPBFGS') .or. (Mode_Old == 'GRADDFP ')) then
  if (abs(S(2)) < Thr) then
    S(2) = One
    !S(2) = One+S(1)/Thr
  else
    S(2) = One+S(1)*S(2)
  end if
  !write(u6,*) 'S(2)(1+...)=',S(2)
  !write(u6,*)
  T(2) = S(1)*S(3)
  T(1) = S(2)*T(2)-S(1)*S(4)
  !write(u6,*) 'T(:)=',T(:)
else if ((Mode_Old == 'DISPDFP ') .or. (Mode_Old == 'GRADBFGS')) then
  if (abs(S(2)) < Thr) then
    S(2) = Zero
    S(2) = One/Thr
  else
    S(2) = One/S(2)
  end if
  T(1) = S(1)*S(3)
  T(2) = S(2)*S(4)
end if

! update the vector w

#ifdef _DEBUGPRINT_
write(u6,*) '(T(i),i=1,2)=',(T(i),i=1,2)
call Check_Vec(W,lVec,'W(2), again')
W(:) = W(:)+T(1)*SCF_V(ipdel)%A(:)
call Check_Vec(SCF_V(ipdel)%A,lvec,'delta_{n-1}')
call Check_Vec(W,lVec,'W(3)')
W(:) = W(:)-T(2)*SCF_V(ipdel)%A(:)
call Check_Vec(SCF_V(ipynm1)%A,lvec,'y_{n-1}')
!call RecPrt('The final W array',' ',W,1,lVec)
call NrmClc(W,lVec,'SOrUpV','The final W array')
call Check_Vec(W,lVec,'W(final)')
#else
W(:) = W(:)+T(1)*SCF_V(ipdel)%A(:)-T(2)*SCF_V(ipynm1)%A(:)
#endif

! so, we've made it, let's clean up workbench

call mma_deallocate(SOGrd)
call mma_deallocate(SODel)
call mma_deallocate(SOScr)
!
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(7) = TimFld(7)+(Cpu2-Cpu1)

return

contains

! Error handling
subroutine Error_handling()

  ! Hmmm, no entry found in LList, that's strange
  write(u6,*) 'SOrUpV: no entry found in LList'
  call Abend()

end subroutine Error_handling

#ifdef _DEBUGPRINT_
subroutine Check_Vec(Vec,nVec,Label)

  integer(kind=iwp), intent(in) :: nVec
  real(kind=wp), intent(in) :: Vec(nVec)
  character(len=*), intent(in) :: Label

  write(u6,*) 'Norm of ',Label,' is : ',sqrt(DDot_(nVec,Vec,1,Vec,1))

end subroutine Check_Vec
#endif

end subroutine Sorupv
