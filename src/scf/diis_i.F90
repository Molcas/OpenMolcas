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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2016,2022, Roland Lindh                                *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _NEW_CODE_
subroutine DIIS_i(CInter,nCI,TrDh,TrDP,TrDD,nTr,nD,iOpt_DIIS,Ind)
!***********************************************************************
!                                                                      *
!     purpose: density matrix optimization                             *
!                                                                      *
!              EDIIS optimization:                                     *
!                   K. N. Kudin, G. E. Scuseria, and E. Cances         *
!                   JCP, 116, , 8255 (2002)                            *
!                   doi:10.1063/1.1470195                              *
!                                                                      *
!     input:                                                           *
!       TrDh    : Traces of D(i)*h of size (nTr)                       *
!       TrDP    : Traces of D(i)*P(j) of size (nTr,nTr)                *
!       TrDD    : Traces of D(i)*D(j) of size (nTr,nTr)                *
!                                                                      *
!     output:                                                          *
!       CInter  : Interpolation coefficients of length nCI             *
!                                                                      *
!***********************************************************************

use SpinAV, only: Do_SpinAV
use InfSCF, only: AccCon, Elst, EmConv, Iter, kOptim, MxIter, MxOptm, TimFld, WarnPOcc
use Constants, only: Zero, One, Half, Quart
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCI, nTr, nD, iOpt_DIIS
real(kind=wp), intent(inout) :: CInter(nCI,nD)
real(kind=wp), intent(in) :: TrDh(nTr,nTr,nD), TrDP(nTr,nTr,nD), TrDD(nTr,nTr,nD)
integer(kind=iwp), intent(out) :: Ind(MxOptm)
integer(kind=iwp) :: i, iD, ii, j, jj, kk, n1, n_Min, nn
real(kind=wp) :: Big, BigOne, BigTwo, CPU1, CPU2, CSUM, DD(MxOptm**2,2), DiFi, DiFj, DiFn, DjFi, DjFj, DnFj, DnFn, E_Min, E_n, &
                 E_n1, E_Pred, E_Tot, Eline(MxOptm,2), EPred(MxIter+1) = Zero, Equad(MxOptm**2,2), h = 0.35_wp, R2, r_SO, Tim1, &
                 Tim2, Tim3, Tmp_A, Tmp_B

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Timing(Cpu1,Tim1,Tim2,Tim3)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'E_pred=',(EPred(i),i=1,iter)
if (nD == 1) then
  write(u6,*) 'E_actu=',(Elst(i,1),i=1,iter)
else
  write(u6,*) 'E_actu=',(Elst(i,1)+Elst(i,2),i=1,iter)
end if
#endif
if (kOptim == 1) then
  AccCon = 'None     '
  return
end if
if (iOpt_DIIS == 1) then
  AccCon = 'EDIIS    '
else
  AccCon = 'ADIIS    '
end if

#ifdef _NEW_CODE_
Ind(1:kOptim) = 0
do i=kOptim,1,-1

  tmp0 = Zero
  do j=1,iter

    if (any(Ind(i+1:kOptim) == j)) cycle

    tmp = sum(Elst(j,1:nD))

    if (tmp < tmp0) then
      tmp = tmp0
      Ind(i) = j
    end if

  end do
end do
#else
do i=1,kOptim
  Ind(i) = iter-kOptim+i
end do
#endif

do iD=1,nD

  ! EDIIS optimization, doi:10.1063/1.1470195, Eq. (8)

  ! Noticed the change in sign - in optim the quadratic terms,
  ! however, are added in the evaluation of Eq. (8). Additionally,
  ! the paper is written in a spin-orbital notation. A factor of
  ! 1/2 for the two-electron term is assumed and does not need to
  ! be included.

  if (iOpt_DIIS == 1) then

    do i=1,kOptim
      ii = Ind(i)
      Eline(i,iD) = Elst(ii,iD)
      do j=1,kOptim
        jj = Ind(j)
        DiFi = TrDh(ii,ii,iD)+TrDP(ii,ii,iD)
        DiFj = TrDh(ii,ii,iD)+TrDP(ii,jj,iD)
        DjFi = TrDh(jj,jj,iD)+TrDP(jj,ii,iD)
        DjFj = TrDh(jj,jj,iD)+TrDP(jj,jj,iD)
        Equad(kOptim*(i-1)+j,iD) = -Half*(DiFi-DiFj-DjFi+DjFj)
        DD(kOptim*(i-1)+j,iD) = TrDD(ii,jj,iD)
      end do
    end do

  else

    ! ADIIS optimization (Eq. 7). We arbitrarily set E(Dn) to zero.
    ! The option is not tested yet.

    kk = iter ! This needs to be properly set!
    do i=1,kOptim
      ii = Ind(i)
      DiFn = TrDh(ii,ii,iD)+TrDP(ii,kk,iD)
      DnFn = TrDh(kk,kk,iD)+TrDP(kk,kk,iD)
      Eline(i,iD) = DiFn-DnFn
      do j=1,kOptim
        jj = Ind(j)
        DiFj = TrDh(ii,ii,iD)+TrDp(ii,jj,iD)
        DnFj = TrDh(kk,kk,iD)+TrDp(kk,jj,iD)
        Equad(kOptim*(i-1)+j,iD) = DiFj-DiFn-DnFj+DnFn
        DD(kOptim*(i-1)+j,iD) = TrDD(ii,jj,iD)
      end do
    end do

  end if

end do

! Tweak for UHF

if (nD == 2) then
  do i=1,kOptim
    tmp_a = Eline(i,1)
    tmp_b = Eline(i,2)
    Eline(i,1) = tmp_a+tmp_b
    Eline(i,2) = tmp_a
    do j=1,kOptim
      tmp_a = Equad(kOptim*(i-1)+j,1)
      tmp_b = Equad(kOptim*(i-1)+j,2)
      Equad(kOptim*(i-1)+j,1) = tmp_a+tmp_b
      Equad(kOptim*(i-1)+j,2) = tmp_a
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (kOptim >= 3) then
  BigOne = Zero
  BigTwo = Zero
  do iD=1,nD
    do i=2,kOptim
      do j=1,i-1
        BigOne = max(BigOne,abs(Eline(i,iD)-Eline(j,iD)))
      end do
    end do
    do i=1,kOptim
      do j=1,kOptim
        BigTwo = max(BigTwo,abs(Equad(i+(j-1)*kOptim,iD)))
      end do
    end do
  end do
  Big = max(BigOne,BigTwo)
  if (Big < 1.0e-8_wp) then
    EmConv = .true.
    WarnPocc = .true.
  else
    EmConv = .false.
    WarnPocc = .false.
  end if
  if (Do_SpinAV) WarnPocc = .false. ! it is not a diagnostic in this case
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find interpolation coefficients with the constaints
!  Sum_i c_i = 1
!  0 =< c_i =< 1

call Optim(E_Pred,Eline,Equad,CInter(1:kOptim,1),kOptim)
EPred(iter+1) = E_Pred

#ifdef _DEBUGPRINT_
write(u6,*) ' Interpolation coefficients:'
write(u6,'(5f16.8)') (CInter(i,1),i=1,kOptim)
#endif

! Temporary fix for UHF

if (nD == 2) CInter(:,2) = CInter(:,1)

#ifdef _DEBUGPRINT_
write(u6,*) ' Interpolation coefficients:'
write(u6,'(5f16.8)') (CInter(i,1),i=1,kOptim)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Update the trust radius, h, by checking to which degree previous
! interpolations estimated the next energy correctly.
!
! The first two energies have not been predicted!

if (iter > 3) then
  E_Pred = EPred(iter)
  E_n1 = sum(Elst(iter,1:nD))
  E_n = sum(Elst(iter-1,1:nD))

# ifdef _DEBUGPRINT_
  write(u6,*) 'iter=',iter
  write(u6,*) 'Energy of iter  =',E_n1
  write(u6,*) 'E_pred of iter  =',E_Pred
  write(u6,*) 'Energy of iter-1=',E_n
  write(u6,*)
# endif

  ! In some cases the DIIS will come out with a set to
  ! coefficient which corresponds to the last density. In this
  ! case the E_Pred(i+1)=E(i). We add a small number to avoid
  ! dividing with zero.

  r_SO = (E_n1-E_n)/(E_Pred-E_n+1.0e-12_wp)
else
  r_SO = One
end if

! Update the trust radius according to this ad hoc scheme.

if (r_SO >= 0.75_wp) then
  h = 1.2_wp*h
else if (r_SO < Quart) then
  h = 0.7_wp*h
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the reference density for the density with the lowest energy.

E_Min = Zero
n_min = 0
n1 = 0
do i=1,kOptim
  E_tot = sum(Elst(Ind(i),1:nD))
  if (E_tot < E_Min) then
    n1 = n_min
    n_min = i
    E_Min = E_tot
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the square of the trace of the difference between the
! reference density and the new interpolated density.

if (.false.) then
  nn = Ind(n_min)
  r2 = Zero
  do i=1,kOptim
    ii = Ind(i)
    do j=1,kOptim
      jj = Ind(j)
      r2 = r2+sum(CInter(i,1:nD)*CInter(j,1:nD)*(TrDD(ii,jj,1:nD)-TrDD(nn,jj,1:nD)-TrDD(ii,nn,1:nD)+TrDD(nn,nn,1:nD)))
    end do
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Perform a RS-DIIS interpolation if required

  if (sqrt(r2) > h) then
    write(u6,*) 'Apply optimization with step restriction'
    write(u6,*) 'r,h =',sqrt(r2),h
    call Abend()
    call Optim2(E_Pred,Eline,Equad,DD,CInter(1:kOptim,1),kOptim,n_min,n1,r2)
    EPred(iter+1) = E_Pred

    ! Temporary fix for UHF

    if (nD == 2) CInter(:,2) = CInter(:,1)

  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Check that the coefficients sum up correctly.

do iD=1,nD
  CSum = sum(CInter(1:kOptim,iD))
  if (abs(CSum-One) > 1.0e-5_wp) then
    write(u6,*) 'diis_i: Abs(CSum - One) > 1.0e-5'
    write(u6,*) 'CSum=',CSum
    call Abend()
  end if

  ! If the coefficient for the last density is zero we are in problem.

  if ((CInter(kOptim,iD) == Zero) .and. (kOptim > 2) .and. (CInter(kOptim-1,iD) == One)) then
    write(u6,*) 'DIIS_I optimization failed!'
    !write(u6,*) 'iD=',iD
    !write(u6,*) (CInter(i,iD),i=1,kOptim)
    CInter(kOptim,iD) = 1.0e-6_wp
    CInter(kOptim,iD-1) = One-1.0e-6_wp
    !EmConv = .true.
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(6) = TimFld(6)+(Cpu2-Cpu1)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine DIIS_i
