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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine ExpHinvv(rdia,v,u,alpha,beta)
! Preconditioning of the state transfer part
! of the  electronic hessian with an subunit
! described with the Explicit hessian and
! the rest with the diagonal
!
!                            -1
! |u> = alpha|u> + beta  (H -E ) |v>
!                          0  0

use Exp, only: H0S, H0F, SBIDT, nExp
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: nConf1

implicit none
real*8 alpha, beta
real*8 v(*), u(*), rdia(*)
real*8, allocatable :: Tmp1(:), Tmp4(:)
integer i, j, iRC

if (nExp /= 0) then
  call mma_allocate(Tmp1,nExp,Label='Tmp1')
  call mma_allocate(Tmp4,nExp,Label='Tmp4')

  do i=1,nExp
    j = SBIDT(i)
    Tmp1(i) = v(j)
    Tmp4(i) = u(j)
  end do

  irc = 0
  call dgetrs_('N',NEXP,1,H0S,nexp,H0F,Tmp1,nexp,irc)
  if (irc /= 0) then
    write(6,*) 'Error in DGETRS called from exphinvv'
    call Abend()
  end if

  if ((alpha == 0.0d0) .and. (beta == 1.0d0)) then
    call DVEM(nConf1,v,1,rdia,1,u,1)
  else if (alpha == 0.0d0) then
    do i=1,nConf1
      u(i) = beta*rDia(i)*v(i)
    end do
  else if (alpha == 1.0d0) then
    do i=1,nConf1
      u(i) = u(i)+beta*rDia(i)*v(i)
    end do
  else
    do i=1,nConf1
      u(i) = alpha*u(i)+beta*rDia(i)*v(i)
    end do
  end if

  do i=1,nExp
    j = SBIDT(i)
    u(j) = alpha*Tmp4(i)+beta*Tmp1(i)
  end do
  call mma_deallocate(Tmp1)
  call mma_deallocate(Tmp4)

else
  if ((alpha == 0.0d0) .and. (beta == 1.0d0)) then
    call DVEM(nConf1,v,1,rdia,1,u,1)
  else if (alpha == 0.0d0) then
    do i=1,nConf1
      u(i) = beta*rDia(i)*v(i)
    end do
  else if (alpha == 1.0d0) then
    do i=1,nConf1
      u(i) = u(i)+beta*rDia(i)*v(i)
    end do
  else
    do i=1,nConf1
      u(i) = alpha*u(i)+beta*rDia(i)*v(i)
    end do
  end if

end if

end subroutine ExpHinvv
