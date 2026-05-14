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

use MCLR_Data, only: H0F, H0S, nConf1, nExp, SBIDT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: rdia(nConf1), v(nConf1), alpha, beta
real(kind=wp), intent(inout) :: u(nConf1)
integer(kind=iwp) :: i, iRC, j
real(kind=wp), allocatable :: Tmp1(:), Tmp4(:)

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
    write(u6,*) 'Error in DGETRS called from exphinvv'
    call Abend()
  end if
end if

if ((alpha == Zero) .and. (beta == One)) then
  u(:) = rDia(:)*v(:)
else if (alpha == Zero) then
  u(:) = beta*rDia(:)*v(:)
else if (alpha == One) then
  u(:) = u(:)+beta*rDia(:)*v(:)
else
  u(:) = alpha*u(:)+beta*rDia(:)*v(:)
end if

if (nExp /= 0) then
  do i=1,nExp
    j = SBIDT(i)
    u(j) = alpha*Tmp4(i)+beta*Tmp1(i)
  end do
  call mma_deallocate(Tmp1)
  call mma_deallocate(Tmp4)
end if

end subroutine ExpHinvv
