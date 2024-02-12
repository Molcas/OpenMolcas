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

subroutine newCF(H,n,A,B,C,Bstev)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: n
complex(kind=8), intent(in) :: H(n,n)
complex(kind=8), intent(out) :: A((n-1),-(n-1):(n-1))
real(kind=8), intent(out) :: B(n,0:n), C(n,0:n)
real(kind=8), intent(out) :: Bstev(n,-n:n)
! local variables:
integer :: ik, iq
real(kind=8) :: rfact, cr, mfact, C0
complex(kind=8) :: trace, cfact
complex(kind=8), allocatable :: Cp(:,:), Cm(:,:)
complex(kind=8) :: mf
real(kind=8) :: knm(12,0:12)
external :: trace
logical :: dbg

!-------------------------------------------
if (n < 1) return
!-------------------------------------------
dbg = .false.
call mma_allocate(Cp,n,n,'operator O')
call mma_allocate(Cm,n,n,'operator W')
!-------------------------------------------
! n=2*J+1;  or   n=2*S+1
Bstev(1:n,-n:n) = 0.0_wp
B(1:n,0:n) = 0.0_wp
C(1:n,0:n) = 0.0_wp
A(1:(n-1),-(n-1):(n-1)) = (0.0_wp,0.0_wp)
call set_knm(knm)

do ik=1,n-1
  do iq=0,ik
    cr = 0.0_wp
    mfact = 0.0_wp
    rfact = 0.0_wp
    cfact = (0.0_wp,0.0_wp)
    C0 = 0.0_wp
    ! generate the operator matrix K=ik, Q=iq, dimension = n
    call ITO(n,ik,iq,C0,Cp,Cm)
    call coeff_redus_sub(n,ik,cr)

    mfact = dble((-1)**iq)
    rfact = C0*C0*dble(2*ik+1)/dble(n)
    cfact = cmplx(mfact*rfact,0.0_wp,wp)

    !-------------------------------------------
    ! Naoya's C/C0 operators:
    a(ik,-iq) = cfact*trace(n,H,Cp)
    a(ik,iq) = cfact*trace(n,H,Cm)

    !-------------------------------------------
    ! make real combinations of CF parameters:
    ! Liviu's ITO operators:
    if (iq == 0) then
      !b(ik, iq)=dble( (0.5_wp,0.0_wp)*(A(ik,iq)+A(ik,-iq)) )
      b(ik,iq) = dble(A(ik,iq))
    else
      mf = cmplx((-1)**iq,0.0_wp,wp)
      b(ik,iq) = dble(A(ik,-iq)+mf*A(ik,iq))
      c(ik,iq) = dble((A(ik,-iq)-mf*A(ik,iq))*(0.0_wp,-1.0_wp))
    end if
    ! scale with the correct ratio:
    b(ik,iq) = b(ik,iq)/(cr*C0)
    c(ik,iq) = c(ik,iq)/(cr*C0)

    !-------------------------------------------
    ! parameters to be used in connection with ESO as in MATLAB
    ! EasySpin program
    ! The ESO operaors are not implemented for k>12 in EasySpin
    ! therefore we do not provide these parameters as well.
    if ((ik <= 12) .and. (iq <= 12)) then
      ! scale with the correct ratio:
      if (iq == 0) then
        bstev(ik,iq) = b(ik,iq)*knm(ik,iq)
      else
        bstev(ik,iq) = b(ik,iq)*knm(ik,iq)
        bstev(ik,-iq) = c(ik,iq)*knm(ik,iq)
      end if
    end if

    if (dbg) write(6,'(A,2I3,5(ES20.13,1x))') 'k,q, b(k,q), c(k,q)',ik,iq,b(ik,iq),c(ik,iq)
  end do
end do

call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine newCF
