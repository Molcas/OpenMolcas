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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne, Onei
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: H(n,n)
complex(kind=wp), intent(out) :: A(n-1,-(n-1):n-1)
real(kind=wp), intent(out) :: B(n,0:n), C(n,0:n), Bstev(n,-n:n)
integer(kind=iwp) :: ik, iq
real(kind=wp) :: C0, cr, knm(12,0:12), mfact, rfact
complex(kind=wp) :: cfact, mf
complex(kind=wp), allocatable :: Cm(:,:), Cp(:,:)
complex(kind=wp), external :: trace

!-------------------------------------------
if (n < 1) return
!-------------------------------------------
call mma_allocate(Cp,n,n,'operator O')
call mma_allocate(Cm,n,n,'operator W')
!-------------------------------------------
! n=2*J+1;  or   n=2*S+1
call set_knm(knm)

Bstev(:,:) = Zero
B(:,:) = Zero
C(:,:) = Zero
A(:,:) = cZero
do ik=1,n-1
  do iq=0,ik
    ! generate the operator matrix K=ik, Q=iq, dimension = n
    call ITO(n,ik,iq,C0,Cp,Cm)
    call coeff_redus_sub(n,ik,cr)

    mfact = real((-1)**iq,kind=wp)
    rfact = C0*C0*real(2*ik+1,kind=wp)/real(n,kind=wp)
    cfact = mfact*rfact*cOne

    !-------------------------------------------
    ! Naoya's C/C0 operators:
    a(ik,-iq) = cfact*trace(n,H,Cp)
    a(ik,iq) = cfact*trace(n,H,Cm)

    !-------------------------------------------
    ! make real combinations of CF parameters:
    ! Liviu's ITO operators:
    if (iq == 0) then
      !b(ik,iq) = real(Half*(A(ik,iq)+A(ik,-iq)))
      b(ik,iq) = real(A(ik,iq))
    else
      mf = real((-1)**iq,kind=wp)*cOne
      b(ik,iq) = real(A(ik,-iq)+mf*A(ik,iq))
      c(ik,iq) = real(-Onei*(A(ik,-iq)-mf*A(ik,iq)))

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

#   ifdef _DEBUGPRINT_
    write(u6,'(A,2I3,5(ES20.13,1x))') 'k,q, b(k,q), c(k,q)',ik,iq,b(ik,iq),c(ik,iq)
#   endif
  end do
end do

call mma_deallocate(Cp)
call mma_deallocate(Cm)

return

end subroutine newCF
