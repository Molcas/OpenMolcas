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

subroutine chngcont(coeffs,coeffst1,coeffst1a,coeffst2,coeffst2a,ncont,nprims,evec,type1,type2,work,work2,work3,MxprimL,rootOVLP, &
                    OVLPinv,exponents)
!#######################################################################
!bs purpose: makes out of old contraction coefficients(in normalized functions)
!bs new coefficients including the kinematical factors
!bs using the diagonal matrices on type1 and type2 (see subroutine kinemat)
!bs coeffst1a and coeffst2a additionally include the exponents alpha
!bs (that is why ....a). So the exponents in the integrals are moved
!bs to the contraction coefficients and not in some way into the primitive
!bs integrals.
!bs
!bs the different cases for contracted integrals differ later on in the
!bs choice of different sets of contraction coefficients.
!#######################################################################
!coeffs    : original contraction coefficients
!coeffst1  : A * cont coeff
!coeffst1a : A * alpha*cont coeff
!coeffst2a : c*A/(E+m) * cont coeff
!coeffst2  : c*A/(E+m) * alpha *cont coeff

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ncont, nprims, MxprimL
real(kind=wp), intent(in) :: coeffs(nprims,ncont), evec(nprims,nprims), type1(*), type2(*), rootOVLP(MxprimL,*), &
                             OVLPinv(MxprimL,*), exponents(*)
real(kind=wp), intent(out) :: coeffst1(nprims,ncont), coeffst1a(nprims,ncont), coeffst2(nprims,ncont), coeffst2a(nprims,ncont), &
                              work(nprims,nprims), work2(nprims,nprims), work3(nprims,nprims)
integer(kind=iwp) :: K

!bs first new coefficients for type1 (A)
!bs generate a transformation matrix on work

!bs build up the transformation matrix
do K=1,nprims
  work2(:,K) = type1(K)*evec(:,K)
end do
call dgemm_('N','T',nprims,nprims,nprims,One,work2,nprims,evec,nprims,Zero,work,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,work,nprims,rootOVLP,MxprimL,Zero,work2,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,rootOVLP,MxprimL,work2,nprims,Zero,work3,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,OVLPinv,MxprimL,work3,nprims,Zero,work,nprims)
!bs now transform the vectors
call dgemm_('N','N',nprims,ncont,nprims,One,work,nprims,coeffs,nprims,Zero,coeffst1,nprims)

!bs now with exponent

do K=1,ncont
  coeffst1a(:,K) = exponents(1:nprims)*coeffst1(:,K)
end do

!bs and now the same for the other type  A/(E+m)

!bs build up the transformation matrix
do K=1,nprims
  work2(:,K) = type2(K)*evec(:,K)
end do
call dgemm_('N','T',nprims,nprims,nprims,One,work2,nprims,evec,nprims,Zero,work,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,work,nprims,rootOVLP,MxprimL,Zero,work2,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,rootOVLP,MxprimL,work2,nprims,Zero,work3,nprims)
call dgemm_('N','N',nprims,nprims,nprims,One,OVLPinv,MxprimL,work3,nprims,Zero,work,nprims)
!bs now transform the vectors
call dgemm_('N','N',nprims,ncont,nprims,One,work,nprims,coeffs,nprims,Zero,coeffst2,nprims)

!bs now with exponent

do K=1,ncont
  coeffst2a(:,K) = exponents(1:nprims)*coeffst2(:,K)
end do

return

end subroutine chngcont
