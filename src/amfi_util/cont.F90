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

subroutine cont(L,breit,ifinite,TKIN,evec,eval,Energy,type1,type2,scratch)
!bs ####################################################################
!bs   cont prepares all required contraction coefficients for functions
!bs   with angular momentum L
!bs ####################################################################

use AMFI_global, only: cntscrtch, contrarray, exponents, MxcontL, MxprimL, ncontrac, normovlp, nprimit, OVLPinv, rootOVLP, &
                       rootOVLPinv
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, ifinite
logical(kind=iwp), intent(in) :: breit
real(kind=wp), intent(out) :: TKIN(MxprimL,MxprimL), evec(MxprimL,MxprimL), eval(MxprimL), Energy(MxprimL), type1(MxprimL), &
                              type2(MxprimL), scratch(MxprimL,MxprimL,3)
logical(kind=iwp), parameter :: breit_finite = .true.

!bs transcon transfers and normalizes contracted functions
!bs ore more precisely the coefficients
call transcon(cntscrtch(:,:,L),MxprimL,MxcontL,normovlp(:,:,L),contrarray(:,0,L),nprimit(L),ncontrac(L))
!bs gentkin generates the matrix of kinetic energy TKIN
call gentkin(L,TKIN,nprimit(L),exponents(:,L),rootOVLPinv(:,:,L))
!bs kindiag diagonalizes TKIN
!bs for finite nucleus
if ((ifinite == 2) .and. (L == 0)) then
  call kindiag(TKIN,nprimit(L),evec,eval,breit_finite)
else
  call kindiag(TKIN,nprimit(L),evec,eval,breit)
end if
!bs kinemat generates kinematic factors in
!bs the basis of eigenvectors
call kinemat(nprimit(L),eval,type1,type2,Energy)
!bs chngcont= changecont generates the contraction coeffs
!bs including kinematic factors and even exponents as factors
call chngcont(contrarray(:,0,L),contrarray(:,1,L),contrarray(:,2,L),contrarray(:,3,L),contrarray(:,4,L),ncontrac(L),nprimit(L), &
              evec,type1,type2,scratch(:,:,1),scratch(:,:,2),scratch(:,:,3),MxprimL,rootOVLP(:,:,L),OVLPinv(:,:,L),exponents(:,L))

return

end subroutine cont
