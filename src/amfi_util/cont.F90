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

subroutine cont(L,breit,ifinite)
!bs ####################################################################
!bs   cont prepares all required contraction coefficients for functions
!bs   with angular momentum L
!bs ####################################################################

implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
dimension tkintria((MxprimL*MxprimL+MxprimL)/2)
logical breit, breit_finite

breit_finite = .true.
!bs transcon transfers and normalizes contracted functions
!bs ore more precizely the coefficients
call transcon(cntscrtch(1,1,L),MxprimL,MxcontL,normovlp(1,1,L),contrarray(iaddori(L)),nprimit(L),ncontrac(L))
!bs gentkin generates the matrix of kinetic energy  TKIN
call gentkin(L,TKIN,nprimit(L),exponents(1,L),rootOVLPinv(1,1,L))
!bs kindiag diagonalizes TKIN
!bs for finite nucleus
if ((ifinite == 2) .and. (L == 0)) then
  call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit_finite)
else
  call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit)
end if
!bs kinemat generates kinematic factors in
!bs the basis of eigenvectors
call kinemat(nprimit(L),eval,type1,type2,Energy)
!bs chngcont= changecont generates the contraction coeffs
!bs including kinematic factors and even exponents as factors
call chngcont(contrarray(iaddori(L)),contrarray(iaddtyp1(L)),contrarray(iaddtyp2(L)),contrarray(iaddtyp3(L)), &
              contrarray(iaddtyp4(L)),ncontrac(L),nprimit(L),evec,type1,type2,scratch4,scratch4(nprimit(L)*nprimit(L)+1), &
              scratch4(2*nprimit(L)*nprimit(L)+1),MxprimL,rootOVLP(1,1,L),OVLPinv(1,1,L),exponents(1,L))

return

end subroutine cont
