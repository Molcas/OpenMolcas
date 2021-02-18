************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine cont(L,breit,ifinite)
cbs###########################################################################
cbs   cont prepares all required contraction coefficients for functions
cbs   with angular momentum L
cbs###########################################################################
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
      dimension tkintria((MxprimL*MxprimL+MxprimL)/2)
      logical breit,breit_finite
      breit_finite=.true.
cbs   transcon transfers and normalizes contracted functions
cbs   ore more precizely the coefficients
      call transcon(cntscrtch(1,1,L),MxprimL,
     *MxcontL,normovlp(1,1,L),
     *contrarray(iaddori(L)),nprimit(L),ncontrac(L))
cbs   gentkin generates the matrix of kinetic energy  TKIN
      call gentkin(L,TKIN,nprimit(L),exponents(1,L),rootOVLPinv(1,1,L))
cbs   kindiag diagonalizes TKIN
cbs   for finite nucleus
      if (ifinite.eq.2.and.L.eq.0) then
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit_finite)
      else
      call kindiag(TKIN,TKINTRIA,nprimit(L),evec,eval,breit)
      endif
cbs   kinemat generates kinematic factors in
cbs   the basis of eigenvectors
      call kinemat(L,nprimit(L),eval,type1,type2,Energy)
cbs   chngcont= changecont generates the contraction coeffs
cbs   including kinematic factors and even exponents as factors
      call chngcont(
     *contrarray(iaddori(L)),
     *contrarray(iaddtyp1(L)),
     *contrarray(iaddtyp2(L)),
     *contrarray(iaddtyp3(L)),
     *contrarray(iaddtyp4(L)),
     *ncontrac(L),nprimit(L),evec,
     *type1,type2,scratch4,scratch4(nprimit(L)*nprimit(L)+1),
     *scratch4(2*nprimit(L)*nprimit(L)+1),MxprimL,
     *rootOVLP(1,1,L),OVLPinv(1,1,L),
     *exponents(1,L))
      return
      end
