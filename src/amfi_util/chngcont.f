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
      subroutine chngcont(coeffs,coeffst1,coeffst1a,coeffst2,
     *coeffst2a,ncont,nprims,evec,
     *type1,type2,work,work2,work3,MxprimL,
     *rootOVLP,OVLPinv,exponents)
c###############################################################################
cbs   purpose: makes out of old contraction coefficients(in normalized functions)
cbs   new coefficients including the kinematical factors
cbs   using the diagonal matrices on type1 and type2 (see subroutine kinemat)
cbs   coeffst1a and coeffst2a additionally include the exponents alpha
cbs   (that is why ....a). So the exponents in the integrals are moved
cbs   to the contraction coefficients and not in some way into the primitive
cbs   integrals.
cbs
cbs   the different cases for contracted integrals differ later on in the
cbs   choice of different sets of contraction coefficients.
cbs
c###############################################################################
      implicit real*8 (a-h,o-z)
      dimension coeffs(nprims,ncont),! original contraction coefficients
     *coeffst1(nprims,ncont),        ! A * cont coeff
     *coeffst1a(nprims,ncont),       ! A * alpha*cont coeff
     *coeffst2a(nprims,ncont),       ! c*A/(E+m) * cont coeff
     *coeffst2(nprims,ncont),        ! c*A/(E+m) * alpha *cont coeff
     *evec(nprims,nprims),
     *work(nprims,nprims) ,
     *work2(nprims,nprims) ,
     *work3(nprims,nprims) ,
     *rootOVLP(MxprimL,*),
     *OVLPinv(MxprimL,*),
     *type1(*),type2(*),
     *exponents(*)
cbs
cbs   first new coefficients for type1 (A)
cbs   generate a transformation matrix on work
cbs
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type1(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst1(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst1(I,K)=coeffst1(I,K)+work(J,I)*coeffs(J,K)
      enddo
      enddo
      enddo
cbs
cbs   now with exponent
cbs
      do K=1,ncont
      do I=1,nprims
      coeffst1a(I,K)=exponents(I)*coeffst1(I,K)
      enddo
      enddo
cbs
cbs   and now the same for the other type  A/(E+m)
cbs
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      work2(I,J)=0d0
      work3(I,J)=0d0
      enddo
      enddo
cbs   build up the transformation matrix
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(I,J)=work(I,J)+evec(I,K)*type2(K)*evec(J,K)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work2(I,J)=work2(I,J)+work(I,K)*rootOVLP(K,J)
      enddo
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work3(I,J)=work3(I,J)+rootOVLP(I,K)*work2(K,J)
      enddo
      enddo
      enddo
      do J=1,nprims
      do I=1,nprims
      work(I,J)=0d0
      enddo
      enddo
      do K=1,nprims
      do J=1,nprims
      do I=1,nprims
      work(J,I)=work(J,I)+OVLPinv(I,K)*work3(K,J)
      enddo
      enddo
      enddo
      do K=1,ncont
      do I=1,nprims
      coeffst2(I,K)=0d0
      enddo
      enddo
cbs   now transform the vectors
      do K=1,ncont
      do J=1,nprims
      do I=1,nprims
      coeffst2(I,K)=coeffst2(I,K)+work(J,I)*coeffs(J,K)
      enddo
      enddo
      enddo
cbs
cbs   now with exponent
cbs
      do K=1,ncont
      do I=1,nprims
      coeffst2a(I,K)=exponents(I)*coeffst2(I,K)
      enddo
      enddo
      return
      end
