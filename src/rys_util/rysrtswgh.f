************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Walter Gautschi                                  *
*               1994, Gene H. Golub                                    *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
*
* Compute Rys roots and weights from scratch:
* - For high t values use the asymptotic scaled Hermite quadrature
* - For lower t, get first alpha and beta from the auxiliary Legendre
*   quadrature, then compute roots and weights
*
      Subroutine RysRtsWgh(TValues,nT,Roots,Weights,Order)
      Use Leg_RW
      Use vRys_RW
      Implicit None
      Integer, Intent(In) :: nT, Order
      Real*8, Intent(In) :: TValues(nT)
      Real*8, Intent(Out) :: Roots(Order,nT), Weights(Order, nT)
#include "stdalloc.fh"
#include "real.fh"
#include "FMM.fh"
      Integer :: i, j, iquad, Err
      Real*8, Dimension(:), Allocatable :: a, b
      Real*8 :: Alpha(Order), Beta(Order)
      Real*8, Parameter :: eps=1.0d-16
      Real*8, External :: TAsymp
      Integer, External :: WhichQuad
*
      Do i = 1, nT
        If (TValues(i).gt.TAsymp(Order) .or. asymptotic_Rys) Then
          Do j = 1, Order
            Roots(j,i)=HerR2(iHerR2(Order)+j-1)/TValues(i)
            Weights(j,i)=HerW2(iHerW2(Order)+j-1)/Sqrt(TValues(i))
          End Do
        Else
          iquad = WhichQuad(Order)
          Call mma_allocate(a,naux(iquad))
          Call mma_allocate(b,naux(iquad))
          Do j = 1, naux(iquad)
            a(j)=Leg_r(j,iquad)
            b(j)=Leg_w(j,iquad)*Exp(-TValues(i)*a(j))
          End Do
          Call Lanczos(Order,naux(iquad),a,b,Alpha,Beta,Err)
          If (Err.ne.0) Then
            write(6,*) Err
            Call WarningMessage(2,'Error in Lanczos')
            Call AbEnd()
          End If
          Call GaussQuad(Order,Alpha,Beta,eps,
     &                   Roots(1,i),Weights(1,i),Err)
          If (Err.ne.0) Then
            write(6,*) Err
            Call WarningMessage(2,'Error in GaussQuad 2')
            Call AbEnd()
          End If
          Call mma_deallocate(a)
          Call mma_deallocate(b)
        End If
      End Do
*
      End Subroutine RysRtsWgh
*
* This function returns the asymptotic limit for the t parameter,
* for values larger than this the scaled Hermite quadrature is
* accurate enough. These values are quite conservative, with
* estimated errors below 1e-16.
*
      Function TAsymp(Order)
      Implicit None
      Real*8 :: TAsymp
      Integer, Intent(In) :: Order
      Select Case(Order)
        Case (1)
          TAsymp=39.0D0
        Case (2)
          TAsymp=47.0D0
        Case (3)
          TAsymp=54.0D0
        Case (4)
          TAsymp=60.0D0
        Case (5)
          TAsymp=66.0D0
        Case (6)
          TAsymp=72.0D0
        Case (7)
          TAsymp=78.0D0
        Case (8)
          TAsymp=83.0D0
        Case (9)
          TAsymp=89.0D0
        Case (10)
          TAsymp=94.0D0
        Case (11)
          TAsymp=99.0D0
        Case (12)
          TAsymp=104.0D0
        Case (13)
          TAsymp=109.0D0
        Case (14)
          TAsymp=115.0D0
        Case (15)
          TAsymp=120.0D0
        Case (16)
          TAsymp=125.0D0
        Case (17)
          TAsymp=130.0D0
        Case (18)
          TAsymp=134.0D0
        Case (19)
          TAsymp=139.0D0
        Case (20)
          TAsymp=144.0D0
        Case Default
* Rough fit
          TAsymp=50.0D0+5*Order
      End Select
      End Function TAsymp
*
* This function returns the number of points to use in the auxiliary
* Legendre quadrature.
*
      Function WhichQuad(Order)
      Use Leg_RW
      Implicit None
      Integer :: WhichQuad
      Integer, Intent(In) :: Order
      Select Case(Order)
        Case (1)
          WhichQuad = 1 !24
        Case (2)
          WhichQuad = 1 !27
        Case (3)
          WhichQuad = 1 !30
        Case (4)
          WhichQuad = 2 !34
        Case (5)
          WhichQuad = 3 !37
        Case (6)
          WhichQuad = 3 !39
        Case (7)
          WhichQuad = 4 !42
        Case (8)
          WhichQuad = 4 !45
        Case (9)
          WhichQuad = 5 !46
        Case (10)
          WhichQuad = 5 !50
        Case (11)
          WhichQuad = 6 !51
        Case (12)
          WhichQuad = 6 !54
        Case (13)
          WhichQuad = 7 !56
        Case (14)
          WhichQuad = 7 !59
        Case (15)
          WhichQuad = 8 !61
        Case (16)
          WhichQuad = 8 !63
        Case (17)
          WhichQuad = 9 !66
        Case (18)
          WhichQuad = 9 !68
        Case (19)
          WhichQuad = 9 !70
        Case (20)
          WhichQuad = 10 !73
        Case Default
* Maximum naux
          WhichQuad = 11 !300
      End Select
      End Function WhichQuad
*
************************************************************************
* Routines GaussQuad and Lanczos adapted from:
*
* Algorithm 726: ORTHPOL -- A package of Routines for Generating
* Orthogonal Polynomials and Gauss-Type Quadrature Rules
*   Walter Gautschi. ACM Trans. Math. Softw. 20 (1994) 21-62
*   doi:10.1145/174603.174605
************************************************************************
*
      Subroutine GaussQuad(n,alpha,beta,eps,roots,weights,ierr)
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  roots(k)=x(k) and the weights as
c weights(k)=w(k), k=1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
c dlambda. The routine computes the nodes as eigenvalues, and the
c weights in term of the first component of the respective normalized
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature
c                      formula; type integer
c                alpha,beta - arrays of dimension  n  to be filled
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: roots - array of dimension  n  containing the Gaussian
c                      nodes (in increasing order)  roots(k)=x(k), k=1,2,
c                      ...,n
c                weights - array of dimension  n  containing the
c                      Gaussian weights  weights(k)=w(k), k=1,2,...,n
c                ierr - an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is
c                      negative.
      Implicit None
      Integer :: n,ierr,i,ii,j,k,l,m,mml
      Integer, Parameter :: maxcyc=30
      Real*8 :: alpha(n),beta(n),roots(n),weights(n),eps,e(n),
     &          b,c,f,g,p,s,r
#include "real.fh"
      ierr=0
      if (n.lt.1) then
        ierr=-1
        return
      end if
c
c Initialization
c
      do k=1,n
        roots(k)=alpha(k)
        if (beta(k).lt.zero) then
          ierr=-2
          return
        end if
        weights(k)=zero
        if (k.gt.1) e(k-1)=sqrt(beta(k))
      end do
      if (n.eq.1) then
        weights(1)=beta(1)
        return
      else
        weights(1)=one
        e(n) = zero
      endif
c
c Loop over roots
c
      do l=1,n
        do j=1,maxcyc
c
c Look for a small subdiagonal element.
c
          do m=l,n
            if (m.eq.n) exit
            if (abs(e(m)).le.eps*(abs(roots(m))+abs(roots(m+1)))) exit
          end do
          p=roots(l)
          if (m.eq.l) exit
c
c Form shift.
c
          g=(roots(l+1)-p)/(two*e(l))
          r=sqrt(g*g+one)
          g=roots(m)-p+e(l)/(g+sign(r,g))
          s=one
          c=one
          p=zero
          mml=m-l
c
c For i=m-1 step -1 until l do ...
c
          do ii=1,mml
            i=m-ii
            f=s*e(i)
            b=c*e(i)
            if (abs(f).lt.abs(g)) then
              s=f/g
              r=sqrt(s*s+one)
              e(i+1)=g*r
              c=one/r
              s=s*c
            else
              c=g/f
              r=sqrt(c*c+one)
              e(i+1)=f*r
              s=one/r
              c=c*s
            endif
            g=roots(i+1)-p
            r=(roots(i)-g)*s+two*c*b
            p=s*r
            roots(i+1)=g+p
            g=c*r-b
c
c Form first component of vector.
c
            f=weights(i+1)
            weights(i+1)=s*weights(i)+c*f
            weights(i)=c*weights(i)-s*f
          end do
          roots(l)=roots(l)-p
          e(l)=g
          e(m)=zero
        end do
c
c Set error - no convergence to an eigenvalue after maxcyc iterations.
c
        if (j.gt.maxcyc) then
          ierr=l
          return
        end if
      end do
c
c Order eigenvalues and eigenvectors.
c
      do ii=2,n
        i=ii-1
        k=i
        p=roots(i)
        do j=ii,n
          if (roots(j).lt.p) then
            k=j
            p=roots(j)
          end if
        end do
        if (k.ne.i) then
          roots(k)=roots(i)
          roots(i)=p
          p=weights(i)
          weights(i)=weights(k)
          weights(k)=p
        end if
      end do
      do k=1,n
        weights(k)=beta(1)*weights(k)*weights(k)
      end do
      return
      end
*
      Subroutine Lanczos(n,ncap,x,w,alpha,beta,ierr)
c This routine carries out the same task as the routine  sti, but
c uses the more stable Lanczos method. The meaning of the input
c and output parameters is the same as in the routine  sti. (This
c routine is adapted from the routine RKPW in W.B. Gragg and
c W.J. Harrod, "The numerically stable reconstruction of Jacobi
c matrices from spectral data", Numer. Math. 44, 1984, 317-335.)
c
c Routine sti:
c
c This routine applies "Stieltjes's procedure" (cf. Section 2.1 of
c W. Gautschi, "On generating orthogonal polynomials", SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow)
c or  +k  (in the case of overflow), where  k  is the recursion index
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close
c to  ncap, since there is a distinct possibility of numerical
c instability developing. (See W. Gautschi, "Is the recurrence relation
c for orthogonal polynomials always stable?", BIT, 1993, to appear.)
c In that case, the routine  lancz  should be used.
      Implicit None
      Integer :: n,ncap,ierr,i,k
      Real*8 :: x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap),
     &          gam,pj,rho,sig,t,tk,tmp,tsig,xlam
#include "real.fh"
      ierr=0
      if (n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      do i=1,ncap
        p0(i)=x(i)
        p1(i)=zero
      end do
      p1(1)=w(1)
      do i=1,ncap-1
        pj=w(i+1)
        gam=one
        sig=zero
        t=zero
        xlam=x(i+1)
        do k=1,i+1
          rho=p1(k)+pj
          tmp=gam*rho
          tsig=sig
          if (rho.le.zero) then
            gam=one
            sig=zero
          else
            gam=p1(k)/rho
            sig=pj/rho
          end if
          tk=sig*(p0(k)-xlam)-gam*t
          p0(k)=p0(k)-(tk-t)
          t=tk
          if (sig.le.zero) then
            pj=tsig*p1(k)
          else
            pj=(t**2)/sig
          end if
          tsig=sig
          p1(k)=tmp
        end do
      end do
      do k=1,n
        alpha(k)=p0(k)
        beta(k)=p1(k)
      end do
      return
      End
