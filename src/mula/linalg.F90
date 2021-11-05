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
! Copyright (C) 1983, Per Ake Malmqvist                                *
!               1994,1995, Niclas Forsberg                             *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
!!
!!  Contains:
!!    Cholesky    : Performs a Cholesky factorization of a positive
!!                  definite matrix A.
!!    Dool_MULA   : Solves the system of linear equations A*X = B.
!!    SolveSecEq  : Solves the seqular equation A*C = S*C*D.
!!    Polfit      : Fit a polynomial of requested dimension to the input
!!                  data.
!!    calcNorm    : Calculates norm of a vector.
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
!!-----------------------------------------------------------------------!
!!

!!
      Subroutine Dool_MULA(A,LA1,LA2,B,LB1,LB2,det)
!!
!!  Purpose:
!!    Solve A*X = B
!!
!!  Input:
!!    A      : Real*8 two dimensional array
!!    B      : Real*8 two dimensional array
!!
!!  Output:
!!    B      : Real*8 two dimensional array - contains
!!             the solution X.
!!  Written by:
!!    Per-AAke Malmquist
!!    Dept. of Theoretical Chemistry, Lund University, 1983.
!!
!!  Modified by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
      Implicit Real*8 ( a-h,o-z )
      Real*8 A(LA1,LA2)
      Real*8 B(LB1,LB2)
#include "WrkSpc.fh"
!!
!!---- Initialize.
      ip=-9999999
      jp=-9999999
      n=LA2
      m=LB2
      Call GetMem('Buf','Allo','Real',ipBuf,n)
      Call GetMem('iPiv','Allo','Inte',ipiPiv,n)
      Call GetMem('jPiv','Allo','Inte',ipjPiv,n)
!!
      Do i = 1,n
      iWork(ipiPiv+i-1) = i
      iWork(ipjPiv+i-1) = i
      End Do
      det = 1.0d0
      Do i = 1,n
!!---- Now find better pivot element.
      Amax = -1.0d0
      Do k = i,n
      Do l = i,n
      Am = abs(A(iWork(ipiPiv+k-1),iWork(ipjPiv+l-1)))
      If ( Amax.le.Am ) Then
      Amax = Am
      ip = k
      jp = l
      End If
      End Do
      End Do
      If ( ip.ne.i ) Then
      det = -det
      iTemp = iWork(ipiPiv+i-1)
      iWork(ipiPiv+i-1) = iWork(ipiPiv+ip-1)
      iWork(ipiPiv+ip-1) = iTemp
      End If
      If ( jp.ne.i ) Then
      det = -det
      jTemp = iWork(ipjPiv+i-1)
      iWork(ipjPiv+i-1) = iWork(ipjPiv+jp-1)
      iWork(ipjPiv+jp-1) = jTemp
      End If
      ip = iWork(ipiPiv+i-1)
      jp = iWork(ipjPiv+i-1)
      diag = A(ip,jp)
!          Buf(i) = diag
      Work(ipBuf+i-1) = diag
      det = det*diag
      Do k = i+1,n
      kp = iWork(ipiPiv+k-1)
      c = A(kp,jp)/diag
      A(kp,jp) = c
      Do l = i+1,n
      lp = iWork(ipjPiv+l-1)
      A(kp,lp) = A(kp,lp)-c*A(ip,lp)
      End Do
      End Do
      End Do
!!
!!---- First resubstitution step.
      Do j = 1,m
      Do i = 2,n
      ip = iWork(ipiPiv+i-1)
      sum = B(ip,j)
      Do k = 1,i-1
      sum = sum-A(ip,iWork(ipjPiv+k-1))*                                &
     &          B(iWork(ipiPiv+k-1),j)
      End Do
      B(ip,j) = sum
      End Do
      End Do
!!
!!---- Second resubstitution step.
      Do j = 1,m
      Do i = n,1,-1
      ip = iWork(ipiPiv+i-1)
      sum = B(ip,j)
      Do k = i+1,n
      sum = sum-A(ip,iWork(ipjPiv+k-1))*                                &
     &          B(iWork(ipiPiv+k-1),j)
      End Do
      B(ip,j) = sum/Work(ipBuf+i-1)
      End Do
      End Do
!!
!!---- Reorganization part.
      Do j = 1,m
      Do i = 1,n
      Work(ipBuf+i-1) = B(iWork(ipiPiv+i-1),j)
      End Do
      Do i = 1,n
      B(iWork(ipjPiv+i-1),j) = Work(ipBuf+i-1)
      End Do
      End Do
!!
      Call GetMem('Buf','Free','Real',ipBuf,n)
      Call GetMem('iPiv','Free','Real',ipiPiv,n)
      Call GetMem('jPiv','Free','Real',ipjPiv,n)
!!
      End
!!
!!-----------------------------------------------------------------------!

!!-----------------------------------------------------------------------!
!!
      Subroutine SolveSecEq(A,n,C,S,D)
!!
!!  Purpose:
!!    Solve the secular equation SAC = CD.
!!
!!  Input:
!!    A        : Real*8 two dimensional array
!!    S        : Real*8 two dimensional array
!!
!!  Output:
!!    C        : Real*8 two dimensional array
!!    D        : Real*8 array
!!
!!  Calls:
!!    Jacob
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1994.
!!
      Implicit Real*8 ( a-h,o-z )
      Integer n
      Real*8 A(n,n),S(n,n),C(n,n), D(n)
!!---- Local declarations.
#include "WrkSpc.fh"

!!
!!---- Initialize.
      nSqr = n**2
      nSqrTri = n*(n+1)/2
!!D Write(6,*)'SolveSecEq test prints.'
!!D Write(6,*)'Matrix A:'
!!D do i=1,n
!!D Write(6,'(1x,5F16.8)')(A(i,j),j=1,n)
!!D end do
!!
!!---- get memory for temporary matrices.
      Call GetMem('Scr','Allo','Real',ipScr,nSqrTri)
      Call GetMem('T','Allo','Real',ipT,nSqr)
      Call GetMem('Temp','Allo','Real',ipTemp,nSqr)
      Call GetMem('Asymm','Allo','Real',ipAsymm,nSqr)

!!
!!---- Transform S to lower packed storage in Scratch.
      k = 1
      Do i = 1,n
      Do j = 1,i
      Work(ipScr+k-1) = S(i,j)
      k = k+1
      End Do
      End Do
!!
!!---- Turn T into a unit matrix.
!vv       T = 0.0d0
      call dcopy_(nSqr,[0.0D0],0,Work(ipT),1)
      Do i = 1,n
      Work(ipT+i+n*(i-1)-1) = 1.0d0
      End Do
!!
!!---- Diagonalize Scratch and scale each column of T with the square
!!     root of the corresponding eigenvalue.
      Call Jacob(Work(ipScr),Work(ipT),n,n)
      Do j = 1,n
      jj = j*(j+1)/2
      Scale = Sqrt(Work(ipScr+jj-1))
      Do i = 1,n
      Work(ipT+i+n*(j-1)-1) = Work(ipT+i+n*(j-1)-1)*Scale
      End Do
      End Do
!!
!!---- Make A symmetric and transform it to lower packed storage in
!!     Scratch.
      Call DGEMM_('N','N',                                              &
     &            n,n,n,                                                &
     &            1.0d0,A,n,                                            &
     &            Work(ipT),n,                                          &
     &            0.0d0,Work(ipTemp),n)
      Call DGEMM_('T','N',                                              &
     &            n,n,n,                                                &
     &            1.0d0,Work(ipT),n,                                    &
     &            Work(ipTemp),n,                                       &
     &            0.0d0,Work(ipAsymm),n)
      k = 1
      Do i = 1,n
      Do j = 1,i
      Work(ipScr+k-1) = Work(ipAsymm+i+n*(j-1)-1)
      k = k+1
      End Do
      End Do
!!
!!---- Diagonalize Scratch.
      Call Jacob(Work(ipScr),Work(ipT),n,n)
      Call JacOrd(Work(ipScr),Work(ipT),n,n)
!!
!!---- Store the eigenvalues in array D.
      Do i = 1,n
      ii = i*(i+1)/2
      D(i) = Work(ipScr+ii-1)
      End Do
!vv       C = T
      call dcopy_(nSqr,Work(ipT),1,C,1)
!!
!!---- Free memory of temporary matrices.
      Call GetMem('Scr','Free','Real',ipScr,nSqrTri)
      Call GetMem('T','Free','Real',ipT,nSqr)
      Call GetMem('Temp','Free','Real',ipTemp,nSqr)
      Call GetMem('Asymm','Free','Real',ipAsymm,nSqr)
!!
      End
!!
!!-----------------------------------------------------------------------!
      Subroutine PolFit(ipow,nvar,var,yin,ndata,                        &
     &       coef,nterm,stand_dev,max_err,diff_vec,                     &
     &       use_weight)
!!
!!  Purpose:
!!    Fit a polynomial to yin.
!!
!!  Written by:
!!    P-AA Malmquist
!!
!!  Modified by:
!!    Niclas Forsberg
!!
      Implicit Real*8 ( a-h,o-z )
      Parameter (mxdeg = 6)
      Integer ipow(nvar,nterm)

      Real*8 var(ndata,nvar)
      Real*8 yin( ndata )
      Real*8 coef( nterm,1 )
      Real*8 stand_dev,max_err
      Real*8 yfit( ndata )
      Real*8 vpow ( 0:mxdeg,nvar )
      Real*8 equmat( nterm,nterm )
      Real*8 rhs( nterm ),term( nterm )
!!
      Real*8 diff_vec(ndata)
      Logical  use_weight
#include "WrkSpc.fh"
!!
!!---- Initialize.
      NrOfVar   = nvar
      nPolyTerm  = nterm
      myCoef2= 1
!vv       rhs    = 0.0d0
      call dcopy_(nterm,[0.0D0],0,rhs,1)
!vv       equmat = 0.0d0
      call dcopy_(nterm*nterm,[0.0D0],0,equmat,1)

!!
!!---- Set up weight vector.
      Call GetMem('weight','Allo','Real',ipweight,ndata)
      If ( use_weight ) Then
      e_min = yin(1)
      e_max = yin(1)
      Do idata = 2,ndata
      If ( yin(idata).lt.e_min ) e_min = yin(idata)
      If ( yin(idata).gt.e_max ) e_max = yin(idata)
      End Do
      e_range = e_max-e_min
      Do idata = 1,ndata
      Work(ipweight+idata-1) = 1.0d0/(1.0d0+1000.0d0*                   &
     &       ((yin(idata)-e_min)/e_range))
      End Do
      Else
!vv          weight = 1.0d0
      call dcopy_(ndata,[0.1D0],0,Work(ipweight),1)
      End If
!!
!!---- Accumulate equation matrix and right-hand-side.
      Do idata = 1,ndata
      !---- Calculate powers of individual variable values.
      Do ivar = 1,NrOfVar
      pow = 1.0d0
      vpow(0,ivar) = 1.0d0
      Do i = 1,mxdeg
      pow = pow*var(idata,ivar)
      vpow(i,ivar) = pow
      End Do
      End Do
!!---- Calculate value of each polynomial term at this point.
      Do iterm = 1,nPolyTerm
      ip = ipow(iterm,1)
      t = vpow(ip,1)
      Do ivar = 2,NrOfVar
      ip = ipow(iterm,ivar)
      t = t*vpow(ip,ivar)
      End Do
      term(iterm) = t
      End Do
!!---- Accumulate equmat and rhs.
      Do iterm = 1,nPolyTerm
      rhs(iterm) = rhs(iterm)+yin(idata)*                               &
     &       term(iterm)*Work(ipweight+idata-1)
      Do jterm = 1,iterm
      equmat(iterm,jterm) = equmat(iterm,jterm)+                        &
     &          term(iterm)*term(jterm)*Work(ipweight+idata-1)
      End Do
      End Do
      End Do
!!
!!---- Set upper triangle of equmat by symmetry.
      Do iterm = 1,nPolyTerm-1
      Do jterm = iterm+1,nPolyTerm
      equmat(iterm,jterm) = equmat(jterm,iterm)
      End Do
      End Do
!!
!!---- Solve the resulting equation system.
      Do i=1,nterm
      coef(i,1) = rhs(i)
      End do
      Call Dool_MULA(equmat,nPolyTerm,nPolyTerm,coef,                   &
     &  nPolyTerm,MyCoef2,det)
      If ( abs(det).eq.0.0 ) Write(6,*)                                 &
     &       'WARNING!! Determinant=0 in PolFit'
!!
!!---- Calculate fitted result.
      Do idata = 1,ndata
!!---- Calculate powers of individual variable values.
      Do ivar = 1,NrOfVar
      pow = 1.0d0
      vpow(0,ivar) = 1.0d0
      Do i = 1,mxdeg
      pow = pow*var(idata,ivar)
      vpow(i,ivar) = pow
      End Do
      End Do
!!---- Calculate value of each polynomial term. Add to yfit.
      pol = 0.0d0
      Do iterm = 1,nPolyTerm
      ip = ipow(iterm,1)
      t = vpow(ip,1)
      Do ivar = 2,NrOfVar
      ip = ipow(iterm,ivar)
      t = t*vpow(ip,ivar)
      End Do
      pol = pol+coef(iterm,1)*t
      End Do
      yfit(idata) = pol
      End Do
!!
!!---- Calculate standard deviation and maximum error.
      sum = 0.0d0
      Do idata = 1,ndata
      diff = abs(yin(idata)-yfit(idata))
      diff_vec(idata) = diff
      If ( idata.eq.1 ) Then
      max_err = diff
      Else
      If ( diff.gt.max_err ) max_err = diff
      End If
      sum = sum+diff**2
      End Do
      stand_dev = sqrt(sum/ndata)
!!
      Call GetMem('weight','Free','Real',ipweight,ndata)
!!
      End
!!
!!-----------------------------------------------------------------------!

!!-----------------------------------------------------------------------!
!!
      Subroutine factor(exponent,nder,rfactor)
!!
!!  Purpose:
!!    Calculate coefficient for n'th derivative.
!!
      Implicit none
      Integer exponent,nder
      Integer i,isum,ntot
      Real*8 rfactor
!!
      isum = 1
      ntot = nder+exponent
      If ( nder.gt.0 ) Then
      Do i = ntot,ntot-nder+1,-1
      isum = isum*i
      End Do
      Else If ( nder.lt.0 ) Then
      isum = 0
      End If
      rfactor = dble(isum)
!!
      End
!!
!!-----------------------------------------------------------------------!


!!-----------------------------------------------------------------------!
!!



!!
      Subroutine Cholesky(A,L,nd)
!!
!!  Purpose:
!!    Perform a Cholesky factorization of the matrix A:
!!
!!                       A = L~ * L
!!
!!    where L is a lower triangular matrix and L~ is L transposed.
!!
!!  Input:
!!    A      : real*8 two dimensional array - positive
!!             definite matrix.
!!
!!  Output:
!!    L      : real*8 two dimensional array - lower
!!             triangular matrix.
!!
!!  Calls:
!!
!       Implicit None
!VV: all calls use nxn
      Real*8 A(nd,nd), L(nd,nd),dd
!!
      Integer n,j
      Integer iRow,jRow
#include "WrkSpc.fh"
!!
!!---- Initialize.
      n=nd
      n2= n*n
      Call GetMem('D','Allo','Real',ipD,n)
      call dcopy_(n2,A,1,L,1)
!vv       L = A
!!
!!---- Take care of the n-1 last rows.
      If ( n.gt.1 ) Then
      Do iRow = n,2,-1
      work(ipd+iRow-1) = L(iRow,iRow)
      do j=1,n
      L(iRow,j)=L(iRow,j)/work(ipd+iRow-1)
      End do
      Do jRow = iRow-1,1,-1
      do j=1,n
      L(jRow,j)=L(jRow,j)-L(jRow,iRow)*L(iRow,j)
      End do
      End Do
      End Do
      End If
!!
!!---- Take care of the first row.
      work(ipd) = L(1,1)
      L(1,1) = 1.0d0
!!
!!---- Multiply Llow with the square root of d.
      Do iRow = 1,n
      If ( work(ipd+iRow-1).lt.0.0 ) Then
      Write(6,*)                                                        &
     &       'Error in Cholesky!!! Matrix not positive definite.'
      Call Abend
      End If
      dd=sqrt(work(ipd+iRow-1))
      do j=1,n
      L(iRow,j) = L(iRow,j)*dd
      End do
      End Do
      Call GetMem('D','Free','Real',ipD,n)
!!
!       call dcopy_(n2,L,1,Llow,1)
!       Llow = L
!!
      End
