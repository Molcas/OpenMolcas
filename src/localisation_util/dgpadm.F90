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

      subroutine DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

      use constants, only: Zero, One, Two
      use definitions, only: iwp,wp,u6
      implicit none
      integer(kind=iwp), intent(in):: ideg, m, ldh, lwsp
      integer(kind=iwp), intent(out):: iexph, ns, iflag, ipiv(m)
      real(kind=wp), intent(in) :: t, H(ldh,m)
      real(kind=wp), intent(out) ::wsp(lwsp)

!-----Purpose----------------------------------------------------------|
!
!     Computes exp(t*H), the matrix exponential of a general matrix in
!     full, using the irreducible rational Pade approximation to the
!     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
!     combined with scaling-and-squaring.
!
!-----Arguments--------------------------------------------------------|
!
!     ideg      : (input) the degre of the diagonal Pade to be used.
!                 a value of 6 is generally satisfactory.
!
!     m         : (input) order of H.
!
!     H(ldh,m)  : (input) argument matrix.
!
!     t         : (input) time-scale (can be < 0).
!
!     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
!
!     ipiv(m)   : (workspace)
!
!>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
!                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
!                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!                 NOTE: if the routine was called with wsp(iptr),
!                       then exp(tH) will start at wsp(iptr+iexph-1).
!
!     ns        : (output) number of scaling-squaring used.
!
!     iflag     : (output) exit flag.
!                      0 - no problem
!                     <0 - problem
!
!----------------------------------------------------------------------|
!     Derived from original code by:
!     Roger B. Sidje (rbs@maths.uq.edu.au)
!     EXPOKIT: Software Package for Computing Matrix Exponentials.
!     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
!----------------------------------------------------------------------|
!
      integer(kind=iwp) mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      real(kind=wp) hnorm,scale,scale2,cp,cq

!---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh<m ) iflag = -1
      if ( lwsp<4*mm+ideg+1 ) iflag = -2
      if ( iflag/=0 ) Then
         Write (u6,*)  'bad sizes (in input of DGPADM)'
         Call Abend()
      endif
!
!---  initialise pointers ...
!
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
!
!---  scaling: seek ns such that ||t*H/2^ns|| < 1/2
!     and set scale = t/2^ns ...
!
      wsp(1:m) = Zero
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         enddo
      enddo
      hnorm = Zero
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      enddo
      hnorm = ABS( t*hnorm )
      if ( hnorm==Zero ) then
         write (u6,*)  'Error - null H in input of DGPADM.'
         Call abend()
      endif
      ns = MAX( 0,INT(LOG(hnorm)/LOG(Two))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
!
!---  compute Pade coefficients ...
!
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = One
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      enddo
!
!---  H2 = scale2*H*H ...
!
      call DGEMM_( 'n','n',m,m,m,scale2,H,ldh,H,ldh,Zero,wsp(ih2),m )
!
!---  initialize p (numerator) and q (denominator) ...
!
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = Zero
            wsp(iq + (j-1)*m + i-1) = Zero
         enddo
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      enddo
!
!---  Apply Horner rule ...
!
      iodd = 1
      k = ideg - 1
      Do
        iused = iodd*iq + (1-iodd)*ip
        call DGEMM_( 'n','n',m,m,m, One,wsp(iused),m,wsp(ih2),m, Zero,wsp(ifree),m )
        do j = 1,m
           wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
        enddo
        ip = (1-iodd)*ifree + iodd*ip
        iq = iodd*ifree + (1-iodd)*iq
        ifree = iused
        iodd = 1-iodd
        k = k-1
        if (k==0) Exit
      End Do
!
!---  Obtain (+/-)(I + 2*(p\q)) ...
!
      if ( iodd == 1 ) then
         call DGEMM_( 'n','n',m,m,m, scale,wsp(iq),m,H,ldh, Zero,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM_( 'n','n',m,m,m, scale,wsp(ip),m,H,ldh, Zero,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY_( mm, -One,wsp(ip),1, wsp(iq),1 )
      call DGESV_( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag/=0 ) Then
         Write (u6,*) 'Problem in DGESV (within DGPADM)'
         Call abend()
      End If
      call DSCAL_( mm, Two, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + One
      enddo
      iput = ip
      if ( ns==0 .and. iodd==1 ) then
         call DSCAL_( mm, -One, wsp(ip), 1 )
      else
!
!--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
!
         iodd = 1
         do k = 1,ns
            iget = iodd*ip + (1-iodd)*iq
            iput = (1-iodd)*ip + iodd*iq
            call DGEMM_( 'n','n',m,m,m,One,wsp(iget),m, wsp(iget),m,Zero,wsp(iput),m )
            iodd = 1-iodd
         enddo
      end if
      iexph = iput

      END
