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
      Complex*16 Function trace( n, A, B )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n,n) and B(n,n)
      Integer, intent(in)          :: n
      Complex(kind=8), intent(in) :: A(n,n), B(n,n)
      ! local variables
      Integer :: i, k

      trace=(0.0_wp,0.0_wp)
      Do i = 1, n
        Do k = 1, n
          trace = trace + A(i,k) * B(k,i)
        End Do
      End Do

      End Function trace



      Complex*16 Function trace2( n, A, B )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n,n) and B(n,n)
      Integer, intent(in)          :: n
      Complex(kind=8), intent(in) :: A(n,n), B(n,n)
      ! local variables
      Integer :: i, k

      trace2=(0.0_wp,0.0_wp)
      Do i = 1, n
        Do k = 1, n
          trace2 = trace2 + A(k,i) * B(i,k)
        End Do
      End Do

      End Function trace2


      Real*8 Function real_1_trace2( n, A )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n,n)
      Integer, intent(in)          :: n
      Real(kind=8), intent(in)    :: A(n,n)
      ! local variables
      Integer :: i

      real_1_trace2=0.0_wp
      Do i = 1, n
        real_1_trace2 = real_1_trace2 + A(i,i)/dble(n)
      End Do

      End Function real_1_trace2


      Complex*16 Function complex_1_trace2( n, A )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n,n)
      Integer, intent(in)          :: n
      Complex(kind=8), intent(in)    :: A(n,n)

      ! local variables
      Integer :: i
      Complex(kind=8) :: FACT

      complex_1_trace2=0.0_wp
      FACT=cmplx(dble(n),0.0_wp,wp)
      Do i = 1, n
        complex_1_trace2 = complex_1_trace2 + A(i,i)/FACT
      End Do

      End Function complex_1_trace2


      Complex*16 Function trace_exch( n1, n2, A, B )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n1,n1,n2,n2) and B(n1,n1,n2,n2)
      Integer, intent(in)          :: n1, n2
      Complex(kind=8), intent(in) :: A(n1,n1,n2,n2), B(n1,n1,n2,n2)
      ! local variables
      Integer :: i1, i2, k1, k2
      Complex(kind=8) :: tr

      trace_exch=(0.0_wp,0.0_wp)
      tr=(0.0_wp,0.0_wp)
      Do i1 = 1, n1
        Do k1 = 1, n1
          Do i2 = 1, n2
            Do k2 = 1, n2
              tr = tr + A(i1,k1, i2,k2) * B( k1,i1, k2,i2)
            End Do
          End Do
        End Do
      End Do
      trace_exch=tr
      End Function trace_exch


      Complex*16 Function trace_exch2( n1, n2, A, O1, O2 )
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
!     size of the square matrices A(n1,n1,n2,n2) and B(n1,n1,n2,n2)
      Integer, intent(in)          :: n1, n2
      Complex(kind=8), intent(in) :: A(n1,n1,n2,n2)
      Complex(kind=8), intent(in) :: O1(n1,n1), O2(n2,n2)
      ! local variables
      Integer :: i1, i2, k1, k2
      Complex(kind=8) :: tr

      trace_exch2=(0.0_wp,0.0_wp)
      tr=(0.0_wp,0.0_wp)
      Do i1 = 1, n1
        Do k1 = 1, n1
          Do i2 = 1, n2
            Do k2 = 1, n2
              tr = tr + A(i1,k1, i2,k2) * O1(k1,i1)*O2(k2,i2)
            End Do
          End Do
        End Do
      End Do
      trace_exch2=tr
      End Function trace_exch2
