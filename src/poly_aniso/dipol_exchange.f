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
      Subroutine Dipol_Exchange( N1, N2, vec, dist, M1, M2,  HDIP )
c this Subroutine computes the dipolar coupling between the two moments
      Implicit None
      Integer, parameter            :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)           :: N1, N2
      Real(kind=wp), intent(in)     :: vec(3), dist
      Complex(kind=wp), intent(in)  :: M1(3,N1,N1)
      Complex(kind=wp), intent(in)  :: M2(3,N2,N2)
      Complex(kind=wp), intent(out) :: HDIP(N1,N1,N2,N2)
c local variables
      Integer          :: m,i1,j1,i2,j2
      Real(kind=wp)    :: MB2
      Complex(kind=wp) :: p2a, p2b, p1, HL, d3, mb2c, threeC, vec1(3)

      Call qEnter('PA_Dipol_Exchange')

!     * in cm-1*T-1   -- the value of (mu_Bohr*mu_Bohr)/(Angstrom^3)  in cm-1
      MB2=0.4329701512063995_wp
c
      If( (N1<=0).OR.(N2<=0) ) Return
      Call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,HDIP,1)
      If(dist.eq.0.0_wp) Then
        Write(6,'(A)') 'DIPOL_EXCHANGE::  dist = 0'
        Write(6,'(A)') 'this is not normal. Stop.'
        Return
      End If
      ! switch to Complex arithmetic:
      ! kind=8, complex double precision
      d3    =(0.0_wp,0.0_wp)
      mb2c  =(0.0_wp,0.0_wp)
      threeC=(3.0_wp,0.0_wp)
      d3    =cmplx(dist*dist*dist,0.0_wp, wp)
      mb2c  =cmplx(MB2,0.0_wp, wp)
      Do m=1,3
        vec1(m)=cmplx(vec(m),0.0_wp, wp)
      End Do
      ! calculate the dipolar coupling:
      Do i1=1,N1
        Do j1=1,N1
          Do i2=1,N2
            Do j2=1,N2
              p2a=(0.0_wp,0.0_wp)
              p2b=(0.0_wp,0.0_wp)
              p1 =(0.0_wp,0.0_wp)
              HL =(0.0_wp,0.0_wp)
              Do m=1,3
                p2a=p2a+M1(m,i1,j1)*vec1(m)
                p2b=p2b+M2(m,i2,j2)*vec1(m)
                p1 =p1 + M1(m,i1,j1) * M2(m,i2,j2)
              End Do
              HL = p1 - threeC * p2a * p2b
              HDIP(i1,j1,i2,j2)=mb2c * HL / d3
            End Do
          End Do
        End Do
      End Do

      Call qExit('PA_Dipol_Exchange')

      Return
      End Subroutine Dipol_Exchange
