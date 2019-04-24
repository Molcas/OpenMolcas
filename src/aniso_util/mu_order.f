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
      Subroutine mu_order( dim, MS, MM, gtens, order, HCF2, AMM,
     &                     AMS, Z, iprint)

C This Subroutine receives the moment matrix dipso(3,dim,dim) and Returns the matrix re-builted using only the 1-st order operators.

      Implicit None
      Integer, parameter       :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)      :: dim, order, iprint
      Real(kind=wp), intent(out) :: gtens(3)
!     initial magnetic moment
      Complex(kind=wp), intent(in)  ::  MM(3,dim,dim)
!     initial spin moment
      Complex(kind=wp), intent(in)  ::  MS(3,dim,dim)
!     transformed magnetic moment
      Complex(kind=wp), intent(out) :: AMM(3,dim,dim)
!     transformed spin moment
      Complex(kind=wp), intent(out) :: AMS(3,dim,dim)
      Complex(kind=wp), intent(out) :: Z(dim,dim)
      Complex(kind=wp), intent(out) :: HCF2(dim,3,dim,dim)

      ! local variables:
      Integer             :: i,j,l,i1,i2,m,n
      Real(kind=wp)       :: maxes(3,3)
      Complex(kind=wp)    :: DIP_O(dim,dim), DIP_W(dim,dim),
     &                       B(3,dim,-dim:dim), BNMC(3,dim,0:dim),
     &                       BNMS(3,dim,0:dim), SP_MOW, SP_DIPO(3),
     &                       SP_DIPW(3), O1, O2, m_fact, trace
      External            :: trace
!------------------------------------------------------------
      Call qEnter('mu_order')
      Z=(0.0_wp,0.0_wp)
      ! get the local pseudospin:
      Call pseudospin( MM, dim, Z, 3,1, iprint)
      ! re-write MM and MS to the new pseudospin basis:
      AMS=(0.0_wp,0.0_wp)
      AMM=(0.0_wp,0.0_wp)
      Do l=1,3
        Do i=1,dim
          Do j=1,dim
            Do i1=1,dim
              Do i2=1,dim
         AMM(l,i,j) = AMM(l,i,j) + MM(l,i1,i2) * conjg(Z(i1,i))*Z(i2,j)
         AMS(l,i,j) = AMS(l,i,j) + MS(l,i1,i2) * conjg(Z(i1,i))*Z(i2,j)
              End Do
            End Do
          End Do
        End Do
      End Do
      If(iprint.gt.2) Then
        Call prMom('MU_ORDER:   AMM(l,i,j):',AMM,dim)
        Call prMom('MU_ORDER:   AMS(l,i,j):',AMS,dim)
      End If
      !project the moment on ITO:
      !obtain the b3m and c3m coefficients:
      B(1:3,1:dim,-dim:dim)=(0.0_wp,0.0_wp)
      Do N=1,dim-1
        Do M=0,N
          DIP_O = (0.0_wp,0.0_wp)
          DIP_W = (0.0_wp,0.0_wp)
          Call Stewens_matrixel(N,M,dim,DIP_O,DIP_W,iprint)
          If(iprint.gt.5) Then
            Write(6,*)
            Write(6,'( 5x,a,i2,a,i3)') 'DIP_O, n = ',N,', m =',m
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F10.6,2x))') (DIP_O(i,j), j=1,dim)
            End Do
            Write(6,*)
            Write(6,'( 5x,a,i2,a,i3)')  'DIP_W, n = ',N,', m =',m
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F10.6,2x))') (DIP_W(i,j), j=1,dim)
            End Do
          End If

          SP_DIPO = (0.0_wp,0.0_wp)
          SP_DIPW = (0.0_wp,0.0_wp)
          SP_MOW  = (0.0_wp,0.0_wp)
          SP_MOW  = trace(dim,DIP_O,DIP_W)
          Do l=1,3
            SP_DIPO(l)=trace(dim,AMS(l,1:dim,1:dim),DIP_O(1:dim,1:dim))
            SP_DIPW(l)=trace(dim,AMS(l,1:dim,1:dim),DIP_W(1:dim,1:dim))

            B(l,n,-m)=SP_DIPO(l)/SP_MOW
            B(l,n, m)=SP_DIPW(l)/SP_MOW
          End Do ! l
        End Do !m
      End Do !n

      BNMC=(0.0_wp,0.0_wp)
      BNMS=(0.0_wp,0.0_wp)
      Do n=1,dim-1
        Do m=0,N
          Do l=1,3
            If(M.eq.0) Then
              BNMC(l,n,m)=(0.5_wp,0.0_wp)*(B(l,n,m)+B(l,n,-m))
            Else
              m_fact=cmplx((-1)**M,0,wp)
              BNMC(l,n,m)=  B(l,n,m)+m_fact*B(l,n,-m)
              BNMS(l,n,m)=( B(l,n,m)-m_fact*B(l,n,-m) )*(0.0_wp,-1.0_wp)
            End If
          End Do
        End Do
      End Do !n

      If(iprint.gt.2) Then
        Write(6,'(100A)') ('-',i=1,47),'|'
        Write(6,'(A)') '  n  |  m  |   |       B       |'//
     &                                 '       C       |'
        Do N=1,dim-1
          Write(6,'(A)') '-----|-----|---|---------------|'//
     &                                   '---------------|'
           Do M=0,N
             If (M.ne.0) Write(6,'(A)')
     &           '     |-----|---|---------------|---------------|'
             Write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))')
     &             N,'|',M,'|','X','|',DBLE(BNMC(1,N,M)),'|',
     &                                 DBLE(BNMS(1,N,M)),'|'
             Write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))')
     &             N,'|',M,'|','Y','|',DBLE(BNMC(2,N,M)),'|',
     &                                 DBLE(BNMS(2,N,M)),'|'
             Write(6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))')
     &             N,'|',M,'|','Z','|',DBLE(BNMC(3,N,M)),'|',
     &                                 DBLE(BNMS(3,N,M)),'|'
           End Do
        End Do
        Write(6,'(100A)') ('-',i=1,47),'|'
      End If

      HCF2=(0.0_wp,0.0_wp)
      Do N=1,dim-1
        Do M=0,N
          DIP_O = (0.0_wp,0.0_wp)
          DIP_W = (0.0_wp,0.0_wp)

          Call Stewens_matrixel( N, M, dim, DIP_O, DIP_W, iprint )

          Do l=1,3
            Do i=1,dim
              Do j=1,dim
                If(M.eq.0) Then
                  HCF2(N,l,i,j) = HCF2(N,l,i,j) + BNMC(l,N,M)*DIP_O(i,j)
                Else
                m_fact=(0.0_wp, 0.0_wp)
                    O1=(0.0_wp, 0.0_wp)
                    O2=(0.0_wp, 0.0_wp)

                  m_fact=cmplx((-1)**M,0,wp)
                    O1=(0.5_wp, 0.0_wp)*( m_fact*DIP_W(i,j)+DIP_O(i,j) )
                    O2=(0.0_wp,-0.5_wp)*( m_fact*DIP_W(i,j)-DIP_O(i,j) )

                  HCF2(N,l,i,j) = HCF2(N,l,i,j) + BNMC(l,N,M) * O1
     &                                          + BNMS(l,N,M) * O2
                End If
              End Do
            End Do
          End Do
        End Do
      End Do !n

      If(iprint.gt.2) Then
        Do N=1,dim-1,2
          Write(6,*)
          Write(6,'( 5x,a,I2,a)')  'HCF2(',N,',l,i,j)'
          Do l=1,3
            Write(6,*)
            Write(6,'(a,i3)') 'PROJECTION =' , l
            Do i=1,dim
              Write(6,'(20(2F12.8,2x))') (HCF2(N,l,i,j), j=1,dim)
            End Do
          End Do
          Write(6,*)
        End Do
      End If

      gtens=0.0_wp
      maxes=0.0_wp
      Call ATENS( HCF2(order,:,:,:), dim, gtens, maxes, 1)
      Call qExit('mu_order')

      Return
      End
