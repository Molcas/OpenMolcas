      Subroutine cZeroMatrix(a,n)
* fills all elements of a square Complex matrix of size n by Complex zero
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          ::  n
      Complex(kind=wp)             :: a(n,n)
      !local
      Integer :: i, j
      Do j=1,n
         Do i=1,n
         a(i,j)=(0.0_wp,0.0_wp)
         End Do
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine rZeroMatrix(a,n)
* fills all elements of a square Real matrix of size n by Complex zero
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Real(kind=wp)                :: a(n,n)
      !local
      Integer :: i, j
      Do j=1,n
         Do i=1,n
            a(i,j)=0.0_wp
         End Do
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine cZeroVector(a,n)
* fills all elements of a Complex vector of size n by Complex zero
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer                      :: n
      Complex(kind=wp)             :: a(n)
      !local
      Integer :: i
      Do i=1,n
         a(i)=(0.0_wp,0.0_wp)
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine rZeroVector(a,n)
* fills all elements of a Complex vector of size n by Complex zero
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Real(kind=wp)                :: a(n)
      !local
      Integer :: i
      Do i=1,n
         a(i)=0.0_wp
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine cZeroMoment(a,n)
* fills all elements of a Complex matrix of size a(3,n,n)  by Complex zero
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Complex(kind=wp)             :: a(3,n,n)
      !local
      Integer :: i, j, l
      Do i=1,n
        Do j=1,n
          Do l=1,3
          a(l,j,i)=(0.0_wp,0.0_wp)
          End Do
        End Do
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine prMom(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Complex(kind=wp)             :: M(3,n,n)
      !local
      Integer           :: i, j, l
      Character(len=1)  :: proj(3)
      Character(len=*)  :: a
      Character(len=50) :: fmtline
      proj(1)='X'
      proj(2)='Y'
      proj(3)='Z'
      Write(6,*)
      Write(6,'(2a)') 'print: ',a
      Write(fmtline,'(a,i2,a)') '(',n,'(2f9.4,1x))'
      Do l=1,3
      Write(6,'(2a)') 'projection: ',proj(l)
          Do i=1,n
             Write(6,fmtline) (M(l,i,j),j=1,n)
          End Do
      Write(6,*)
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine prMom_herm(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Complex(kind=wp)             :: M(3,n,n)
      !local
      Integer           :: i, j, l
      Real(kind=wp)     :: R
      Character(len=*)  :: a
      Write(6,*)
      Write(6,'(2a)') 'print: ',a
      Do i=1,n
          Do j=1,i
          R=0.0_wp
          R=( ABS(M(1,i,j)) + ABS(M(2,i,j)) + ABS(M(3,i,j)) ) / 3.0_wp
             Write(6,'(A,2I3,A,3(2F16.7,2x), 2F20.7)')
     &                 'i j: ',i,j,' <i|O|j>=',(M(l,i,j),l=1,3), R
          End Do
      Write(6,*)
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine pa_prMat(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Complex(kind=wp)             :: M(n,n)
      !local
      Integer        :: i, j
      Character*(*) a
      Character*(50) fmtline
      Write(6,*)
      Write(6,'(2a)') 'print: ',a
      Write(fmtline,'(a,i2,a)') '(',n,'(2f12.4,1x))'
      Do i=1,n
         Write(6,fmtline) (M(i,j),j=1,n)
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine pa_prMatR(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Real(kind=wp)                :: M(n,n)
      !local
      Integer        :: i, j
      Character*(*) a
      Character*(50) fmtline
      Write(6,*)
      Write(6,'(2a)') 'print: ',a
      Write(fmtline,'(a,i2,a)') '(',n,'(f19.14,1x))'
      Do i=1,n
         Write(6,fmtline) (M(i,j),j=1,n)
      End Do
      Return
      End
c------------------------------------------------------------------------
      Subroutine print_ZFS(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer                      :: i,j,n,k,jEnd
      Complex(kind=wp), intent(in) :: M(n,n)
      Character*(*)                :: a


      Write(6,'(/)')
      Write(6,'(100A)') ('-',i=1,87)
      Write(6,'(A)') a
      Do j=1,n,4
         jEnd=MIN(n,J+3)
         If (MOD(n,2)==0) Then
         ! code for odd N
            Write(6,'(150A)') ('-',i=1,10),(('-',i=1,24),k=j,jEnd),'|'
            Write(6,'(10x,A,50(8x,A,I3,A,7x,A))') '|',
     &                         ('|',2*i-n-1,'/2 >','|',i=j,jEnd)
            Write(6,'(150A)')  ('-',i=1,10),'|',
     &                         ('---- Real ----- Imag --|',k=j,jEnd)
            ! print the matrix
            Do i=1,n
               Write(6,'(1x,A,I3,A,1x,A,50(2F11.5,1x,A))')
     &                                    '<',2*i-n-1,'/2','| |',
     &                                    (M(k,i),'|' ,k=j,jEnd)
            End Do
            Write(6,'(150A)') ('-',i=1,10),(('-',i=1,24),k=j,jEnd),'|'

         Else
         ! code for odd N

            Write(6,'(150A)') ('-',i=1,8),(('-',i=1,24),k=j,jEnd),'|'
            Write(6,'(8x,A,50(8x,A,I3,A,9x,A))') '|',
     &                         ('|',-(n-1)/2-1+i,' >','|',i=j,jEnd)
            Write(6,'(150A)')  ('-',i=1,8),'|',
     &                         ('---- Real ----- Imag --|',k=j,jEnd)
            Do i=1,n
               Write(6,'(1x,A,I3,1x,A,50(2F11.5,1x,A))')
     &                                '<',-(n-1)/2-1+i,'| |',
     &                                  (M(k,i),'|' ,k=j,jEnd)
            End Do
            Write(6,'(150A)') ('-',i=1,8),(('-',i=1,24),k=j,jEnd),'|'
         End If
      End Do !j

      Return
      End
c------------------------------------------------------------------------
      Subroutine print_ZFS_eigenvectors(a,m,n)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer                      :: i,j,n,k,jEnd
      Complex(kind=wp), intent(in) :: M(n,n) ! eigenvectors
      Character*(1)                :: a


      Write(6,'(/)')
      Do j=1,n,2
        jEnd=MIN(n,j+1) ! '|  > |'
        Write(6,'(A,6A)') '--------|',
     &                ('-----------------------------|',i=j,jEnd)
        Write(6,'(3A,6(6x,a,i3,5x,a))') ' | ',a,'M > |',
     &                         ('ab initio state',i,'|',i=j,jEnd)
        Write(6,'(A,6A)') '--------|',
     &                ('-- Real ---- Imag --|-Weight-|',i=j,jEnd)

        Do i=1,n
          If(MOD(n,2)==1) Then
            Write(6,'(1x,A,1x,i2,A,   6(2(E22.14,1x),a,F6.1,1x,a))')
     &                '|',-(n-1)/2-(1-i),' > |',
     &         (DBLE(M(i,k)), AIMAG(M(i,k)), '|',
     &         100.0_wp*(DBLE(M(i,k))**2 + AIMAG(M(i,k))**2),
     &         '%|',k=j,jEnd)
          Else
            Write(6,'(A,i3,a,a,       6(2(E22.14,1x),a,F6.1,1x,a))')
     &                '|',-(n-1)+2*(i-1),'/2> ','|',
     &         (DBLE(M(i,k)), AIMAG(M(i,k)), '|',
     &         100.0_wp*(DBLE(M(i,k))**2 + AIMAG(M(i,k))**2),
     &         '%|',k=j,jEnd)
          End If
        End Do  !i


      Write(6,'(A,6A)') '--------|',
     &                ('-----------------------------|',i=j,jEnd)
      End Do ! j


      Return
      End
c------------------------------------------------------------------------
      Subroutine print_ZFS_naoya(A,M,N)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: N ! dimension of the pseudospin
      Complex(kind=wp), intent(in) :: M(n,n) ! complex parameters to print
      Character*(1)                :: a
      ! local variables:
      Integer       :: k,i,j,jEnd
      Real(kind=wp) :: Mr(n), Mi(n), Weight(n)
      Character(1)  :: cRsign(n), cIsign(n)

      Write(6,'(/)')
      Do j=1,n,2
        jEnd=MIN(n,j+1) ! '|  > |'
        If(j==1) Write(6,'(150A)') '--------|',
     &                                     (('-',k=1,58),'|',i=j,jEnd)
        Write(6,'(3A,6(16x,a,i3,24x,a))') ' | ',a,'M > |',
     &                              ('ab initio state',i,'|',i=j,jEnd)
        Write(6,'(A,6A)') '--------|',
     & ('-------  Real  -------|------  Imaginary  -------|-Weight-|',
     &                                                       i=j,jEnd)

        Do i=1,n
            ! fix the sign:
            Do k=j,jEnd
               Mr(k)=0.0_wp
               Mi(k)=0.0_wp
               Weight(k)=0.0_wp

               Mr(k)=DBLE( M(i,k))
               Mi(k)=AIMAG(M(i,k))
               Weight(k)=100.0_wp*( Mr(k)*Mr(k) + Mi(k)*Mi(k) )

               If (Mr(k) >= 0.0_wp) Then
                  cRsign(k)='+'
               Else
                  cRsign(k)='-'
               End If
               If (Mi(k) >= 0.0_wp) Then
                  cIsign(k)='+'
               Else
                  cIsign(k)='-'
               End If
            End Do

          ! print it
          If(MOD(n,2)==1) Then
            Write(6,'(1x,A,1x,i2,A, 2(2(1x,A,E20.14,1x),a,F6.1,1x,a))')
     &                '|',-(n-1)/2-(1-i),' > |',
     &         (cRsign(k),ABS(Mr(k)), cIsign(k),ABS(Mi(k)), '*I |',
     &                Weight(k),'%|',k=j,jEnd)
          Else
            Write(6,'(A,i3,a,a,     2(2(1x,A,E20.14,1x),a,F6.1,1x,a))')
     &                '|',-(n-1)+2*(i-1),'/2> ','|',
     &         (cRsign(k),ABS(Mr(k)), cIsign(k),ABS(Mi(k)), '*I |',
     &                Weight(k),'%|',k=j,jEnd)
          End If
        End Do  !i



        Write(6,'(150A)') '--------|',(('-',k=1,58),'|',i=j,jEnd)
      End Do ! j

      Return
      End Subroutine print_ZFS_naoya






      Subroutine print_CFP_alpha(nlanth,n,B,C)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n        ! dimension of the pseudospin
      Integer, intent(in)          :: nlanth   ! number of the lanthanide
      Real(kind=wp), intent(in)    :: B(n,0:n), C(n,0:n) ! real and imaginary CF paraemters
      !Logical, intent(in), optional:: print_all
      ! local variables:
      Integer       :: k,q,i
      Real(kind=wp) :: a(6)


      a=0.0_wp
      Call set_an(nlanth,a)

!          O(m1,m2)=(0.5_wp,0.0_wp) * ( Om + mQ*Op )
!          W(m1,m2)=(0.0_wp,0.5_wp) * ( Om - mQ*Op )
      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') 'The Crystal-Field Hamiltonian:'
      Write(6,'(A)') '   Hcf = SUM_{k,q} alpha(k) * [ B(k,q) '//
     &               '* O(k,q) +  C(k,q) * W(k,q) ];'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(k,q) =  0.5 * ( (-1)**q * Y(k,+q)'//
     &               ' + Y(k,-q) );'
      Write(6,'(A)') '   W(k,q) = -0.5 * ( (-1)**q * Y(k,+q)'//
     &               ' - Y(k,-q) ) * I;   (I = imaginary unit)'
      Write(6,'(A)') '   k - the rank of the ITO, = 2, 4, 6;'
      Write(6,'(A)') '   q - the component of the ITO, = 0, 1, 2, '//
     &               '... k;'
      Write(6,'(A)') '   alpha(k) - Stevens coefficients;'
      Write(6,'(A)') 'These operators have been defined in: '
      Write(6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., '//
     &               '137, 064112 (2012).'

      Write(6,'(100A)') ('-',i=1,76),'|'
      Write(6,'(A)') '  k  |  q  |    1/alpha(k)  |'//
     &               '         B(k,q)        |         C(k,q)        |'
      Do k=2,6,2
        If (A(k).ne.0.0_wp) Then
        Write(6,'(A)') '-----|-----|----------------|'//
     &                 '-----------------------|'//
     &                 '-----------------------|'
        Do q=0,k
          If ( q==k/2) Then
            Write(6,'(2(2x,I1,2x,A),F14.5,2x,A,2(E22.14,1x,A))')
     &                k,'|',q,'|',1.0_wp/a(k),'|',
     &                B(k,q)/a(k),'|',  C(k,q)/a(k),'|'
          Else
            Write(6,'(2(2x,I1,2x,A),16x,A,2(E22.14,1x,A))')
     &                k,'|',q,'|','|',
     &                B(k,q)/a(k),'|',  C(k,q)/a(k),'|'
          End If
        End Do
        End If
      End Do
      Write(6,'(100A)') ('-',i=1,76),'|'

      Return
      End Subroutine print_CFP_alpha







      Subroutine print_CFP_LCLU(n,B,C,print_all)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n                  ! dimension of the pseudospin
      Real(kind=wp), intent(in)    :: B(n,0:n), C(n,0:n) ! real and imaginary CF paraemters
      Logical, intent(in), optional:: print_all
      ! local variables:
      Integer       :: k,q,i

!          O(m1,m2)=(0.5_wp,0.0_wp) * ( Om + mQ*Op )
!          W(m1,m2)=(0.0_wp,0.5_wp) * ( Om - mQ*Op )
      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') 'The Crystal-Field Hamiltonian:'
      Write(6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) +  '//
     &                                      'C(k,q) * W(k,q) ];'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(k,q) =  0.5 * ( (-1)**q * Y(k,+q)'//
     &               ' + Y(k,-q) );'
      Write(6,'(A)') '   W(k,q) = -0.5 * ( (-1)**q * Y(k,+q)'//
     &               ' - Y(k,-q) ) * I;   (I = imaginary unit)'
      Write(6,'(A)') '   k - the rank of the ITO, = 2, 4, 6;'
      Write(6,'(A)') '   q - the component of the ITO, = 0, 1, 2, '//
     &               '... k;'
      Write(6,'(A)') 'These operators have been defined in: '
      Write(6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., '//
     &               '137, 064112 (2012).'

      Write(6,'(100A)') ('-',i=1,59),'|'
      Write(6,'(A)') '  k  |  q  |         B(k,q)        |'//
     &                           '         C(k,q)        |'
      If(print_all) Then
         Do k=2,n-1
           Write(6,'(A)') '-----|-----|-----------------------|'//
     &                              '-----------------------|'
           Do q=0,k
             Write(6,'(2(1x,I2,2x,A),2(E22.14,1x,A))')
     &                k,'|',q,'|',B(k,q),'|',C(k,q),'|'
           End Do
         End Do
      Else
         Do k=2,n-1,2
           Write(6,'(A)') '-----|-----|-----------------------|'//
     &                              '-----------------------|'
           Do q=0,k
             Write(6,'(2(1x,I2,2x,A),2(E22.14,1x,A))')
     &                k,'|',q,'|',B(k,q),'|',C(k,q),'|'
           End Do
         End Do
      End If
      Write(6,'(100A)') ('-',i=1,59),'|'
      Return
      End Subroutine print_CFP_LCLU



      Subroutine print_CFP_stev(n,B,print_all)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n                  ! dimension of the pseudospin
      Real(kind=wp), intent(in)    :: B(n,-n:n)          ! real and imaginary CF paraemters
      Logical, intent(in), optional:: print_all
      ! local variables:
      Integer       :: k,q,i,kmax,iq
      Real(kind=wp) :: knm(12,0:12), f

      Call set_knm( knm )
      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') 'The Crystal-Field Hamiltonian:'
      Write(6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO)'//
     &               'as defined in:'
      Write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State '//
     &                   'Phys.,18(1985) 1415-1430.'
      Write(6,'(10x,A)') '2. Implemented in the "EasySpin" function '//
     &                   'in MATLAB, www.easyspin.org.'
      Write(6,'(A)') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
      Write(6,'(A)') '   q - the component of the ITO, = -k, -k+1, '//
     &               '... 0, 1, ... k;'
      If((n-1)>12) Then
        Write(6,'(A)') 'k = 12 may not be the highest rank of '//
     &                 'the ITO for this case, but it '
        Write(6,'(A)') 'is the maximal k implemented in the'//
     &                 ' "EasySpin" function in MATLAB.'
      End If
      Write(6,'(A)') 'Knm are proportionality coefficients between'//
     &               ' the ESO and operators defined in '
      Write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
      Write(6,'(100A)') ('-',i=1,48),'|'
      Write(6,'(A)') '  k |  q  |    (K)^2    |'//
     &               '         B(k,q)        |'
      If((n-1)>12) Then
        kmax=12
      Else
        kmax=n-1
      End If

      If(print_all) Then
         Do k=2,kmax
           Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'
           Do q=-k,k
             iq=abs(q)
             f=knm(k,iq)*knm(k,iq)
               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,'//
     &                   '2(E22.14,1x,A))')
     &                   k,'|',q,'|',f,'|',B(k,q),'|'
           End Do
         End Do
      Else
         Do k=2,kmax,2
           Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'
           Do q=-k,k
             iq=abs(q)
             f=knm(k,iq)*knm(k,iq)
               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,'//
     &                   '2(E22.14,1x,A))')
     &                   k,'|',q,'|',f,'|',B(k,q),'|'
           End Do
         End Do
      End If
      Write(6,'(90A)') ('-',i=1,48),'|'
      Return
      End Subroutine print_CFP_stev







      Subroutine print_CFP_naoya(N,A,print_all)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: N ! dimension of the pseudospin
      Complex(kind=wp), intent(in) :: A( (n-1), -(n-1):(n-1) ) ! complex parameters to print
      Logical, intent(in),optional :: print_all
      ! local variables:
      Integer       :: k,q,i
      Real(kind=wp) :: Ar, Ai
      Character(1)  :: cRsign, cIsign

      Call qEnter('SA_PRCF')

      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') 'The Crystal-Field Hamiltonian:'
      Write(6,'(A)')
      Write(6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
      Write(6,'(A)')
      Write(6,'(A)') ' where:                                  '
      Write(6,'(A)') '   O(k,q) =  Irreducible Tensor Operators'
      Write(6,'(A)') '             defined as follows:         '
      Write(6,'(A)')
      Write(6,'(A)') '          Y(k,q)             CG(J,M2,k,q,J,M1)'
      Write(6,'(A)') ' < J,M1 | ------ | J,M2 >  = -----------------'
      Write(6,'(A)') '          Y(k,0)              CG(J,J,k,0,J,J) '
      Write(6,'(A)')
      Write(6,'(A)') '  CG - Clebsh-Gordan Coefficient:'
      Write(6,'(A)') '                            c,gm     '
      Write(6,'(A)') '      CG(a,al,b,bt,c,gm) = C         '
      Write(6,'(A)') '                            a,al,b,bt'
      Write(6,'(A)')
      Write(6,'(A)') '   k - the rank of the ITO, = 2, 4, 6, 8, 10, 12.'
      Write(6,'(A)') '   q - the component of the ITO, = -k, -k+1, '//
     &               '... 0, 1, ... k;'
      ! thee table:
      Write(6,'(100A)') ('-',i=1,59),'|'
      Write(6,'(A,11x,A,12x,A)') '  k |  q  |',
     &           'Complex parameter  A(k,q)','|'
      Write(6,'(A)') '----|-----|'//
     &           '--------  Real  ------|-----  Imaginary  -------|'

      If(print_all) Then
         Do k=2,N-1
            Do q=-k,k
               Ar=DBLE( A(k,q))
               Ai=AIMAG(A(k,q))
               If (Ar.ge.0.0_wp) Then
                  cRsign='+'
               Else
                  cRsign='-'
               End If
               If (Ai.ge.0.0_wp) Then
                  cIsign='+'
               Else
                  cIsign='-'
               End If

               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),'//
     &                   'A,ES20.14,1x,A,ES20.14,A)')
     &                   k,'|',q,'| ', cRsign,ABS(Ar), cIsign,ABS(Ai),
     &                   ' *I |'

            End Do
            If(k.ne.(N-1-mod(N-1,2)))
     &      Write(6,'(A)') '----|-----|'//
     &              '----------------------|-------------------------|'
         End Do
      Else
         Do k=2,N-1,2
            Do q=-k,k
               Ar=DBLE( A(k,q))
               Ai=AIMAG(A(k,q))
               If (Ar.ge.0.0_wp) Then
                  cRsign='+'
               Else
                  cRsign='-'
               End If
               If (Ai.ge.0.0_wp) Then
                  cIsign='+'
               Else
                  cIsign='-'
               End If

               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),'//
     &                   'A,ES20.14,1x,A,ES20.14,A)')
     &                   k,'|',q,'| ', cRsign,ABS(Ar), cIsign,ABS(Ai),
     &                   ' *I |'

            End Do
            If(k.ne.(N-1-mod(N-1,2)))
     &      Write(6,'(A)') '----|-----|'//
     &              '----------------------|-------------------------|'
         End Do
      End If

      Write(6,'(100A)') ('-',i=1,59),'|'
      Call qExit('SA_PRCF')
      Return
      End Subroutine print_cfp_naoya





      Subroutine print_MOM_ITO_stev(n,B,print_all)
      Implicit None
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n                  ! dimension of the pseudospin
      Real(kind=wp), intent(in)    :: B(3,n,-n:n)        ! real and imaginary CF parameters, for x,y,z
      Logical, intent(in), optional:: print_all
      ! local variables:
      Integer       :: k,q,i,kmax,iq
      Real(kind=wp) :: knm(12,0:12), f

      Call set_knm( knm )
      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') 'The magnetic moment is decomposed in Stev ITO:'
      Write(6,'(A)') '   Hcf = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO)'//
     &               'as defined in:'
      Write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State '//
     &                   'Phys.,18(1985) 1415-1430.'
      Write(6,'(10x,A)') '2. Implemented in the "EasySpin" function '//
     &                   'in MATLAB, www.easyspin.org.'
      Write(6,'(A)') '   k - the rank of the ITO, = 1, 3, 5, 7, 9, 11.'
      Write(6,'(A)') '   q - the component of the ITO, = -k, -k+1, '//
     &               '... 0, 1, ... k;'
      If((n-1)>12) Then
        Write(6,'(A)') 'k = 12 may not be the highest rank of '//
     &                 'the ITO for this case, but it '
        Write(6,'(A)') 'is the maximal k implemented in the'//
     &                 ' "EasySpin" function in MATLAB.'
      End If
      Write(6,'(A)') 'Knm are proportionality coefficients between'//
     &               ' the ESO and operators defined in '
      Write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
      Write(6,'(100A)') ('-',i=1,96),'|'
      Write(6,'(A)') '  k |  q  |    (K)^2    |'//
     &               '        B(k,q) - X     |'//
     &               '        B(k,q) - Y     |'//
     &               '        B(k,q) - Z     |'
      If((n-1)>12) Then
        kmax=12
      Else
        kmax=n-1
      End If

      If(print_all) Then

         Do k=1,kmax
           Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'//
     &                    '-----------------------|'//
     &                    '-----------------------|'
           Do q=-k,k
             iq=abs(q)
             f=knm(k,iq)*knm(k,iq)
               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,'//
     &                   '3(E22.14,1x,A))')
     &                   k,'|',q,'|',f,'|',
     &                   B(1,k,q),'|',B(2,k,q),'|',B(3,k,q),'|'
           End Do
         End Do

      Else

         Do k=1,kmax,2
           Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'//
     &                    '-----------------------|'//
     &                    '-----------------------|'
           Do q=-k,k
             iq=abs(q)
             f=knm(k,iq)*knm(k,iq)
               Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),F11.2,2x,A,'//
     &                   '3(E22.14,1x,A))')
     &                   k,'|',q,'|',f,'|',
     &                   B(1,k,q),'|',B(2,k,q),'|',B(3,k,q),'|'
           End Do
         End Do

      End If
      Write(6,'(100A)') ('-',i=1,96),'|'
      Return
      End Subroutine print_MOM_ITO_stev
