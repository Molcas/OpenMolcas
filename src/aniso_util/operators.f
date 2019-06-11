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
      Subroutine ITO(N,k,q,C0,Cp,Cm)
! This soubroutine implements the core ITO generation
! using the following ITO convention:
!                               S,M2
!                             CG
!                               S,M1,k,q
!  < S,M1 | O_k_q | S,M2 > = -------------
!                               S,S
!                             CG
!                               S,S,k,0
!
!  The denominator is computed using a simpler formula:
!  Varshalovich, p252, formula (42):
!
!  N = 2S+1, dimension   (input)
!  k = rank of ITO       (input)
!  q = projection of ITO (input)
!  Cp = O+ operator (output), complex
!  Cm = O- operator (output), complex
!  C0 = CG0  (output), real number, positive

      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n,k,q
      Real(kind=wp), intent(out)   :: C0
      Complex(kind=wp), intent(out):: Cp(n,n), Cm(n,n)
      ! local
      Integer       :: m1, m2
      Real(kind=wp) :: rm1, rm2, rS, rK, rQ, CGp, CGm, fct
      External      :: fct

      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

      rS=dble(n-1)/2.0_wp
      rK=dble(k)
      rQ=dble(q)
      C0=fct(n-1)*sqrt(dble(n)/(fct(n-k-1)*fct(n+k)))
      ! M1 and M2 go from max value to min value
      Do m1=1,n
       Do m2=1,n
        rm1=rS-dble(m1-1)
        rm2=rS-dble(m2-1)
        Call Clebsh_Gordan( rS,rm2,  rK, rQ,  rS,rm1,  CGp )
        Call Clebsh_Gordan( rS,rm2,  rK,-rQ,  rS,rm1,  CGm )
        Cp(m1,m2)=cmplx( CGp/C0, 0.0_wp, wp )
        Cm(m1,m2)=cmplx( CGm/C0, 0.0_wp, wp )
       End Do
      End Do

      Return
      End subroutine ITO


      Subroutine Liviu_ITO(n,k,q,O,W,redME)
! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators
      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n,k,q
      Complex(kind=wp), intent(out):: O(n,n), W(n,n), redME
      ! local
      Real(kind=wp) :: CR, C0

      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
      CR=0.0_wp
      C0=0.0_wp
      redME=(0.0_wp,0.0_wp)
      Call coeff_redus_sub(n,k,CR)
      Call ITO(n,k,q,C0,O,W)
      redME=cmplx(C0*CR,0.0_wp,wp)
      Call zscal_(n*n,redME,W,1)
      Call zscal_(n*n,redME,O,1)
      Return
      End subroutine Liviu_ITO



      Subroutine Stev_ITO(n,k,q,O,W,redME)
! generate O, W as in previous Stewens_matrixel function:
! redME =  ratio between Naoya's and Liviu operators
! scaled as for Stevens Operators by knm
      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n,k,q
      Complex(kind=wp), intent(out):: O(n,n), W(n,n), redME
      ! local
      Real(kind=wp) :: F, knm(12,0:12)

      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
      redME=(0.0_wp,0.0_wp)
      If((k>12).or.(q>12)) Return

      Call set_knm( knm )
      Call Liviu_ITO(n,k,q,O,W,redME)
      F=1.0_wp/knm(k,q)
      redME=cmplx(F,0.0_wp,wp)
      Call zscal_(n*n,redME,W,1)
      Call zscal_(n*n,redME,O,1)
      Return
      End subroutine Stev_ITO


      Subroutine ESO(n,k,q,O,W,redME)
!  generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
!  only for q >= 0
      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)          :: n,k,q
      Complex(kind=wp), intent(out):: O(n,n), W(n,n), redME
      ! local
      Integer                      :: m1, m2
      Real(kind=wp)                :: CR, C0, knm(12,0:12), F
      Complex(kind=wp)             :: mQ, HALF_R, FALF_I
      Complex(kind=wp), allocatable:: Cp(:,:), Cm(:,:)

      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,O,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,W,1)
      redME=(0.0_wp,0.0_wp)
      If((k>12).OR.(q>12)) Return ! not available in MATLAB

      Call mma_allocate(Cp,n,n,'Cp')
      Call mma_allocate(Cm,n,n,'Cm')
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

      Call set_knm( knm )
      Call coeff_redus_sub(n,k,CR)
      Call ITO(n,k,q,C0,Cp,Cm)

      F=C0*CR/knm(k,q)
      redME=cmplx( F, 0.0_wp, wp )
      mQ   =cmplx( (-1)**q, 0.0_wp, wp )
      HALF_R=(0.5_wp,0.0_wp)
      FALF_I=(0.0_wp,0.5_wp)
      Do m1=1,n
        Do m2=1,n
          O(m1,m2)=HALF_R * redME*( Cm(m1,m2) + mQ*Cp(m1,m2) )
          W(m1,m2)=FALF_I * redME*( Cm(m1,m2) - mQ*Cp(m1,m2) )
        End Do
      End Do
      Call mma_deallocate(Cp)
      Call mma_deallocate(Cm)
      Return
      End subroutine ESO


      Subroutine Liviu_ESO(n,k,q,O,W,redME)
!  generate Hermitian ESO operators, as in MATLAB EasySpin(stev) function
!  only for q >= 0
!  Om = Ok-q
!  Op = Ok+q
! the only difference with MATLAB's ESO are scaling factor Knm
      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)          :: n,k,q
      Complex(kind=wp), intent(out):: O(n,n), W(n,n), redME
      ! local
      Integer                      :: m1, m2
      Real(kind=wp)                :: CR, C0
      Complex(kind=wp)             :: Om, Op, mQ
      Complex(kind=wp), allocatable:: Cp(:,:), Cm(:,:)

      Call mma_allocate(Cp,n,n,'Cp')
      Call mma_allocate(Cm,n,n,'Cm')
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cp,1)
      Call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,Cm,1)

      Call coeff_redus_sub(n,k,CR)
      Call ITO(n,k,q,C0,Cp,Cm)
      redME=cmplx( C0*CR, 0.0_wp, wp )
      mQ   =cmplx( dble((-1)**q), 0.0_wp, wp )
      Do m1=1,n
        Do m2=1,n
          Om = Cm(m1,m2) * redME
          Op = Cp(m1,m2) * redME
          O(m1,m2)=(0.5_wp,0.0_wp) * ( Om + mQ*Op )
          W(m1,m2)=(0.0_wp,0.5_wp) * ( Om - mQ*Op )
        End Do
      End Do
      Call mma_deallocate(Cp)
      Call mma_deallocate(Cm)
      Return
      End subroutine Liviu_ESO



      Subroutine Stewens_matrixel( N, M, dim, ITO_O, ITO_W, IPRINT )
C
C  This Subroutine calculates the matrix elements of the ITO  (On)
C  on the basis of the eigenfunctions of the effective spin.
C
C
C    N -- the rank of the ITO (On)
C    dim -- the multiplicity of the effective spin
C    Dip_Stewens(N, L, dim, dim) -- the matrix elements of the ITO tensor
C                  operators in the basis of effective spin eigenfunctions
C
      Implicit None
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer,intent(in)            :: N, M, dim, IPRINT
      Complex(kind=wp), intent(out) :: ITO_O(dim,dim),
     &                                 ITO_W(dim,dim)
      Integer                       :: npar,i,j,ms1,ms2
      Real(kind=wp)                 :: a,al,b,bt,c,gm, COEFF_REDUS,
     &                                 coeffCG
      Complex(kind=wp)              ::
     &                             ITO_PLUS(-dim:dim,-dim:dim),
     &                            ITO_MINUS(-dim:dim,-dim:dim)
!***********************************************************************
      Call qEnter('Stewens_m')

      NPAR=MOD(dim,2)
      COEFF_REDUS=0.0_wp
      ITO_PLUS( -dim:dim,-dim:dim)=(0.0_wp,0.0_wp)
      ITO_MINUS(-dim:dim,-dim:dim)=(0.0_wp,0.0_wp)
      ITO_O(1:dim,1:dim)=(0.0_wp,0.0_wp)
      ITO_W(1:dim,1:dim)=(0.0_wp,0.0_wp)


      Call COEFF_REDUS_SUB(dim,N,COEFF_REDUS)
      !COEFF_REDUS=1.0_wp
      a = DBLE(N)
      al= DBLE(M)
      c = (DBLE(dim)-1.0_wp)/2.0_wp
      b = (DBLE(dim)-1.0_wp)/2.0_wp
      Do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
      If((ms1==0).AND.(NPAR==0)) Go To 160
        If(NPAR==0) Then
          If(ms1<0) Then
             gm=DBLE(ms1)+0.5_wp
          Else
             gm=DBLE(ms1)-0.5_wp
          End If
        Else
          gm=DBLE(ms1)
        End If

        Do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
          If((ms2==0).AND.(NPAR==0)) Go To 150
          If(NPAR==0) Then
            If(ms2<0) Then
              bt=DBLE(ms2)+0.5_wp
            Else
              bt=DBLE(ms2)-0.5_wp
            End If
          Else
            bt=DBLE(ms2)
          End If
          coeffCG=0.0_wp

          Call Clebsh_Gordan(a,al,b,bt,c,gm, coeffCG)

          ITO_PLUS(ms1,ms2)=CMPLX(coeffCG*COEFF_REDUS,0.0_wp,wp)

          If(iprint>5) Then
            Write(6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))')
     &              'ms1=',ms1,'ms2=',ms2,
     &              'a=',a,'al=',al,
     &              'b=',b,'bt=',bt,
     &              'c=',c,'gm=',gm,
     &              'coeffCG=',coeffCG, 'coeffCG^2=',coeffCG**2
          End If
  150     Continue
        End Do !ms2
  160   Continue
      End Do !ms1

      al=-DBLE(M)

      Do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
        If((ms1==0).AND.(NPAR==0)) Go To 161
        If(NPAR==0) Then
          If(ms1<0) Then
            gm=DBLE(ms1)+0.5_wp
          Else
            gm=DBLE(ms1)-0.5_wp
          End If
        Else
          gm=DBLE(ms1)
        End If

        Do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
          If((ms2==0).AND.(NPAR==0)) Go To 151

          If(NPAR==0) Then
            If(ms2<0) Then
              bt=DBLE(ms2)+0.5_wp
            Else
              bt=DBLE(ms2)-0.5_wp
            End If
          Else
            bt=DBLE(ms2)
          End If
          coeffCG=0.0_wp

          Call Clebsh_Gordan(a,al,b,bt,c,gm, coeffCG)

          ITO_MINUS(ms1,ms2)=CMPLX(coeffCG*COEFF_REDUS,0.0_wp,wp)

          If(iprint>5) Then
            Write(6,'(5x,2(a,i3,2x),6(a,f4.1,2x),2(a,f14.10,2x))')
     &              'ms1=',ms1,'ms2=',ms2,
     &              'a=',a,'al=',al,
     &              'b=',b,'bt=',bt,
     &              'c=',c,'gm=',gm,
     &              'coeffCG=',coeffCG, 'coeffCG^2=',coeffCG**2
          End If
  151     Continue
        End Do
  161   Continue
      End Do

      i=0
      Do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
        If((ms1==0).AND.(NPAR==0)) Go To 180
        i=i+1
        j=0

        Do ms2=-(dim-NPAR)/2,(dim-NPAR)/2
          If((ms2==0).AND.(NPAR==0)) Go To 170
          j=j+1
          ITO_O(i,j)=ITO_PLUS( ms1,ms2)
          ITO_W(i,j)=ITO_MINUS(ms1,ms2)
  170     Continue
        End Do
  180   Continue
      End Do

      If(IPRINT>3) Then
        Write(6,'(/)')
        Write(6,'(5X,A,2I3)') 'Operator ITO_PLUS',N,M
        Write(6,*)
        Do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
          If((ms1==0).AND.(NPAR==0)) Go To 158
          If(NPAR==1)
     &      Write(6,'(16(2X,2E12.3))') (ITO_PLUS(ms1,ms2),
     &                            ms2=-(dim-NPAR)/2,(dim-NPAR)/2)
          If(NPAR==0)
     &      Write(6,'(16(2X,2E12.3))')
     &                    (ITO_PLUS(ms1,ms2),ms2=-dim/2,    -1),
     &                    (ITO_PLUS(ms1,ms2),ms2=     1, dim/2)
  158     Continue
        End Do   ! ms1

        Write(6,'(/)')
        Write(6,'(5X,A,2I3)') 'Operator ITO_MINUS',N,M
        Write(6,*)
        Do ms1=-(dim-NPAR)/2,(dim-NPAR)/2
          If((ms1==0).AND.(NPAR==0)) Go To 159
          If(NPAR==1)
     &      Write(6,'(16(2X,2E12.3))') (ITO_MINUS(ms1,ms2),
     &                            ms2=-(dim-NPAR)/2,(dim-NPAR)/2)
          If(NPAR==0)
     &      Write(6,'(16(2X,2E12.3))')
     &                    (ITO_MINUS(ms1,ms2),ms2=-dim/2,   -1),
     &                    (ITO_MINUS(ms1,ms2),ms2=     1,dim/2)
  159     Continue
        End Do   ! ms1


        Write(6,'(/////)')
        Write(6,'(5X,A,2I3)') 'ITO_O',N,M
        Write(6,*)
        Do i=1,dim
          Write(6,'(16(2X,2E12.3))') (ITO_O(i,j),j=1,dim)
        End Do   ! i
        Write(6,'(/)')
        Write(6,'(5X,A,2I3)') 'ITO_W',N,M
        Write(6,*)
        Do i=1,dim
          Write(6,'(16(2X,2E12.3))') (ITO_W(i,j),j=1,dim)
        End Do   ! i
      End If ! iPrint

      Call qExit('Stewens_m')

      Return
      End


      Real*8 Function fct(n)
      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: n
      Integer             :: i
      Real(kind=wp)       :: xct
      ! this function provides correct answer till n=169 only

      xct=1.0_wp
      fct=1.0_wp
      If ( n<0 ) Then
        Write(6,'(A,i0)') 'FCT:  N<0 !'
        Write(6,'(A,i0)') 'N = ', N
        Write(6,'(A   )') 'It is an impossible case.'
        fct=-9.d99
        Return

      Else If ( n==0 ) Then
        Return

      Else If ( (n>0) .AND. (n.le.169) ) Then
        Do i=1,n
          xct=xct*DBLE(i)
        End Do

      Else
        Write(6,'(A,i0)') 'FCT:   N = ',N
        Write(6,'(A)') 'Factorial of N>169 overflows on x86_64'
        Write(6,'(A)') 'Use higher numerical precision, ' //
     &                 'or rethink your algorithm.'
      End If

      fct=xct

      Return
      End function fct

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine COEFF_REDUS_SUB(dim,N,COEFF_REDUS)
!
!    THIS Subroutine ReturnS THE VALUE OF THE REDUCED MATRIX ELEMENT  <S|| On ||S>
!   WHERE S-EFFECTIVE SPIN, On -- THE HIGHER ORDER SPIN OPERATORS
!
!     The convention is to follow the paper:
!     C. Rudowicz and C.Y. Chung
!     J. Phys.: Condens. Matter, 2004, 16, pp. 5825
!     DOI: https://doi.org/10.1088/0953-8984/16/32/018
!
!     with a minor addition of the Norm(N) factor which depends on the operator rank
!     This norm makes the matrix elements of ESO identical to those of the
!     operators from EasySpin (stev).

      Implicit None
      Integer, Parameter         :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)        :: N, dim
      Real(kind=wp), intent(out) :: COEFF_REDUS
      Integer                    :: i
      Real(kind=wp)              :: SZ, FCT, Norm(100), s1, s2
      external                   :: fct

      COEFF_REDUS=0.0_wp
      SZ=0.0_wp
      Do i=1,100
        Norm(i)=0.0_wp
      End Do
      !SZ=DBLE(dim-1)/2.0_wp
      Norm(1) = 1.0_wp
      Norm(2) = 2.0_wp
      Norm(3) = 2.0_wp
      Norm(4) = 8.0_wp
      Norm(5) = 8.0_wp
      Norm(6) = 16.0_wp
      Norm(7) = 16.0_wp
      Norm(8) = 128.0_wp
      Norm(9) = 128.0_wp
      Norm(10)= 256.0_wp
      Norm(11)= 256.0_wp
      Norm(12)= 1024.0_wp
      Norm(13)= 1024.0_wp
      Norm(14)= 2048.0_wp
      Norm(15)= 2048.0_wp
      Norm(16)= 32768.0_wp
      Norm(17)= 32768.0_wp
      Norm(18)= 65536.0_wp
      Norm(19)= 65536.0_wp
      Norm(20)= 262144.0_wp
      Norm(21)= 262144.0_wp
      Norm(22)= 524288.0_wp
      Norm(23)= 524288.0_wp
      Norm(24)= 4194304.0_wp
      Norm(25)= 4194304.0_wp
      Norm(26)= 8388608.0_wp
      Norm(27)= 8388608.0_wp
      Norm(28)= 33554432.0_wp
      Norm(29)= 33554432.0_wp
      Norm(30)= 67108867.0_wp
      Norm(31)= 67108867.0_wp
      Norm(32)= 2147483648.0_wp
      Norm(33)= 2147483648.0_wp
      Norm(34)= 4294967296.0_wp
      Norm(35)= 4294967296.0_wp
      Norm(36)= 17179869184.0_wp
      Norm(37)= 17179869184.0_wp
      Norm(38)= 34359738368.0_wp
      Norm(39)= 34359738368.0_wp
      Norm(40)= 274877906944.0_wp
      Norm(41)= 274877906944.0_wp
      Norm(42)= 549755813888.0_wp
      Norm(43)= 549755813888.0_wp
      Norm(44)= 2199023255552.0_wp
      Norm(45)= 2199023255552.0_wp
      Norm(46)= 4398046511104.0_wp
      Norm(47)= 4398046511104.0_wp
      Norm(48)= 70368744177664.0_wp
      Norm(49)= 70368744177664.0_wp
      Norm(50)= 140737488355328.0_wp
      Norm(51)= 140737488355328.0_wp
      Norm(52)= 562949953421312.0_wp
      Norm(53)= 562949953421312.0_wp
      Norm(54)= 1125899906842624.0_wp
      Norm(55)= 1125899906842624.0_wp
      Norm(56)= 9007199254740992.0_wp
      Norm(57)= 9007199254740992.0_wp
      Norm(58)= 18014398509481984.0_wp
      Norm(59)= 18014398509481984.0_wp
      Norm(60)= 72057594037927936.0_wp
      Norm(61)= 72057594037927936.0_wp
      Norm(62)= 144115188075855872.0_wp
      Norm(63)= 144115188075855872.0_wp
      Norm(64)= 9223372036854775808.0_wp
      Norm(65)= 9223372036854775808.0_wp
      Norm(66)= 18446744073709551616.0_wp
      Norm(67)= 18446744073709551616.0_wp
      Norm(68)= 73786976294838206464.0_wp
      Norm(69)= 73786976294838206464.0_wp
      Norm(70)= 147573952589676412928.0_wp
      Norm(71)= 147573952589676412928.0_wp
      Norm(72)= 1180591620717411303424.0_wp
      Norm(73)= 1180591620717411303424.0_wp
      Norm(74)= 2361183241434822606848.0_wp
      Norm(75)= 2361183241434822606848.0_wp
      Norm(76)= 9444732965739290427392.0_wp
      Norm(77)= 9444732965739290427392.0_wp
      Norm(78)= 18889465931478580854784.0_wp
      Norm(79)= 18889465931478580854784.0_wp
      Norm(80)= 302231454903657293676544.0_wp
      Norm(81)= 302231454903657293676544.0_wp
C  new method
      s1=0.0_wp
      s2=0.0_wp
      s1=SQRT( fct(dim+N)/fct(dim-N-1) )
      s2=dble(2**N)*SQRT(DBLE(dim))

      COEFF_REDUS= Norm(N) * s1 / s2

      Return
      End



      Subroutine Clebsh_Gordan(a,al,b,bt,c,gm,coeffCG)
      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Real(kind=wp),intent(in)  :: a, al, b, bt, c, gm
      Real(kind=wp),intent(out) :: coeffCG
      Real(kind=wp)             :: u, fct, s1, s2
      Integer                   :: lb1, lb2, i
      Logical                   :: check_triangle
      External                  :: check_triangle, fct
c exclude the cases for which CG coefficients are exactly zero
      coeffCG=0.0_wp
      If((al+bt).ne.gm) Return
      If(a<0.0_wp) Return
      If(b<0.0_wp) Return
      If(c<0.0_wp) Return
      If(ABS(al)>a) Return
      If(ABS(bt)>b) Return
      If(ABS(gm)>c) Return
      If((ABS(a-b)>c).or.((a+b)<c)) Return
      If((ABS(b-c)>a).or.((b+c)<a)) Return
      If((ABS(c-a)>b).or.((c+a)<b)) Return
      If(MOD(nint(2.0_wp*a),2) .ne. MOD(nint(2.0_wp*ABS(al)),2)) Return
      If(MOD(nint(2.0_wp*b),2) .ne. MOD(nint(2.0_wp*ABS(bt)),2)) Return
      If(MOD(nint(2.0_wp*c),2) .ne. MOD(nint(2.0_wp*ABS(gm)),2)) Return
      u=0.0_wp
      lb1=INT(MIN(c-b+al,c-a-bt))
      lb2=INT(MIN(a+b-c,a-al,b+bt))
      If(lb1<0) Then
        If(-lb1>lb2) Then
          coeffCG=0.0_wp
          Return
        Else

          Do i=-lb1,lb2
            u = u+DBLE( (-1)**i )/
     &            DBLE(fct(i)*fct(NINT( a-al -i ))
     &                       *fct(NINT( b+bt -i ))
     &                       *fct(NINT( a+b-c-i ))
     &                       *fct(NINT( c-b+al+i))
     &                       *fct(NINT( c-a-bt+i)) )
          End Do
        End If
      Else
        Do i =0,lb2
          u = u+DBLE( (-1)**i )/
     &          DBLE(fct(i)*fct(NINT( a-al -i ))
     &                     *fct(NINT( b+bt -i ))
     &                     *fct(NINT( a+b-c-i ))
     &                     *fct(NINT( c-b+al+i))
     &                     *fct(NINT( c-a-bt+i)) )
        End Do
      End If

      s1=0.0_wp
      s2=0.0_wp
      s1=SQRT( DBLE( fct(NINT( a+b-c  ))
     &              *fct(NINT( a-b+c  ))
     &              *fct(NINT(-a+b+c  )) )
     &        /DBLE( fct(NINT( a+b+c+1)) )  )

      s2=SQRT( DBLE( fct(NINT( a+al ))
     &              *fct(NINT( a-al ))
     &              *fct(NINT( b+bt ))
     &              *fct(NINT( b-bt ))
     &              *fct(NINT( c+gm ))
     &              *fct(NINT( c-gm ))
     &              *(2*c+1)             )  )

      coeffCG = u*s1*s2

      Return
      End


      Real*8 function W9J(a,b,c, d,e,f, g,h,j)
c Calculates a Wigner 9-j symbol. Argument a-j are Integer and are
c twice the true value of the 9-j's arguments, in the form
c { a b c }
c { d e f }
c { g h j }
c this is the implementation the formula 10.2.4. (20) from:
c   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.
c
      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: a,b,c,d,e,f,g,h,j
      Integer             :: n,nlow,nhig
      Real(kind=wp)       :: dlt,fct,W6J
      Logical             :: check_triangle
      External            :: fct,dlt,W6J,check_triangle

      W9j =0.0_wp
      If(MOD(a+b,2) .ne. MOD(c,2)) Return
      If(MOD(d+e,2) .ne. MOD(f,2)) Return
      If(MOD(g+h,2) .ne. MOD(j,2)) Return
      If(MOD(a+d,2) .ne. MOD(g,2)) Return
      If(MOD(b+e,2) .ne. MOD(h,2)) Return
      If(MOD(c+f,2) .ne. MOD(j,2)) Return
      If((ABS(a-b) > c) .or. (a+b < c)) Return
      If((ABS(d-e) > f) .or. (d+e < f)) Return
      If((ABS(g-h) > j) .or. (g+h < j)) Return
      If((ABS(a-d) > g) .or. (a+d < g)) Return
      If((ABS(b-e) > h) .or. (b+e < h)) Return
      If((ABS(c-f) > j) .or. (c+f < j)) Return
      If(check_triangle(a,b,c) .eqv. .false.) Return
      If(check_triangle(d,e,f) .eqv. .false.) Return
      If(check_triangle(g,h,j) .eqv. .false.) Return
      If(check_triangle(a,d,g) .eqv. .false.) Return
      If(check_triangle(b,e,h) .eqv. .false.) Return
      If(check_triangle(c,f,j) .eqv. .false.) Return

      nlow=MAX( ABS(a-j)/2, ABS(d-h)/2, ABS(b-f)/2  )
      nhig=MIN(    (a+j)/2,    (d+h)/2,    (b+f)/2  )

      Do n=nlow,nhig
      W9j=W9j + DBLE(2*n+1)*DBLE((-1)**(2*n))
     &        * W6J(a,b,c,  f,  j,2*n )
     &        * W6J(d,e,f,  b,2*n,  h )
     &        * W6J(g,h,j,2*n,  a,  d )
      End Do
      Return
      End function W9J


      Real*8 function W6J(a,b,c,d,e,f)
c Calculates a Wigner 6-j symbol. Argument a-f are positive Integer
c and are twice the true value of the 6-j's arguments, in the form
c { a b c }
c { d e f }
c
c this is the implementation the formula 9.2.1. (1) from:
c   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: a,b,c,d,e,f
      Integer             :: n,nlow,nhig
      Real(kind=wp)       :: dlt,sum,fct,isum
      Logical             :: check_triangle
      External            :: fct,dlt,check_triangle

      W6J=0.0_wp
      If(MOD(a+b,2) .ne. MOD(c,2)) Return
      If(MOD(c+d,2) .ne. MOD(e,2)) Return
      If(MOD(a+e,2) .ne. MOD(f,2)) Return
      If(MOD(b+d,2) .ne. MOD(f,2)) Return
      If((ABS(a-b) > c) .or. (a+b < c)) Return
      If((ABS(c-d) > e) .or. (c+d < e)) Return
      If((ABS(a-e) > f) .or. (a+e < f)) Return
      If((ABS(b-d) > f) .or. (b+d < f)) Return

      If(check_triangle(a,b,c).eqv. .false.) Return
      If(check_triangle(c,d,e).eqv. .false.) Return
      If(check_triangle(a,e,f).eqv. .false.) Return
      If(check_triangle(b,d,f).eqv. .false.) Return

      nlow=0
      nhig=0
      nlow = MAX( (a+b+c)/2, (c+d+e)/2, (b+d+f)/2, (a+e+f)/2 )
      nhig = MIN( (a+b+d+e)/2, (b+c+e+f)/2, (a+c+d+f)/2)

      sum =0.0_wp
      Do n=nlow,nhig
      isum=DBLE((-1)**n)*fct(n+1)
     & /fct(  ( a+c+d+f)/2-n)
     & /fct(  ( b+c+e+f)/2-n)
     & /fct(n-( a+b+c  )/2  )
     & /fct(n-( c+d+e  )/2  )
     & /fct(n-( a+e+f  )/2  )
     & /fct(n-( b+d+f  )/2  )
     & /fct(  ( a+b+d+e)/2-n)
      sum=sum+isum
      End Do
      W6J=dlt(a,b,c)*dlt(c,d,e)*dlt(a,e,f)*dlt(b,d,f)*sum
      Return
      End function W6J


      Real*8 function W3J(j1,j2,j3,m1,m2,m3)
c Calculates a Wigner 3-j symbol. Argument j1,j2,j3 are positive Integer
c and are twice the true value of the 3-j's arguments, in the form
c { j1 j2 j3 }
c { m1 m2 m3 }
c
c this is the implementation the formula 8.1.2. (11) from:
c   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: j1, j2, j3, m1, m2, m3
      Real(kind=wp)       :: fct, dlt, coeffCG
      Logical             :: check_triangle
      External            :: check_triangle, fct, dlt

      W3J=0.0_wp
      coeffCG=0.0_wp
      Call Clebsh_Gordan(DBLE(j1)/2.0_wp, DBLE(m1)/2.0_wp,
     &                   DBLE(j2)/2.0_wp, DBLE(m2)/2.0_wp,
     &                   DBLE(j3)/2.0_wp,-DBLE(m3)/2.0_wp, coeffCG )
      If(coeffCG==0.0_wp) Return
      W3J=DBLE((-1)**((j1-j2-m3)/2))*coeffCG/SQRT(DBLE(j3+1))
      Return
      End function W3J


      Real*8 function WCG(a, al, b, bt, c, gm)
c Calculates a Clebsch-Gordan Coefficient. Argument a, al, b, bt, c, gm are Integer,
c Double their actual value.
c (   c/2, gm/2            )
c ( C                      )
c (   a/2, al/2, b/2, bt/2 )
c this is the implementation the formula 8.2.1. (3) from:
c   D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c   "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: a, al, b, bt, c, gm
      Integer             :: lb1,lb2,i
      Real(kind=wp)       :: u,fct,dlt
      External            :: fct,dlt

      WCG=0.0_wp

      If((al+bt).ne.gm) Return
      If(a<0) Return
      If(b<0) Return
      If(c<0) Return
      If(ABS(al)>a) Return
      If(ABS(bt)>b) Return
      If(ABS(gm)>c) Return
      If((ABS(a-b)>c).or.((a+b)<c)) Return
      If((ABS(b-c)>a).or.((b+c)<a)) Return
      If((ABS(c-a)>b).or.((c+a)<b)) Return
      If(MOD(a,2) .ne. MOD(ABS(al),2)) Return
      If(MOD(b,2) .ne. MOD(ABS(bt),2)) Return
      If(MOD(c,2) .ne. MOD(ABS(gm),2)) Return
      u=0.0_wp
      lb1=MIN( (c-b+al)/2,(c-a-bt)/2 )
      lb2=MIN( (a+b-c)/2 , (a-al)/2 , (b+bt)/2 )
      If(lb1<0) Then
        If(-lb1>lb2) Then
          WCG=0.0_wp
          Return
        Else
          Do i=-lb1,lb2
            u = u + DBLE( (-1)**i )/( fct(i) * fct( (a+b-c -2*i)/2)
     &                                       * fct( (c-b+al+2*i)/2)
     &                                       * fct( (c-a-bt+2*i)/2)
     &                                       * fct( (a-al  -2*i)/2)
     &                                       * fct( (b+bt  -2*i)/2) )
          End Do
        End If
      Else
        Do i=0,lb2
          u = u + DBLE( (-1)**i )/( fct(i) * fct( (a+b-c -2*i)/2)
     &                                     * fct( (c-b+al+2*i)/2)
     &                                     * fct( (c-a-bt+2*i)/2)
     &                                     * fct( (a-al  -2*i)/2)
     &                                     * fct( (b+bt  -2*i)/2)  )
        End Do
      End If
      WCG=u*dlt(a,b,c)*SQRT( fct(( a+al)/2) * fct(( a-al)/2)
     &                     * fct(( b+bt)/2) * fct(( b-bt)/2)
     &                     * fct(( c+gm)/2) * fct(( c-gm)/2)
     &                     * (c+1)
     &                     )

c      Write(6,'(A)') 'a,  al,  b,  bt,  c,  gm'
c      Write(6,'(6(F4.1,2x),F20.14)')
c     & DBLE(a)/2.0_wp, DBLE(al)/2.0_wp,DBLE(b)/2.0_wp,
c     & DBLE(bt)/2.0_wp, DBLE(c)/2.0_wp, DBLE(gm)/2.0_wp, WCG
      Return
      End function WCG



      Real*8 function dlt(a,b,c)
c  calculates the delta(a,b,c) function using the formula 8.2.1. from:
c    D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c    "Quantum Theory of Angular Momentum", World ScientIfic, 1988.
c
c  a,b,c are positive Integer numbers,
c  their values are DoUBLE than their original value
      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: a,b,c
      Real(kind=wp)       :: fct
      Logical             :: check_triangle
      External            :: check_triangle, fct

      dlt=0.0_wp
      If((ABS(a-b)>c).or.(a+b<c)) Return
      If((ABS(b-c)>a).or.(b+c<a)) Return
      If((ABS(c-a)>b).or.(c+a<b)) Return
      If(MOD(( a+b-c),2) == 1) Return
      If(MOD(( a-b+c),2) == 1) Return
      If(MOD((-a+b+c),2) == 1) Return
      If(MOD(( a+b+c),2) == 1) Return
      If(check_triangle(a,b,c).eqv. .false.) Return
      ! special cases:
      If(a==0) Then
        dlt=1.0_wp/SQRT(DBLE(b+1))
      End If
      If(b==0) Then
        dlt=1.0_wp/SQRT(DBLE(a+1))
      End If
      If(c==0) Then
        dlt=1.0_wp/SQRT(DBLE(a+1))
      End If

      dlt=SQRT(fct((a+b-c)/2)*fct((a-b+c)/2)*fct((-a+b+c)/2)
     &         /fct((a+b+c)/2+1))
      Return
      End function dlt



      Logical function check_triangle(a,b,c)
c  checks If the values a,b,c comply with the triangle rule
      Implicit None
      Integer, intent(in) :: a, b, c

      check_triangle=.false.

      If ( (a<=0) .OR. (b<=0) .OR. (c<=0) ) then
        Write(6,'(A)') 'a=',a
        Write(6,'(A)') 'b=',b
        Write(6,'(A)') 'c=',c
        Write(6,'(A)') 'The rule is: a>0, b>0 and c>0!'
        Write(6,'(A)') 'Please check this issue, or report a bug!'
        Return
      End If

      If( ((a+b)>=c) .and. ((b+c)>=a) .and. ((c+a)>=b) ) then
        check_triangle=.true.
      End If
      Return
      End function check_triangle



      Complex*16 function WignerD(J,M1,M2,al,bt,gm)
c the function Returns the Wigner-D function specIfying the rotation
c of the |J,M1,M2> around three  angles(alpha,beta,gamma).
c The rotation is either active (ik=1) or passive (ik=2).
c   J, M1, M2 are specIfied as Integer numbers, with a value DoUBLE than their actual size;
c   i.e. J = 2*J (Real); M1= 2*M1(Real); M2=2*M2(Real);
c   alpha, beta, gamma are specIfied as Double precision (Real(kind=wp) ::). These values must be defined in
c   radians ( i.e. in units of Pi)
c
c
c  This is the implementation of the formula 4.3.(1) from
c    D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c    "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)       :: J, M1, M2
      Real(kind=wp), intent(in) :: al, bt, gm
      Complex(kind=wp)          :: m1_fact, m2_fact, wig_fac
      Real(kind=wp)             :: wigner_d
      External                  :: wigner_d

c  check correctness
      WignerD=(0.0_wp,0.0_wp)
      wig_fac=(0.0_wp,0.0_wp)
      m1_fact=(0.0_wp,0.0_wp)
      m2_fact=(0.0_wp,0.0_wp)
      If(ABS(M1)>J) Return
      If(ABS(M2)>J) Return
      If(ABS(J )<0) Return

      m1_fact=EXP( (0.0_wp,-1.0_wp)**(DBLE(al*M1)/2.0_wp) )
      m2_fact=EXP( (0.0_wp,-1.0_wp)**(DBLE(gm*M2)/2.0_wp) )
      wig_fac=CMPLX( wigner_d(J,M1,M2,bt), 0.0_wp, wp )
      WignerD=m1_fact*m2_fact*wig_fac
      Return
      End function WignerD

      Real*8 function wigner_d(J,M1,M2,bt)
c  This is the implementation of the formula 4.3.1 (2) from
c    D.A. Varshalovich, A.N. Moskalev, V.K. Khersonskii,
c    "Quantum Theory of Angular Momentum", World ScientIfic, 1988.

      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)       :: J, M1, M2
      Real(kind=wp), intent(in) :: bt
      Real(kind=wp)             :: ksum, fct
      Integer                   :: kmin, kmax, i
      External                  :: fct

      wigner_d=0.0_wp
      ksum=0.0_wp
      kmax=MIN((J -M1)/2,(J -M2)/2)
      kmin=MAX(0,-(M1+M2)/2)
      Do i=kmin,kmax
        ksum=DBLE((-1)**(i))*(COS(bt/2.0_wp)**DBLE(  M1/2+M2/2+2*i))
     &                      *(SIN(bt/2.0_wp)**DBLE(J-M1/2-M2/2-2*i))
     &                     /fct(            i )
     &                     /fct( (J -M1)/2 -i )
     &                     /fct( (J -M2)/2 -i )
     &                     /fct( (M1+M2)/2 +i )
        wigner_d=wigner_d+ksum
      End Do
      wigner_d=wigner_d*DBLE( (-1)**((J-M2)/2) )*
     &         SQRT( fct((J+M1)/2)*fct((J-M1)/2)*
     &               fct((J+M2)/2)*fct((J-M2)/2)  )
      Return
      End function wigner_d

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Real*8 function RedME( La,Sa,  LaP,SaP,  L,S )
c function evaluates the reduced matrix elements of the ground atomic J multipet
c
c  the function evaluates the formula:
c <L_A, S_A || Operator ( iL_T, iS_T)  || L_A, S_A >>   formula (S.20) derived by Naoya Iwahara
c
c iL, iS, iL_T, iS_T, Jp, Sp and Lp  are defined as Integers with a value DoUBLE of their true value
c
c the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only

      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in) :: La, Sa, LaP, SaP, L, S
      Integer             :: Ja, JaP, jm, js, s_orb, l_orb
      Real(kind=wp)       :: WCG, W9J, temp, factor
      external            :: WCG, W9J

      RedME=0.0_wp
      temp=0.0_wp
      l_orb=6  ! Double of true value l_orb = 3
      s_orb=1  ! Double of true value s_orb = 1/2
      JaP  = LaP + SaP
      Ja   = La  + Sa
      If(WCG(La, La, L, 0, La, La)==0.0_wp) Return
      If(WCG(Sa, Sa, S, 0, Sa, Sa)==0.0_wp) Return
      If(WCG(La, La, l_orb, LaP-La, LaP, LaP)==0.0_wp) Return
      If(WCG(Sa, Sa, s_orb, SaP-Sa, SaP, SaP)==0.0_wp) Return

      factor=SQRT( DBLE((La+1)*(Sa+1)) ) / WCG(La, La, L, 0, La, La)
     &                                   / WCG(Sa, Sa, S, 0, Sa, Sa)

      Do jm = -l_orb, l_orb
        Do js = -s_orb, s_orb
          temp = temp + DBLE( (-1)**((l_orb+jm+s_orb+js)/2) )
     &                 * WCG(l_orb, -jm, l_orb, jm, L, 0)
     &                 * WCG(s_orb, -js, s_orb, js, S, 0)
     &                 * WCG(LaP, La+jm, SaP, Sa+js, JaP, La+jm+Sa+js)
     &                 * WCG(LaP, La+jm, SaP, Sa+js, JaP, La+jm+Sa+js)
     &                 * WCG(La, La, l_orb, jm, LaP, La+jm)
     &                 * WCG(La, La, l_orb, jm, LaP, La+jm)
     &                 * WCG(Sa, Sa, s_orb, js, SaP, Sa+js)
     &                 * WCG(Sa, Sa, s_orb, js, SaP, Sa+js)
     &                 / WCG(La, La, l_orb, LaP-La, LaP, LaP)
     &                 / WCG(La, La, l_orb, LaP-La, LaP, LaP)
     &                 / WCG(Sa, Sa, s_orb, SaP-Sa, SaP, SaP)
     &                 / WCG(Sa, Sa, s_orb, SaP-Sa, SaP, SaP)
        End Do !is
      End Do !im

      RedME=temp*factor

      Return
      End function RedME


      Real*8 function jot1( t, L, ML, S, MS,  La, Sa,  LaP, SaP )
c function evaluates the J1 exchange matrix element ( formula S.19)
c
c iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
c the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
c    Substitutions:

      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)       :: L,ML, S,MS, La,Sa, LaP,SaP
      Integer                   :: Ja, JaP, l_orb
      Real(kind=wp), intent(in) :: t
      Real(kind=wp)             :: W9J, WCG, RedME, txt
      External                  :: W9J, WCG, RedME

      jot1 = 0.0_wp
      txt=0.0_wp
      l_orb=6  ! Double of true value l_orb = 3
c      s_orb=1  ! Double of true value s_orb = 1/2
      If (MOD(L,2)==0) Then
        If (ML == 0) Then
          txt =        DBLE((-1)**(l_orb/2))*WCG(l_orb,-4,l_orb,4,L,0)
        Else If (ML == 8) Then
          txt = -0.5_wp*DBLE((-1)**(l_orb/2))*WCG(l_orb, 4,l_orb,4,L,8)
        Else If (ML ==-8) Then
          txt = -0.5_wp*DBLE((-1)**(l_orb/2))*WCG(l_orb, 4,l_orb,4,L,8)
        End If
      End If
      Ja  = La  + Sa
      JaP = LaP + SaP

      jot1= t * txt * SQRT( DBLE( (Ja+1)*(L+s+1) ))
     &              * W9J(Ja,La,Sa, Ja,La,Sa,L+s,L,2)
     &              * WCG(L,ML,2,MS,L+s,ML+MS)
     &              * WCG(Ja,Ja,L+s,0,Ja,Ja)
     &              * RedME( La, Sa, LaP, SaP, L, 2  )/SQRT(2.0_wp)

      Return
      End function jot1




      Real*8 function jot0( t,  L, ML,   La, Sa,  LaP, SaP )
c function evaluates the J1 exchange matrix element ( formula S.19)
c
c iL_T,iL_TM,iS_T,iS_TM,iL,iS,Sp,Lp are defined as Integers with a value DoUBLE of their true value
c the formula is valid for Tb, Dy, Ho, Er, Tm and Yb only
c    Substitutions:

      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)       :: L, ML, La, Sa, LaP, SaP
      Integer                   :: Ja, l_orb
      Real(kind=wp), intent(in) :: t
      Real(kind=wp)             :: W6J, WCG, RedME, txt, W9Jl
      External                  :: W6J, WCG, RedME

      jot0 = 0.0_wp
      l_orb=6  ! Double of true value l_orb = 3
c     s_orb=1  ! Double of true value s_orb = 1/2
      txt=0.0_wp
      If (MOD(L,4)==0) Then
              If (ML == 0) Then
              txt = DBLE((-1)**(l_orb/2))*WCG(l_orb,-4,l_orb,4,L,0)
          Else If (ML == 8) Then
              txt = -0.5_wp*DBLE((-1)**(l_orb/2))*
     &                                       WCG(l_orb,4,l_orb,4,L,8)
          Else If (ML ==-8) Then
              txt = -0.5_wp*DBLE((-1)**(l_orb/2))*
     &                                       WCG(l_orb,4,l_orb,4,L,8)
           End If
      End If
       Ja  = La  + Sa
      W9Jl=0.0_wp
      W9Jl=DBLE((-1)**( (La+Ja+Sa+L)/2 ))
     & *SQRT(DBLE( (Sa+1)*(L+1) ))
     & *W6J(Ja,La,Sa, La,Ja, L )

      jot0= -t * txt * SQRT( DBLE( (Ja+1)*(L+1) ))
     & * W9Jl
     & * WCG(Ja,Ja,L,0,Ja,Ja)
     & * RedME( La,Sa,  LaP,SaP,  L,0 )/SQRT(2.0_wp)

      Return
      End function jot0


      Subroutine verify_CG(N)
      Implicit none
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      integer :: n,k,q,m1,m2
      real(wp) :: rJ,rK,rQ,mf,rM1,rM2,rfE,rfG
      real(wp) :: CG_A,CG_B,CG_C,CG_D,CG_E,CG_F,CG_G,CG_H
      logical  :: prn

      !2J=n-1
      rJ=dble(n-1)/2.0_wp

      Do k=1,n-1
       Do q=0,k
        rK=dble(k)
        rQ=dble(q)

        Do m1=1,n
         Do m2=1,n
          ! real value for the projections:
          rM1=-rJ+dble(m1-1)
          rM2=-rJ+dble(m2-1)

          CG_A=0.0_wp
          CG_B=0.0_wp
          CG_C=0.0_wp
          CG_D=0.0_wp
          CG_E=0.0_wp
          CG_F=0.0_wp
          CG_G=0.0_wp
          CG_H=0.0_wp
          mf=(-1)**nint(rK)
                            !(a , alpha, b, beta,    c,  gamma)
          Call Clebsh_Gordan(rJ,  rM2,   rK,    rQ,   rJ,  rM1, CG_A)
          Call Clebsh_Gordan(rK,   rQ,   rJ,   rM2,   rJ,  rM1, CG_B)
          Call Clebsh_Gordan(rJ, -rM2,   rK,   -rQ,   rJ, -rM1, CG_C)
          Call Clebsh_Gordan(rK,  -rQ,   rJ,  -rM2,   rJ, -rM1, CG_D)


          rfE=((-1)**(rJ-rM2))*(sqrt( dble(n)/(2.0_wp*rK+1.0_wp) ))
          Call Clebsh_Gordan(rJ,  rM2,   rJ,  -rM1,   rK,  -rQ, CG_E)
          Call Clebsh_Gordan(rJ,  rM1,   rJ,  -rM2,   rK,   rQ, CG_F)

          rfG=(-1)**(rK+rQ)
          Call Clebsh_Gordan(rJ, -rM1,   rK,    rQ,   rJ, -rM2, CG_G)
          Call Clebsh_Gordan(rK,  -rQ,   rJ,   rM1,   rJ,  rM2, CG_H)

          prn=(CG_A.ne.0.0_wp).OR.(CG_B.ne.0.0_wp).OR.(CG_C.ne.0.0_wp)
     &    .OR.(CG_D.ne.0.0_wp).OR.(CG_E.ne.0.0_wp).OR.(CG_F.ne.0.0_wp)
     &    .OR.(CG_G.ne.0.0_wp).OR.(CG_H.ne.0.0_wp)

         If(prn)
     &     print '(A,1x,8F12.6)','n,k,q,CG:',CG_A,mf*CG_B,mf*CG_C,CG_D,
     &            rfE*CG_E, rfE*CG_F, rfG*CG_G, rfG*CG_H

         End Do
        End Do

       End Do
      End Do
      Return
      End subroutine verify_CG
