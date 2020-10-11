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
      Subroutine g_high( esom, GRAD, s_som, dipsom, imltpl, dim,
     &                   Do_structure_abc, cryst, coord, gtens, maxes,
     &                   iprint )

c     this routine calculates the g-tensor and d-tensor in the basis of the any effective spin,
c     (coming from 1 molecular term)
c
      Implicit None
      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer, intent(in):: imltpl, dim, iprint
      logical, intent(in):: Do_structure_abc, GRAD

      Real(kind=8), intent(in) :: esom(dim), cryst(6), coord(3)
      Real(kind=8), intent(out):: gtens(3), maxes(3,3)
      Complex(kind=8),intent(in) :: dipsom(3,dim,dim),
     &                                s_som(3,dim,dim)
      ! local variables:
      Integer :: i

C intializations
      If(Iprint>2) Then
        CALL prMom('G_HIGH:  DIPSOM(l,i,j):',dipsom,dim)
        CALL prMom('G_HIGH:   S_SOM(l,i,j):', s_som,dim)
      End If
c--------------------------------------------------------------------------
      Write(6,'(/)')
      Write(6,'(100A)') ('%',i=1,95)
      If(MOD(dim,2)==0) Then
        Write(6,'(5X,A,I2,A,I2,A)') 'CALCULATION OF PSEUDOSPIN '//
     &          'HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL,
     &          ' ( effective S = ',dim-1,'/2)'
      Else
        Write(6,'(5X,A,I2,A,I1,A)') 'CALCULATION OF PSEUDOSPIN '//
     &          'HAMILTONIAN TENSORS FOR THE MULTIPLET',iMLTPL,
     &          ' ( effective S = ',(dim-1)/2,')'
      End If
      Write(6,'(100A)') ('%',i=1,95)
      Write(6,'(A)') 'The pseudospin is defined in the basis of the '//
     &               'following spin-orbit states:'
      Do i=1,dim
        If(dim>9) Then
          Write(6,'(a,i2,a,i2,a,f11.3,a)') 'spin-orbit state',i,
     &                      '; energy(',i,') = ',ESOM(i),' cm-1.'
        Else
          Write(6,'(a,i1,a,i1,a,f11.3,a)') 'spin-orbit state ',i,
     &                      '. energy(',i,') = ',ESOM(i),' cm-1.'
        End If
      End Do
      If (dim==2) Write(6,'(a,f17.10,a)') 'Tunnelling splitting:',
     &                                     ESOM(2)-ESOM(1),' cm-1.'

       CALL G_HIGH_1( iMLTPL, dim, ESOM, GRAD, S_SOM, dipsom,
     &                Do_structure_abc, cryst, coord, gtens, maxes,
     &                iprint)

      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Subroutine G_HIGH_1( iMLTPL, dim, ESOM, GRAD, S_SOM, dipsom,
     &                     Do_structure_abc, cryst, coord, gtens, maxes,
     &                     iprint)
C
C     This routine calculates the g-tensor and D-tensor in the basis of the any effective spin,
C     (COMING FROM 1 MOLECULAR TERM)
C
C     dim ---  the multiplicity of the effective spin
C
      Implicit None
      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "barrier.fh"
#include "stdalloc.fh"
      Integer, intent(in)         :: dim, iMLTPL, iprint
      Real (kind=8), intent(in)  :: ESOM(dim), cryst(6), coord(3)
      Real (kind=8), intent(out) :: gtens(3), maxes(3,3)
      Complex (kind=8),intent(in):: s_som(3,dim,dim), dipsom(3,dim,dim)
      Logical, intent(in)         :: Do_structure_abc, GRAD
      ! local variables:
      Integer           :: I, L, M, J, N, I1, I2, nmax, IsFreeUnit,
     &                     LuDgrad, rc
      Real (kind=8)    :: ESUM, E0, CHECK_SGN2
      Real (kind=8)    :: knm(12,0:12)
      Real (kind=8), allocatable :: ELOC(:)
      Real (kind=8), allocatable :: axes_in_abc(:,:)

      Complex (kind=8), allocatable ::
     &                   DIP_O(:,:), DIP_W(:,:), MUX(:,:), MUY(:,:),
     &                   MUZ(:,:), MUXZ(:,:), MUZX(:,:), HZFS(:,:),
     &                   DIP_MOW(:,:), HZFS_MONM(:,:), HZFS_MWNM(:,:),
     &                   ZOUT(:,:), AMS(:,:,:), AMSSPIN(:,:,:),
     &                   DIPSO2(:,:,:), S_SO2(:,:,:),
     &                   HCF2(:,:,:,:),           !dim,3,dim,dim),
     &                   SP_DIPO(:), SP_DIPW(:)  !3), SP_DIPW(3)
      Complex (kind=8) :: B(3,dim,-dim:dim),
     &                     C(  dim,-dim:dim),
     &                  BNMC(3,dim,0:dim),
     &                  BNMS(3,dim,0:dim),
     &                  CNMC(  dim,0:dim),
     &                  CNMS(  dim,0:dim), ES(0:2), FS(0:2),
     &               SP_HZFSO, SP_HZFSW, SP_MOW, CHECK_SGN, m_fact,
     &                  trace
      External trace, IsFreeUnit
!----------------------------------------------------------------------

      Call mma_allocate(ELOC,dim,'ELOC')
      Call mma_allocate(axes_in_abc,3,3,'axes_in_abc')

      Call mma_allocate(DIP_O,dim,dim,'DIP_O')
      Call mma_allocate(DIP_W,dim,dim,'DIP_W')
      Call mma_allocate(MUX,dim,dim,'MUX')
      Call mma_allocate(MUY,dim,dim,'MUY')
      Call mma_allocate(MUZ,dim,dim,'MUZ')
      Call mma_allocate(MUXZ,dim,dim,'MUXZ')
      Call mma_allocate(MUZX,dim,dim,'MUZX')
      Call mma_allocate(HZFS,dim,dim,'HZFS')
      Call mma_allocate(DIP_MOW,dim,dim,'DIP_MOW')
      Call mma_allocate(HZFS_MONM,dim,dim,'HZFS_MONM')
      Call mma_allocate(HZFS_MWNM,dim,dim,'HZFS_MWNM')
      Call mma_allocate(ZOUT,dim,dim,'ZOUT')
      Call mma_allocate(AMS,3,dim,dim,'AMS')
      Call mma_allocate(AMSSPIN,3,dim,dim,'AMSSPIN')
      Call mma_allocate(DIPSO2,3,dim,dim,'DIPSO2')
      Call mma_allocate(S_SO2,3,dim,dim,'S_SO2')
      Call mma_allocate(HCF2,dim,3,dim,dim,'HCF2')
      Call mma_allocate(SP_DIPO,3,'SP_DIPO')
      Call mma_allocate(SP_DIPW,3,'SP_DIPW')

!----------------------------------------------------------------------
      CALL atens(dipsom, dim, gtens, maxes, 2 )
c  save data for construction of the blocking barriers
      Do i=1,3
        Do j=1,3
          axes(iMLTPL,i,j)=0.0_wp
          axes(iMLTPL,i,j)=maxes(i,j)
        End Do
      End Do
c compute the magnetic axes in the crystalographic coordinate system, If requested:
      If(do_structure_abc) Then
        axes_in_abc=0.0_wp
        If(iprint>4) Then
          Write(6,'(A, 6F12.6)') 'cryst = ', (cryst(i),i=1,6)
          Write(6,'(A, 3F12.6)') 'coord = ', (coord(i),i=1,3)
        End If

        rc = 0
        CALL abc_axes(cryst, coord, maxes, axes_in_abc, 1, rc)

        Write(6,'(19x,32a,3x,a)') '|',('-',i=1,4),'|',
     &      ('-',i=1,5),' a ',('-',i=1,7),' b ',('-',i=1,7),' c ',
     &      ('-',i=1,3),'|','a , b , c  -- crystallographic axes'
        Write(6,'(A,F12.9,A,3F10.6,1x,A,16x,a)') ' gX = ',gtens(1),
     &      ' | Xm |',(axes_in_abc(j,1),j=1,3),'|',
     &      '(defined in the input)'
        Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gY = ',gtens(2),
     &      ' | Ym |',(axes_in_abc(j,2),j=1,3),'|'
        Write(6,'(A,F12.9,A,3F10.6,1x,A)') ' gZ = ',gtens(3),
     &      ' | Zm |',(axes_in_abc(j,3),j=1,3),'|'
        Write(6,'(83a)') ('-',i=1,56),'|'
      End If ! do_structure_abc
C Compute the matrix elements of the magnetic moment in the coordinate system
C of magnetic axes.  ==> I.e. ROTATE the matrix DipSO to the coordinate system of magnetic axes
c      maxes=0.0_wp
c      maxes(1,1)=1.0_wp
c      maxes(2,2)=1.0_wp
c      maxes(3,3)=1.0_wp
      CALL rotmom2( dipsom, dim, maxes, dipso2 )
      CALL rotmom2(  s_som, dim, maxes,  s_so2 )

      If(iprint>2) Then
        CALL prMom('G_HIGH_1:  DIPSO2(l,i,j):',dipso2,dim)
        CALL prMom('G_HIGH_1:   S_SO2(l,i,j):', s_so2,dim)
      End If

      Call zcopy_(  dim*dim,[(0.0_wp,0.0_wp)],0,ZOUT,1)
      Call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,AMS,1)
      Call zcopy_(3*dim*dim,[(0.0_wp,0.0_wp)],0,AMSSPIN,1)
      Call zcopy_(dim*3*dim*dim,[(0.0_wp,0.0_wp)],0,HCF2,1)

      CALL mu_order( dim,s_so2,dipso2,gtens,1,HCF2,AMS,AMSSPIN,ZOUT,
     &               iprint)

      check_sgn  =(0.0_wp,0.0_wp)
      check_sgn2 = 0.0_wp
      mux = (0.0_wp,0.0_wp)
      muy = (0.0_wp,0.0_wp)
      muz = (0.0_wp,0.0_wp)
      muxz= (0.0_wp,0.0_wp)
      muzx= (0.0_wp,0.0_wp)
      Do i=1,dim
        Do j=1,dim
          mux(i,j) = HCF2(1,1,i,j)
          muy(i,j) = HCF2(1,2,i,j)
          muz(i,j) = HCF2(1,3,i,j)
        End Do
      End Do
      CALL ZGEMM_( 'N','N', dim,dim,dim,(1.0_wp,0.0_wp),mux,dim,muz,dim,
     &            (0.0_wp,0.0_wp),muxz,dim )
      CALL ZGEMM_( 'N','N', dim,dim,dim,(1.0_wp,0.0_wp),muz,dim,mux,dim,
     &            (0.0_wp,0.0_wp),muzx,dim )

      If(ABS(muy(1,2)) > 1.d-25 ) Then
        check_sgn=(0.0_wp,-1.0_wp)*(muxz(1,2)-muzx(1,2))/muy(1,2)
        check_sgn2=DBLE(check_sgn)
      Else
        Write(6,'(A)') 'Is it an Ising Doublet?'
        Write(6,'(A)') 'For an Ising Doublet gX=gY=0, therefore, '//
     &                 'the product gX * gY * gZ is also zero'
      End If

      If (iprint>2) Then
        Write(6,'(5x,A,2F20.14)') 'check_sgn  = ', check_sgn
        Write(6,'(5x,A,F20.14)') 'check_sgn2 = ', check_sgn2
      End If

      Write(6,'(A,F11.6)') 'CHECK-SIGN parameter = ',check_sgn2
      If (check_sgn2<0.0_wp) Then
        Write(6,'(A,i2,a)') 'The sign of the product gX * gY * gZ '//
     &                      'for multiplet',iMLTPL,': < 0.'
      Else If(check_sgn2>0.0_wp) Then
        Write(6,'(A,i2,a,F9.6)') 'The sign of the product gX * gY'//
     &                           ' * gZ for multiplet',iMLTPL,': > 0.'
      End If
C  Obtain the b3m and c3m coefficients:
      B(1:3,1:dim,-dim:dim)=(0.0_wp,0.0_wp)
      Do N=1,dim-1,2
        Do M=0,N
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_O,1)
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_W,1)

          CALL Stewens_matrixel(N,M,dim,DIP_O,DIP_W,IPRINT)

          If(iprint>3) Then
            Write(6,*)
            Write(6,'( 5x,a,i2,a,i3)') 'DIP_O, N = ',N,', M =',m
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F10.6,2x))') (DIP_O(i,j), j=1,dim)
            End Do

            Write(6,*)
            Write(6,'( 5x,a,i2,a,i3)')  'DIP_W, N = ',N,', M =',m
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
            SP_DIPO(l)=trace(dim,AMS(l,1:dim,1:dim),DIP_O)
            SP_DIPW(l)=trace(dim,AMS(l,1:dim,1:dim),DIP_W)

            B(l,n,-m)=SP_DIPO(l)/SP_MOW
            B(l,n, m)=SP_DIPW(l)/SP_MOW
          End Do ! l
        End Do !m
      End Do !n

      BNMC(1:3,1:dim,0:dim)=(0.0_wp,0.0_wp)
      BNMS(1:3,1:dim,0:dim)=(0.0_wp,0.0_wp)
      Do n=1,dim-1,2
        Do m=0,N
          Do l=1,3
            If(M==0) Then
              BNMC(l,n,m)=(0.5_wp,0.0_wp)*(B(l,n,m)+B(l,n,-m))
            Else
              m_fact=dcmplx((-1)**M,0.0)
              BNMC(l,n,m)=   B(l,n,m) + m_fact*B(l,n,-m)
              BNMS(l,n,m)= ( B(l,n,m) - m_fact*B(l,n,-m) )*
     &                     (0.0_wp,-1.0_wp)
            End If
          End Do
        End Do
      End Do !n

      Write(6,*)
      Write(6,'(100A)') ('-',i=1,80)
      Write(6,'(A)') 'DECOMPOSITION OF THE MAGNETIC MOMENT Mu_i IN'//
     &               ' IRREDUCIBLE TENSOR OPERATORS (ITO):'
      Write(6,'(100A)') ('-',i=1,80)
      Write(6,*)
      Write(6,'(A)') 'The quantization axis is the main magnetic'//
     &               ' axis of this multiplet (Zm).'
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') '   Mu_i = SUM_{n,m}: [ B(i,n,m) * O(n,m)'//
     &               ' +  C(i,n,m) * W(n,m) ]'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m)'//
     &               ' + Y(n,-m) );'
      Write(6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m)'//
     &               ' - Y(n,-m) ) * I;   (I = imaginary unit)'
      Write(6,'(A)') '   n - the rank of the ITO, = 1, 3, 5, ...'//
     &               ' 2*spin;'
      Write(6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
      Write(6,'(A)') '   i - the Cartesian projection of the '//
     &               'magnetic moment, i = x,y,z;'
      Write(6,'(A)') 'These operators have been defined in: '
      Write(6,'(A)') '  L. F. Chibotaru, L.Ungur, J. Chem. Phys., '//
     &               '137, 064112 (2012).'
      Write(6,'(100A)') ('-',i=1,63),'|'
      Write(6,'(A)') '  n  |  m  | i |'//
     &               '        B(i,n,m)       |        C(i,n,m)       |'
      Do N=1,dim-1,2
        Write(6,'(A)') '-----|-----|---|-----------------------|'//
     &                 '-----------------------|'
        Do M=0,N
          If (M.ne.0)
     &      Write(6,'(A)') '     |-----|---|-----------------------|'//
     &                     '-----------------------|'
            Write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(E22.14,1x,A))')
     &                N,'|',M,'|','X','|', DBLE(BNMC(1,N,M)),'|',
     &                                     DBLE(BNMS(1,N,M)),'|'
            Write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(E22.14,1x,A))')
     &                N,'|',M,'|','Y','|', DBLE(BNMC(2,N,M)),'|',
     &                                     DBLE(BNMS(2,N,M)),'|'
            Write(6,'(2(1x,I2,2x,A),1x,A,1x,A,2(E22.14,1x,A))')
     &                N,'|',M,'|','Z','|', DBLE(BNMC(3,N,M)),'|',
     &                                     DBLE(BNMS(3,N,M)),'|'
        End Do
      End Do
      Write(6,'(100A)') ('-',i=1,63),'|'
c decomposition of the magnetic moment in Extended Stevens Operators

      CALL Set_knm( knm )

      Write(6,'(/)')
      Write(6,'(100A)') ('*',i=1,80)
      Write(6,'(A)') '   Mu_i = SUM_{k,q} * [ B(i,k,q) * O(k,q) ];'
      Write(6,'(A)') 'where:'
      Write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO)'//
     &               'as defined in:'
      Write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State '//
     &                   'Phys.,18(1985) 1415-1430.'
      Write(6,'(10x,A)') '2. Implemented in the "EasySpin" function '//
     &                   'in MATLAB, www.easyspin.org.'
      Write(6,'(A    )') '   k - the rank of the ITO, = 1, 3, 5, 7, '//
     &                   '9, 11;'
      Write(6,'(A    )') '   q - the component of the ITO, = -k, '//
     &                   '-k+1, ... 0, 1, ... k;'
      If((dim-1)>11) Then
        Write(6,'(A)') 'k = 11 may not be the highest rank of '//
     &                 'the ITO for this case, but it '
        Write(6,'(A)') 'is the maximal k implemented in the'//
     &                 ' "EasySpin" function in MATLAB.'
      End If
      Write(6,'(A)') 'Knm are proportionality coefficients between'//
     &               ' the ESO and operators defined in '
      Write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
      Write(6,'(100A)') ('-',i=1,51),'|'
      Write(6,'(A)') '  k |  q  | i |   (Knm)^2  |'//
     &               '         B(k,q)        |'

      If((dim-1)>11) Then
        Nmax=11
      Else
        nmax=dim-1
      End If

      Do N=1,nmax,2
        Write(6,'(A)') '----|-----|---|------------|'//
     &                 '-----------------------|'
        Do M=-N,N
          If (M<0) Then
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,'//
     &              '2(E22.14,1x,A))')
     &               N,'|',M,'|','X','|','|',
     &               DBLE(BNMS(1,N,ABS(M)))*knm(n,ABS(m)),'|'
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,'//
     &              'F10.3,1x,A,2(E22.14,1x,A))')
     &               N,'|',M,'|','Y','|',
     &               knm(n,ABS(m))*knm(n,ABS(m)),'|',
     &               DBLE(BNMS(2,N,ABS(M)))*knm(n,ABS(m)),'|'
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,'//
     &              '2(E22.14,1x,A))')
     &               N,'|',M,'|','Z','|','|',
     &               DBLE(BNMS(3,N,ABS(M)))*knm(n,ABS(m)),'|'
            Write(6,'(A)') '    |-----|---|------------|'//
     &                     '-----------------------|'
          Else
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,'//
     &              '2(E22.14,1x,A))')
     &               N,'|',M,'|','X','|','|',
     &               DBLE(BNMC(1,N,ABS(M)))*knm(n,ABS(m)),'|'
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,1x,'//
     &              'F10.3,1x,A,2(E22.14,1x,A))')
     &               N,'|',M,'|','Y','|',
     &               knm(n,ABS(m))*knm(n,ABS(m)),'|',
     &               DBLE(BNMC(2,N,ABS(M)))*knm(n,ABS(m)),'|'
            Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),1x,A,1x,A,12x,A,'//
     &              '2(E22.14,1x,A))')
     &               N,'|',M,'|','Z','|','|',
     &               DBLE(BNMC(3,N,ABS(M)))*knm(n,ABS(m)),'|'
             If (M.ne.N)
     &            Write(6,'(A)') '    |-----|---|------------|'//
     &                           '-----------------------|'
          End If !M<0
        End Do !M
      End Do !N
      Write(6,'(100A)') ('-',i=1,51),'|'

C  Calculation of the ZFS tensors and the coefficients of the higher order spin-operators Enm and Fnm
      If(dim>2) Then
        ESUM=0.0_wp
        Do I=1,dim
          ELOC(i)=0.0_wp
          ESUM=ESUM+ESOM(I)
        End Do
        E0=ESUM/DBLE(dim)
        Do I=1,dim
          ELOC(I)=ESOM(I)-E0
        End Do

        Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS,1)
        Do i1=1,dim
          Do i2=1,dim
            Do i=1,dim
              HZFS(i1,i2)=HZFS(i1,i2) + ELOC(i)
     &                                * conjg( ZOUT(i,I1) )
     &                                *        ZOUT(i,I2)
            End Do
          End Do
        End Do

        If(iprint>2) Then
          Write(6,*)
          Write(6,'(5X,A)') 'SPIN-ORBIT ENERGIES OF THE FIRST '//
     &                      'MOLECULAR TERM SHIfTED TO THE MASS CENTER'
          Write(6,*)
          Write(6,'(5X,A,F10.4)') 'E0 = ',E0
          Write(6,'(15X,A,11X,A)') 'ESOM', 'ESO_LOC'

          Do I=1,dim
            Write(6,'(5X,F15.6,2X,F15.6)') ESOM(I), ELOC(I)
          End Do
        End If

        Write(6,*)
        Write(6,'(100A)') ('-',i=1,87)
        Write(6,'(A)') 'DECOMPOSITION OF THE ZERO-FIELD SPLITTING '//
     &                 '(ZFS) IN IRREDUCIBLE TENSOR OPERATORS (ITO):'
        Write(6,'(100A)') ('-',i=1,87)
        Write(6,*)
        Write(6,'(A)') 'Ab Initio Calculated Zero-Field Splitting '//
     &                 'Matrix written in the basis of Pseudospin '//
     &                 'Eigenfunctions'
        If (MOD(dim,2)==0) Then
          Write(6,'(950A)') ('-',i=1,10),(('-',i=1,24),j=1,dim),'|'
          Write(6,'(10x,A,50(8x,A,I3,A,7x,A))') '|',
     &      ('|',2*i-dim-1,'/2 >','|',i=1,dim)
          Write(6,'(950A)') ('-',i=1,10),'|',(('-',i=1,23),'|',j=1,dim)
          Do i=1,dim
            Write(6,'(1x,A,I3,A,1x,A,50(2F11.5,1x,A))')
     &              '<',2*i-dim-1,'/2','| |', ( HZFS(j,i),'|' ,j=1,dim)
          End Do
          Write(6,'(950A)') ('-',i=1,10),(('-',i=1,24),j=1,dim),'|'
        Else
          Write(6,'(950A)') ('-',i=1,8),(('-',i=1,24),j=1,dim),'|'
          Write(6,'(8x,A,50(8x,A,I3,A,9x,A))')
     &                     '|',('|',-(dim-1)/2-1+i,' >','|',i=1,dim)
          Write(6,'(950A)') ('-',i=1,8),'|',(('-',i=1,23),'|',j=1,dim)
            Do I=1,dim
              Write(6,'(1x,A,I3,1x,A,50(2F11.5,1x,A))')
     &             '<',-(dim-1)/2-1+i,'| |',  ( HZFS(j,i),'|' ,j=1,dim)
            End Do
          Write(6,'(950A)') ('-',i=1,8),(('-',i=1,24),j=1,dim),'|'
        End If

      C(1:dim,-dim:dim)=(0.0_wp,0.0_wp)
      Do N=2,dim-1,2
        Do M=0,N
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_O,1)
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_W,1)
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS_MONM,1)
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,HZFS_MWNM,1)
          Call zcopy_(dim*dim,[(0.0_wp,0.0_wp)],0,DIP_MOW,1)

          CALL Stewens_matrixel(N,M,dim,DIP_O,DIP_W,IPRINT)

          If(IPRINT>5) Then
            Write(6,'(/)')
            Write(6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_O  N = ',N,'M =',M
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2X,2E20.10))') (DIP_O(i,j), j=1,dim)
            End Do
            Write(6,'(/)')
            Write(6,'(5x,a,i3,3x,A,I3)') 'DIP_STEWENS_W  N = ',N,'M =',M
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2X,2E20.10))') (DIP_W(i,j), j=1,dim)
            End Do
          End If

          SP_HZFSO =(0.0_wp,0.0_wp)
          SP_HZFSW =(0.0_wp,0.0_wp)
          SP_MOW   =(0.0_wp,0.0_wp)
          SP_MOW   = trace(dim,DIP_O,DIP_W)
          SP_HZFSO = trace(dim, HZFS,DIP_O)
          SP_HZFSW = trace(dim, HZFS,DIP_W)

          C(N,-M)=SP_HZFSO/SP_MOW
          C(N, M)=SP_HZFSW/SP_MOW

          If(IPRINT>5) Then
            Write(6,'(/)')
            Write(6,'( 5x,a)')  'HZFS_MONM(i,j)'
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F18.10,2x))') (HZFS_MONM(i,j), j=1,dim)
            End Do
            Write(6,*)
            Write(6,'(5X,a,2F18.10)') 'SP_HZFSO = ', SP_HZFSO
            Write(6,'(/)')
            Write(6,'( 5x,a)')  'HZFS_MWNM(i,j)'
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F18.10,2x))') (HZFS_MWNM(i,j), j=1,dim)
            End Do
            Write(6,*)
            Write(6,'(5X,a,2F18.10)') 'SP_HZFSW = ', SP_HZFSW
            Write(6,'(/)')
            Write(6,'( 5x,a)')  'HZFS_MOW(i,j)(i,j)'
            Write(6,*)
            Do i=1,dim
              Write(6,'(20(2F18.10,2x))') (DIP_MOW(i,j), j=1,dim)
            End Do
            Write(6,*)
            Write(6,'(5X,a,2F18.10)') 'SP_MOW = ', SP_MOW
          End If
        End Do !M
      End Do !N

      CNMC(1:dim,0:dim)=(0.0_wp,0.0_wp)
      CNMS(1:dim,0:dim)=(0.0_wp,0.0_wp)
      Do N=2, dim-1,2
        Do M=0,N
          If(M==0) Then
            CNMC(N,M)=(0.5_wp,0.0_wp)*(C(N,M)+C(N,-M))
          Else
            m_fact=cmplx((-1)**M,0,wp)
            CNMC(N,M)=  C(N,M) + m_fact * C(N,-M)
            CNMS(N,M)=( C(N,M) - m_fact * C(N,-M) )*(0.0_wp,-1.0_wp)
          End If
        End Do
      End Do

        Write(6,'(A)') 'The ZFS Hamiltonian:'
        Write(6,'(A)') '   ZFS = SUM_{n,m}: [ E(n,m) * O(n,m)'//
     &                 ' +  F(n,m) * W(n,m) ]'
        Write(6,'(A)') 'where:'
        Write(6,'(A)') '   O(n,m) =  0.5 * ( (-1)**m * Y(n,+m)'//
     &                 ' + Y(n,-m) );'
        Write(6,'(A)') '   W(n,m) = -0.5 * ( (-1)**m * Y(n,+m)'//
     &                 ' - Y(n,-m) ) * I;    (I = imaginary unit)'
        Write(6,'(A)') '   n - the rank of the ITO, = 2, 4, 6, ... '//
     &                 '2*spin;'
        Write(6,'(A)') '   m - the component of the ITO, = 0, 1, ... n;'
        Write(6,'(A)') 'The quantization axis is the main magnetic'//
     &                 ' axis of this multiplet (Zm).'
        Write(6,'(100A)') ('-',i=1,59),'|'
        Write(6,'(A)') '  n  |  m  |         E(n,m)        |'//
     &                 '         F(n,m)        |'
        Do N=2, dim-1,2
          Write(6,'(A)') '-----|-----|-----------------------|'//
     &                   '-----------------------|'
          Do M=0,N
            Write(6,'(2(1x,I2,2x,A),2(E22.14,1x,A))') N,'|',M,'|',
     &                    DBLE(CNMC(N,M)),'|',DBLE(CNMS(N,M)),'|'
          End Do
        End Do
        Write(6,'(100A)') ('-',i=1,59),'|'
c
c decomposition of the ZFS matrix in ExtEnded Stevens Operators
        Write(6,'(100A)') ('*',i=1,80)
        Write(6,'(A)') 'The ZFS Hamiltonian:'
        Write(6,'(A)') '   ZFS = SUM_{k,q} * [ B(k,q) * O(k,q) ];'
        Write(6,'(A)') 'where:'
        Write(6,'(A)') '   O(k,q) =  Extended Stevens Operators (ESO)'//
     &                 'as defined in:'
        Write(6,'(10x,A)') '1. Rudowicz, C.; J.Phys.C: Solid State '//
     &                     'Phys.,18(1985) 1415-1430.'
        Write(6,'(10x,A)') '2. Implemented in the "EasySpin" '//
     &                     'function in MATLAB, www.easyspin.org.'
        Write(6,'(A    )') '   k - the rank of the ITO, = 2, 4, 6, '//
     &                     '8, 10, 12.'
        Write(6,'(A)')     '   q - the component of the ITO, = -k, '//
     &                     '-k+1, ... 0, 1, ... k;'

        If((dim-1)>12) Then
          Write(6,'(A)') 'k = 12 may not be the highest rank of '//
     &                   'the ITO for this case, but it '
          Write(6,'(A)') 'is the maximal k implemented in the'//
     &                   ' "EasySpin" function in MATLAB.'
        End If

        Write(6,'(A)') 'Knm are proportionality coefficients between'//
     &                 ' the ESO and operators defined in '
        Write(6,'(A)') 'J. Chem. Phys., 137, 064112 (2012).'
        Write(6,'(100A)') ('-',i=1,48),'|'
        Write(6,'(A)') '  k |  q  |    (Knm)^2  |'//
     &                 '         B(k,q)        |'
        If((dim-1)>12) Then
          Nmax=12
        Else
          Nmax=dim-1
        End If
        Do N=2,Nmax,2
          Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'
          Do M=-N,N
            If (M<0) Then
              Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),'//
     &                                 'F11.2,2x,A,2(E22.14,1x,A))')
     &                 N,'|',M,'|',knm(n,ABS(m))*knm(n,ABS(m)),'|',
     &                 DBLE(CNMS(N,ABS(M)))*knm(n,ABS(m)),'|'
            Else
              Write(6,'((1x,I2,1x,A),(1x,I3,1x,A),'//
     &                                 'F11.2,2x,A,2(E22.14,1x,A))')
     &                 N,'|',M,'|',knm(n,ABS(m))*knm(n,ABS(m)),'|',
     &                 DBLE(CNMC(N,ABS(M)))*Knm(n,ABS(m)),'|'
            End If
          End Do
        End Do
        Write(6,'(100A)') ('-',i=1,48),'|'
      !-----------------------------
      ! for the interface related to CF gradient calculation:
      If (GRAD) Then
         LuDgrad=IsFreeUnit(81)
         Call molcas_open(LuDgrad,'DMAT')
         Do N=2,Nmax,2
           Write(6,'(A)') '----|-----|-------------|'//
     &                    '-----------------------|'
           Do M=-N,N
             If (M.lt.0) Then
               Write(LuDgrad,'(I4,I4,1x,2(E25.15))')
     &                   N,M, DBLE(CNMS(n,ABS(n)))*knm(n,ABS(m))
             Else
               Write(LuDgrad,'(I4,I4,1x,2(E25.15))')
     &                   N,M, DBLE(CNMC(n,ABS(m)))*Knm(n,ABS(m))
             End If
           End Do
         End Do
         Close(LuDgrad)
      End If
      !-----------------------------
c
        Do i=0,2
          ES(i)=0.0_wp
          FS(i)=0.0_wp
          ES(i)=CNMC(2,i)
          FS(i)=CNMS(2,i)
        End Do
        CALL DMATRIX(Es,Fs,maxes,2)
      End If ! decomposition of ZFS in higher-order ITO operators

      Write(6,*)
      Write(6,'(A,I2,A,I3,A)') 'ANGULAR MOMENTS ALONG THE MAIN'//
     &                         ' MAGNETIC AXES'
      CALL moments(dim,s_so2,dipso2,iprint)


!----------------------------------------------------------------------
      Call mma_deallocate(ELOC)
      Call mma_deallocate(axes_in_abc)
      Call mma_deallocate(DIP_O)
      Call mma_deallocate(DIP_W)
      Call mma_deallocate(MUX)
      Call mma_deallocate(MUY)
      Call mma_deallocate(MUZ)
      Call mma_deallocate(MUXZ)
      Call mma_deallocate(MUZX)
      Call mma_deallocate(HZFS)
      Call mma_deallocate(DIP_MOW)
      Call mma_deallocate(HZFS_MONM)
      Call mma_deallocate(HZFS_MWNM)
      Call mma_deallocate(ZOUT)
      Call mma_deallocate(AMS)
      Call mma_deallocate(AMSSPIN)
      Call mma_deallocate(DIPSO2)
      Call mma_deallocate(S_SO2)
      Call mma_deallocate(HCF2)
      Call mma_deallocate(SP_DIPO)
      Call mma_deallocate(SP_DIPW)



      Return
      End
