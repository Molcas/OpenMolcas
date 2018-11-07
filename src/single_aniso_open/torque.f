      Subroutine torque(Nss,NM,AngPoints,EM,eso,dipm,sm,zJ,thrs,mem,
     &                  m_paranoid,smagn,H_torq,T_torq,ma,dbg)

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "cntrl.fh"
#include "mgrid.fh"
#include "mvect.fh"
#include "stdalloc.fh"
      Integer, intent(in)          :: nss, nm
      Integer, intent(in)          :: AngPoints
      Integer, intent(in)          :: mem
      Logical, intent(in)          :: smagn
      Logical, intent(in)          :: m_paranoid
c ab initio data:
      Real(kind=wp), intent(in)    :: ESO(nss)  ! exchange energies printed out in the previous part
      Real(kind=wp), intent(in)    :: EM, zJ, thrs
      Real(kind=wp), intent(in)    :: H_torq, T_torq
      Real(kind=wp), intent(in)    :: ma(3,3) ! main magnetic axes
      Complex(kind=wp), intent(in) :: DIPM(3,nss,nss)
      Complex(kind=wp), intent(in) ::   SM(3,nss,nss)

! local data:
c magnetic field strength and orientation data:
      Integer       :: nPlanes
      Parameter (nPlanes=1)
!      Real(kind=wp) :: dlth
      Real(kind=wp), allocatable :: W(:)    ! W(NM) ! Zeeman exchange energies
      Real(kind=wp), allocatable :: ZT    ! ZT(nTempMagn) ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST(:) ! ST(3) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT(:) ! MT(3) ! total magnetisation
      Real(kind=wp), allocatable :: dX(:) ! dX(AngPoints)
      Real(kind=wp), allocatable :: dY(:) ! dY(AngPoints)
      Real(kind=wp), allocatable :: dZ(:) ! dZ(AngPoints)
      Real(kind=wp), allocatable :: Ang(:)  ! Ang(AngPoints)
      Complex(kind=wp), allocatable :: M(:,:,:)
      Complex(kind=wp), allocatable :: S(:,:,:)
c magnetic torque
      !Real(kind=wp), allocatable :: tx(:,:) ! tx(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
      Real(kind=wp), allocatable :: ty(:) ! ty(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
      !Real(kind=wp), allocatable :: tz(:,:) ! tz(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque
c magnetic and spin moments (i.e. the BIG matrices):
      Character(len=99):: STLNE1, STLNE2
      Real(kind=wp)    :: cm3tomB, g(3),mg(3,3)!,ma_inv(3,3)!,det
      Real(kind=wp)    :: AngStep,AngRad,pi
      Logical          :: DBG
      Integer          :: IM,I,J,l,mem_local,RtoB,CtoB,nH_torq,nT_torq

c      Boltz_k=0.6950356000_wp   !   in cm^-1*K-1
c      mu_Bohr=0.4668643740_wp   !   in cm-1*T-1
      cm3tomB=0.5584938904_wp   !   in cm3 * mol-1 * T
      pi=3.1415926535897932384626433832795028841971_wp

      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,'(20X,A)')
     & 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Write(6,'(2X,A,i3,A)') 'Magnetization torque is '//
c     & 'calculated for the ',NH,' field points, in the field domain:'
c      Write(6,'(2X,F4.1,1x,a,1X,F4.1,a,10(F6.3,a))') HMIN,'--',HMAX,
c     & ' T., at the following temperatures:',
c     & (TempMagn(l),' K.;',l=1,nTempMagn)
c      Do i=11,nTempMagn,10
c        j=MIN(nTempMagn,i+9)
c        Write(6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
c      End Do
      Write(6,'(2X,A,F10.5,A,F10.5)') 'Magnetization torque is '//
     & 'calculated for one field point,',H_torq,
     & ' T., at the following temperature:', T_torq

      Write(6,'(2x,A,i3,A)')
     & 'Angular dependence of the magnetization torque '//
     & 'is computed for ',AngPoints,' angular points distributed'
      Write(6,'(2x,A)') 'in the domain 0-180 deg.'
      Write(6,'(2X,A)') 'The cut-off energy for the '//
     & 'exact diagonalization of the Zeeman Hamiltonian is:'
      Write(6,'(2x,a,F15.9,A)') 'E = ',EM ,' cm(-1).'
      If(NM.lt.10) Then
        Write(6,'(2X,A,i2,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.10) .AND. (NM.lt.100) ) Then
        Write(6,'(2X,A,i3,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.100) .AND. (NM.lt.1000) ) Then
        Write(6,'(2X,A,i4,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.1000) .AND. (NM.lt.10000) ) Then
        Write(6,'(2X,A,i5,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      End If
      nH_torq=1
      nT_torq=1

c      If (dbg) Then
        Do i=1,3
           Write(6,*) (ma(j,i),j=1,3)
        End Do
c      End If
!-----------------------------------------------------------------------
! Allocate memory for this computation:
      mem_local=0
      RtoB=8
      CtoB=16
      If(nM>0) Then
         ! Zeeman exchange energy spectrum
         Call mma_allocate(W,nM,'W')
         Call dcopy_(nM,0.0_wp,0,W,1)
         mem_local=mem_local+nM*RtoB
      End If

      Call mma_allocate(ST,3,'ST')
      Call dcopy_(3,0.0_wp,0,ST,1)
      mem_local=mem_local+3*RtoB

      Call mma_allocate(MT,3,'MT')
      Call dcopy_(3,0.0_wp,0,MT,1)
      mem_local=mem_local+3*RtoB

      If(AngPoints>0) Then
         Call mma_allocate(dX,AngPoints,'dX')
         Call mma_allocate(dY,AngPoints,'dY')
         Call mma_allocate(dZ,AngPoints,'dZ')
         Call dcopy_(AngPoints,0.0_wp,0,dX,1)
         Call dcopy_(AngPoints,0.0_wp,0,dY,1)
         Call dcopy_(AngPoints,0.0_wp,0,dZ,1)
         mem_local=mem_local+3*AngPoints*RtoB

         Call mma_allocate(Ang,AngPoints,'Ang')
         Call dcopy_(AngPoints,0.0_wp,0,Ang,1)
         mem_local=mem_local+AngPoints*RtoB
      End If

      If( (nPlanes>0).and.(AngPoints>0) ) Then
         Call mma_allocate(ty,AngPoints,'ty')
         Call dcopy_(AngPoints,0.0_wp,0,ty,1)
         mem_local=mem_local+AngPoints*RtoB
      End If

      If (nss>0) Then
         Call mma_allocate(M,3,nss,nss,'Mrot')
         Call mma_allocate(S,3,nss,nss,'Srot')
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,M,1)
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,S,1)
         mem_local=mem_local+2*3*nss*nss*CtoB
      End If


      If(dbg) Write(6,*) 'TORQ:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'TORQ:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local
!-----------------------------------------------------------------------
          ! rotate the moments to the coordiante system of
          ! the ground state
         ! ma_inv=0.0_wp
         ! Call REVERSE(ma,ma_inv,DET)
          Call rotmom2( DIPM, nss, ma, M )
          Call rotmom2(   SM, nss, ma, S )

          g=0.0_wp
          mg=0.0_wp
          Call atens( M(1:3,1:2,1:2), 2, g, mg, 2)
!-----------------------------------------------------------------------
      Call dcopy_(AngPoints,0.0_wp,0,dX,1)
      Call dcopy_(AngPoints,0.0_wp,0,dY,1)
      Call dcopy_(AngPoints,0.0_wp,0,dZ,1)
      Call dcopy_(AngPoints,0.0_wp,0,Ang,1)
      !Call hdir2(AngPoints,2,dX(:),dY(:),dZ(:),Ang,2)
      AngStep=0.0_wp
      AngRad =0.0_wp
      AngStep=360.0_wp/dble(AngPoints-1)
      dX(1)=1.0_wp
      dZ(1)=0.0_wp
      Do i=1,AngPoints
        AngRad=dble(i-1)*AngStep*Pi/180.0_wp!+122.625_wp*Pi/180.0_wp
        Ang(i)=dble(i-1)*AngStep
         dX(i)=cos(AngRad)
         dZ(i)=sin(AngRad)
      End Do
      If(dbg) Then
        Write(6,'(A,I5)') 'Angular grid for Magnetization Torque, '//
     &                    'Cartesian Component =',L
        Write(6,'(2x,A,4x,A,5x,3(10X,A,10x))') 'Nr.','Angle','X','Y','Z'
        Do i=1,AngPoints
          Write(6,'(I4,F10.3,3x,3F21.14)') i,Ang(i),dX(i),dY(i),dZ(i)
        End Do
      End If
!-----------------------------------------------------------------------


      Do IM=1,AngPoints
         WRITE(STLNE1,'(A   )') 'SINGLE_ANISO:  torque:'
         WRITE(STLNE2,'(A,I3)') 'Magnetization at point ',IM
         Call StatusLine(STLNE1,STLNE2)
         ZT=0.0_wp
         Call dcopy_(nM,0.0_wp,0,W,1)
         Call dcopy_(3,0.0_wp,0,MT,1)
         Call dcopy_(3,0.0_wp,0,ST,1)
c exchange magnetization:
         Call MAGN( nss, NM, dX(iM), dY(iM), dZ(iM),
     &              H_torq,
     &              eso, zJ, THRS,
     &              M, !DIPM,
     &              S, !  SM,
     &              nT_torq, T_torq, smagn,
     &              W,
     &              ZT,
     &              ST,
     &              MT, m_paranoid , DBG )
        If(dbg) Write(6,'(A,3F18.10)') 'TORQ: MT=',MT(1),MT(2),MT(3)

       !  tx(iPl,iM) = ( MT(2)*dZ(iM) - MT(3)*dY(iM) ) *H_torq
         ty(iM) = ( MT(3)*dX(iM) - MT(1)*dZ(iM) ) * H_torq
       !  tz(iPl,iM) = ( MT(1)*dY(iM) - MT(2)*dX(iM) ) *H_torq
      End Do ! iM
c -------------------------------------------------------------------
C   WRITING SOME OF THE OUTPUT....
C -------------------------------------------------------------------
      Write(6,*)
      Write(6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION '//
     &                   'TORQUE'
      Write(6,'(30X,A)') '(Units of torque: [energy, cm-1])'
      Write(6,*)

        Write(6,'(5x,A)') 'Orientation of the applied magnetic'//
     & ' field employed:'
        Write(6,'(10A)') '--------|','---------------------------|'

        Write(6,'(2x,A,10x,A)') 'Angle |',
     &                             'rotation in the XZ plane          |'
        Write(6,'(10A)') '--------|',
     &             '--- proj X ---|','--- proj Y ---|','--- proj Z ---|'
         Do iM=1,AngPoints
           Write(6,'(F7.3,1x,A,F13.10,1x,A,F13.10,1x,A,F13.10,1x,'//
     &               'A,F20.14,1x,A)')
     &        Ang(iM),'|',dX(iM),' ',dY(iM),' ',dZ(iM),'|',ty(iM),'|'
           End Do
        Write(6,'(10A)') '--------|','---------------------------|'

        Write(6,*)
        Write(6,'(10A)') '--------|','---------------------------|'
        Write(6,'(A,F9.4,A)') 'Magnetic field strength = ',H_torq,
     &                        ' Tesla'
        Write(6,'(12x,A,F9.4,A)') 'Temperature = ',T_torq,' Kelvin'
        Write(6,'(10A)') '--------|','---------------------------|'

        Write(6,'(2x,A,3(10x,A))') 'Angle |',
     &                             'rotation in the XZ plane          |'
        Write(6,'(10A)') '--------|', '-- torque along Y --|'
         Do iM=1,AngPoints
            Write(6,'(F7.3,1x,A,E18.10)') Ang(iM),'|',ty(iM)
         End Do
        Write(6,'(10A)') '--------|','---------------------------|'

!-----------------------------------------------------------------------
! deallocate memory for this computation:
      If(nM>0) Then
         Call mma_deallocate(W)
      End If

      Call mma_deallocate(ST)
      Call mma_deallocate(MT)

      If(AngPoints>0) Then
         Call mma_deallocate(dX)
         Call mma_deallocate(dY)
         Call mma_deallocate(dZ)
         Call mma_deallocate(Ang)
      End If

      If( (nPlanes>0).and.(AngPoints>0) ) Then
         Call mma_deallocate(ty)
      End If
      If(nss>0) Then
         Call mma_deallocate(m)
         Call mma_deallocate(s)
      End If

      If(dbg) Write(6,*) 'TORQ: allocated memory was sucessfully '//
     &                   'deallocated'

      Return
      End
