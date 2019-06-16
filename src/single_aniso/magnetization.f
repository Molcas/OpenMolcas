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
      Subroutine magnetization( nss, nM, nTempMagn, nDirTot, nDir,
     &                          nDirZee, nH, iPrint, LUZee, mem,
     &                          compute_Mdir_vector, zeeman_energy,
     &                          hinput, m_paranoid, smagn, doplot,
     &                          TempMagn, eso, dirX, dirY, dirZ,
     &                          dir_weight, hexp, magn_exp, zJ, hmin,
     &                          hmax, EM, thrs,
     &                          dipm, sm, dbg )
************************************************************************
*                                                                      *
*     MAGNETIZATION control section                                    *
*                                                                      *
*     calling arguments:                                               *
*     NSS     : number of spin-orbit states (total)                    *
*               scalar integer                                         *
*     NM      : size of the Zeeman Hamiltonian matrix                  *
*               scalar integer                                         *
*     EM      : cut-off energy (energy of the last s-o state which is  *
*               included in the Zeeman matrix                          *
*               scalar real*8                                          *
*     EM      : cut-off energy (energy of the last s-o state which is  *
*               included in the Zeeman matrix                          *
*               scalar real*8                                          *
*     IFINAL  : integer                                                *
*               termination flag                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     Liviu Ungur                                                      *
*     University of Leuven, Belgium, 2008-2017                         *
*                                                                      *
*----------------------------------------------------------------------*
*     Hystory:                                                         *
*     Liviu Ungur, 2008-2017 various modifications                     *
************************************************************************

      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"
#include "stdalloc.fh"
!----------------------------------------------------------------
      ! input data (30)
      Integer, intent(in)          :: nDir
      Integer, intent(in)          :: nDirZee
      Integer, intent(in)          :: nDirTot
      Integer, intent(in)          :: nss
      Integer, intent(in)          :: nM
      Integer, intent(in)          :: nH
      Integer, intent(in)          :: iprint
      Integer, intent(in)          :: nTempMagn
      Integer, intent(in)          :: mem
      Integer, intent(in)          :: LUZee(nDirZee)

      Logical, intent(in)          :: compute_Mdir_vector
      Logical, intent(in)          :: zeeman_energy
      Logical, intent(in)          :: DoPlot
      Logical, intent(in)          :: hinput
      Logical, intent(in)          :: smagn
      Logical, intent(in)          :: m_paranoid
      Logical, intent(in)          :: dbg

      Real(kind=wp), intent(in)    :: dirX(nDir), dirY(nDir), dirZ(nDir)
      Real(kind=wp), intent(in)    :: dir_weight(nDirZee,3)
      Real(kind=wp), intent(in)    :: hmin, hmax
      Real(kind=wp), intent(in)    :: zj, thrs
      Real(kind=wp), intent(in)    :: eso(nss)
      Real(kind=wp), intent(in)    :: EM
      Real(kind=wp), intent(in)    :: TempMagn(nTempMagn)
      Real(kind=wp), intent(in)    :: hexp(nH)
      Real(kind=wp), intent(in)    :: magn_exp(nH,nTempMagn)

      Complex(kind=wp), intent(in) ::   sm(3,nss,nss)
      Complex(kind=wp), intent(in) :: dipm(3,nss,nss)
!----------------------------------------------------------------
!     local variables
      Integer       :: nP,iTEnd,iT,IM,I,L,J,IC,IDIR,IH,iTemp
      Real(kind=wp) :: DLTH, mv, sv, dev, Boltz_k,mu_Bohr
      Character(len=99):: STLNE1, STLNE2

      Real(kind=wp), allocatable :: WM(:)         ! WM(nm)
      Real(kind=wp), allocatable :: MT(:,:,:)     ! MT(3,nH,nTempMagn)
      Real(kind=wp), allocatable :: ST(:,:,:)     ! ST(3,nH,nTempMagn)
      ! magnetization and spin vectors
      Real(kind=wp), allocatable :: MVEC(:,:,:,:)
!                                   MVEC(nDirTot,nH,nTempMagn,3)
      Real(kind=wp), allocatable :: SVEC(:,:,:,:)
!                                   SVEC(nDirTot,nH,nTempMagn,3)
      Real(kind=wp), allocatable :: H(:)          ! H(nH)
      ! average powder M and S:
      Real(kind=wp), allocatable :: MAV(:,:) ! MAV(nH,nTempMagn)
      Real(kind=wp), allocatable :: SAV(:,:) ! SAV(nH,nTempMagn)

      Real(kind=wp), allocatable :: ZT(:,:)  ! ZT(nH,nTempMagn)
      Real(kind=wp), allocatable :: STDEV(:) ! STDEV(nTempMagn)
      Real(kind=wp), allocatable :: dHX(:)   ! dHX(nDirTot)
      Real(kind=wp), allocatable :: dHY(:)   ! dHY(nDirTot)
      Real(kind=wp), allocatable :: dHZ(:)   ! dHZ(nDirTot)
      Real(kind=wp), allocatable :: dHW(:)   ! dHW(nDirTot)

      Integer          :: mem_local, RtoB

      External      :: dev

      Call qEnter('SA_magn')
      Boltz_k=0.6950356_wp                    !   in cm-1*K-1
      mu_Bohr=0.466864374_wp                  !   in cm-1*T-1
!-----------------------------------------------------------------------
! Allocate necessary memory
      mem_local=0
      RtoB=8

      If(nM>=0) Then
         ! Zeeman exchange energy spectrum
         Call mma_allocate(WM,nM,'W')
         Call dcopy_(nM,[0.0_wp],0,WM,1)
         mem_local=mem_local+nM*RtoB
      End If

      If((nH>=0).and.(nTempMagn>=0)) Then
         Call mma_allocate(MT,3,nH,nTempMagn,'MT')
         Call dcopy_(3*nH*nTempMagn,[0.0_wp],0,MT,1)
         mem_local=mem_local+3*nH*nTempMagn*RtoB

         Call mma_allocate(ST,3,nH,nTempMagn,'ST')
         Call dcopy_(3*nH*nTempMagn,[0.0_wp],0,ST,1)
         mem_local=mem_local+3*nH*nTempMagn*RtoB

         Call mma_allocate(MAV,nH,nTempMagn,'MAV')
         Call dcopy_(nH*nTempMagn,[0.0_wp],0,MAV,1)
         mem_local=mem_local+nH*nTempMagn*RtoB

         Call mma_allocate(SAV,nH,nTempMagn,'SAV')
         Call dcopy_(nH*nTempMagn,[0.0_wp],0,SAV,1)
         mem_local=mem_local+nH*nTempMagn*RtoB

         Call mma_allocate(ZT,nH,nTempMagn,'ZT')
         Call dcopy_(nH*nTempMagn,[0.0_wp],0,ZT,1)
         mem_local=mem_local+nH*nTempMagn*RtoB

         If(nDirTot>=0) Then
            Call mma_allocate(MVEC,nDirTot,nH,nTempMagn,3,'MVEC')
            Call mma_allocate(SVEC,nDirTot,nH,nTempMagn,3,'SVEC')
            Call dcopy_(3*nDirTot*nH*nTempMagn,[0.0_wp],0,MVEC,1)
            Call dcopy_(3*nDirTot*nH*nTempMagn,[0.0_wp],0,SVEC,1)
            mem_local=mem_local+6*nDirTot*nH*nTempMagn*RtoB
         End If
      End If

      If(nH>=0) Then
         Call mma_allocate(H,nH,'H')
         Call dcopy_(nH,[0.0_wp],0,H,1)
         mem_local=mem_local+nH*RtoB
      End If

      If((nTempMagn>=0).and.hinput) Then
         Call mma_allocate(STDEV,nTempMagn,'H')
         Call dcopy_(nTempMagn,[0.0_wp],0,STDEV,1)
         mem_local=mem_local+nTempMagn*RtoB
      End If

      If(nDirTot>=0) Then
            Call mma_allocate(dHX,nDirTot,'dHX')
            Call mma_allocate(dHY,nDirTot,'dHY')
            Call mma_allocate(dHZ,nDirTot,'dHZ')
            Call mma_allocate(dHW,nDirTot,'dHW')
            Call dcopy_(nDirTot,[0.0_wp],0,dHX,1)
            Call dcopy_(nDirTot,[0.0_wp],0,dHY,1)
            Call dcopy_(nDirTot,[0.0_wp],0,dHZ,1)
            Call dcopy_(nDirTot,[0.0_wp],0,dHW,1)
            mem_local=mem_local+4*nDirTot*RtoB
      End If
      If(dbg) Write(6,*) 'MAGNETIZATION:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'MAGNETIZATION:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local
!-----------------------------------------------------------------------
      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      If(DBG.or.iprint>3) Then
        Write(6,'(A, I5)')  'nH                  = ',nH
        Write(6,'(A, I5)')  'nM                  = ',nM
        Write(6,'(A, I5)')  'nsymm               = ',nsymm
        Write(6,'(A, I5)')  'ngrid               = ',ngrid
        Write(6,'(A, I5)')  'nDir                = ',nDir
        Write(6,'(A, I5)')  'nDirZee             = ',nDirZee
        Write(6,'(A, I5)')  'nDirTot             = ',nDirTot
        Write(6,'(A, F9.5)')'HMIN                = ',hmin
        Write(6,'(A, F9.5)')'HMAX                = ',hmax
        Write(6,'(A, F9.5)')'zJ                  = ',zJ
        Write(6,*)          'compute_Mdir_vector = ',compute_Mdir_vector
        Write(6,*)          'hinput              = ',hinput
        Write(6,*)          'zeeman_energy       = ',zeeman_energy
        Write(6,*)          'hinput              = ',hinput
        Write(6,*)          'smagn               = ',smagn
        Write(6,'(A)') 'dir_weight'
        Do i=1,nDirZee
          Write(6,'(3F10.6)') (dir_weight(i,j),j=1,3)
        End Do
        Write(6,'(A)') 'nDir'
        Do i=1,nDir
          Write(6,'(3F10.6)') dirX(i),dirY(i),dirZ(i)
        End Do
        Write(6,'(30(F6.3,a))')
     &                       (TempMagn(iTemp),' K.;',iTemp=1,nTempMagn)
        If (zeeman_energy) Then
          Write(6,'(A)') 'dir_weight'
          Do i=1,nDirZee
            Write(6,'(3F10.6)') (dir_weight(i,j),j=1,3)
          End Do
        End If
      End If

      nP=get_nP(nsymm,ngrid)


      Call hdir( nDir,nDirZee, dirX,dirY,dirZ,dir_weight,
     &           nP,nsymm,ngrid,nDirTot,dHX,dHY,dHZ,dHW)


      Write(6,'(2X,A,i3,A)') 'Molar magnetization will be '//
     &                       'calculated in ',NH,
     &                       ' points, equally distributed in '//
     &                       'magnetic field range'
      Write(6,'(2X,F4.1,1x,a,1X,F4.1,a,10(F6.3,a))') HMIN,'--',HMAX,
     &                       ' T., at the following temperatures:'
      Do i=1,nTempMagn,10
        j=MIN(nTempMagn,i+9)
        Write(6,'(10(F8.4,A))') (TempMagn(l),' K.;',l=i,j)
      End Do
      Write(6,'(2X,A,I4,A)') 'Powder molar magnetization will be '//
     &                       'averaged on ',nP,' directions of the '//
     &                       'applied magnetic field.'
      Write(6,'(2x,10A)') ('--------', i=1,10)
      If(nsymm.eq.1) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a hemisphere:'
        Write(6,'(38x,A)') 'z >= 0;'
      Else If(nsymm.eq.2) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a 4th-of-a-sphere:'
        Write(6,'(34x,A)') 'x >= 0; z >= 0;'
      Else If(nsymm.eq.3) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a 8th-of-a-sphere:'
        Write(6,'(30x,A)') 'x >= 0; y >= 0; z >= 0;'
      End If
      Write(6,'(2x,10A)') ('--------', i=1,10)
      Write(6,'(2x,A,12x,A,2(18x,A),16x,A)') 'Nr.','x','y','z',
     &                                         'weight'
      Do i=1,nP
        Write(6,'(i4,2x,4(F18.12,1x))')
     &                      i,dHX(i+nDir+nDirZee),dHY(i+nDir+nDirZee),
     &                        dHZ(i+nDir+nDirZee),dHW(i+nDir+nDirZee)
      End Do

      If( nDir>0 ) Then
        Write(6,'(2x,10A)') ('--------', i=1,10)
        Write(6,'(23x,A)') ' Magnetization vector will be computed'
        Write(6,'(2x,10A)') ('--------', i=1,10)
        Write(6,'(2x,A,12x,A,2(18x,A))') 'Nr.','x','y','z'
        Do i=1,nDir
          Write(6,'(i4,2x,4(F18.12,1x))') i, dHX(i), dHY(i), dHZ(i)
        End Do
      End If

      If (zeeman_energy) Then
        Write(6,'(2x,10A)') ('--------', i=1,10)
        Write(6,'(23x,A)') 'Zeeman Energy Splitting will be computed'
        Write(6,'(2x,10A)') ('--------', i=1,10)
        Write(6,'(2x,A,12x,A,2(18x,A))') 'Nr.','x','y','z'
        Do i=1,nDirZee
           j=i+nDir
          Write(6,'(i4,2x,4(F18.12,1x))') j, dHX(j), dHY(j), dHZ(j)
        End Do
      End If

      Write(6,'(2x,10A)') ('--------', i=1,10)
      Write(6,'(2X,A)') 'The cut-off energy for the exact '//
     &                  'diagonalization of the Zeeman Hamiltonian is:'
      Write(6,'(2x,a,F15.9,A)') 'E = ',EM ,' cm(-1).'
      If(NM.lt.10) Then
        Write(6,'(2X,A,i2,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.10) .AND. (NM.lt.100) ) Then
        Write(6,'(2X,A,i3,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.100) .AND. (NM.lt.1000) ) Then
        Write(6,'(2X,A,i4,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.1000) .AND. (NM.lt.10000) ) Then
        Write(6,'(2X,A,i5,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      End If
      If(compute_Mdir_vector) Then
        Write(6,'(2x,A,i2,a)') 'The magnetization vector for ',nDir,
     &                         ' directions of the applied '//
     &                         'magnetic field will be calculated.'
      Else
        Write(6,'(2X,A)') 'The magnetization vector was not calculated.'
      End If
      If(zeeman_energy) Then
        Write(6,'(2x,A,i2,a)') 'The Zeeman splitting for ',nDirZee,
     &                         ' directions of the applied '//
     &                         'magnetic field will be calculated.'
        Write(6,'(2x,a     )') 'The Zeeman energies for each '//
     &                         'direction of the applied magnetic '//
     &                         'field are written in files '//
     &                         '"zeeman_energy_xxx.txt".'
      Else
        Write(6,'(2X,A)') 'Computation of the Zeeman splitting was '//
     &                    'not requested.'
      End If
!      smagn=.false.
!      m_paranoid=.true.
!      THRS=1.d-10 !threshold for convergence of average spin, in case (zJ .ne. 0)
C /// opening the loop over the field points
      Do iH=1,nH
c /// ---------------------------------------------------------------
        If (HINPUT) Then
           H(iH)=HEXP(iH)
           If(H(iH).eq.0.0_wp) Then
             H(iH)=0.0001_wp
           End If
        Else
          DLTH=(HMAX-HMIN)/DBLE(NH-1)
          If (iH.eq.1) Then
            H(IH)=HMIN+0.0001_wp
          Else
            H(IH)=HMIN+DLTH*DBLE(IH-1)
          End If
          If(H(iH).eq.0.0_wp) Then
            H(iH)=0.0001_wp
          End If
        End If

         !    If(DBG)
         If(DBG) Write(6,'(A,i0,A,F10.5,A,L2,A,L2)')
     &          'MAGNETIZATION::  H(',iH,') = ', H(iH),
     &          'smagn=',smagn,' m_paranoid=',m_paranoid

c ///  opening the loop over dIfferent directions of the magnetic field
        Do iM=1,NDIRTOT
          MT(1:3,iH,1:nTempMagn)=0.0_wp
          ST(1:3,iH,1:nTempMagn)=0.0_wp
          Call dcopy_(nM,[0.0_wp],0,WM,1)
c         Entry into monitor: Status line
          WRITE(STLNE1,'(A)') 'SINGLE_ANISO:  powder magnetization:'
          WRITE(STLNE2,'(A,I4,A,I4,A,I4,A,I4)')
     &    'Field:',IH,'from ',nH,' at direction ',IM,'from ', NDIRTOT
          Call StatusLine(STLNE1,STLNE2)
c  actual calculation of the MT and ST, ZT
          Call MAGN( NSS, NM, dHX(iM),dHY(iM),dHZ(iM),H(iH),
     &               ESO, zJ, THRS,
     &               DIPM,
     &                 SM,
     &               nTempMagn,TempMagn,smagn,
     &               WM,
     &               ZT(iH,1:nTempMagn),
     &               ST(1:3,iH,1:nTempMagn),
     &               MT(1:3,iH,1:nTempMagn), m_paranoid, DBG )
          If (DBG .and. (iH.eq.nH) .and. (iM.eq.23)) Then
            Write(6,'(A,3ES16.8)') 'iM:',dHX(iM),dHY(iM),dHZ(iM)
            Write(6,'(2(A,3ES16.8,1x),A,ES16.8)')
     &                                'MT:',(MT(l,iH,1),l=1,3),
     &                                'ST:',(ST(l,iH,1),l=1,3),
     &                             'ZSTAT:',ZT(iH,1)
            Write(6,'(A,3ES16.8)') 'WM:',(WM(l),l=1,nM)
          End If
C  ---------------------------------------------------------------------
          If (zeeman_energy) Then
             If((iH.eq.1).and.(iM.eq.nDir+1))
     &         Write(6,'(A)') 'Energies of the Zeeman Hamiltonian '//
     &                        'for the following directions of the '//
     &                        'applied field:'
             If ( (iH.eq.1).and.(iM.gt.nDir) .and.
     &                          (iM.le.(nDir+nDirZee)) ) Then
               Write(6,'(A,I3,A,3F10.6,3x,5A)')
     &                 'direction Nr.',iM-nDir, ' : ',
     &                  dHX(iM), dHY(iM), dHZ(iM),
     &                 'written in file "zeeman_energy_',
     &                  CHAR(48+mod(int((iM-nDir)/100),10)),
     &                  CHAR(48+mod(int((iM-nDir)/10 ),10)),
     &                  CHAR(48+mod(    (iM-nDir)     ,10)),
     &                  '.txt".'

               Write(LUZee(iM-nDir),'(A,3F24.15)')
     &                  '# direction of the applied magnetic field:',
     &                  dHX(iM), dHY(iM), dHZ(iM)
               Write(LUZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',
     &                  ' State =>',(i,i=1,nm)
             End If

             If ( (iM.gt.nDir) .and. (iM.le.(nDir+nDirZee)) ) Then
               Write(LUZee(iM-nDir),'(F8.4,1000F10.3)')
     &                               H(IH),(WM(I),I=1,NM)
             End If
          End If !zeeman_energy
C  ---------------------------------------------------------------------
!         computing the AVERAGE MOMENTS calculated at different temperatures
!         (TempMagn(i))
          Do iT=1,nTempMagn
            Do ic=1,3
              MVEC(iM,iH,iT,ic)=MT(ic,iH,iT)
              SVEC(iM,iH,iT,ic)=ST(ic,iH,iT)
            End Do  ! ic

            If(iM.gt.nDir+nDirZee) Then
c              MAV(iH,iTemp) = MAV(iH,iTemp)+MT(1,iH,iT)*dHX(iM)*dHW(iM)
c     &                                     +MT(2,iH,iT)*dHY(iM)*dHW(iM)
c     &                                     +MT(3,iH,iT)*dHZ(iM)*dHW(iM)
c
c              SAV(iH,iTemp) = SAV(iH,iTemp)+ST(1,iH,iT)*dHX(iM)*dHW(iM)
c     &                                     +ST(2,iH,iT)*dHY(iM)*dHW(iM)
c     &                                     +ST(3,iH,iT)*dHZ(iM)*dHW(iM)
              ! accumulate contributions:
              Call daxpy_(1,dHX(iM)*dHW(iM),MT(1,iH,iT),1,MAV(iH,iT),1)
              Call daxpy_(1,dHY(iM)*dHW(iM),MT(2,iH,iT),1,MAV(iH,iT),1)
              Call daxpy_(1,dHZ(iM)*dHW(iM),MT(3,iH,iT),1,MAV(iH,iT),1)
              Call daxpy_(1,dHX(iM)*dHW(iM),ST(1,iH,iT),1,SAV(iH,iT),1)
              Call daxpy_(1,dHY(iM)*dHW(iM),ST(2,iH,iT),1,SAV(iH,iT),1)
              Call daxpy_(1,dHZ(iM)*dHW(iM),ST(3,iH,iT),1,SAV(iH,iT),1)
            End If
            If (iprint > 2) Then
               If((iM==1).AND.(iH==1))
     &            Write(6,'(2x,A,1x,A,4x,A,7x,A,7x,A)')
     &                  'iH','iM','iT',
     &                  'moment(iM,iT)',
     &                  'spin(iM,iT)'
                  Write(6,'(3i4,3(F21.15,1x))')
     &                   iH,iM,iT,
     &                (  MT(1,iH,iT)*dHX(iM)
     &                 + MT(2,iH,iT)*dHY(iM)
     &                 + MT(3,iH,iT)*dHZ(iM) ),

     &                (  ST(1,iH,iT)*dHX(iM)
     &                 + ST(2,iH,iT)*dHY(iM)
     &                 + ST(3,iH,iT)*dHZ(iM) )
            End If
          End Do !iT
C  ---------------------------------------------------------------------
c ///  closing the loop over directions of magnetic field
        End Do   !  iM
c /// -------------------------------------------------------------------
c ///  closing the loop over the points of magnetic field:
      End Do ! IH
c /// -------------------------------------------------------------------
c Close Zeeman files, if opened
      If (Zeeman_Energy) Then
        Do i=1,nDirZee
          Close( LUZee(i) )
        End Do
      End If

c /// -------------------------------------------------------------------
      If(nDir .gt. 0) Then

      If (smagn) Then
        Do iT=1,nTempMagn
          Write(6,*)
          Do iDir=1,nDir
            Write(6,*)
            Write(6,'(A,A,1x,A)') '--------|',
     &                            '------------------------------'//
     &                            '------------------------------|',
     &                            '|------------------------------'//
     &                            '------------------------------|'
            Write(6,'(A,i3,26x,A,1x,A,60x,A)')
     &                   'Direction of the applied magnetic field:',
     &                    iDir,'|','|','|'
            Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                   'proj X=',dHX(iDIR),'|','|','|'
            Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                   'proj Y=',dHY(iDir),'|','|','|'
            Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                   'proj Z=',dHZ(iDir),'|','|','|'
            Write(6,'(A,F7.4,A,41x,A,1x,A,60x,A)')
     &                   'Temperature = ',TempMagn(iT),
     &                   ' Kelvin','|','|','|'
            Write(6,'(A,A,1x,A)') '--------|',
     &                            '------------------------------'//
     &                            '------------------------------|',
     &                            '|------------------------------'//
     &                            '------------------------------|'
            Write(6,'(2x,A,12x,2A,1x,A,10x,2A)') 'Field |',
     &                   'Magnetization Vector            |',
     &                   '   Total Magn. |','|',
     &                   'Spin Magnetization Vector         |',
     &                   '   Total Magn. |'
            Write(6,'(5A,1x,5A)') '--------|',
     &                   '--- proj X ---|','--- proj Y ---|',
     &                   '--- proj Z ---|','- in this dir.-|',
     &                   '|--- proj X ---|','--- proj Y ---|',
     &                   '--- proj Z ---|','- in this dir.-|'
            Do iH=1,nH
              mv=0.0_wp
              sv=0.0_wp
              mv = MVEC(iDir,iH,iT,1)*dHX(iDir)
     &           + MVEC(iDir,iH,iT,2)*dHY(iDir)
     &           + MVEC(iDir,iH,iT,3)*dHZ(iDir)
              sv = SVEC(iDir,iH,iT,1)*dHX(iDir)
     &           + SVEC(iDir,iH,iT,2)*dHY(iDir)
     &           + SVEC(iDir,iH,iT,3)*dHZ(iDir)

              Write(6,'(F7.3,1x,A, 3(E13.6,1x,A),E14.7,1x,A,1x,A, '//
     &                            '3(E13.6,1x,A),E14.7,1x,A)')
     &                   H(iH),'|', MVEC(iDir,iH,iT,1),' ',
     &                              MVEC(iDir,iH,iT,2),' ',
     &                              MVEC(iDir,iH,iT,3),'|', mv,'|','|',
     &                              SVEC(iDir,iH,iT,1),' ',
     &                              SVEC(iDir,iH,iT,2),' ',
     &                              SVEC(iDir,iH,iT,3),'|', sv,'|'
            End Do
            Write(6,'(A,A,1x,A)') '--------|',
     &                            '------------------------------'//
     &                            '------------------------------|',
     &                            '|------------------------------'//
     &                            '------------------------------|'
          End Do !iDir
        End Do !iT

      Else
      ! smagn == .false.

        Do iT=1,nTempMagn
          Write(6,*)
          Do iDir=1,nDir
            Write(6,*)
            Write(6,'(2A)') '--------|',
     &                      '------------------------------'//
     &                      '------------------------------|'
            Write(6,'(A,i3,26x,A)')
     &                      'Direction of the applied magnetic field:',
     &                       iDir,'|'
            Write(6,'(A,F18.14,44x,A)') 'proj X=',dHX(iDIR),'|'
            Write(6,'(A,F18.14,44x,A)') 'proj Y=',dHY(iDir),'|'
            Write(6,'(A,F18.14,44x,A)') 'proj Z=',dHZ(iDir),'|'
            Write(6,'(A,F7.4,A,41x,A)') 'Temperature = ',TempMagn(iT),
     &                                  ' Kelvin','|'
            Write(6,'(2A)') '--------|',
     &                      '------------------------------'//
     &                      '------------------------------|'
            Write(6,'(2x,A,12x,2A)') 'Field |',
     &                      'Magnetization Vector            |',
     &                      '   Total Magn. |'
            Write(6,'(5A)') '--------|',
     &                      '--- proj X ---|','--- proj Y ---|',
     &                      '--- proj Z ---|','- in this dir.-|'
            Do iH=1,nH
              mv=0.0_wp
              mv = MVEC(iDir,iH,iT,1)*dHX(iDir)
     &           + MVEC(iDir,iH,iT,2)*dHY(iDir)
     &           + MVEC(iDir,iH,iT,3)*dHZ(iDir)
              Write(6,'(F7.3,1x,A,3(E13.6,1x,A),E14.7,1x,A)') H(iH),'|',
     &                MVEC(iDir,iH,iT,1),' ',
     &                MVEC(iDir,iH,iT,2),' ',
     &                MVEC(iDir,iH,iT,3),'|', mv,'|'
            End Do
            Write(6,'(2A)') '--------|',
     &                      '------------------------------'//
     &                      '------------------------------|'
          End Do !iDir
        End Do !iT
        End If !(smagn)
      End If !(nDir>0)
c /// -------------------------------------------------------------------
c     COMPUTING THE STANDARD DEVIATION OF THE MAGNETIZATION
      If(HINPUT) Then
        Do iT=1,nTempMagn
          STDEV(iT) = 0.0_wp
          STDEV(iT) = dev(nH,MAV(:,iT),magn_exp(:,iT))
        End Do
      End If
c /// -------------------------------------------------------------------

      Write(6,*)
      Write(6,'(25X,A)') 'HIGH-FIELD POWDER MAGNETIZATION'
      Write(6,'(30X,A)') '(Units: Bohr magneton)'
      Write(6,*)
      Do iT=1,nTempMagn,5
        iTEnd=min(nTempMagn,iT+4)

        Write(6,'(A,11A)') '-----------|',
     &                    ('---------------|',i=iT,iTEnd+1)
        Write(6,'(A,10(F10.3,A))') '   H(T)    |STATISTICAL SUM|',
     &                    (TempMagn(i),' K.  |',i=iT,iTEnd)
        Write(6,'(A,11A)') '-----------|',
     &                    ('---------------|',i=iT,iTEnd+1)

        Do iH=1,nH
          Write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))')
     &          H(iH),'|',ZT(iH,1),'|', (MAV(iH,i),'|',i=iT,iTEnd)
        End Do

        Write(6,'(A,11A)') '-----------|',
     &                    ('---------------|',i=iT,iTEnd+1)
        If(HINPUT) Then
          Write(6,'(A,15x,11(f14.10,1x,A) )')
     &                       'ST.DEV.M   |',
     &                      ( STDEV(i),'|',i=1,nTempMagn)
          Write(6,'(A,11A)') '-----------|',
     &                      ('---------------|',i=iT,iTEnd+1)
        End If
      End Do


      If(smagn) Then
         Write(6,*)
         Write(6,'(15X,A)') 'HIGH-FIELD POWDER SPIN MAGNETIZATION'
         Write(6,'(20X,A)') '(Units: Bohr magneton)'
         Write(6,*)
         Do iT=1,nTempMagn,5
           iTEnd=min(nTempMagn,iT+4)
           Write(6,'(A,11A)') '---------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
           Write(6,'(A,11(F10.3,A))') '   H(T)  |STATISTICAL SUM|',
     &                               (TempMagn(i),' K.  |',i=iT,iTEnd)
           Write(6,'(A,11A)') '---------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
           Do iH=1,nH
            Write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))')
     &              H(iH),'|',ZT(iH,1),'|', (SAV(iH,i),'|',i=iT,iTEnd)
           End Do
           Write(6,'(A,11A)') '---------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
         End Do
      End If!smagn


      If(DoPlot) Then
        IF ( hinput ) THEN
           Call plot_MH_with_Exp( nH, H, nTempMagn, TempMagn, MAV,
     &                            magn_exp, zJ )
        ELSE
           Call plot_MH_no_Exp( nH, H, nTempMagn, TempMagn, MAV, zJ )
        END IF
!        IF ( zeeman_energy ) THEN
!           Call plot_zeeman( nH, nM, nDirZee, H, LuZee )
!        END IF
      End If


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Call Add_Info('MAGN_AVERAGED',MAV(1:nH,1:nTempMagn),
     &               nH*nTempMagn,5)
      If(compute_Mdir_vector) Then
        Call Add_Info('MAGN_VECT_X(2)     ',MVEC(1, 2,1,1),1,5)
        Call Add_Info('MAGN_VECT_X(nH/2)  ',
     &                     MVEC(1,(NH-1)/2,1,1),1,5)
        Call Add_Info('MAGN_VECT_X(nH)    ',MVEC(1,NH,1,1),1,5)
        Call Add_Info('MAGN_VECT_Y(2)     ',MVEC(1, 2,1,2),1,5)
        Call Add_Info('MAGN_VECT_Y(nH/2)  ',
     &                     MVEC(1,(NH-1)/2,1,2),1,5)
        Call Add_Info('MAGN_VECT_Y(nH)    ',MVEC(1,NH,1,2),1,5)
        Call Add_Info('MAGN_VECT_Z(2)     ',MVEC(1, 2,1,3),1,5)
        Call Add_Info('MAGN_VECT_Z(nH/2)  ',
     &                     MVEC(1,(NH-1)/2,1,3),1,5)
        Call Add_Info('MAGN_VECT_Z(nH)    ',MVEC(1,NH,1,3),1,5)
      End If


!-----------------------------------------------------------------------
! Deallocate necessary memory
      If(nM>=0) Then
         Call mma_deallocate(WM)
      End If

      If((nH>=0).and.(nTempMagn>=0)) Then
         Call mma_deallocate(MT)
         Call mma_deallocate(ST)
         Call mma_deallocate(MAV)
         Call mma_deallocate(SAV)
         Call mma_deallocate(ZT)
         If(nDirTot>=0) Then
            Call mma_deallocate(MVEC)
            Call mma_deallocate(SVEC)
         End If
      End If

      If(nH>=0) Then
         Call mma_deallocate(H)
      End If

      If((nTempMagn>=0).and.hinput) Then
         Call mma_deallocate(STDEV)
      End If

      If(nDirTot>=0) Then
            Call mma_deallocate(dHX)
            Call mma_deallocate(dHY)
            Call mma_deallocate(dHZ)
            Call mma_deallocate(dHW)
      End If

      Call qExit('SA_magn')
      Return
      End




