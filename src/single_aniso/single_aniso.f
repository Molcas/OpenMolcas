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
      Subroutine SINGLE_ANISO_OPEN(IReturn)

      Implicit None

      Integer ::  iReturn,NSS,NSTATE
      Integer ::  i, input_to_read, nH, nT, nTempMagn
      Integer ::  nDir, nDirZee, nMult
      Logical ::  ifrestart,GRAD
      Character(180) :: input_file_name
      Logical :: dbg

C----------------------------------------------------------------------
C*  initializations
      ireturn  =0
      ifrestart=.false.
      GRAD     =.false.
      nT       =0
      nH       =0
      nTempMagn=0
      nDir     =0
      nDirZee  =0
      nMult    =0
      nss      =0
      nstate   =0
      dbg = .false.

      ! check for the "restart" option:
      If(dbg) Write(6,*) 'Enter restart_check'
      Call restart_check( ifrestart, input_to_read, input_file_name,
     &                    nT, nH, nTempMagn, nDir, nDirZee, nMult,
     &                    GRAD )
      If(dbg) Write(6,*) 'Exit restart_check'

      If ( ifrestart ) Then
         Call restart_sa( input_to_read, input_file_name,
     &                    nss, nstate )
      Else ! not a restart job -- take all data form RUNFILE
         Call fetch_data_RunFile_init( nss, nstate )
      End If !Ifrestart

      Write(6,'(120A)') ('@',i=1,95)
      Write(6,'(A)') '   SINGLE_ANISO (OPEN)'
      Write(6,'(A)') '(last updated on 12-March-2018)'
      Write(6,'(A)') '   New features: '
      Write(6,*)
      Write(6,'(A)') '1.  Calculation of the SIGN of the product '//
     &               'gX * gY * gZ for any moment;'
      Write(6,'(A)') '2.  Calculation of the parameters of the '//
     &               'Crystal-Field for lanthanides (CRYS).'
      Write(6,'(A)') '3.  Automatic generation of various plot: (PLOT)'
      Write(6,'(A)') '     -- Powder magnetic susceptibilty:  XT=f(T)'
      Write(6,'(A)') '     -- Powder molar magnetization:  M=f(H,T)'
      Write(6,'(A)') '4.  Support for various restart options: (REST)'
      Write(6,'(A)') '     -- from  $Project.rassi.h5 file.'
      Write(6,'(A)') '     -- from  $Project.aniso (binary) file.'
      Write(6,'(A)') '     -- from  ANISOINPUT (ascii) file.'
      Write(6,'(A)') '     -- from  $Project.RunFile file.'
      Write(6,'(A)') '5.  RASSI was adjusted to provide more '//
     &               'accurate tunnelling gaps between near-'//
     &               'degenerate states'
      Write(6,'(A)') 'Check the MOLCAS manual for details and input '//
     &               'examples.'
      Write(6,*)
      Write(6,'(120A)') ('@',i=1,95)
      Call xFlush(6)

      Call SINGLE_ANISO2(nH,nT,nTempMagn,nDir,nDirZee,nss,nstate,
     &                   nMult,input_file_name,
     &                   Ifrestart,IReturn,GRAD)

      Return
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine SINGLE_ANISO2(nH,nT,nTempMagn,nDir,nDirZee,nss,nstate,
     &                         nMult,input_file_name,
     &                         Ifrestart,IReturn,GRAD)

      Implicit None
      Integer, Parameter :: wp=selected_real_kind(p=15,r=307)
#include "cntrl.fh"
#include "stdalloc.fh"
      Integer                   :: mem,RtoB,CtoB,ItoB
      Integer, intent(in)       :: nss,nstate
      Integer                   :: i,j,iReturn
      Integer                   :: idim,Ifunct,imltpl,iprint
      Integer                   :: axisoption
      Integer                   :: ndimcf,ldimcf
      Integer                   :: input_to_read,nK,MG,nm
      Integer                   :: encut_definition
      Integer                   :: ncut,nDirTot
      Integer                   :: nBlock
      Integer                   :: imanIfold
      Integer                   :: i1,i2,ldim2,lDIM
      Integer                   :: AngPoints
      Integer                   :: nlanth
      Real(kind=wp)             :: zmagn(3,3)
      Real(kind=wp)             :: rdiff(nstate)
      Real(kind=wp)             :: a(6),HMIN,HMAX
      Real(kind=wp)             :: cryst(6),coord(3)
      Real(kind=wp)             :: encut_rate,em,Boltz_k,mu_Bohr
      Real(kind=wp)             :: zJ, thrs
      Integer, allocatable      :: multiplicity(:)
      !---g-tens-------------------
      Integer, intent(in)        :: nMult
      Integer, allocatable       :: ndim(:)
      Real(kind=wp), allocatable :: gtens(:,:), maxes(:,:,:)
      !---M-------------------
      Integer       :: nH, nTempMagn
      !---MVEC and ZEEM-------------------
      Integer       :: nDir, nDirZee
      Integer, allocatable :: LuZee(:)
      Real(kind=wp), allocatable :: dir_weight(:,:)
      Real(kind=wp), allocatable :: dirX(:), dirY(:), dirZ(:)
      !---XT-------------------
      Integer       :: nT
      Real(kind=wp), allocatable :: Texp(:)
      Real(kind=wp), allocatable :: chit_exp(:)
      Real(kind=wp) :: xfield
      Real(kind=wp) :: tmin, tmax
      !---Oscillator strength----------
c      Real(kind=wp) :: F, Fx,Fy,Fz, AT, Ax,Ay,Az, AF, dnrm, dE
      Real(kind=wp) :: DZNRM2
      External      :: DZNRM2
      Real(kind=wp) :: H_torq, T_torq
      !----BIG ARRAYS------------------
      Real(kind=wp), allocatable :: eso(:)
      Real(kind=wp), allocatable :: esfs(:)
      Real(kind=wp), allocatable :: t(:)
      Real(kind=wp), allocatable :: XTexp(:)
      Real(kind=wp), allocatable :: XT_no_field(:)
      Real(kind=wp), allocatable :: hexp(:)
      Real(kind=wp), allocatable :: magn_exp(:,:)
      Real(kind=wp), allocatable :: angmom(:,:,:)
      Real(kind=wp), allocatable ::  eDmom(:,:,:)
      Real(kind=wp), allocatable ::   amfi(:,:,:)
      Real(kind=wp), allocatable :: TempMagn(:)
      Complex(kind=wp), allocatable :: MM(:,:,:), MS(:,:,:), HSO(:,:),
     &                                 ML(:,:,:), DM(:,:,:), U(:,:)
      Character(180), intent(in) :: input_file_name

      Logical :: poly_file
      Logical :: ifrestart
      Logical :: Do_structure_abc
      Logical :: compute_magnetization
      Logical :: m_paranoid
      Logical :: compute_Mdir_vector
      Logical :: zeeman_energy
      Logical :: compute_torque
      Logical :: smagn
      Logical :: compute_cf
      Logical :: compute_g_tensors
      Logical :: compute_barrier
      Logical :: tinput, hinput
      Logical :: GRAD
      Logical :: DoPlot
      Logical :: DBG
      Integer :: l

      Call qEnter('SA_main1')
      DBG=.false.

      Boltz_k=0.6950356_wp                    !   in cm-1*K-1
      mu_Bohr=0.466864374_wp                  !   in cm-1*T-1

c---------------------------------------------------------------------
      ! Allocate memory for all arrays:
c---------------------------------------------------------------------
      If(dbg) Write(6,*) 'S_A2:: Memory Allocation Parameters'
      If(dbg) Write(6,*) 'S_A2:: nss             =',nss
      If(dbg) Write(6,*) 'S_A2:: nstate          =',nstate
      If(dbg) Write(6,*) 'S_A2:: nH              =',nH
      If(dbg) Write(6,*) 'S_A2:: nT              =',nT
      If(dbg) Write(6,*) 'S_A2:: nTempMagn       =',nTempMagn
      If(dbg) Write(6,*) 'S_A2:: nMult           =',nMult
      If(dbg) Write(6,*) 'S_A2:: nDir            =',nDir
      If(dbg) Write(6,*) 'S_A2:: nDirZee         =',nDirZee
      If(dbg) Write(6,*) 'S_A2:: input_file_name =',input_file_name
      If(dbg) Write(6,*) 'S_A2:: GRAD            =',GRAD

      mem=0
      RtoB=8
      CtoB=16
      ItoB=8

      If(nstate>0) Then
         ! spin free energies
         Call mma_allocate(esfs,nstate,'esfs')
         Call dcopy_(nstate,0.0_wp,0,ESFS,1)
         mem=mem+nstate*RtoB
         ! angular momentum
         Call mma_allocate(ANGMOM,3,nstate,nstate,'angmom')
         Call dcopy_(3*nstate*nstate,0.0_wp,0,ANGMOM,1)
         mem=mem+3*nstate*nstate*RtoB
         ! electric dipole moment
         Call mma_allocate( EDMOM,3,nstate,nstate,'edmom')
         Call dcopy_(3*nstate*nstate,0.0_wp,0, EDMOM,1)
         mem=mem+3*nstate*nstate*RtoB
         ! amfi integrals
         Call mma_allocate( AMFI,3,nstate,nstate,'amfi')
         Call dcopy_(3*nstate*nstate,0.0_wp,0, AMFI,1)
         mem=mem+3*nstate*nstate*RtoB
         ! multiplicity of each state
         Call mma_allocate(multiplicity,nstate,'multiplicity')
         Call icopy(nstate,0,0,multiplicity,1)
         mem=mem+nstate*ItoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 1 =',mem
      End If

      If(nss>0) Then
         ! spin orbit energies
         Call mma_allocate(eso,nss,'eso')
         Call dcopy_(nss,0.0_wp,0,eso,1)
         mem=mem+nss*RtoB
         ! spin orbit eigenstates
         Call mma_allocate(U,nss,nss,'U')
         Call zcopy_(nss*nss,(0.0_wp,0.0_wp),0,U,1)
         mem=mem+nss*nss*CtoB
         ! spin orbit hamiltonian
         Call mma_allocate(HSO,nss,nss,'HSO')
         Call zcopy_(nss*nss,(0.0_wp,0.0_wp),0,HSO,1)
         mem=mem+nss*nss*CtoB
         ! magnetic moment
         Call mma_allocate(MM,3,nss,nss,'MM')
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,MM,1)
         mem=mem+3*nss*nss*CtoB
         ! spin moment
         Call mma_allocate(MS,3,nss,nss,'MS')
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,MS,1)
         mem=mem+3*nss*nss*CtoB
         ! orbital mooment
         Call mma_allocate(ML,3,nss,nss,'ML')
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,ML,1)
         mem=mem+3*nss*nss*CtoB
         ! electric dipole moment
         Call mma_allocate(DM,3,nss,nss,'DM')
         Call zcopy_(3*nss*nss,(0.0_wp,0.0_wp),0,DM,1)
         mem=mem+3*nss*nss*CtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 2 =',mem
      End If

      If( (nH>0).and.(nTempMagn>0) ) Then
         ! experimental magnetic field points
         Call mma_allocate(Hexp,nH,'Hexp')
         Call dcopy_(nH,0.0_wp,0,Hexp,1)
         mem=mem+nH*RtoB
         ! experiemental magnetization
         Call mma_allocate(magn_exp,nH,nTempMagn,'magn_exp')
         Call dcopy_(nH*nTempMagn,0.0_wp,0,magn_exp,1)
         mem=mem+nH*nTempMagn*RtoB
         ! temperature points for magnetization
         Call mma_allocate(TempMagn,nTempMagn,'TempMagn')
         Call dcopy_(nTempMagn,0.0_wp,0,TempMagn,1)
         mem=mem+nTempMagn*RtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 3 =',mem
      End If

      If(nMult>0) Then
         ! dimensions of pseudospins
         Call mma_allocate(ndim,nMult,'ndim')
         Call icopy(nMult,0,0,ndim,1)
         mem=mem+nMult*ItoB
         ! temperature points for magnetization
         Call mma_allocate(gtens,nMult,3,'gtens')
         Call dcopy_(3*nMult,0.0_wp,0,gtens,1)
         mem=mem+3*nMult*RtoB
         ! temperature points for magnetization
         Call mma_allocate(maxes,nMult,3,3,'maxes')
         Call dcopy_(3*3*nMult,0.0_wp,0,maxes,1)
         mem=mem+3*3*nMult*RtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 4 =',mem
      End If

      If((nT+nTempMagn)>0) Then
         Call mma_allocate(T,(nT+nTempMagn),'Temperature')
         Call dcopy_((nT+nTempMagn),0.0_wp,0,T,1)
         mem=mem+(nT+nTempMagn)*RtoB
         Call mma_allocate(XTexp,(nT+nTempMagn),'XTexp')
         Call dcopy_((nT+nTempMagn),0.0_wp,0,XTexp,1)
         mem=mem+(nT+nTempMagn)*RtoB
         Call mma_allocate(XT_no_field,(nT+nTempMagn),'XT_no_field')
         Call dcopy_((nT+nTempMagn),0.0_wp,0,XT_no_field,1)
         mem=mem+(nT+nTempMagn)*RtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 5 =',mem
      End If

      If(nDirZee>0) Then
         ! unit numbers for the files with Zeeman energies
         Call mma_allocate(LuZee,nDirZee,'LUZee')
         Call icopy(nDirZee,0,0,LuZee,1)
         mem=mem+nDirZee*ItoB
         ! directions for applied field for Zeeman states
         Call mma_allocate(dir_weight,nDirZee,3,'dir_weight')
         Call dcopy_(3*nDirZee,0.0_wp,0,dir_weight,1)
         mem=mem+3*nDirZee*RtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 6 =',mem
      End If

      If(nDir>0) Then
         ! magnetization vectors
         Call mma_allocate(dirX,nDir,'dirX')
         Call mma_allocate(dirY,nDir,'dirY')
         Call mma_allocate(dirZ,nDir,'dirZ')
         mem=mem+3*nDir*RtoB
         Call dcopy_(nDir,0.0_wp,0,dirX,1)
         Call dcopy_(nDir,0.0_wp,0,dirY,1)
         Call dcopy_(nDir,0.0_wp,0,dirZ,1)
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 7 =',mem
      End If

      If(nT>0) Then
        ! T expeirimental given by user in the input
        Call mma_allocate(Texp,nT,'Texp')
        Call dcopy_(nT,0.0_wp,0,Texp,1)
        mem=mem+nT*RtoB
        ! XT expeirimental given by user in the input
        Call mma_allocate(chit_exp,nT,'chit_exp')
        Call dcopy_(nT,0.0_wp,0,chit_exp,1)
        mem=mem+nT*RtoB
         ! allocated memory counter
         If(dbg) Write(6,'(A,I16)') 'mem 8 =',mem
      End If

      Write(6,'(A,I16,A)') 'The code allocated initially:',mem,
     &                     ' bytes of memory for this run.'
      Call xFlush(6)
c---------------------------------------------------------------------
      IReturn=0
      IPRINT=2
      lDIM=1
      iDIM=1
      NM=0
      EM=0.0_wp
      POLY_FILE=.FALSE.
      compute_CF=.FALSE.
      axisoption=1
      nDIMcf=1
      lDIMcf=1
      a(:)=0.0_wp
      H_torq=0.1_wp ! in Tesla
      T_torq=2.0_wp ! in K
C  read the input
      IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter readin_single'
      Call readin_single(iprint,nmult,ndim,ldim,ndimcf,ldimcf,
     & nlanth,axisoption,poly_file,Ifrestart,input_to_read,nk,mg,
     & zmagn,Do_structure_abc,cryst,coord,encut_definition,
     & compute_g_tensors,compute_CF,nDirTot,nss,nstate,
     & compute_magnetization, compute_torque, smagn, tinput, hinput,
     & compute_Mdir_vector, zeeman_energy, LUZee, doplot,
     & encut_rate,ncut,nTempMagn,TempMagn, m_paranoid,
     & compute_barrier,nBlock,AngPoints,input_file_name,
     & nT,nH,texp,chit_exp,zJ,hexp,magn_exp,hmin,hmax,
     & nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,thrs,
     & H_torq,T_torq )
      IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit readin_single'

      If ( ifrestart ) Then
         ! if restart, fetch "big data" from the input file
         If ( input_to_read .eq. 1 ) Then
            Call read_binary_aniso( nss, nstate, multiplicity, eso,
     &                              esfs, U, MM, MS, ML, DM, ANGMOM,
     &                              EDMOM )

         Else If ( input_to_read .eq. 2 ) Then
            ! get the information from formatted aniso.input file:
            IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter '//
     &                         'read_formatted_aniso'
            IF(DBG) Write(6,*) 'SA:',input_file_name
            Call read_formatted_aniso( input_file_name, nss, nstate,
     &                                 multiplicity, eso, esfs, U, MM,
     &                                 MS, ML, DM, ANGMOM, EDMOM )
            IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit  '//
     &                         'read_formatted_aniso'

         Else If ( input_to_read .eq. 3 ) Then
            ! get the information from RASSI-HDF5 file:
            IF(DBG) Write(6,*) 'SA:',input_file_name
#ifdef _HDF5_
            Call read_hdf5_all( input_file_name, nss, nstate,
     &                         multiplicity, eso, esfs,  U, MM,
     &                         MS, ML, DM, ANGMOM, EDMOM, AMFI, HSO )
#else
            Call WarningMessage(2,'File '//trim(input_file_name)//
     &                          ' cannot be opened. Molcas was'//
     &                          ' compiled without HDF5 option.')
            Call Quit_OnUserError()
#endif

         Else If ( input_to_read .eq. 4 ) Then
            IF(DBG) Write(6,*) 'SA:',input_file_name
            ! get the information from formatted aniso.input file:
            Call read_formatted_aniso_old( input_file_name, nss, nstate,
     &                                multiplicity, eso, MM, MS, ML )
         End If ! input_to_read

      Else
         ! ifrestart = .false., i.e. usual S-A calculation
         Call fetch_data_RunFile_all( nss, nstate, multiplicity, eso,
     &                                esfs, U, MM, MS, ML, DM,
     &                                ANGMOM, EDMOM, AMFI, HSO )
         If (DBG) Then
            Write(6,'(A)') 'SA: ANGMOM(x,y,z)'
            Do i=1,nstate
               Do j=1,nstate
                  Write(6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (ANGMOM(l,i,j),l=1,3)
               End Do
            End Do
            Write(6,'(/)')
            Write(6,'(A)') 'SA: EDMOM(x,y,z)'
            Do i=1,nstate
               Do j=1,nstate
                  Write(6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (EDMOM(l,i,j),l=1,3)
               End Do
            End Do
            Write(6,'(/)')
            Write(6,'(A)') 'SA: AMFI(x,y,z)'
            Do i=1,nstate
               Do j=1,nstate
                  Write(6,'(2i4,A,3ES24.14)') i,j,' |',
     &                                       (AMFI(l,i,j),l=1,3)
               End Do
            End Do
            Write(6,'(/)')
            Write(6,'(A)') 'SA: HSO(i,j)'
            Do i=1,nss
               Do j=1,nss
                  Write(6,'(2i4,A,4ES24.14)') i,j,' |',HSO(i,j),HSO(j,i)
               End Do
            End Do
         End If !DBG

      End If ! Ifrestart

      ! print some input data in the beginning of the output:
      ! so that the user knows which ws the input ...
          Write(6,'(A)') 'LOW-LYING SPIN-ORBIT ENERGIES:'
          Do i=1,nss
            Write(6,'(A,I4,A,F25.14)')
     &               'ENERGY OF THE SPIN-ORBIT STATE (',i,') =',ESO(i)
          End Do
          Write(6,'(A)') 'LOW-LYING SPIN-FREE ENERGIES:'
          Do i=1,nstate
            Write(6,'(A,I4,A,F25.14)')
     &               'ENERGY OF THE SPIN-FREE STATE  (',i,') =',ESFS(i)
          End Do
          rdiff=0.0_wp
          lDIM2=1 ! the same as the defult for lDIM
          Do i=2,nstate
            rdiff(i)=ABS(ESFS(i)-ESFS(i-1))
            If( rdiff(i) < 20.0_wp ) Then
              lDIM2=lDIM2+1
            Else
            goto 107
            End If
          End Do
107       Continue
          ! set the lDIM:
          If (lDIM2>lDIM) lDIM=lDIM2


!----- input processing finished -----
! save some important data, regardless of the following execution
      ! -- binary $Project.aniso
      Call write_binary_aniso( nss, nstate, multiplicity, eso,
     &                         esfs, U, MM, MS, DM, ANGMOM, EDMOM )
      ! ASCII -- anisoinput:
      Call write_formatted_aniso( nss, nstate, multiplicity, eso,
     &                            esfs, U, MM, MS, DM, ANGMOM, EDMOM )
!
!----- compute various properties ----------|
! calculation of magnetic Hamiltonians:
!      IDIM=0
      IFUNCT=0
      If ( compute_g_tensors .AND. nMult>0 ) Then
         Do IMLTPL=1,NMULT
           IReturn=0
           If(ndim(imltpl).eq.1) Then
             Write(6,'(5X,A,I2,A)') 'THE DIMENSION OF THE ',IMLTPL,
     &                              ' MULTIPLET IS 1.'
             Write(6,'(5X,A)') 'THERE IS NO G TENSOR FOR '//
     &                         'EFFECTIVE S = 0.'
             Go To 10
           End If

           i1=1           +Ifunct
           i2=ndim(imltpl)+Ifunct

           If ( i2<=nss) then
              IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter g_high',IMLTPL
              Call g_high( eso(i1:i2), GRAD,
     &                      MS(1:3, i1:i2, i1:i2 ),
     &                      MM(1:3, i1:i2, i1:i2 ),
     &                     imltpl,
     &                     ndim(imltpl),
     &                     Do_structure_abc,
     &                     cryst,
     &                     coord,
     &                     gtens(imltpl,1:3), maxes(imltpl,1:3,1:3),
     &                     iprint )
              IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit g_high',IMLTPL
           End If
10         Continue

           IFUNCT=IFUNCT+NDIM(IMLTPL)
           If(IReturn.ne.0) Call ABEnd()
         End Do
      End If

!----------------------------------------------------------------------|
      !   >> AB INITIO CRYSTAL-FIELD <<

      If(compute_CF .AND. nDIMCF>0 ) Then
        If(axisoption==1 .AND. nMult>0 .AND. nDIM(1)>1 ) Then
          iDIM=NDIM(1)
          lDIM=NDIM(1)
        Else
          iDIM=nDIMCF
          lDIM=lDIMCF
        End If

        ! compute the CF of the ground |J,MJ> multiplet:
        If((nDIMCF>1).AND.(nDIMCF<=nss)) Then
          IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter CF',nDIMCF
          Call CRYSTALFIELD( ESO(1:nDIMCF),
     &                        MM(1:3,1:nDIMCF,1:nDIMCF),
     &                        MS(1:3,1:nDIMCF,1:nDIMCF),
     &                       nDIMcf, iDIM, nlanth,
     &                       zmagn, axisoption, GRAD,
     &                       iPrint)
          IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit CF',nDIMCF
        End If
      End If

      IF(DBG) Write(6,*) 'SINGLE_ANISO2:: nlanth=', nlanth
      IF(DBG) Write(6,*) 'SINGLE_ANISO2:: lDIMCF=', lDIMCF
      IF(DBG) Write(6,*) 'SINGLE_ANISO2:: nstate=', nstate


      If(compute_CF .AND. ( (lDIMCF>0).AND.(lDIMCF<=nstate) )) Then
        If(axisoption==1 .AND. nMult>0 .AND. nDIM(1)>1 ) Then
          iDIM=NDIM(1)
          lDIM=NDIM(1)
        Else
          iDIM=nDIMCF
          lDIM=lDIMCF
        End If
        ! compute the CF of the ground |L,ML> term:
        If( .not.(ifrestart.and.(input_to_read.eq.4)) ) Then
             IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter t-CF',lDIMCF
             Call termCF( angmom(1:3,1:lDIMCF,1:lDIMCF),
     &                    AMFI(1:3,1:lDIMCF,1:lDIMCF),
     &                    esfs(1:lDIMCF),
     &                    lDIMCF, lDIM,
     &                    zmagn, axisoption, nlanth,
     &                    iPrint )
             IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit t-CF',lDIMCF
        End If !ifrestart
      End If !compute_CF



!----------------------------------------------------------------------|
      !   >> AB INITIO BLOCKING BARRIER FOR SMMs <<
      If ( compute_barrier ) Then
        If ( nBlock.ne.0 ) Then
          imanifold=1

          IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter barrier',nBlock
          Call  BARRIER( nBlock, MM(1:3,1:nBlock,1:nBlock),
     &                   eso(1:nBlock), imanIfold, nMult, nDim,
     &                   iPrint )
          IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit barrier',nBlock

        Else
           Write(6,'(A)') 'nBlock parameter is not defined. '
           Write(6,'(A)') 'Did you specify the MLTP keyword in the '//
     &                    'input?'
           Write(6,'(A)') 'If the problem persists,please, '//
     &                    'submit a bug report.'
        End If
      End If


!----------------------------------------------------------------------|

      Call set_T( nT, nTempMagn, TINPUT, TempMagn, Tmin, Tmax,
     &            chit_exp, Texp, T, XTexp )

      If( compute_torque .OR. compute_magnetization .OR.
     &    (Xfield .ne. 0.0_wp) ) Then
         ! set the NM- number of states to be exactly diagonalized
         ! in Zeeman Interaction
         IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter set_nm'
         Call set_nm( nss, ncut, encut_definition, nk, mg,
     &                nTempMagn, hmax, ESO, encut_rate, TempMagn,
     &                nM, EM, dbg )
         IF(DBG) Write(6,*) 'SINGLE_ANISO2::  exit set_nm'

      End If


!----------------------------------------------------------------------|
      !   >> MAGNETIC SUSCEPTIBILITY <<
      iReturn=0
      IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter SUSCEPTIBILITY'
      Call SUSCEPTIBILITY( NSS, ESO, MS, MM, nT, nTempMagn, T, tmin,
     &                     tmax, XTexp, zJ, tinput, XT_no_field,
     &                     doplot, iPrint, mem)
      IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit SUSCEPTIBILITY'

      If (Xfield .ne. 0.0_wp ) Then
         IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter XT_dMoverdH_single'
         Call XT_dMoverdH_single( nss, nTempMagn, nT, nM, Tmin, Tmax,
     &                            XTexp, ESO, T, zJ, Xfield, EM, MM, MS,
     &                            XT_no_field, tinput, smagn, mem )
         IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit XT_dMoverdH_single'
      End If
!----------------------------------------------------------------------|
      !           >>  TORQUE  <<
      If( compute_torque ) Then

        IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter TORQUE'
        Call torque( Nss, NM, AngPoints, EM, eso, mm, ms, zJ, thrs, mem,
     &                  m_paranoid, smagn, H_torq, T_torq, zmagn, dbg)
        IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit TORQUE'

      End If
!----------------------------------------------------------------------|
      !           >>  MAGNETIZATION <<

      If( compute_magnetization ) Then

        IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Enter magnetization'
      Call magnetization( nss, nM, nTempMagn, nDirTot, nDir, nDirZee,
     &                    nH, iPrint, LUZee, mem, compute_Mdir_vector,
     &                    zeeman_energy, hinput, m_paranoid, smagn,
     &                    doplot, TempMagn, eso, dirX, dirY, dirZ,
     &                    dir_weight, hexp, magn_exp, zJ, hmin, hmax,
     &                    EM, thrs,  mm, ms, dbg )
        IF(DBG) Write(6,*) 'SINGLE_ANISO2::  Exit magnetization'

      Else
        Write(6,*)
        Write(6,'(5X,A)') 'ON USER REQUEST, THE MAGNETIZATION VECTOR '//
     &                    'AND THE MEAN MAGNETIZATION WAS NOT '//
     &                    'CALCULATED'
      End If


c---------------------------------------------------------------------
      ! Deallocate memory for all arrays:
c---------------------------------------------------------------------
      If(nstate>0) Then
         Call mma_deallocate(esfs)
         Call mma_deallocate(ANGMOM)
         Call mma_deallocate( EDMOM)
         Call mma_deallocate(  AMFI)
         Call mma_deallocate(multiplicity)
      End If

      If(nss>0) Then
         Call mma_deallocate(eso)
         Call mma_deallocate(U)
         Call mma_deallocate(HSO)
         Call mma_deallocate(MM)
         Call mma_deallocate(MS)
         Call mma_deallocate(ML)
         Call mma_deallocate(DM)
      End If

      If( (nH>0).and.(nTempMagn>0) ) Then
         Call mma_deallocate(Hexp)
         Call mma_deallocate(magn_exp)
         Call mma_deallocate(TempMagn)
      End If

      If(nMult>0) Then
         Call mma_deallocate(ndim)
         Call mma_deallocate(gtens)
         Call mma_deallocate(maxes)
      End If

      If((nT+nTempMagn)>0) Then
         Call mma_deallocate(T)
         Call mma_deallocate(XTexp)
         Call mma_deallocate(XT_no_field)
      End If

      If(nDirZee>0) Then
         Call mma_deallocate(LuZee)
         Call mma_deallocate(dir_weight)
      End If

      If(nDir>0) Then
         Call mma_deallocate(dirX)
         Call mma_deallocate(dirY)
         Call mma_deallocate(dirZ)
      End If

      If(nT>0) Then
         Call mma_deallocate(Texp)
         Call mma_deallocate(chit_exp)
      End If

      Call qExit('SA_main1')

      Return
      End
