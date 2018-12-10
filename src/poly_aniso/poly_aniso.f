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
      Subroutine POLY_ANISO_1( nneq, neqv, nmax, exch, nLoc,
     &                         nCenter, nT, nH, nTempMagn, nDir,
     &                         nDirZee, nMult, nPair, MxRank1, MxRank2,
     &                         iReturn )

      Implicit None
#include "warnings.fh"
#include "stdalloc.fh"
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
!======================================================================
c  definition of the cluster:

      Integer                       :: mem,RtoB,CtoB,ItoB
      Integer                       :: nneq      ! number of non-equivalent sites
      Integer                       :: neqv      ! number of equivalent sites of each type, neqv = MAXVAL(neq(:))
      Integer                       :: nCenter
      Integer, allocatable          :: neq(:)
c  definition of the local metal sites
      Integer                       :: nLoc                           ! number of spin-orbit states. nLoc = MAXVAL(nss(:));
      Integer, allocatable          :: nss(:), nsfs(:)
c      Integer       :: nsfs(nneq), multiplicity(nneq,nLoc)
      Real(kind=wp), allocatable    :: gtens_input(:,:)  !gtens_input(3,nneq)
      Real(kind=wp), allocatable    :: D_fact(:), EoverD_fact(:)
      Real(kind=wp), allocatable    :: riso(:,:,:)       ! definition of the main local axes in general coord system
      Real(kind=wp), allocatable    :: R_LG( :,:,:,:)    !R_LG( nneq,neqv,3,3)
      Real(kind=wp), allocatable    :: R_ROT(:,:,:,:)    !R_ROT(nneq,neqv,3,3)
      Real(kind=wp), allocatable    :: eso(:,:)          ! spin orbit energies on individual metal sites
      Complex(kind=wp), allocatable :: dipso(:,:,:,:)
      Complex(kind=wp), allocatable ::  s_so(:,:,:,:)
      Character(1)                  :: itype(nneq)
      Character(180)                :: namefile_aniso(nneq)
      Logical                       :: ifHDF
      Logical                       :: DoPlot
c  definition of the exchange:
      Integer                       :: exch                           ! total number of exchange states
      Integer                       :: nPair                          ! number of metal pairs (number of interactions)
      Integer                       :: nmax                           ! exchange basis, nmax = MAXVAL(nexch(:))
      Integer                       :: MxRank1, MxRank2
      Integer, allocatable          :: nexch(:)
      Integer, allocatable          :: i_pair(:,:)            ! index of the metal site in a given interacting pair
      Integer, allocatable          :: imaxrank(:,:)          ! index of rank ITOs for JITO exchange definition
      Logical                       :: AnisoLines1
      Logical                       :: AnisoLines3
      Logical                       :: AnisoLines9
      Logical                       :: Dipol
      Logical                       :: DM_exchange
      Logical                       :: decompose_exchange
      Logical                       :: JITO_exchange
      Real(kind=wp), allocatable    :: Jex(:)            ! Lines exchange    ( 1 parameter / interacting pair)
      Real(kind=wp), allocatable    :: JAex(:,:)         ! Anisotropic Lines ( 3 parameter / interacting pair)
      Real(kind=wp), allocatable    :: JAex9(:,:,:)      ! Anisotropic Lines full ( 9 parameters / interacting pair)
      Real(kind=wp), allocatable    :: JDMex(:,:)         ! Anisotropic Lines ( 3 parameter / interacting pair)
      Real(kind=wp), allocatable    :: JITOexR(:,:,:,:,:)
      Real(kind=wp), allocatable    :: JITOexI(:,:,:,:,:)
      Real(kind=wp), allocatable    :: W(:)              ! exchange spectrum
      Complex(kind=wp), allocatable :: Z(:,:)         ! exchange eigenstates
      Complex(kind=wp), allocatable :: dipexch(:,:,:) ! exchange magnetic moment
      Complex(kind=wp), allocatable ::  s_exch(:,:,:) ! exchange spin moment
      ! options used in connection with KE
      Integer                       :: lant, multLn, KEOPT
      Logical                       :: KE
      Real(kind=wp)                 :: tpar, upar
      ! options used in connection with Dipol-Dipol interaction
      Real(kind=wp), allocatable    :: MagnCoords(:,:) !MagnCoords(nneq,3)

      Integer                       :: nTempMagn
c  definition of g and D tensors
      Logical                       :: compute_g_tensors
      Integer                       :: nMult
      Integer, allocatable          :: nDim(:) ! multiplicity of each multiplet
      Real(kind=wp), allocatable    :: gtens(:,:) ! gtensor of each multiplet
      Real(kind=wp), allocatable    :: maxes(:,:,:) ! main axes of each multiplet
c  definition of data for susceptibility
      Integer                       :: nT
      Logical                       :: compute_susceptibility
      Logical                       :: tinput
      Real(kind=wp)                 :: tmin, tmax, dltT0

      Real(kind=wp), allocatable    ::     T(:)
      Real(kind=wp), allocatable    :: XTexp(:)
      Real(kind=wp), allocatable    :: XT_no_field(:)
      Real(kind=wp), allocatable    :: chit_exp(:)
      Real(kind=wp), allocatable    :: Texp(:)

      Real(kind=wp), allocatable    :: XLM(:,:,:,:)
      Real(kind=wp), allocatable    :: ZLM(:,:)
      Real(kind=wp), allocatable    :: XRM(:,:,:,:)
      Real(kind=wp), allocatable    :: ZRM(:,:)
      ! options related to XT_MoverH
      Real(kind=wp)                 :: Xfield

c  definition of data for magnetization:
      Integer                       :: nH, nM
      Integer                       :: iopt
      Real(kind=wp), allocatable    :: TempMagn(:)
      Real(kind=wp), allocatable    :: Hexp(:), Mexp(:,:)
      Real(kind=wp)                 :: thrs, em
      Real(kind=wp)                 :: hmin, hmax, dltH0
      Logical                       :: hinput
      Logical                       :: compute_magnetization
      Logical                       :: compute_Mdir_vector
      Logical                       :: zeeman_energy
      Logical                       :: m_paranoid
      Logical                       :: m_accurate
      Logical                       :: smagn
      ! options used to set up nM and EM
      Integer                       :: encut_definition
      Integer                       :: nK, mG ! encut_definition=1;
      Integer                       :: ncut   ! encut_definition=2;
      Real(kind=wp)                 :: encut_rate ! encut_definition=3;

c  magnetization torque
      Logical                       :: compute_torque
      Integer                       :: AngPoints, nP
c  Zeeman energy and M vector
      Integer                       :: nDir, nDirZee, nDirTot
      Integer, allocatable          :: LuZee(:)
      Real(kind=wp), allocatable    :: dirX(:), dirY(:), dirZ(:)
      Real(kind=wp), allocatable    :: dir_weight(:,:)
c  definition of mean field parameter
      Real(kind=wp)                 :: zJ
c  definintion of the crystal axes:
      Logical                       :: Do_structure_abc
      Real(kind=wp)                 :: cryst(6) ! a, b, c, alpha, beta, gamma
      Real(kind=wp)                 :: coord(3) ! Cartesian coordinates of the main metal site, or center
c  definitions for blocking barrier
      Integer                       :: nBlock
      Logical                       :: compute_barrier
c  options for automatic fitting of parameters:
      Logical                       :: fitCHI !-- not used so far
      Logical                       :: fitM !-- not used so far
c  fundamental constants:
      Real(kind=wp)                 :: boltz_k
      Real(kind=wp)                 :: mu_bohr

      Integer                       :: iPrint
      Integer                       :: idim
      Integer                       :: imltpl
      Integer                       :: Ifunct
      Integer                       :: iReturn
      Integer                       :: i,i1,i2,j,l
      Integer                       :: l1(2),l2(2),l3(2),l4(2),l5(2)
      Integer                       :: imanIfold
      Integer                       :: ibuf
c      Integer                      :: nsta
c      Integer                      :: icase, nmagmult
      Logical                       :: check_title
      Character(180)                :: Title
      Logical                       :: GRAD
      Logical                       :: lvM
      Logical                       :: dbg
      Character(len=180)            :: fname

      Call qEnter('PA_1')
      dbg=.false.
c---------------------------------------------------------------------
      ! Constants:
      boltz_k=0.6950356_wp                    !   in cm^-1*K-1
      mu_bohr=0.466864374_wp                  !   in cm-1*T-1
      GRAD=.false.
c---------------------------------------------------------------------
      ! Allocate memory for all arrays:
c---------------------------------------------------------------------
      If(dbg) Write(6,*) '      exch = ',exch
      If(dbg) Write(6,*) '     nPair = ',nPair
      If(dbg) Write(6,*) '     nMult = ',nMult
      If(dbg) Write(6,*) '      nneq = ',nneq
      If(dbg) Write(6,*) '      neqv = ',neqv
      If(dbg) Write(6,*) '      nLoc = ',nLoc
      If(dbg) Write(6,*) '   nDirZee = ',nDirZee
      If(dbg) Write(6,*) '      nDir = ',nDir
      If(dbg) Write(6,*) '        nH = ',nH
      If(dbg) Write(6,*) '        nT = ',nT
      If(dbg) Write(6,*) ' nTempMagn = ',nTempMagn

      mem=0
      RtoB=8
      CtoB=16
      ItoB=8

      If(exch>0) Then
        ! exchange energy spectrum
        Call mma_allocate(W,exch,'W')
        Call dcopy_(exch,0.0_wp,0,W,1)
        mem=exch*RtoB
        ! exchange eigenvectors:
        Call mma_allocate(Z,exch,exch,'Z')
        Call zcopy_(exch*exch,(0.0_wp,0.0_wp),0,Z,1)
        mem=mem+exch*exch*CtoB
        ! magnetic moment
        Call mma_allocate(dipexch,3,exch,exch,'dipexch')
        Call zcopy_(3*exch*exch,(0.0_wp,0.0_wp),0,dipexch,1)
        mem=mem+3*exch*exch*CtoB
        ! spin moment
        Call mma_allocate( s_exch,3,exch,exch,'s_exch')
        Call zcopy_(3*exch*exch,(0.0_wp,0.0_wp),0, s_exch,1)
        mem=mem+3*exch*exch*CtoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 1 =',mem
      End If

      If(nPair>0) Then
        ! index of metal site for each interacting pair:
        Call mma_allocate(i_pair,nPair,2,'i_pair')
        Call icopy(nPair*2,0,0,i_pair,1)
        mem=mem+nPair*2*ItoB
        ! exchange Lines-1 parameter
        Call mma_allocate(Jex,nPair,'Jex')
        Call dcopy_(nPair,0.0_wp,0,Jex,1)
        mem=mem+nPair*RtoB
        ! exchange Lines-3 parameter
        Call mma_allocate(JAex,nPair,3,'Jex')
        Call dcopy_(nPair*3,0.0_wp,0,JAex,1)
        mem=mem+nPair*3*RtoB
        ! exchange Lines-9 parameter
        Call mma_allocate(JAex9,nPair,3,3,'Jex')
        Call dcopy_(nPair*3*3,0.0_wp,0,JAex9,1)
        mem=mem+nPair*3*3*RtoB
        ! exchange Dzyaloshinsky-Morya parameter
        Call mma_allocate(JDMex,nPair,3,'Jex')
        Call dcopy_(nPair*3,0.0_wp,0,JDMex,1)
        mem=mem+nPair*3*RtoB
        ! index of ITO ranks for for each interacting pair:
        Call mma_allocate(imaxrank,nPair,2,'imaxrank')
        Call icopy(nPair*2,0,0,imaxrank,1)
        mem=mem+nPair*2*ItoB
        ! exchange JITO
        If((MxRank1>0).AND.(MxRank2>0)) Then
           l1(1)=1;         l1(2)=nPair
           l2(1)=1;         l2(2)=MxRank1
           l3(1)=-MxRank1;  l3(2)=MxRank1
           l4(1)=1;         l4(2)=MxRank2
           l5(1)=-MxRank2;  l5(2)=MxRank2
           ibuf=nPair*MxRank1*MxRank2*(2*MxRank1+1)*(2*MxRank2+1)
           If(dbg) Write(6,'(A,I10)') 'nPair  =',nPair
           If(dbg) Write(6,'(A,I10)') 'MxRank1=',MxRank1
           If(dbg) Write(6,'(A,I10)') 'MxRank2=',MxRank2
           If(dbg) Write(6,'(A,I10)') 'ibuf   =',ibuf
           Call xFlush(6)
           Call mma_allocate(JITOexR,l1,l2,l3,l4,l5,'JITOexR')
           Call mma_allocate(JITOexI,l1,l2,l3,l4,l5,'JITOexI')
           Call dcopy_(ibuf,0.0_wp,0,JITOexR,1)
           Call dcopy_(ibuf,0.0_wp,0,JITOexI,1)
           mem=mem+2*ibuf*RtoB
        End If
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 2 =',mem
      End If

      If(nMult>0) Then
        ! index of metal site for each interacting pair:
        Call mma_allocate(nDim,nMult,'nDim')
        Call icopy(nMult,0,0,nDim,1)
        mem=mem+nMult*ItoB
        ! g_tensor for each multiplet
        Call mma_allocate(gtens,nMult,3,'gtens')
        Call dcopy_(nMult*3,0.0_wp,0,gtens,1)
        mem=mem+nMult*3*RtoB
        ! main axes of the g tensor for each multiplet
        Call mma_allocate(maxes,nMult,3,3,'maxes')
        Call dcopy_(nMult*3*3,0.0_wp,0,maxes,1)
        mem=mem+nMult*3*3*RtoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 3 =',mem
      End If

      If(nneq>0) Then
        ! number of equivalent centers, per type
        Call mma_allocate(neq,nneq,'neq')
        Call icopy(nneq,0,0,neq,1)
        mem=mem+nneq*ItoB
        ! local number of spin orbit states:
        Call mma_allocate(nss,nneq,'nss')
        Call icopy(nneq,0,0,nss,1)
        mem=mem+nneq*ItoB
        ! local number of spin free states:
        Call mma_allocate(nsfs,nneq,'nsfs')
        Call icopy(nneq,0,0,nsfs,1)
        mem=mem+nneq*ItoB
        ! local basis for exchange:
        Call mma_allocate(nexch,nneq,'nexch')
        Call icopy(nneq,0,0,nexch,1)
        mem=mem+nneq*ItoB
        ! gtens_input(:,:)
        Call mma_allocate(gtens_input,3,nneq,'gtens_input')
        Call dcopy_(3*nneq,0.0_wp,0,gtens_input,1)
        mem=mem+3*nneq*RtoB
        ! D_fact
        Call mma_allocate(D_fact,nneq,'D_fact')
        Call dcopy_(nneq,0.0_wp,0,D_fact,1)
        mem=mem+nneq*RtoB
        ! EoverD_fact
        Call mma_allocate(EoverD_fact,nneq,'EoverD_fact')
        Call dcopy_(nneq,0.0_wp,0,EoverD_fact,1)
        mem=mem+nneq*RtoB
        ! MagnCoords()
        Call mma_allocate(MagnCoords,nneq,3,'MagnCoords')
        Call dcopy_(nneq*3,0.0_wp,0,MagnCoords,1)
        mem=mem+nneq*3*RtoB
        ! riso
        Call mma_allocate(riso,nneq,3,3,'riso')
        Call dcopy_(3*3*nneq,0.0_wp,0,riso,1)
        mem=mem+9*nneq*RtoB


        If(neqv>0) Then
          ! R_LG
          Call mma_allocate(r_lg,nneq,neqv,3,3,'r_lg')
          Call dcopy_(nneq*neqv*3*3,0.0_wp,0,r_lg,1)
          mem=mem+nneq*neqv*3*3*RtoB
          ! R_ROT
          Call mma_allocate(r_rot,nneq,neqv,3,3,'r_rot')
          Call dcopy_(nneq*neqv*3*3,0.0_wp,0,r_rot,1)
          mem=mem+nneq*neqv*3*3*RtoB
        End If

        If(nLoc>0) Then
          ! local spin orbit states
          Call mma_allocate(eso,nneq,nLoc,'eso')
          Call dcopy_( nneq*nLoc,0.0_wp,0,eso,1)
          mem=mem+nneq*nLoc*RtoB
          ! local magnetic moment
          Call mma_allocate(dipso,nneq,3,nLoc,nLoc,'dipso')
          Call zcopy_(3*nneq*nLoc*nLoc,(0.0_wp,0.0_wp),0,dipso,1)
          mem=mem+3*nneq*nLoc*nLoc*CtoB
          ! local spin moment
          Call mma_allocate( s_so,nneq,3,nLoc,nLoc,'s_so')
          Call zcopy_(3*nneq*nLoc*nLoc,(0.0_wp,0.0_wp),0, s_so,1)
          mem=mem+3*nneq*nLoc*nLoc*CtoB
        End If

        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 4 =',mem
      End If

      If(nDirZee>0) Then
        ! unit numbers of the files containing Zeeman states
        Call mma_allocate(LuZee,nDirZee,'LuZee')
        Call icopy(nDirZee,0,0,LuZee,1)
        mem=mem+nDirZee*ItoB
        ! directional vectors for the magnetic field, for computing Zeeman states
        Call mma_allocate(dir_weight,nDirZee,3,'dir_weight')
        Call dcopy_(3*nDirZee,0.0_wp,0,dir_weight,1)
        mem=mem+3*nDirZee*ItoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 5 =',mem
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
        If(dbg) Write(6,'(A,I16)') 'mem 6 =',mem
      End If

      If((nH>0).and.(nTempMagn>0)) Then
        ! experimental field points:
        Call mma_allocate(Hexp,nH,'Hexp')
        Call dcopy_(nH,0.0_wp,0,Hexp,1)
        mem=mem+nH*RtoB
        ! experimental magnetisation:
        Call mma_allocate(Mexp,nH,nTempMagn,'Mexp')
        Call dcopy_(nH*nTempMagn,0.0_wp,0,Mexp,1)
        mem=mem+nH*nTempMagn*RtoB
        ! temperature points for magnetization:
        Call mma_allocate(TempMagn,nTempMagn,'TempMagn')
        Call dcopy_(nTempMagn,0.0_wp,0,TempMagn,1)
        mem=mem+nTempMagn*RtoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 7 =',mem
      End If

      If((nCenter>0).and.(nTempMagn>0)) Then
        ! XT for local centers, all states
        Call mma_allocate(XLM,nCenter,nTempMagn,3,3,'XLM')
        Call dcopy_(nCenter*nTempMagn*3*3,0.0_wp,0,XLM,1)
        mem=mem+nCenter*nTempMagn*3*3*RtoB
        ! Statistical sum for local centers, all states
        Call mma_allocate(ZLM,nCenter,nTempMagn,'ZLM')
        Call dcopy_(nCenter*nTempMagn,0.0_wp,0,ZLM,1)
        mem=mem+nCenter*nTempMagn*RtoB
        ! XT for local centers, exchange states
        Call mma_allocate(XRM,nCenter,nTempMagn,3,3,'XRM')
        Call dcopy_(nCenter*nTempMagn*3*3,0.0_wp,0,XRM,1)
        mem=mem+nCenter*nTempMagn*3*3*RtoB
        ! Statistical sum for local centers, exchange states
        Call mma_allocate(ZRM,nCenter,nTempMagn,'ZRM')
        Call dcopy_(nCenter*nTempMagn,0.0_wp,0,ZRM,1)
        mem=mem+nCenter*nTempMagn*RtoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 8 =',mem
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
      End If

      If((nT+nTempMagn)>0) Then
        ! -- add nTempMagn points, so that all measurables are computed at once...
        ! temperature points for which XT will be computed
        Call mma_allocate(T,(nTempMagn+nT),'Temperature')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,T,1)
        mem=mem+(nT+nTempMagn)*RtoB
        ! XT experimental tmagn XTexp+Tmagn
        Call mma_allocate(XTexp,(nTempMagn+nT),'XTexp')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,XTexp,1)
        mem=mem+(nT+nTempMagn)*RtoB
        ! XT in the absence of the magnetic field
        Call mma_allocate(XT_no_field,(nTempMagn+nT),'XT_no_field')
        Call dcopy_((nT+nTempMagn),0.0_wp,0,XT_no_field,1)
        mem=mem+(nT+nTempMagn)*RtoB
        ! allocated memory counter
        If(dbg) Write(6,'(A,I16)') 'mem 9 =',mem
      End If

      Write(6,'(A,I16,A)') 'The code allocated at least:',mem,
     &                     ' bytes of memory for this run.'
      Call xFlush(6)
c---------------------------------------------------------------------
      ! set default values of the main variables and arrays:
      If(dbg) Write(6,*) 'Enter set_defaults'
        Call set_defaults( nneq, nTempMagn, nDir, nDirZee, nMult,
     &                     neq, nexch, nK, mG, ncut, nP,
     &                     AngPoints, nBlock, encut_definition,
     &                     iopt, iPrint,
     &                     dltT0, dltH0, zJ, tmin, tmax, hmin, hmax,
     &                     XField, thrs, TempMagn,
     &                     cryst, coord, encut_rate, gtens_input,
     &                     D_fact, EoverD_fact, riso,
     &                     decompose_exchange, AnisoLines1,
     &                     AnisoLines3, AnisoLines9, DM_exchange,
     &                     Dipol, KE, JITO_exchange, fitCHI, fitM,
     &                     Do_structure_abc, doplot,
     &                     compute_g_tensors, tinput,
     &                     compute_susceptibility,
     &                     compute_torque, compute_barrier,
     &                     compute_magnetization, hinput,
     &                     compute_Mdir_vector, zeeman_energy,
     &                     m_paranoid, m_accurate, smagn,
     &                     itype )
      If(dbg) Write(6,*) 'Exit set_defaults'
      If(dbg.and.(neqv.gt.1)) Write(6,*) 'Enter fetch_neq'
      If(neqv.gt.1) Call fetch_neq(nneq,neq(1:nneq),nexch(1:nneq))
      If(dbg.and.(neqv.gt.1)) Write(6,*) 'Exit fetch_neq'

c---------------------------------------------------------------------
c      ! read the Standard Input:
      If(dbg) Write(6,*) 'Enter Readin_poly'
      Call Readin_poly(nneq, neq, neqv, exch, nCenter,
     &                 nT, nH, nTempMagn, nDir, nDirZee,
     &                 nMult, nPair, nexch, nDim, i_pair,
     &                 lant, multLn, iPrint, keopt,
     &                 encut_definition, nK, mG, iopt, nP,
     &                 AngPoints, ncut, LUZee, MxRank1, MxRank2,
     &                 imaxrank,
     &                 TempMagn, R_LG, R_ROT, Jex, JAex, JAex9, JDMex,
     &                 JITOexR, JITOexI,
     &                 tpar, upar, cryst, coord,
     &                 Xfield, gtens_input, D_fact, EoverD_fact, riso,
     &                 MagnCoords, thrs, tmin, tmax,
     &                 hmin, hmax, Texp(1:nT), chit_exp(1:nT),
     &                 Hexp, Mexp, encut_rate,
     &                 zJ, dirX, dirY, dirZ,
     &                 dir_weight,
     &                 Title, itype,
     &                 ifHDF,
     &                 compute_g_tensors, compute_magnetization,
     &                 TINPUT, HINPUT, Do_structure_abc, doplot,
     &                 compute_Mdir_vector, zeeman_energy,
     &                 m_paranoid, m_accurate, smagn,
     &                 compute_susceptibility, decompose_exchange,
     &                 KE, fitCHI, fitM, compute_torque,
     &                 compute_barrier, Dipol, check_title,
     &                 AnisoLines1, AnisoLines3, AnisoLines9,
     &                 DM_exchange, JITO_exchange )
      If(dbg) Write(6,*) 'Exit Readin_poly'
      Call xFlush(6)
c---------------------------------------------------------------------
c     ! fetch the data from aniso_x.input files: (formatted ANISOINPUT)
      Do i=1,nneq
         If(dbg) Write(6,'(A,A)') 'itype(i)=',itype(i)
         If( (itype(i).eq.'B').OR.(itype(i).eq.'C') ) Then

            If(dbg) Write(6,*) 'Enter generate_isotrop_site'
            Call generate_isotrop_site( nss(i), nsfs(i), nexch(i),
     &                                  nLoc, gtens_input(1:3,i),
     &                                  riso(i,1:3,1:3), D_fact(i),
     &                                  EoverD_fact(i),
     &                                  eso(i,1:nexch(i)),
     &                             dipso(i,1:3,1:nexch(i),1:nexch(i)),
     &                              s_so(i,1:3,1:nexch(i),1:nexch(i)) )
            If(dbg) Write(6,*) 'Exit generate_isotrop_site'
            If(dbg) Call xFlush(6)

         Else If (itype(i).eq.'A') Then

            If(ifHDF) Then
               ! set their names:
               If(i<10) Then
                  Write(namefile_aniso(i),'(4A)') 'aniso_hdf_',
     &                      CHAR(48+mod( int( i     ),10)),'.input'
                  If(dbg) Write(6,'(A,i2,A,A)') 'namefile_aniso(',i,
     &              ')=', namefile_aniso(i)
               Else If(i>=10 .and. i<=99) Then
                  Write(namefile_aniso(i),'(4A)') 'aniso_hdf_',
     &                      CHAR(48+mod( int((i)/10 ),10)),
     &                      CHAR(48+mod( int( i     ),10)),'.input'
               End If
               If(dbg) Write(6,*) 'PA:  namefile_aniso(i)=',
     &                             trim(namefile_aniso(i))

#ifdef _HDF5_
               If(dbg) Write(6,*) 'Enter read_hdf5_poly'
               Call read_hdf5_init(namefile_aniso(i),nsfs(i),nss(i))
               If(dbg) Write(6,*) 'i=', i,'nsfs(i),nss(i)=',
     &                                     nsfs(i),nss(i)
               Call read_hdf5_poly( namefile_aniso(i),
     &                              nss(i),
     &                              nsfs(i),
     &                              eso(i,1:nss(i)),
     &                             dipso(i,1:3,1:nss(i),1:nss(i)),
     &                              s_so(i,1:3,1:nss(i),1:nss(i)),
     &                              iReturn )
               If(dbg) Write(6,*) 'Exit read_hdf5_poly'
               If(dbg) Write(6,*) 'ESO(i)=',(ESO(i,j),j=1,nss(i))
#else
               Call WarningMessage(2,'File '//trim(namefile_aniso(i))//
     &                              ' cannot be opened. Molcas was'//
     &                              ' compiled without HDF5 option.')
               Call Quit_OnUserError()
#endif
            Else
               ! set their names:
               If(i<10) Then
                  Write(namefile_aniso(i),'(4A)') 'aniso_',
     &                      CHAR(48+mod( int( i     ),10)),'.input'
               Else If(i>=10 .and. i<=99) Then
                  Write(namefile_aniso(i),'(4A)') 'aniso_',
     &                      CHAR(48+mod( int((i)/10 ),10)),
     &                      CHAR(48+mod( int( i     ),10)),'.input'
               End If
               ! get the information from formatted aniso.input file:
               If(dbg) Write(6,*) 'Enter read_formatted_aniso_poly'
               Call read_formatted_aniso_poly(
     &                              namefile_aniso(i),
     &                              nss(i),
     &                              nsfs(i),
     &                              nLoc,
     &                              eso(i,1:nLoc),
     &                             dipso(i,1:3,1:nLoc,1:nLoc),
     &                              s_so(i,1:3,1:nLoc,1:nLoc),
     &                              iReturn )
               If(dbg) Write(6,*) 'Exit read_formatted_aniso_poly'
            End if ! ifHDF
         Else
            Call quit(_RC_INTERNAL_ERROR_)
         End If
      End Do

      If(dbg) Write(6,*) 'Enter input_process'
      Call input_process( nneq, neq, neqv, nmax, nCenter, nexch,
     &                    nDir, nDirZee, nDirTot, nH, nT, exch,
     &                    nBlock, nTempMagn, iopt, nMult, nDim,
     &                    nPair, i_pair, nP, AngPoints, lant,
     &                    multLn, KEOPT, encut_definition, nK,
     &                    mG, ncut, nsfs, nss, nLoc,
     &                    MxRank1, MxRank2, imaxrank, iPrint,

     &                    R_LG, gtens_input, dirX, dirY, dirZ,
     &                    dir_weight, zJ, cryst, coord, hmin,
     &                    hmax, TempMagn, thrs, Hexp, Mexp,
     &                    tmin, tmax, chit_exp, Texp, Xfield,
     &                    Jex, JAex, JAex9, JDMex, tpar, upar,
     &                    MagnCoords, encut_rate, eso,
     &                    JITOexR, JITOexI,

     &                    Title, itype, namefile_aniso,

     &                    Do_structure_abc,
     &                    compute_barrier, fitCHI, fitM, hinput,
     &                    tinput, compute_magnetization,
     &                    compute_Mdir_vector, zeeman_energy,
     &                    m_paranoid, m_accurate, smagn,
     &                    compute_g_tensors,compute_torque,
     &                    compute_susceptibility, AnisoLines1,
     &                    AnisoLines3, AnisoLines9, Dipol,
     &                    check_title, KE, DM_exchange, JITO_exchange )
      If(dbg) Write(6,*) 'Exit input_process'

c at this moment all input values are set. proceed to compute various
c properties:
c---------------------------------------------------------------------
c in case J coupling parameters are not known, proceed to get them from
c a minimization scheme between caclulated and measured magnetic susceptibility
c
c the function below optimizes the N+1 parameters:
c N - exchange couplings J, and one parameter for total shIft of all experimental points;
c      If (fitCHI) Then
c      Call fitCHI()
c      End If
cc at this point all J parameters are known; we can proceed to compute
c the energy spectrum and the resulting properties
c---------------------------------------------------------------------
c      Lines=.true.
c      If(AnisoLines.eqv..true.) Then
c      Lines=.false.
c      End If
c      Dipol=dipole_included
c      KE=exch_long


      If(dbg) Write(6,*) 'Dipol         = ', Dipol
      If(dbg) Write(6,*) 'AnisoLines1   = ', AnisoLines1
      If(dbg) Write(6,*) 'AnisoLines3   = ', AnisoLines3
      If(dbg) Write(6,*) 'AnisoLines9   = ', AnisoLines9
      If(dbg) Write(6,*) 'DM_exchange   = ', DM_exchange
      If(dbg) Write(6,*) 'JITO_exchange = ', JITO_exchange
      If(dbg) Write(6,*) 'nmax          = ', nmax

      Call exchctl( exch, nneq, neqv, neq, nexch,
     &              nmax, nCenter, npair, i_pair,
     &              MxRank1, MxRank2, imaxrank,
     &              Jex, JAex, JAex9, JDMex, JITOexR, JITOexI,
     &              eso(1:nneq,1:nmax),
     &              s_so(1:nneq,1:3,1:nmax,1:nmax),
     &              dipso(1:nneq,1:3,1:nmax,1:nmax),
     &             MagnCoords, R_ROT, R_LG, riso, tpar, upar, lant,
     &              itype, Dipol, AnisoLines1, AnisoLines3,
     &              AnisoLines9, KE, KEOPT, DM_exchange, JITO_exchange,
     &               W, Z, S_EXCH, DIPEXCH, iPrint, mem )
cc---------------------------------------------------------------------
!     generate an 'aniso -like file' for testign exchange proj in SA
      fname='ANISOINPUT_POLY'
      Call write_formatted_aniso_poly( fname, exch, W,
     &                                  dipexch, s_exch)
cc---------------------------------------------------------------------
c compute the magnetic moments on individual metal sites, for
c the lowest NSTA exchange states;
C     NMAX = maximal number of local states which participate into the exchange coupling
      Call MOMLOC2( exch, nmax, nneq, neq, neqv,
     &              r_rot, nCenter, nExch, W, Z,
     &              dipexch, s_exch,
     &              dipso(1:nneq,1:3,1:nmax,1:nmax),
     &               s_so(1:nneq,1:3,1:nmax,1:nmax)  )

c---------------------------------------------------------------------
      If(compute_g_tensors) Then
         If(dbg) Write(6,'(A,90I3)') 'nmult= ',nmult
         If(dbg) Write(6,'(A,90I3)') 'ndim() = ',ndim(1:nmult)

         Ifunct=0
         Do imltpl=1,nmult
            iReturn=0
            i1=1           +Ifunct
            i2=ndim(imltpl)+Ifunct

            If (ndim(imltpl)<=1) Go To 10
            If (i2>exch) Go To 20

            If(dbg) Write(6,'(A,90I3)') 'ndim(imltpl)=',ndim(imltpl)
            If(dbg) Call prmom('PA: s_exch:',
     &               s_exch(1:3,i1:i2,i1:i2),ndim(imltpl))
            If(dbg) Call prmom('PA: dip_exch:',
     &              dipexch(1:3,i1:i2,i1:i2),ndim(imltpl))

            Call g_high( w(i1:i2), GRAD,
     &                    s_exch(1:3, i1:i2, i1:i2),
     &                   dipexch(1:3, i1:i2, i1:i2),
     &                   imltpl, ndim(imltpl),
     &                   Do_structure_abc,
     &                   cryst, coord,
     &                   gtens(imltpl,1:3),
     &                   maxes(imltpl,1:3,1:3),
     &                   iprint )

10          Continue
            Ifunct=Ifunct+ndim(imltpl)
            If( iReturn.ne.0) stop
         End Do ! imltpl
20       Continue
      End If

!----------------------------------------------------------------------|
      !   >> AB INITIO BLOCKING BARRIER FOR SMMs <<
      If ( compute_barrier ) Then
        If ( nBlock.ne.0 ) Then
          imanifold=1

          Write(6,'(A)') 'UBAR:: matrix elements of the (input) '//
     &                  'magnetic moment '
          Write(6,'(A)') '                |'//
     &                   '         - projection - X -        |'//
     &                   '         - projection - Y -        |'//
     &                   '         - projection - Z -        |'
          Write(6,'(A)') '----------------|'//
     &                   '----- Real ------|---- Imaginary --|'//
     &                   '----- Real ------|---- Imaginary --|'//
     &                   '----- Real ------|---- Imaginary --|'
          Do i=1,nBlock
            Do j=1,i
            Write(6,'(A,i2,A,i2,A,6ES18.10)') '< ',i,'|M_xyz|',j,' > ',
     &                                        (dipexch(l,i,j),l=1,3)
            End Do
          End Do

          Call  BARRIER( nBlock, dipexch(1:3,1:nBlock,1:nBlock),
     &                   W(1:nBlock), imanIfold, nMult, nDim(1:nMult),
     &                   iprint )

        Else
           Write(6,'(A)') 'nBlock parameter is not defined. '
           Write(6,'(A)') 'Did you specify the MLTP keyword in the '//
     &                    'input?'
           Write(6,'(A)') 'If the problem persists,please, '//
     &                    'submit a bug report.'
        End If
      End If



!---------------------------------------------------------------------
!             MAGNETIC SUSCEPTIBILITY
!---------------------------------------------------------------------
       If (compute_susceptibility .AND. (nT>0) ) Then
          lvM =.false.
          lvM = compute_magnetization .OR. (Xfield.ne.0.0_wp) .OR.
     &          compute_torque

         ! set nT, T(i) and XTexp(i) arrays:
          Call set_T( nT, nTempMagn, TINPUT, TempMagn, Tmin, Tmax,
     &                chit_exp, Texp, T, XTexp )

          ! XT at H=0
          Call susceptibility( exch, nLoc, nCenter, nneq,
     &                         neqv, neq, nss, nexch, nTempMagn,
     &                         nT, Tmin, Tmax, XTexp,
     &                         eso, dipso, s_so, W, dipexch,
     &                         s_exch, T, R_LG, zJ, tinput,
     &                         XLM, ZLM, XRM, ZRM, iopt, XT_no_field,
     &                         DoPlot, mem )

         If( Xfield.ne.0.0_wp ) Then
          ! XT at H>0:

           Call set_nm( exch, ncut, encut_definition, nk, mg,
     &                  nTempMagn, hmax, w, encut_rate, TempMagn,
     &                  nM, EM, dbg )
          Do i=1,nT+nTempMagn
             Write(6,*) i,  T(i)
          End Do
          If(dbg) Write(6,*) 'nm    =',nm
          If(dbg) Write(6,*) 'Tmin  =',Tmin
          If(dbg) Write(6,*) 'Tmax  =',Tmax
          If(dbg) Write(6,*) 'THRS  =',THRS
          If(dbg) Write(6,*) 'smagn =',smagn
          If(dbg) Write(6,*) 'm_par =',m_paranoid
          If(dbg) Write(6,*) 'm_acc =',m_accurate
          If(dbg) Write(6,*) 'tinpu =',tinput

           Call XT_dMoverdH( exch, nLoc, nCenter, nneq, neqv, neq,
     &                       nss, nexch, nTempMagn, nT, NM, iopt, mem,
     &                       Tmin, Tmax, XTexp, eso, w, T, R_ROT,
     &                       zJ, Xfield, EM, THRS, XT_no_field,
     &                       dipso, s_so, dipexch, s_exch,
     &                       tinput, smagn, m_paranoid, m_accurate )
         End If
      Else
         Write(6,'(A)') 'Computation of the magnetic susceptibility'//
     &                  '... skipped by the user'
       End If
!---------------------------------------------------------------------
!             MAGNETIC TORQUE
!---------------------------------------------------------------------

       If (compute_torque) Then
         ! set the NM- number of states to be exactly diagonalized in Zeeman Interaction
         Call set_nm( exch, ncut, encut_definition, nk, mg,
     &                nTempMagn, hmax, w, encut_rate, TempMagn,
     &                nM, EM, dbg )

         Call torque( nneq, nCenter, neq, neqv, nLoc, exch,
     &                nTempMagn, nH, nM, AngPoints, nexch,
     &                iopt, nss, mem,
     &                smagn, m_paranoid, m_accurate,
     &                TempMagn, w, hmin, hmax, dltH0, EM, zJ, THRS,
     &                hexp,
     &                dipexch, s_exch, dipso, s_so, eso,
     &                hinput, r_rot, XLM, ZLM, XRM, ZRM )

      Else
         Write(6,'(A)') 'Computation of the magnetization torque ... '//
     &                  'skipped by the user'
      End If !compute_torque
!---------------------------------------------------------------------
!             MOLAR MAGNETIZATION
!---------------------------------------------------------------------
      If(compute_magnetization .AND. (nH>0) ) Then
         ! set the NM- number of states to be exactly diagonalized in Zeeman Interaction
         Call set_nm( exch, ncut, encut_definition, nk, mg,
     &                nTempMagn, hmax, w, encut_rate, TempMagn,
     &                nM, EM, dbg )

         Call magnetization( exch, nLoc, nM, nH, nneq, neq, neqv,
     &                       nCenter, nTempMagn, nDir, nDirZee,
     &                       nDirTot, nss, nexch, iopt, LUZee,
     &                       TempMagn, hexp, mexp, hmin, hmax, em,
     &                       zJ, thrs, dirX, dirY, dirZ, dir_weight,
     &                       w, dipexch, s_exch, dipso, s_so, eso,
     &                       hinput, r_rot, XLM, ZLM, XRM, ZRM,
     &                       zeeman_energy, compute_Mdir_vector,
     &                       m_paranoid, m_accurate, smagn, mem )
      Else
         Write(6,'(A)') 'Computation of the molar magnetization ... '//
     &                  'skipped by the user'
      End If !compute_magnetization

!---------------------------------------------------------------------
! Deallocate memory for big arrays:
!---------------------------------------------------------------------
      If(exch>0) Then
        Call mma_deallocate(W)
        Call mma_deallocate(Z)
        Call mma_deallocate(dipexch)
        Call mma_deallocate(s_exch)
      End If

      If(nPair>0) Then
        Call mma_deallocate(i_pair)
        Call mma_deallocate(Jex)
        Call mma_deallocate(JAex)
        Call mma_deallocate(JAex9)
        Call mma_deallocate(JDMex)
        Call mma_deallocate(imaxrank)
        If((MxRank1>0).AND.(MxRank2>0)) Then
          Call mma_deallocate(JITOexR)
          Call mma_deallocate(JITOexI)
        End If
      End If

      If(nMult>0) Then
        Call mma_deallocate(nDim)
        Call mma_deallocate(gtens)
        Call mma_deallocate(maxes)
      End If

      If(nneq>0) Then
        Call mma_deallocate(neq)
        Call mma_deallocate(nss)
        Call mma_deallocate(nsfs)
        Call mma_deallocate(nexch)
        Call mma_deallocate(gtens_input)
        Call mma_deallocate(D_fact)
        Call mma_deallocate(EoverD_fact)
        Call mma_deallocate(MagnCoords)
        Call mma_deallocate(riso)
        If(neqv>0) Then
          Call mma_deallocate(r_lg)
          Call mma_deallocate(r_rot)
        End If
        If(nLoc>0) Then
          Call mma_deallocate(eso)
          Call mma_deallocate(dipso)
          Call mma_deallocate(s_so)
        End If
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

      If((nH>0).and.(nTempMagn>0)) Then
        ! experimental field points:
        Call mma_deallocate(Hexp)
        Call mma_deallocate(Mexp)
        Call mma_deallocate(TempMagn)
      End If

      If((nCenter>0).and.(nTempMagn>0)) Then
        Call mma_deallocate(XLM)
        Call mma_deallocate(ZLM)
        Call mma_deallocate(XRM)
        Call mma_deallocate(ZRM)
      End If

      If(nT>0) Then
        Call mma_deallocate(Texp)
        Call mma_deallocate(chit_exp)
      End If

      If((nT+nTempMagn)>0) Then
        Call mma_deallocate(T)
        Call mma_deallocate(XTexp)
        Call mma_deallocate(XT_no_field)
      End If
!---------------------------------------------------------------------
      Write(6,*)
      Write(6,'(10A)') (('-@-#-$-%-&-*-'),idim=1,10)
      Write(6,*)
      Write(6,'(10X,A)') 'HAPPY   LANDING !!!   POLY_ANISO ENDED  OK !'
      Write(6,*)
      Write(6,'(10A)') (('-*-&-%-$-#-@-'),idim=1,10)
      Call xFlush(6)
      Call qExit('PA_1')

      Return
      End
