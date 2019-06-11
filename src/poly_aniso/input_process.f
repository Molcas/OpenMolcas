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
      Subroutine input_process( nneq, neq, neqv, nmax, nCenter, nexch,
     &                          nDir, nDirZee, nDirTot, nH, nT, exch,
     &                          nBlock, nTempMagn, iopt, nMult, nDim,
     &                          nPair, i_pair, nP, AngPoints, lant,
     &                          multLn, KEOPT, encut_definition, nK,
     &                          mG, ncut, nsfs, nss, nLoc,
     &                          MxRank1, MxRank2, imaxrank, iPrint,

     &                          R_LG, gtens_input, dirX, dirY, dirZ,
     &                          dir_weight, zJ, cryst, coord, hmin,
     &                          hmax, TempMagn, thrs, Hexp, Mexp,
     &                          tmin, tmax, chit_exp, Texp, Xfield,
     &                          Jex, JAex, JAex9, JDMex, tpar, upar,
     &                          MagnCoords, encut_rate, eso,
     &                          JITOexR, JITOexI,

     &                          Title, itype, namefile_aniso,

     &                          Do_structure_abc,
     &                          compute_barrier, fitCHI, fitM, hinput,
     &                          tinput, compute_magnetization,
     &                          compute_Mdir_vector, zeeman_energy,
     &                          m_paranoid, m_accurate, smagn,
     &                          compute_g_tensors,compute_torque,
     &                          compute_susceptibility, Lines,
     &                          AnisoLines3, AnisoLines9, Dipol,
     &                          check_title, KE, DM_exchange,
     &                          JITO_exchange )

      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"

      Integer, intent(in)       :: iPrint
!     number of non-equivalent sites
      Integer, intent(in)       :: nneq
!     number of equivalent sites of each type, neqv = MAXVAL(neq(:))
      Integer, intent(in)       :: neq(nneq), neqv
!     number of equivalent sites of each type, neqv = MAXVAL(neq(:))
      Integer, intent(in)       :: nexch(nneq), nmax
      Character(1), intent(in)  :: itype(nneq)
      Integer, intent(in)       :: nCenter
      Integer, intent(in)       :: nLoc
      Real(kind=wp), intent(in) :: R_LG( nneq,neqv,3,3)
      Real(kind=wp), intent(in) :: gtens_input(3,nneq)
      Real(kind=wp), intent(in) :: eso(nneq,nLoc)
      Integer, intent(in)       :: nss(nneq), nsfs(nneq)
      Character(180), intent(in) :: namefile_aniso(nneq)
c  definition of the exchange:
!     total number of exchange states
      Integer, intent(in)       :: exch
!     number of metal pairs (number of interactions)
      Integer, intent(inout)    :: nPair
!     index of the metal site in a given interacting pair
      Integer, intent(inout)    :: i_pair(nPair,2)
!     index of the ITO ranks for each pair
      Integer, intent(in)       :: imaxrank(nPair,2)
      Integer, intent(in)       :: MxRank1, MxRank2
      Logical, intent(in)       :: Lines, AnisoLines3, AnisoLines9
      Logical, intent(in)       :: Dipol, DM_exchange, JITO_exchange
!     Lines exchange    ( 1 parameter / interacting pair)
      Real(kind=wp), intent(in) :: Jex(nPair)
!     Anisotropic Lines ( 3 parameter / interacting pair)
      Real(kind=wp), intent(in) :: JAex(nPair,3)
!     Anisotropic Lines full ( 9 parameters / interacting pair)
      Real(kind=wp), intent(in) :: JAex9(nPair,3,3)
!     Dzyaloshinsky-Morya exchange
      Real(kind=wp), intent(in) :: JDMex(nPair,3)
      Real(kind=wp), intent(in) ::
     &                          JITOexR(nPair,MxRank1,-MxRank1:MxRank1,
     &                                        MxRank2,-MxRank2:MxRank2),
     &                          JITOexI(nPair,MxRank1,-MxRank1:MxRank1,
     &                                        MxRank2,-MxRank2:MxRank2)

      ! options used in connection with KE
      Integer, intent(in)       :: lant, KEOPT, multLn
      Logical, intent(in)       :: KE
      Real(kind=wp), intent(in) :: tpar, upar
      ! options used in connection with Dipol-Dipol interaction
      Real(kind=wp), intent(in) :: MagnCoords(nneq,3)

c  definition of data for susceptibility
      Integer, intent(inout)    :: nT
      Logical, intent(in)       :: tinput
      Logical, intent(inout)    :: compute_susceptibility
      Real(kind=wp), intent(inout) :: tmin, tmax
      Real(kind=wp), intent(in) :: chit_exp(nT), Texp(nT)
      ! options related to XT_MoverH
      Real(kind=wp), intent(in) :: Xfield
      Integer, intent(in)       :: nH
      Integer, intent(in)       :: nTempMagn
      Integer, intent(in)       :: iopt
      Real(kind=wp), intent(inout) :: TempMagn(nTempMagn)
      Real(kind=wp), intent(in) :: Hexp(nH), Mexp(nH,nTempMagn)
      Real(kind=wp), intent(in) :: thrs
      Real(kind=wp), intent(in) :: hmin, hmax
      Logical, intent(in)       :: hinput
      Logical, intent(in)       :: compute_magnetization
      Logical, intent(in)       :: compute_Mdir_vector
      Logical, intent(in)       :: zeeman_energy
      Logical, intent(inout)    :: m_paranoid
      Logical, intent(in)       :: m_accurate
      Logical, intent(in)       :: smagn
      ! options used to set up nM and EM
      Integer, intent(in)       :: encut_definition
      Integer, intent(in)       :: nK, mG ! encut_definition=1;
      Integer, intent(in)       :: ncut   ! encut_definition=2;
      Real(kind=wp), intent(in) :: encut_rate ! encut_definition=3;

c  definition of g and D tensors
      Integer, intent(in)       :: nMult
      Integer, intent(in)       :: nDim(nMult)
      Logical, intent(in)       :: compute_g_tensors
c  magnetization torque
      Integer, intent(in)       :: nP
      Integer, intent(in)       :: AngPoints
      Logical, intent(in)       :: compute_torque
c  Zeeman energy and M vector
      Integer, intent(in)       :: nDir, nDirZee
      Real(kind=wp), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir)
      Real(kind=wp), intent(in) :: dir_weight(nDirZee,3)
c  definition of mean field parameter
      Real(kind=wp), intent(in) :: zJ
c  definintion of the crystal axes:
      Logical, intent(in)       :: Do_structure_abc
!     a, b, c, alpha, beta, gamma
      Real(kind=wp), intent(in) :: cryst(6)
!     Cartesian coordinates of the main metal site, or center
      Real(kind=wp), intent(in) :: coord(3)
c  definitions for blocking barrier
      Integer, intent(inout)    :: nBlock
      Logical, intent(inout)    :: compute_barrier
c  options for automatic fitting of parameters:
      Logical, intent(in)       :: fitCHI !-- not used so far
      Logical, intent(in)       :: fitM !-- not used so far

      Logical, intent(in)       :: check_title
      Character(180),intent(in) :: Title


! local variables
      Character(180) fmtline
      Logical  :: nosym
      Logical  :: ab_initio_all
      Logical  :: DBG
      Integer  :: nDirTot
      Integer  :: i,j,k,l,m,n,iT,iH,jEnd,irank1,irank2,iproj1,iproj2
      Integer  :: icount_B_sites


      Call qEnter('PA_input_process')
      DBG=.false.
c-----------------------------------------------------------------------
c print the data from this Subroutine:
      If(check_title) Then
        Write(6,'(A,A72)')  'TITL :         = ', Title
      End If
! ======================================================================
c     INFORMATION about individual magnetic sites
! ======================================================================
      Write(6,'(A,  I5)')   'PRLV :         = ', iPrint
      Write(6,'(A)') 'Information about individual magnetic sites:'
      Write(6,'(A)')
      Write(6,'(A,  I5)')   'NNEQ :         = ', nneq
      Write(6,'(A,  I5)')   'NEQV :         = ', neqv
      Write(6,'(A,50I5)')   'NEQ() :        = ', (neq(i),i=1,nneq)
      Write(6,'(A,  I5)')   'nCenter :      = ', nCenter

      Write(6,'(A,50I5  )') 'NEXCH():       = ', (nexch(i),i=1,nneq)
      Write(6,'(A,50I5  )') 'NMAX:          = ', nmax
      Write(6,'(A,50I5  )') 'NSTATE():      = ', (nsfs(i),i=1,nneq)
      Write(6,'(A,50I5  )') 'NSS():         = ', (nss(i),i=1,nneq)

      Do i=1, nneq
         If(itype(i).eq.'A') Then
            Write(6,'(i3,1x,A15,1x,3i5,A)') i, namefile_aniso(i),
     &                                   nsfs(i), nss(i), nexch(i),
     &                                   ' -- computed ab initio'
            Write(6,'(A,i2,A )') 'ESO(',i,' ):'
            Write(6,'(10F10.3)') (ESO(i,j),j=1,nss(i))
         Else
            Write(6,'(i3,A15,1x,3i5,A)') i, ' - - - - - - - ',
     &                                   nsfs(i), nss(i), nexch(i),
     &                                   ' -- generated as isotropic'
         End If
      End Do

      nosym=.true.
      icount_B_sites=0
      ab_initio_all=.true.
      Do i=1,nneq
        If (neq(i)>1) nosym=.false.
        If (itype(i).eq.'B') Then
           icount_B_sites=icount_B_sites+1
           ab_initio_all=.false.
        End If
      End Do
      If(DBG) Write(6,'(A,i3)') 'input_process:  icount_B_sites=',
     &                                           icount_B_sites

      If ( nosym.eqv..false. )  Then
         Write(6,'(A,A)') 'SYMM :         = ',' rotation matrices '//
     &                     'for equivalent magnetic centers:'
         Do i=1,nneq
            Write(6,'(8x,A,i3)') 'center type:',i
            Do j=1,neq(i)
               Do m=1,3
                  Write(6,'(20x,3F11.6)') (R_lg(i,j,m,n),n=1,3)
               End Do
            End Do
         End Do
      Else
         Write(6,'(A,A )') 'SYMM :         = ',' rotation matrices '//
     &                     'are "Identity Matrices" for all '//
     &                     'magnetic centers'
      End If

      If(ab_initio_all) Then
         Write(6,'(17x,A,I2,A)') 'all centers have been computed ' //
     &                           'ab initio and ',nneq,
     &                           ' file(s) of "aniso_i.input" '//
     &                           'type exist.'
      Else
         Write(6,'(17x,I2,A,I2,A)') nneq-icount_B_sites,' center(s) ' //
     &                             'have been computed ab initio and',
     &                              nneq-icount_B_sites ,
     &                             ' file(s) of "aniso_i.input" '//
     &                             'type exist.'
         Write(6,'(17x,I2,A,15F7.3)') icount_B_sites,' center(s) ' //
     &                             'have been defined empirical: '//
     &                             'no ZFS and g_tens = '
         Write(6,'(17x,A,99F10.5)') 'gX = ',(gtens_input(1,i),i=1,nneq)
         Write(6,'(17x,A,99F10.5)') 'gY = ',(gtens_input(2,i),i=1,nneq)
         Write(6,'(17x,A,99F10.5)') 'gZ = ',(gtens_input(3,i),i=1,nneq)
      End If
! ======================================================================
c     INFORMATION about exchange
! ======================================================================
      Write(6,'(A)') 'Information about exchange and magnetic couplings'
      If(Lines) Then
         Write(6,'(A,  I3  )') 'PAIR or LIN1 : = ', npair
         Do i=1,npair
            Write(6,'(17x,2I3,2x, F9.5)') i_pair(i,1),i_pair(i,2),
     &                                    Jex(i)
         End Do
      End If

      If(AnisoLines3) Then
         Write(6,'(A,  I3  )') 'LIN3 :         = ', npair
         Do i=1,npair
            Write(6,'(17x,2I3,2x,3F9.5)') i_pair(i,1),i_pair(i,2),
     &                                    (JAex(i,j),j=1,3)
         End Do
      End If

      If(AnisoLines9) Then
         Write(6,'(A,  I3  )') 'LIN9 :         = ', npair
         Do i=1,npair
            Write(6,'(17x,2I3,2x,9F9.5)') i_pair(i,1),i_pair(i,2),
     &                                    (JAex9(i,1,j),j=1,3),
     &                                    (JAex9(i,2,j),j=1,3),
     &                                    (JAex9(i,3,j),j=1,3)
         End Do
      End If

      If(DM_exchange) Then
         Write(6,'(A,  I3  )') 'DMEX :         = ', npair
         Do i=1,npair
            Write(6,'(17x,2I3,2x,3F9.5)') i_pair(i,1),i_pair(i,2),
     &                                    (JDMex(i,j),j=1,3)
         End Do
      End If


      If(Dipol) Then
         Write(6,'(A,A )') 'DIPO :         = ',' dipolar coupling '//
     &                     'for the above exchange pairs is included'
         Write(6,'(A,A )') 'COOR :         = ',' symmetrized '//
     &                     'cartesian coordinates for each '//
     &                     'non-equivalent magnetic centers '
         Do i=1,nneq
            Write(6,'(17x,i3,3F11.6)') i,(MagnCoords(i,l),l=1,3)
         End Do
      Else
         Write(6,'(A,A )') 'DIPO :         = ',' dipolar coupling '//
     &                     'for the above exchange pairs is NOT '//
     &                     'included.'
      End If


      If(KE) Then
         Write(6,'(A,A )') 'LONG :         = ','kinetic exchange '//
     &                     'coupling is requested:'
         Write(6,'(A,2F14.7)') 't     = ', tpar
         Write(6,'(A,2F14.7)') 'U     = ', upar
         Write(6,'(A,i2    )') 'lant  = ', lant
         Write(6,'(A,i2    )') 'multLn= ', multLn
         Write(6,'(A,i2    )') 'keopt = ', KEOPT
      End If

      If(JITO_exchange) Then
         Write(6,'(A,A )') 'ITOJ :         = ',' Full anisotopic ITO'//
     &                     ' exchange for the above exchange pairs'//
     &                     ' is included (in cm-1)'
         Do i=1,npair
           Write(6,'(17x,2I3)') i_pair(i,1),i_pair(i,2)
           Do irank1=1,imaxrank(i,1),2
            Do iproj1=-irank1,irank1
             Do irank2=1,imaxrank(i,2),2
              Do iproj2=-irank2,irank2
               Write(6,'(23x,4I3,2x,2ES23.14)')
     &                        irank1, iproj1, irank2, iproj2,
     &             JITOexR(i, irank1, iproj1, irank2, iproj2),
     &             JITOexI(i, irank1, iproj1, irank2, iproj2)
              End Do
             End Do
            End Do
           End Do
         End Do
      End If
c-----------------------------------------------------------------------
c     print Exchange Hamiltonian
c     print Exchange Eigenstates
c     print decomposition of exchange
c    ...
! ======================================================================
!  Print out of the UBAR:
      If (compute_barrier) Then
         If (compute_g_tensors) Then
            Nblock=0
            Do i=1,nmult
               Nblock=Nblock+ndim(i)
            End Do
         Else
            Write(6,'(A)') 'UBAR keyword needs to know the sizes of'
            Write(6,'(A)') 'the manIfolds forming the barrier'
            Write(6,'(A)') 'Was the MLTP keyword defined?'
            Write(6,'(A)') 'The program will continue, but UBAR is '//
     &                     'deactivated'
            compute_barrier=.false.
         End If

         If (Nblock > exch) Then
            Write(6,'(A)') 'It seems that you have asked for '//
     &                     'properties of inexistent states.'
            Write(6,'(A)') 'Nblock > Nexch!'
            Write(6,'(A,I5)') 'Nblock = ',nBlock
            Write(6,'(A,I5)') 'Nexch  = ',exch
            Write(6,'(A   )') 'Reset   NBLOCK=EXCH and continue'
            nBlock=exch
         End If
      End If


! ======================================================================
!  Print out of the g and D tensors:

      If ( compute_g_tensors .AND. (nMult > 0)) Then
         Write(6,'(A)') 'Information about g and D tensors:'
         Write(6,'(A,I3)')   'MLTP :         = ', NMULT
         If( nMult <= 20 ) Then
            Write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,NMULT)
         Else
            Write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,20)
            Do j=21,NMULT,20
               jEnd=MIN(NMULT,J+19)
            Write(6,'(A,20I3)') '                 ',(NDIM(i),i=j,jEnd)
           End Do
         End If
      Else
         Write(6,'(A)') 'Computation of pseudospin Hamiltonians '//
     &                  'was not requested'
      End If

      If ( compute_g_tensors .AND. (nMult == 0)) Then
         Write(6,'(A,i3)') 'nMult =',nMult
         Write(6,'(A   )') 'Computation of pseudospin Hamiltonians '//
     &                     'is deactivated.'
      End If
! ======================================================================

      If( Do_structure_abc .AND.
     &    ((nMult>0).OR.compute_susceptibility) ) Then
         Write(6,'(2A)') 'ABCC :         = ','the main axes '//
     &                   'for the computed tensors are '//
     &                   'written also'
         Write(6,'(A)') 'in the crystallographic "abc" axes'
         Write(6,'(10x,A,F9.4)') 'a       = ', cryst(1)
         Write(6,'(10x,A,F9.4)') 'b       = ', cryst(2)
         Write(6,'(10x,A,F9.4)') 'c       = ', cryst(3)
         Write(6,'(10x,A,F9.4)') 'alpha   = ', cryst(4)
         Write(6,'(10x,A,F9.4)') 'beta    = ', cryst(5)
         Write(6,'(10x,A,F9.4)') 'gamma   = ', cryst(6)
         Write(6,'(10x,a,3F9.4)') 'coords: = ',(coord(i),i=1,3)
      End If


! ======================================================================
!  Print out of the SUSCEPTIBILITY
      compute_susceptibility=.true.
      If(tmin==0.0_wp) tmin=0.0_wp
      If(tmax==0.0_wp) tmax=0.0_wp
      If(compute_susceptibility) Then
       !-----------------------------------------!
         Write(6,'(A)') 'Magnetic susceptibility will be computed'//
     &                  'in the limit of zero applied magnetic field.'
       !-----------------------------------------!
         Write(6,'(A,A )') 'CHIT :         = ',' molar magnetic '//
     &                     'susceptibility is computed using'
       !-----------------------------------------!
         If (IOPT.eq.1) Then
            Write(6,'(18x,A)') 'newly-derived formula for total '//
     &                         'XT and M'
         Else If (IOPT.eq.2) Then
            Write(6,'(18x,A)') 'additive formula employed for '//
     &                        'total XT and M'
         Else If (IOPT.eq.3) Then
            Write(6,'(18x,A)') 'old formula employed for total XT and M'
         Else
            Write(6,'(18x,A)') 'IOPT parameter out of range.'
            Return
         End If
       !-----------------------------------------!
         If(TINPUT) Then
            Write(6,'(A)') 'TEXP :         = the experimental '//
     &                     'temperature interval is provided:'
            Do iT=1,nT
               Write(6,'(5X,i4,F11.5,F14.7)') iT, Texp(iT), chit_exp(iT)
            End Do
         Else
            Write(6,'(A, I4)')  'TINT :      nT = ', nT
            Write(6,'(A,F7.3)') '          Tmin = ', Tmin
            Write(6,'(A,F7.3)') '          Tmax = ', Tmax
         End If
       !-----------------------------------------!
         If (xField.ne.0.0_wp) Then
            Write(6,'(A)') 'Magnetic susceptibility will be computed '//
     &                     'also by using:'
            Write(6,'(A)') ' X= dM/dH formula and'
            Write(6,'(A)') ' X=  M/H formula (commonly used by '//
     &                     'experimentalists)'
            Write(6,'(A,F9.5,A)') 'FIELD :          ',Xfield,' Telsa'
         End If
       !-----------------------------------------!
      Else
        Write(6,'(A)') 'Computation of magnetic susceptibility '//
     &                 'was not requested.'
      End If

! ======================================================================
!  Print out of the TORQUE
      If(compute_torque) Then
         Write(6,'(A)') 'TORQ :         = ',' magnetic torque '//
     &                  'is computed'
         Write(6,'(A,I4)') 'AngPoints = ', AngPoints
         Write(6,'(A,i4)') 'nP        = ', nP
      Else
         Write(6,'(A)') 'Computation of magnetic torque was not '//
     &                  'requested.'
      End If

! ======================================================================




! ======================================================================
!  Print out of the MAGNETIZATION
      If(compute_magnetization) Then
         Write(6,'(2A )') 'MAGN :         = ',' molar magnetization '//
     &                    'is computed'
       !-----------------------------------------!
        If(nH > 0) Then
           If(HINPUT) Then
              Write(6,'(A)') 'HEXP :         = the experimental '//
     &                       'field interval is provided:'
              Write(fmtline,'(A,i3,A)')
     &                       '(5x,i4,F11.5,',nTempMagn,'F14.7)'
              Do iH=1,nH
                 Write(6,fmtline)
     &               iH, Hexp(iH), (Mexp(iH,iT),iT=1,nTempMagn)
              End Do
           Else
              Write(6,'(A, I4)')  'HINT :      nH = ', nH
              Write(6,'(A,F7.3)') '          Hmin = ', Hmin
              Write(6,'(A,F7.3)') '          Hmax = ', Hmax
           End If
        End If
       !-----------------------------------------!
        If(encut_definition.eq.1) Then
           Write(6,'(A, I4)')  'NCUT :         = ', ncut
        Else If(encut_definition.eq.2) Then
           Write(6,'(A,2I4)')  'ECUT :         = ', nk,mg
        Else If(encut_definition.eq.3) Then
           Write(6,'(A,F7.3)') 'ERAT :         = ', encut_rate
        End If
       !-----------------------------------------!
        Write(6,'(A, I4   )') 'TMAG :         = ',nTempMagn
        ! check for negative T
        Do i=1, nTempMagn
          If(TempMagn(i) <= 0.0_wp) TempMagn(i)=2.0_wp
        End Do

       !-----------------------------------------!
        If (nTempMagn.le.10) Then
           Write(6,'(6x,A,10F7.3)') 'TempMagn = ',
     &                            (TempMagn(i),i=1,nTempMagn)
        Else
           Do i=1,nTempMagn,10
              j=MIN(nTempMagn,i+9)
              Write(6,'(17x,10F7.3)') (TempMagn(k),k=i,j)
           End Do
        End If
       !-----------------------------------------!
        If(zeeman_energy) Then
           Write(6,'(2A,I3,A)') 'ZEEM :         = ',' Zeeman '//
     &                          'splitting for the following ',nDirZee
           Write(6,'(17x,A  )') 'direction(s) of the magnetic '//
     &                          'field is given in the'//
     &                          ' "zeeman_energy_xxx.output" file(s)'
           Do i=1, nDirZee
              Write(6,'(17x,3F11.6)')  (dir_weight(i,l),l=1,3)
           End Do
        End If
       !-----------------------------------------!
        If(compute_Mdir_vector) Then
           Write(6,'(2A,I2,a)') 'MVEC :         = ','Magnetization '//
     &                          'vector(s) for the following ',nDir,
     &                          ' direction(s) of the applied magnetic'
           Write(6,'(17x,A)')   'field will be computed;'
           Do i=1,nDir
              Write(6,'(17x,A,I2,A,3F10.6)') 'Dir :',i,' : ',
     &                                       DirX(i), DirY(i), DirZ(i)
           End Do
        End If

        If(nH > 0) Then
       !-----------------------------------------!
           Write(6,'(A, I4)') 'MAVE :   nsymm = ', nsymm
           Write(6,'(A, I4)') '         ngrid = ', ngrid
           Write(6,'(A, I4)') '       nPoints = ', get_nP(nsymm,ngrid)
       !-----------------------------------------!
           nDirTot = 0
           nDirTot = nDir + nDirZee + get_nP(nsymm,ngrid)
           Write(6,'(A,I4,A)') 'Magnetization will be computed in',
     &                          nDirTot,' directions.'
           Write(6,'(A,I6,A)') 'There will be ',nDirTot*nH,
     &                         ' calculations in total.'

           !-----------------------------------------------------------|
           If (m_accurate) Then
              Write(6,'(2A   )') 'MACC :         = ','Contribution '//
     &                  'to magnetization from local excited states'
              Write(6,'(17x,A)') 'is computed exactly.'
           Else
              Write(6,'(17x,A)') 'is accounted by formula: '//
     &                           ' M_excited = X_excited * H'
           End If
           !-----------------------------------------------------------|

           !-----------------------------------------------------------|
           If((exch.lt.256) .AND. (zJ.ne.0.0_wp)) m_paranoid=.true.

           If ( m_paranoid .AND. (zJ.ne.0.0_wp)) Then
              Write(6,'(2A)') 'MPAR :         = ','Average spin '//
     &                        '< S > of the environment is computed'
              Write(6,'(2A)') '                 ','for each '//
     &                        'temperature point (more accurate)'
           Else
              Write(6,'(2A)') 'MPAR :         = ','Average spin '//
     &                        '< S > of the environment is computed'
              Write(6,'(2A)') '                 ','for the highest '//
     &                        'temperature point only and used for '//
     &                        'all computed temperature points.'
           End If
           !-----------------------------------------------------------|

           !-----------------------------------------------------------|
           If (smagn) Then
              Write(6,'(2A)') 'SMAG :         = ','Spin-only '//
     &                        'magnetisation is computed'
           Else
              Write(6,'(A)') 'Computation of spin magnetization '//
     &                       'was not requested.'
           End If
           !-----------------------------------------------------------|

           Write(6,'(A,ES14.6)') 'THRS = ', thrs
        End If
      Else
           Write(6,'(A)') 'Computation of molar magnetization was '//
     &                    'not requested.'
      End If
!  End Print out of the MAGNETIZATION
! ======================================================================
      If( compute_susceptibility .OR. compute_torque .OR.
     &    compute_magnetization ) Then
         If ( zJ .ne. 0.0_wp ) Then
            Write(6,'(A)') 'Mean field intermolecular interaction '//
     &                     'parameter'
            Write(6,'(A,F12.7)') 'ZJPR :         = ', zJ
         End If
      End If
! ======================================================================
      If( ( fitCHI .OR. fitM ) .AND. ( tinput .OR. hinput ) ) Then
            Write(6,'(A)') 'Automatic fitting of exchange parameters'//
     &                     'is yet in the development'
      End If
! ======================================================================
      Call qExit('PA_input_process')

      Return
      End subroutine input_process
