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
      Subroutine Readin_poly(nneq, neq, neqv, exch, nCenter,
     &                       nT, nH, nTempMagn, nDir, nDirZee,
     &                       nMult, nPair, nexch, nDim,
     &                       i_pair, lant, multLn, iPrint, keopt,
     &                       encut_definition, nK, mG, iopt, nP,
     &                       AngPoints, ncut, LUZee, MxRank1, MxRank2,
     &                       imaxrank,

     &                       TempMagn, R_LG, R_ROT, Jex, JAex, JAex9,
     &                       JDMex, JITOexR, JITOexI,
     &                       tpar, upar, cryst, coord, Xfield,
     &                       gtens_input, D_fact, EoverD_fact, riso,
     &                       MagnCoords, thrs, tmin, tmax, hmin, hmax,
     &                       Texp, chit_exp, Hexp, Mexp, encut_rate,
     &                       zJ, dirX, dirY, dirZ, dir_weight,

     &                       Title, itype,

     &                       ifHDF,
     &                       compute_g_tensors, compute_magnetization,
     &                       TINPUT,HINPUT, Do_structure_abc, DoPlot,
     &                       compute_Mdir_vector, zeeman_energy,
     &                       m_paranoid, m_accurate, smagn,
     &                       compute_susceptibility, decompose_exchange,
     &                       KE, fitCHI, fitM, compute_torque,
     &                       compute_barrier, Dipol, check_title,
     &                       AnisoLines1, AnisoLines3, AnisoLines9,
     &                       DM_exchange, JITO_exchange)
C
C  THIS ROUTINE READS THE standard input.
C
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"
#include "warnings.fh"

c  definition of the cluster:
      Integer       :: nneq, neqv, neq(nneq), nCenter
      Logical       :: ifHDF
c  definition of the local metal sites
      Real(kind=wp) :: R_LG( nneq,neqv,3,3)
      Real(kind=wp) :: R_ROT(nneq,neqv,3,3)
      Real(kind=wp) :: gtens_input(3,nneq)
      Real(kind=wp) :: D_fact(nneq)
      Real(kind=wp) :: EoverD_fact(nneq)
      Real(kind=wp), intent(out) :: riso(nneq,3,3)
      Character(1)  :: itype(nneq)
c  definition of the exchange:
!     total number of exchange states
      Integer       :: exch
!     number of metal pairs (number of interactions)
      Integer       :: nPair
!     exchange basis, nmax= MAX(nexch(:))
      Integer       :: nexch(nneq)
!     index of the metal site in a given interacting pair
      Integer       :: i_pair(nPair,2)
      Logical       :: AnisoLines1, AnisoLines3, AnisoLines9
      Logical       :: Dipol, DM_exchange
!     Lines exchange    ( 1 parameter / interacting pair)
      Real(kind=wp) :: Jex(nPair)
!     Anisotropic Lines ( 3 parameter / interacting pair)
      Real(kind=wp) :: JAex(nPair,3)
!     Anisotropic Lines full ( 9 parameters / interacting pair)
      Real(kind=wp) :: JAex9(nPair,3,3)
      Real(kind=wp) :: JDMex(nPair,3)
      ! options used in connection with ITO exchange:
      Logical             :: JITO_exchange
      Integer, intent(in) :: MxRank1, MxRank2
      Integer, intent(out):: imaxrank(npair,2)
      Real(kind=wp), intent(out) ::
     &                 JITOexR(nPair,MxRank1,-MxRank1:MxRank1,
     &                               MxRank2,-MxRank2:MxRank2)
      Real(kind=wp), intent(out) ::
     &                 JITOexI(nPair,MxRank1,-MxRank1:MxRank1,
     &                               MxRank2,-MxRank2:MxRank2)

      ! options used in connection with KE
      Integer       :: lant, KEOPT, multLn
      Logical       :: KE
      Real(kind=wp) :: tpar, upar
      ! options used in connection with Dipol-Dipol interaction
      Real(kind=wp) :: MagnCoords(nneq,3)

c  definition of g and D tensors
      Integer       :: nMult
      Integer       :: nDim(nMult)
      Logical       :: compute_g_tensors
c  definition of data for susceptibility
      Integer       :: nT
      Logical       :: tinput, compute_susceptibility
      Real(kind=wp) :: tmin, tmax
      Real(kind=wp) :: chit_exp(nT), Texp(nT)
      ! options related to XT_MoverH
      Real(kind=wp) :: Xfield
c  definition of data for magnetization:
      Integer       :: nH
      Integer       :: nTempMagn
      Integer       :: iopt
      Real(kind=wp) :: TempMagn(nTempMagn)
      Real(kind=wp) :: Hexp(nH), Mexp(nH,nTempMagn)
      Real(kind=wp) :: thrs
      Real(kind=wp) :: hmin, hmax
      Logical       :: hinput
      Logical       :: compute_magnetization
      Logical       :: compute_Mdir_vector
      Logical       :: zeeman_energy
      Logical       :: m_paranoid
      Logical       :: m_accurate
      Logical       :: smagn
      ! options used to set up nM and EM
      Integer       :: encut_definition
      Integer       :: nK, mG ! encut_definition=1;
      Integer       :: ncut   ! encut_definition=2;
      Real(kind=wp) :: encut_rate ! encut_definition=3;
c decompose exchange
      Logical       :: decompose_exchange

c  magnetization torque
      Integer       :: nP
      Integer       :: AngPoints
      Logical       :: compute_torque
c  Zeeman energy and M vector
      Integer       :: nDir, nDirZee
      Integer       :: LUZee(nDirZee)
      Real(kind=wp) :: dirX(nDir), dirY(nDir), dirZ(nDir)
      Real(kind=wp) :: dir_weight(nDirZee,3)
c  definition of mean field parameter
      Real(kind=wp) :: zJ
c  definintion of the crystal axes:
      Logical       :: Do_structure_abc
      Real(kind=wp) :: cryst(6) ! a, b, c, alpha, beta, gamma
!     Cartesian coordinates of the main metal site, or center
      Real(kind=wp) :: coord(3)
c  definitions for blocking barrier
      Logical       :: compute_barrier
c  options for automatic fitting of parameters:
      Logical       :: fitCHI !-- not used so far
      Logical       :: fitM !-- not used so far
c  definition of print level
      Integer       :: iPrint
      Logical       :: check_title
      Character(180):: Title
      Logical       :: DoPlot


c--------- LOCAL VARIABLES --------------------
      Integer       :: Input
      Real(kind=wp) :: check_dir_weight(3*nDirZee),tmp,sum
      Integer       :: ll,i,j,l,k,m,n,icount_b_sites,lp,ic,jc
      Integer       :: linenr, inneq, irank1, irank2, iproj1, iproj2
      Integer       :: duplicate_check(nPair), nind(exch,2)
      Integer       :: nst, ASUM, jrank1,jrank2,jproj1,jproj2
      Integer       :: i1,i2,lb1,lb2
!     index of the metal site in a given interacting pair
c      Integer       :: nind(nPair,2)
c      Integer       :: ind_exch(nneq)
c      Real(lind=wp) :: magncoords(2*maxanisofiles,3)
      Real(kind=wp) :: finddetr, detR
      Real(kind=wp) :: tmpR(3,3)
      Logical       :: nosym
      External      :: finddetr

c variables connected to computation of g and d tensors
      Logical       :: check_nneq_presence, readgfact, ab_initio_all
      Logical       :: tcheck, hcheck, encut_check
      Logical       :: check_symm_presence
      Logical       :: checktmag

      Real(kind=wp) :: t2, t1

      Character(2)  :: lanth
c      Character(14) :: namefile_energy(nDirZee)
      Character(21) :: namefile_energy
      Character(288):: Line, ctmp, string

      Integer :: IsFreeUnit
      External :: IsFreeUnit

      Logical :: DBG



      Call qEnter('PA_readinp')

      DBG=.false.

      check_nneq_presence   = .true.
      check_title           = .false.
      icount_B_sites        = 0
      i_pair                = 0
      nosym                 = .true.
      ENCUT_check           = .false.
      TINPUT                = .false.
      HINPUT                = .false.
      TCHECK                = .false.
      HCHECK                = .false.
      Do i=1,nneq
         Do j=1,Neq(i)
            R_rot(i,j,1,1)=1.0_wp
            R_rot(i,j,2,2)=1.0_wp
            R_rot(i,j,3,3)=1.0_wp
             R_lg(i,j,1,1)=1.0_wp
             R_lg(i,j,2,2)=1.0_wp
             R_lg(i,j,3,3)=1.0_wp
         End Do
      End Do
      check_symm_presence=.true.
      If (neqv>1) check_symm_presence=.false.


      If(DBG) Write(6,'(A,  i6)') 'RDIN:      nneq=',nneq
      If(DBG) Write(6,'(A,  i6)') 'RDIN:      neqv=',neqv
      If(DBG) Write(6,'(A,  i6)') 'RDIN:     nPair=',nPair
c      If(DBG) Write(6,'(A,99I6)') 'RDIN:  i_pair(:,1)=',
c     &                                    (i_pair(i,1),i=1,nPair)
c      If(DBG) Write(6,'(A,99I6)') 'RDIN:  i_pair(:,2)=',
c     &                                    (i_pair(i,2),i=1,nPair)
      If(DBG) Write(6,'(A,99i6)') 'RDIN:     neq()=',(neq(i),i=1,nneq)
      If(DBG) Write(6,'(A,99i6)') 'RDIN:   nexch()=',(nexch(i),i=1,nneq)


C=========== End of default settings====================================
      Input = 5
      Rewind(Input)
50    READ(Input,'(A72)',End=998) LINE
      If(DBG) write(6,'(A)') LINE
      Call NORMAL(LINE)
      If(LINE(1:5).ne.'&POLY') Go To 50
      LINENR=0

100   Call xFlush(6)
      READ(Input,'(A72)',End=998) LINE
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If(LINE(1:1).eq.'*')   Go To 100
      If(LINE     .eq.' ')   Go To 100
      If(LINE(1:3).eq.'END') Go To 210 !End


* ------------ TITL ---------------------------------------------------**
      If (LINE(1:4).eq.'TITL') Then

         READ(Input,*,ERR=997) ctmp

         If(DBG) write(6,'(A)') ctmp
         check_title=.true.
         Title = trim(ctmp)
         LINENR=LINENR+1
         Go To 100
      End If



*---  process MLTP command --------------------------------------------**
      If (LINE(1:4).eq.'MLTP') Then

         READ(Input,*,ERR=997) NMULT

         If(DBG) write(6,'(A,i4)') 'NMULT =',NMULT

         READ(Input,*,Err=997) (NDIM(i),i=1,NMULT)

         Do i=1,NMULT
            If( NDIM(i)<0 ) Then
               ctmp=''
               Write(ctmp,'(A,I2,A)') 'MLTP: the dimension of the '//
     &                                'multiplet ',i,' is negative!!!'
               Call WarningMessage(2,ctmp)
               Call quit(_RC_INPUT_ERROR_)
            Else If ( NDIM(i)==0 ) Then
               ctmp=''
               Write(ctmp,'(A,I2,A)') 'MLTP: the dimension of the '//
     &                                'multiplet ',i,' is zero!!!'
               Call WarningMessage(2,ctmp)
               Call quit(_RC_INPUT_ERROR_)
            End If
         End Do

         If(DBG) write(6,'(A,100i4)') 'NDIM: ',(NDIM(i),i=1,NMULT)
         compute_g_tensors=.true.
         LINENR=LINENR+2
         Go To 100
      End If



*---  process TINT command --------------------------------------------**
      If (LINE(1:4).eq.'TINT') Then
         compute_susceptibility=.true.
         If(TINPUT.EQV..FALSE.) Then
            TCHECK=.TRUE.
            t1=0.0_wp
            t2=0.0_wp

            READ(Input,*,ERR=997) t1, t2, nT

            If( (t1<0).OR.(t2<0)) Then
               Call WarningMessage(2,
     &                  'TINT: negative temperature requested! ')
               Call quit(_RC_INPUT_ERROR_)
            End If
            If( (t1-t2)>0.0_wp) Then
               Tmin=t2
               Tmax=t1
            Else If ( (t1-t2)<0.0_wp ) Then
               Tmin=t1
               Tmax=t2
            Else ! t1==t2
               Call WarningMessage(2,
     &                  'TINT: temperature interval == 0! ' )
               Call quit(_RC_INPUT_ERROR_)
            End If

            If(DBG) write(6,'(A,2E15.7,i6)') 'Tmin, Tmax, nT: ',
     &                                    Tmin, Tmax, nT
         Else
            goto 590
         End If
         LINENR=LINENR+1
         Go To 100
      End If



*---  process HINT command --------------------------------------------**
      If (LINE(1:4).eq.'HINT') Then
         If(HINPUT.EQV..FALSE.) Then
            HCHECK=.TRUE.
            compute_magnetization=.true.
            t1=0.0_wp
            t2=0.0_wp

            READ(Input,*,ERR=997) t1, t2, nH

            If ( (t1<0).OR.(t2<0) ) Then
               Call WarningMessage(2,
     &                  'HINT: negative field requested! ')
               Call quit(_RC_INPUT_ERROR_)
            End If
            If ( (t1-t2)>0.0_wp ) Then
               Hmin=t2
               Hmax=t1
            Else If ( (t1-t2)<0.0_wp ) Then
               Hmin=t1
               Hmax=t2
            Else ! t1==t2
               Call WarningMessage(2,
     &                  'HINT: temperature interval == 0! ')
               Call quit(_RC_INPUT_ERROR_)
            End If


            If(DBG) write(6,'(A,2E15.7,i6)') 'Hmin, Hmax, nH: ',
     &                                    Hmin, Hmax, nH
         Else
            Go To 591
         End If
         LINENR=LINENR+1
         Go To 100
      End If



*---  process THRS command --------------------------------------------**
      If (LINE(1:4).eq.'THRS') Then

         READ(Input,*,ERR=997) THRS

         If ( thrs<0.0_wp ) Then
            Call WarningMessage(2,'THRS: negative threshold for '//
     &                            ' average M!!! ')
            Write(6,'(A)') 'Set to default thrs=1d-10'
         End If

         If(DBG) write(6,'(A,E15.7)') 'THRS: ', THRS
         LINENR=LINENR+1
         Go To 100
      End If


*---  process XFIE command --------------------------------------------**
      If (LINE(1:4).eq.'XFIE') Then
         compute_susceptibility=.true.

         READ(Input,*,ERR=997) tmp

         If ( tmp<0.0_wp ) Then
            Call WarningMessage(2,'XFIE: negative value of the'//
     &                            ' applied field !')
            Write(6,'(A)') 'Set to positive !'
            Xfield=abs(tmp)
         Else If (tmp==0.0_wp) Then
            Call WarningMessage(2,'XFIE: zero value of the'//
     &                            ' applied field !')
            Write(6,'(A)') 'Field-applied XT will not be computed!'
         Else
            Xfield=tmp
         End If
         If(DBG) write(6,'(A,E15.7)') 'XFIE: ',xField

         LINENR=LINENR+1
         Go To 100
      End If


*---  process PLOT command --------------------------------------------**
      If (LINE(1:4).eq.'PLOT') Then
         DoPlot=.true.

         If(DBG) write(6,'(A,L2)') 'PLOT: ',DoPlot

         LINENR=LINENR+1
         Go To 100
      End If


*---  process IOPT command --------------------------------------------**
      If (LINE(1:4).eq.'IOPT') Then

         READ(Input,*,ERR=997)  I  !option for computing MSUM and XTSUM

         If ( (i<0).or.(i>3) ) Then
            Call WarningMessage(2,'IOPT: value out of range!!!')
            Write(6,'(A)') 'Set to default !'
            iopt=1
         Else
            iopt=i
         End If

         If(DBG) write(6,'(A,I6)') 'IOPT: ',IOPT
         LINENR=LINENR+1
         Go To 100
      End If


*---  process SMAG command --------------------------------------------**
      If (LINE(1:4).eq.'SMAG') Then
         smagn=.true.
         compute_magnetization=.true.
         If(DBG) write(6,'(A,L2)') 'SMAG: ',smagn
         LINENR=LINENR+1
         Go To 100
      End If


*---  process MACC command --------------------------------------------**
      If (LINE(1:4).eq.'MACC') Then
         compute_magnetization=.true.  ! request for computation of M(H)
         m_accurate=.true.             ! request for computation of M(H)
         If(DBG) write(6,'(A,L2)') 'MACC: ',m_accurate
         LINENR=LINENR+1
         Go To 100
      End If


*---  process FITX command --------------------------------------------**
      If (LINE(1:4).eq.'FITX') Then
         fitCHI=.true.
         LINENR=LINENR+1
         Go To 100
      End If


*---  process FITM command --------------------------------------------**
      If (LINE(1:4).eq.'FITM') Then
         fitM=.true.
         LINENR=LINENR+1
         Go To 100
      End If


*---  process MPAR command --------------------------------------------**
      If (LINE(1:4).eq.'MPAR') Then
         m_paranoid=.true.             ! request for computation of M(H)
         compute_magnetization=.true.
         If(DBG) write(6,'(A,L2)') 'MPAR: ',m_paranoid
         LINENR=LINENR+1
         Go To 100
      End If


*---  process TORQ command --------------------------------------------**
      If (LINE(1:4).eq.'TORQ') Then

         compute_torque=.true.         ! request for computation of M(H)
         READ(Input,*,ERR=997) i        ! number of angular points

         If ( i<=0 ) Then
            Call WarningMessage(2,'TORQ: nP value out of range!!!')
            Write(6,'(A)') 'Set to default !'
            nP=45
         Else
            nP=i
         End If

         If(DBG) write(6,'(A)') 'TORQ: nP=',nP
         AngPoints=nP+1
         LINENR=LINENR+1
         Go To 100
      End If


*---  process MAVE command --------------------------------------------**
      If (LINE(1:4).eq.'MAVE') Then
         compute_magnetization=.true.  ! request for computation of M(H)

         READ(Input,*,ERR=997) i,j  !nsymm, ngrid

         If ( (i<=0).or.(i>=4) ) Then
            Call WarningMessage(2,'MAVE: nSYMM value out of range!!!')
            Write(6,'(A)') 'Set to default !'
            nsymm=1
         Else
            nsymm=i
         End If

         If ( (j<=0).or.(j>=33) ) Then
            Call WarningMessage(2,'MAVE: nGRID value out of range!!!')
            Write(6,'(A)') 'Set to default !'
            ngrid=15
         Else
            ngrid=i
         End If

         If(DBG) write(6,'(2(A,i6))')' nsymm:  ', nsymm,
     &                               ' ngrid:  ', ngrid
         LINENR=LINENR+1
         Go To 100
      End If


*---  process NCUT command --------------------------------------------**
      If (LINE(1:4).eq.'NCUT') Then
         If(ENCUT_check) Then
            Go To 595
         Else
            ENCUT_check=.true.
!           request for computation of M(H)
            compute_magnetization=.true.
            encut_definition=1

            READ(Input,*,ERR=997) i ! NCUT
!           E_cut= exchange_energy(Ncut)

            If ( (i<=0).or.(i>exch) ) Then
               Call WarningMessage(2,'NCUT: value out of range!!!')
               Write(6,'(A)') 'Set to full exchange basis.'
               nCUT=exch
            Else
               nCUT=i
            End If
            If(DBG) write(6,'(A,2i6)') 'ncut:  ', nCut

            LINENR=LINENR+1
            Go To 100
         End If
      End If


*---  process ENCU command --------------------------------------------**
      If (LINE(1:4).eq.'ENCU') Then
         If(ENCUT_check) Then
            Go To 595
         Else
            ENCUT_check=.true.
!           request for computation of M(H)
            compute_magnetization=.true.
            encut_definition=2

            READ(Input,*,ERR=997) NK, MG
!           E_cut = NK * K_Boltz + MG * mu_Bohr
            If(DBG) write(6,'(A,2i6)') 'encu:  nK, mG=',NK, MG

            LINENR=LINENR+1
            Go To 100
         End If
      End If

*---  process ERAT command --------------------------------------------**
      If (LINE(1:4).eq.'ERAT') Then
         If(ENCUT_check) Then
            Go To 595
         Else
            ENCUT_check=.true.
!           request for computation of M(H)
            compute_magnetization=.true.
            encut_definition=3

            READ(Input,*,ERR=997) encut_rate
!           Ncut = INT(nexch*encut_rate); E_cut= E(Ncut)
            If(DBG) write(6,'(A,i6)') 'encut_rate=',encut_rate

            LINENR=LINENR+1
            Go To 100
         End If
      End If

*---  process ZJPR command --------------------------------------------*
      If (LINE(1:4).eq.'ZJPR') Then
         READ(Input,*,ERR=997) ZJ
         If(DBG) write(6,'(A,E18.10)') 'zJ    =',zJ
         LINENR=LINENR+1
         Go To 100
      End If

*---  process PRLV command --------------------------------------------*
      If (LINE(1:4).eq.'PRLV') Then
         READ(Input,*,ERR=997) iPrint
         If(DBG) write(6,'(A,i6)') 'iPrint=',iPrint
         LINENR=LINENR+1
         Go To 100
      End If

*---  process COOR command --------------------------------------------*
      If (LINE(1:4).eq.'COOR') Then
         Dipol=.true.
         If(DBG) write(6,'(A)') 'isite   MagnCoords:'
         Do i=1,nneq
            READ(Input,*,ERR=997) (MagnCoords(i,l),l=1,3)
            If(DBG) write(6,'(i3,5x,3E18.10)') i,(MagnCoords(i,l),l=1,3)
         End Do
         LINENR=LINENR+nneq
         Go To 100
      End If

*---  process MVEC command --------------------------------------------**
      If (LINE(1:4).eq.'MVEC') Then
         compute_magnetization=.true.  ! request for computation of M(H)
         compute_Mdir_vector=.true.

         READ(Input,*,ERR=997) nDir
         If(DBG) write(6,'(A,i3)') 'nDir = ',nDir

         Do i=1,nDir
            READ(Input,*,ERR=997) DirX(i), DirY(i), DirZ(i)
            If(DBG) write(6,'(i3,5x,3E18.10)') i,DirX(i),DirY(i),DirZ(i)
         End Do
c  some processing:
         Do i=1,nDir
            sum = 0.0_wp
            sum = DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
            If ( sum .eq. 0.0_wp ) Then
               Write(6,'(a,i3,a)') 'error: MVEC  vector ',i,
     &                             'has the modulus = 0.0_wp.'
               Write(6,'(a     )') 'the program will stop now.'
               Call quit(_RC_INPUT_ERROR_)
            End If
            If ( ABS(sum-1.0_wp).gt.0.5d-13 ) Then
               Write(6,'(a,i3,a)') 'the vector ',i,'was re-normalized.'
               tmp=dirX(i)/sqrt(sum)
               dirX(i)=tmp
               tmp=dirY(i)/sqrt(sum)
               dirY(i)=tmp
               tmp=dirZ(i)/sqrt(sum)
               dirZ(i)=tmp
            End If
         End Do
c
         LINENR=LINENR+NDIR+1
         Go To 100
      End If

*---  process TEXP command --------------------------------------------*
      If (LINE(1:4).eq.'TEXP') Then
         compute_susceptibility=.true.
         If(TCHECK.EQV..FALSE.) Then
            TINPUT=.TRUE.

            READ(Input,*,ERR=997) NT
            If(DBG) write(6,'(A,i3)') 'nT = ',nT

            Do i=1,NT
               texp(i)=0.0_wp
               chit_exp(i)=0.0_wp

               READ(Input,*,ERR=997) texp(i), chit_exp(i)

               ! check and clean negative values:
               If(    texp(i)<0.0_wp)     texp(i)=abs(    texp(i))
               If(chit_exp(i)<0.0_wp) chit_exp(i)=abs(chit_exp(i))
            End Do
            tmin=texp(1)
            tmax=texp(nT)
         Else
            Go To 590
         End If
         LINENR=LINENR+NT+1
         Go To 100
      End If

*---  process HEXP command --------------------------------------------*
      If (LINE(1:4).eq.'HEXP') Then
         compute_magnetization=.true.
         If(checkTMAG) Then
            Write(6,'(A)') 'The data provided in TMAG will be ignored.'
         End If

         If(HCHECK.EQV..FALSE.) Then
            HINPUT=.TRUE.

            READ(Input,*) nTempMagn, (TempMagn(i),i=1,nTempMagn)
            READ(Input,*) nH
            If(DBG) Write(6,*) 'HEXP: nTempMagn =',nTempMagn
            If(DBG) Write(6,*) 'HEXP: TempMagn()=',
     &                          (TempMagn(i),i=1,nTempMagn)
            If(DBG) Write(6,*) 'HEXP: nH        =',nH

            If ( nH <0 ) nH=abs(nH)
            If ( nH==0 ) Call Quit_OnUserError()

            Do i=1,nH
               hexp(i)=0.0_wp
               Do j=1,nTempMagn
                  Mexp(i,j)=0.0_wp
               End Do
            End Do

            Do i=1,nH
               READ(Input,*,ERR=997) Hexp(i),
     &                              (Mexp(i,j),j=1,nTempMagn)
               ! check and clean negative values:
               If(hexp(i)<0.0_wp) hexp(i)=abs(hexp(i))
               Do j=1,nTempMagn
                  If(Mexp(i,j)<0.0_wp) Mexp(i,j)=abs(Mexp(i,j))
               End Do
            End Do
            hmin=hexp(1)
            hmax=hexp(nH)
         Else
            Go To 591
         End If
         LINENR=LINENR+NH+2
         Go To 100
      End If
*---  process TMAG command --------------------------------------------*
      If (LINE(1:4).eq.'TMAG') Then
         If(HINPUT.EQV..FALSE.) Then
            compute_magnetization=.true.
            checkTMAG=.true.

            READ(Input,*,ERR=997) nTempMagn,(TempMagn(i),i=1,nTempMagn)
            If(DBG) Write(6,*) 'TMAG: nTempMagn =',nTempMagn
            If(DBG) Write(6,*) 'TMAG: TempMagn()=',
     &                               (TempMagn(i),i=1,nTempMagn)

            ! check and clean for zero / negative values:
            Do i=1, nTempMagn
               If ( TempMagn(i) < 0.0_wp) Then
                  Call WarningMessage(2,'TMAG: negative temperature '//
     &                                  'requested! ')
                  WRITE(6,'(A)') 'Set to positive.'
                       TempMagn(i)=abs(TempMagn(i))
               Else If ( TempMagn(i)==0.0_wp ) Then
                  Call WarningMessage(2,'TMAG: zero temperature '//
     &                                  'requested! ')
                  WRITE(6,'(A)') 'Set to 0.0001 K.'
                       TempMagn(i)=0.0001_wp
               End If
            End Do
         Else
            Write(6,'(A)') 'TMAG data is taken from HEXP.'
         End If
         LINENR=LINENR+1
         Go To 100
      End If
*---  process NNEQ command --------------------------------------------*
c this is the most important keyword for Poly_Aniso
      If (LINE(1:4).eq.'NNEQ') Then
!        number of non-equivalent centers; type of all centers
         READ(Input,*,ERR=997)  NNEQ, ab_initio_all, ifHDF
         If(DBG) Write(6,'(A,i4,A,L2,A,L2)') 'NNEQ=', NNEQ,
     &                          ' ab_initio_all=',ab_initio_all,
     &                          ' ifHDF=', ifHDF
!        number of equivalent centers of type "i"
         READ(Input,*,ERR=997) (NEQ(i)  ,i=1,Nneq)
         If(DBG) Write(6,'(A,100I4)') 'NEQ(I)=',(NEQ(i),i=1,nneq)
!        number of RASSI wf for exchange
         READ(Input,*,ERR=997) (Nexch(i),i=1,Nneq)
         If(DBG) Write(6,'(A,100I4)') 'NExch(I)=',(NExch(i),i=1,nneq)

         Do i=1,nneq
            If(Nexch(i).lt.0) Then
               Write(6,'(A,i1,a)') 'The number of functions taken '//
     &                             'in the exchange interaction '//
     &                             'from center of type ',i,
     &                             ' is negative!'
               Write(6,'(A     )') 'The program has to stop, '//
     &                             'since the input is not reasonable.'
               Call quit(_RC_INPUT_ERROR_)
            Else If(Nexch(i).eq.0) Then
               Write(6,'(A,i1,a)') 'The number of functions taken '//
     &                             'in the exchange interaction '//
     &                             'from center of type ',i,
     &                             ' is 0 (zero).'
               Write(6,'(A)') 'The program has to stop, since '//
     &                        'the input is not reasonable.'
               Call quit(_RC_INPUT_ERROR_)
            End If
         End Do
         check_nneq_presence=.true.
         Do i=1,nneq
            If (neq(i).gt.1) Then
               nosym = .false.
            End If
         End Do
c        Write(6,'(A,i5)') 'exch = ', exch
         If(exch.eq.1) Then
            Write(6,'(3/)')
            Write(6,'(100A)') ('#',i=1,100)
            Write(6,'(3/)')
            Write(6,'(A)') 'The size of the exchange matrix is 1. '//
     &                     'Is this really what you intended '//
     &                     'to compute?'
            Write(6,'(A)') 'The program will continue...'
            Write(6,'(3/)')
            Write(6,'(100A)') ('#',i=1,100)
            Write(6,'(3/)')
         End If

c         Do i=1,nCenter
c            ind_exch(i)=i
c         End Do

!        If the EXCH is above the limit => exit with an error
         If(exch > 15000) Then
            Write(6,'(A)') 'The number of exchange states is very large'
            Write(6,'(A)') 'EXCH=',exch
            Write(6,'(A)') 'The calculation will continue, but might '//
     &                     'take a LOT of time'
            Write(6,'(A)') 'We recomend to switch OFF the '//
     &                     'computation of powder magnetization'
         End If

         If(ab_initio_all .eqv. .false.) Then
            readgfact=.true.

!           type of the center:   A -- the information is read from aniso_ion.input
!                                 B -- the center is isotropic with g factor read from the input
!                                 C -- the center is anisotropic with
            READ(Input,*,ERR=997) (itype(i),i=1,Nneq)
            If(DBG) Write(6,'(A,100A3)') 'itype: ',(itype(i),i=1,nneq)
            If(DBG) Call xFlush(6)
            icount_B_sites=0
            Do i=1,nneq
               If ((itype(i).eq.'B').OR.(itype(i).eq.'C')) Then
                  icount_B_sites=icount_B_sites+1
                  READ(Input,*,ERR=997) (gtens_input(l,i),l=1,3),
     &                                   D_fact(i), EoverD_fact(i)
                  If(DBG) Write(6,'(A,i4,A,3E20.10, 2(A,E20.10) )')
     &                   'gtens_input(',i,')=',(gtens_input(l,i),l=1,3),
     &                   ' D = ', D_fact(i),' E/D =', EoverD_fact(i)
                  If( (itype(i).eq.'C').AND.(
     &               (gtens_input(1,i).ne.gtens_input(2,i)).OR.
     &               (gtens_input(1,i).ne.gtens_input(3,i)).OR.
     &               (gtens_input(2,i).ne.gtens_input(3,i))) ) Then
                      Do ic=1,3
                         READ(Input,*,ERR=997) (riso(i,jc,ic),jc=1,3)
                      End Do

                  Else
                      riso(i,1:3,1:3)=0.0_wp
                      Do ic=1,3
                          riso(i,ic,ic)=1.0_wp
                      End Do
                  End if
               End If

               If(abs(D_fact(i)) > 0.0_wp) Then
                   If( abs(EoverD_fact(i)/D_fact(i)) >
     &                                 0.3333333333333333333_wp ) Then
                   Write(string,'(A,i3,A)') 'NNEQ: |E/D| for center',i,
     &              ' > 1/3'
                   Call WarningMessage(2,trim(string))
                   WRITE(6,'(A,i3,A)') '|E/D| for center ',i,
     &                                 ' is set to 1/3'
                   End If
               End If
            End Do
         Else
            Do i=1,nneq
               itype(i)='A'
            End Do
         End If !ab_initio_all
         LINENR=LINENR+3+icount_B_sites
         If(DBG) Call xFlush(6)
         Go To 100
      End If
*---  process LIN9 command --------------------------------------------*
      If (LINE(1:4).eq.'LIN9') Then
         AnisoLines9=.true.

         READ(Input,*,ERR=997) npair
         If(DBG) Write(6,'(A,i6)')'LIN9:  nPair=',nPair

         Do i=1,npair
            Do j=1,3
               Do k=1,3
                  JAex9(i,j,k)=0.0_wp
               End Do
            End Do
            i_pair(i,1)=0
            i_pair(i,2)=0

            ! the convention for 9 exchange interactions is
            ! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
            READ(Input,*,ERR=997) i_pair(i,1),i_pair(i,2),
     &                              (JAex9(i,1,j),j=1,3),
     &                              (JAex9(i,2,j),j=1,3),
     &                              (JAex9(i,3,j),j=1,3)
           If(DBG) Write(6,'(A,2I3,9F10.6)') 'LIN9: ',
     &      i_pair(i,1),i_pair(i,2),
     &                                       (JAex9(i,1,j),j=1,3),
     &                                       (JAex9(i,2,j),j=1,3),
     &                                       (JAex9(i,3,j),j=1,3)
         End Do
         LINENR=LINENR+npair+1
         Go To 100
      End If



*---  process LIN3 command --------------------------------------------*
      If (LINE(1:4).eq.'ALIN' .OR. LINE(1:4).eq.'LIN3') Then
         AnisoLines3=.true.

         READ(Input,*,ERR=997) npair
         If(DBG) Write(6,'(A,i6)')'nPair=',nPair

         Do i=1,npair
            Do j=1,3
               JAex(i,j)=0.0_wp
            End Do
            i_pair(i,1)=0
            i_pair(i,2)=0

            ! the convention for 3 exchange interactions is
            ! Jxx, Jyy, Jzz
            READ(Input,*,ERR=997) i_pair(i,1),i_pair(i,2),
     &                            (JAex(i,j),j=1,3)
            If(DBG) Write(6,'(A,i6)') i_pair(i,1),i_pair(i,2),
     &                            (JAex(i,j),j=1,3)
         End Do
         LINENR=LINENR+npair+1
         Go To 100
      End If


*---  process JITO command --------------------------------------------*
      If (LINE(1:4).eq.'ITOJ') Then
         JITO_exchange=.true.

         READ(Input,*,ERR=997) npair
         If(DBG) Write(6,'(A,i6)') 'nPair=',nPair
         JITOexR=0.0_wp
         JITOexI=0.0_wp

         Do i=1,npair
            i_pair(i,1)=0
            i_pair(i,2)=0
            imaxrank(i,1)=0
            imaxrank(i,2)=0

            ! the convention for exchange interactions is
            !
            READ(Input,*,ERR=997) i_pair(i,1),i_pair(i,2),
     &                          imaxrank(i,1),imaxrank(i,2)
              Do irank1=1,imaxrank(i,1),2
               Do iproj1=-irank1,irank1
                Do irank2=1,imaxrank(i,2),2
                 Do iproj2=-irank2,irank2
                  READ(Input,*,ERR=997) jrank1, jproj1, jrank2, jproj2,
     &               JITOexR(i,jrank1,jproj1,jrank2,jproj2),
     &               JITOexI(i,jrank1,jproj1,jrank2,jproj2)
                 End Do
                End Do
               End Do
              End Do

              If(DBG) Then
               Write(6,'(A,I3)') 'ITO Exchange parameters for pair:',i
               Do jrank1=1,imaxrank(i,1),2
                Do jproj1=-jrank1,jrank1
                 Do jrank2=1,imaxrank(i,2),2
                  Do jproj2=-jrank2,jrank2
                   WRITE(6,'(4I3,2x,2ES21.14)')
     &                          jrank1,jproj1,jrank2,jproj2,
     &                JITOexR(i,jrank1,jproj1,jrank2,jproj2),
     &                JITOexI(i,jrank1,jproj1,jrank2,jproj2)
                  End Do
                 End Do
                End Do
               End Do
              End If
         End Do
         LINENR=LINENR+npair+1
         Go To 100
      End If


*---  process PAIR = LIN1 command --------------------------------------------*
      If (LINE(1:4).eq.'PAIR' .OR. LINE(1:4).eq.'LIN1') Then
        AnisoLines1=.true.

         READ(Input,*,ERR=997) npair
         If(DBG) Write(6,'(A,i6)')'nPair=',nPair

         Do i=1,npair
            Jex(i)=0.0_wp
            i_pair(i,1)=0
            i_pair(i,2)=0

            READ(Input,*,ERR=997) i_pair(i,1), i_pair(i,2), Jex(i)
            If(DBG) Write(6,'(i4,2x,2I4,2x,E18.10)')
     &                         i, i_pair(i,1), i_pair(i,2), Jex(i)

         End Do
         LINENR=LINENR+npair+1
         Go To 100
      End If



*---  process DMEX command --------------------------------------------*
      If (LINE(1:4).eq.'DMEX') Then
         DM_exchange=.true.
         If(DBG) write(6,'(A,L2)') 'DMEX::  DM_exchange=',DM_exchange

         READ(Input,*,ERR=997) npair
         If(DBG) Write(6,'(A,i6)')'nPair=',nPair

         Do i=1,npair
            Do j=1,3
               JDMex(i,j)=0.0_wp
            End Do
            i_pair(i,1)=0
            i_pair(i,2)=0

            ! the convention for 3 exchange interactions is
            ! JDMex(x, y, z)
            READ(Input,*,ERR=997) i_pair(i,1),i_pair(i,2),
     &                            (JDMex(i,j),j=1,3)
         End Do
         LINENR=LINENR+npair+1
         Go To 100
         LINENR=LINENR+1
         Go To 100
      End If



*---  process ZEEM command --------------------------------------------*
      If (LINE(1:4).eq.'ZEEM') Then
         zeeman_energy=.true.
         compute_magnetization=.true.
         LUZEE=0
         Do i=1,3*nDirZee
            check_dir_weight(i)=0.0_wp
         End Do

         READ(Input,*,ERR=997) nDirZee

        Do i=1,nDirZee
!         open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
!         be further written in mangetization() subroutine
          Write(namefile_energy,'(5A)') 'zeeman_energy_',
     &                     CHAR(48+mod( int((i)/100),10)),
     &                     CHAR(48+mod( int((i)/10 ),10)),
     &                     CHAR(48+mod( int( i     ),10)),'.txt'
          If(DBG) Write(6,'(2A)') 'namefile_energy: ', namefile_energy
          LUZee(i)=IsFreeUnit(30+i)
          Call molcas_open(LUZee(i),namefile_energy)

          READ(Input,*,ERR=997) (dir_weight(i,l),l=1,3)

          check_dir_weight(i)=0.0_wp
          check_dir_weight(i)=sqrt( dir_weight(i,1)**2 +
     &                              dir_weight(i,2)**2 +
     &                              dir_weight(i,3)**2 )

          If( (check_dir_weight(i).lt.0.995_wp).OR.
     &        (check_dir_weight(i).gt.1.005_wp) ) Then
            Write(6,'(A)') 'The directions for the magnetic field '//
     &                     'for the computation of the Zeeman '//
     &                     'splitting are wrong.'
            Write(6,'(A)') '( px^2 + py^2 + pz^2 ) must give 1.!'
            Write(6,'(A,I3,2x,A,F9.5)') 'In the present case for '//
     &                                  'direction Nr.', i,
     &                                  ' the dir_weight = px^2 + '//
     &                                  'py^2 + pz^2 = ',
     &                                   check_dir_weight(i)**2
            LINENR=LINENR+2+i
            Go To 997
          End If

        End Do
         LINENR=LINENR+nDirZee+1
         Go To 100
      End If



*---  process MAGN command --------------------------------------------*
c      If (LINE(1:4).eq.'MAGN') Then
c         compute_magnetization=.true.
c         Go To 100
c      End If



*---  process SYMM command --------------------------------------------*
      If (LINE(1:4).eq.'SYMM') Then
         nosym=.false.
         check_symm_presence=.true.
         R_lg=0.0_wp
         R_rot=0.0_wp
         If(DBG) Write(6,'(A,i6)')'SYMM - at init'
         ll=0
         Do i=1,nneq
         If(DBG) Write(6,'(A,i6)')'SYMM:  i=',i

            READ(Input,*,ERR=997) inneq
            If(DBG) Write(6,'(A,i6)')'inneq=',inneq

            Do j=1,Neq(i)
            If(DBG) Write(6,'(A,i6)')'SYMM:  j=',j
               Do m=1,3
                  ll=ll+1
                  READ(Input,*,ERR=997) (R_lg(i,j,m,n),n=1,3)
                  If(DBG) Write(6,'(3E20.12)') (R_lg(i,j,m,n),n=1,3)
               End Do
            End Do

            Do j=1,neq(i)
               detR=0.0_wp
               tmpR=0.0_wp
               Do m=1,3
                  Do n=1,3
                     tmpR(m,n)=R_lg(i,j,m,n)
                  End Do
               End Do

               detR=FindDetR( tmpR(1:3,1:3), 3)
               If(DBG) Write(6,'(A,3E20.12)')'SYMM:  detR=',detR

               If ((abs(detR).lt.0.999).OR.(abs(detR).gt.1.001)) Then
                  Write(6,'(A)') 'The rotation matrices must be '//
     &                           'UNITARY.'
                  Write(6,'(A)') 'and of RIGHT hand system'
                  Write(6,'(A,F11.6)') 'DET = ', detR
                  LINENR=LINENR+ll+i+1
                  Go To 997
               End If
               If(detR.lt.0.0_wp) Then
                  Do m=1,3
                     Do n=1,3
                        R_rot(i,j,m,n)=-R_lg(i,j,m,n)
                     End Do !n
                  End Do !m
               Else If(detR.gt.0.0_wp) Then
                  Do m=1,3
                     Do n=1,3
                        R_rot(i,j,m,n)= R_lg(i,j,m,n)
                     End Do !n
                  End Do !m
               End If
            End Do !neq(i)
         End Do !nneq
         LINENR=LINENR+nneq+ll
         Go To 100
      End If



*---  process EXCH command --------------------------------------------*
      If (LINE(1:4).eq.'EXCH') Then
        decompose_exchange=.true.
        Go To 100
      End If

*---  process EXCH command --------------------------------------------*
c      If (LINE(1:4).eq.'END') Then
c        Go To 100
c      End If



*---  process UBAR command --------------------------------------------*
      If (LINE(1:4).eq.'UBAR') Then
        compute_barrier=.true.
c        READ(Input,*,ERR=997) icase
cc       icase =1  --> magnetic field is applied along the main magnetic
cc                     axis of each Doublet (multiplet)
c        If     ( icase.eq.1 ) Then
c          continue
c        Else If ( icase.eq.2 ) Then
cc       icase =2  --> magnetic field is applied along the main magnetic
cc                     axis of the specIfied Doublet (number) number (NDim)
c          READ(Input,*,ERR=997) NmagMult
c        Else If ( icase.eq.3 ) Then
cc       icase =3  --> a new coordination system is defined by the user,
cc                     and the magnetic field is applied along gZ (third axis)
c          Do i=1,3
c          READ(Input,*,ERR=997) (uBar_Rot(i,j),j=1,3)
c          End Do
c          LINENR=LINENR+3
c        Else
c          Write(6,'(A)') 'Is the UBAR keyword used correctly?'
c          Write(6,'(A)') 'The ICASE parameter is not understood.'
c        End If
        LINENR=LINENR+1
        Go To 100
      End If



C-------------------------------------------
      If (LINE(1:4).eq.'ABCC') Then
      Do_structure_abc = .TRUE.
      READ(Input,*,ERR=997) (cryst(i),i=1,6)
      coord(1)=0.0_wp
      coord(2)=0.0_wp
      coord(3)=0.0_wp
c      READ(Input,*,ERR=997) (coord(i),i=1,3)
        LINENR=LINENR+2
        Go To 100
      End If
c array "cryst" collects the crystallographic data:
c  cryst(1)= a
c  cryst(2)= b
c  cryst(3)= c
c  cryst(4)= alpha
c  cryst(5)= beta
c  cryst(6)= gamma
c  coord(i) =the coordinates of the magnetic center in "abc" axes
c  logical variable 'Do_structure_abc' will make the program compute
c  the magnetic and anisotropy axes in the "abc" coordinate system
C-------------------------------------------




*---  process LONG command --------------------------------------------*
      If (LINE(1:4).eq.'LONG') Then
         KE=.true.

         READ(Input,*,ERR=997) lanth, tpar, upar, KEOPT

         If( (lanth.eq.'GD').OR.(lanth.eq.'gd').OR.
     &       (lanth.eq.'gD').OR.(lanth.eq.'Gd') ) Then
            lant=1
            multLn=8
         Else If( (lanth.eq.'TB').OR.(lanth.eq.'tb').OR.
     &       (lanth.eq.'tB').OR.(lanth.eq.'Tb') ) Then
            lant=2
            multLn=13
         Else If( (lanth.eq.'DY').OR.(lanth.eq.'dy').OR.
     &            (lanth.eq.'dY').OR.(lanth.eq.'Dy') ) Then
            lant=3
            multLn=16
         Else If( (lanth.eq.'HO').OR.(lanth.eq.'ho').OR.
     &            (lanth.eq.'hO').OR.(lanth.eq.'Ho') ) Then
            lant=4
            multLn=17
         Else If( (lanth.eq.'ER').OR.(lanth.eq.'er').OR.
     &            (lanth.eq.'eR').OR.(lanth.eq.'Er') ) Then
            lant=5
            multLn=16
         Else If( (lanth.eq.'TM').OR.(lanth.eq.'tm').OR.
     &            (lanth.eq.'tM').OR.(lanth.eq.'Tm') ) Then
            lant=6
            multLn=13
         Else If( (lanth.eq.'YB').OR.(lanth.eq.'yb').OR.
     &            (lanth.eq.'yB').OR.(lanth.eq.'Yb') ) Then
            lant=7
            multLn=8
         Else
            Write(6,'( A)') 'Error in getting the type of '//
     &                      'the lanthanide!'
            Write(6,'(2A)') 'The program has this info: '//
     &                      'lanth =',lanth
            Write(6,'(26x,A,i2)') 'multLn =',multLn
         End If
         LINENR=LINENR+2
         Go To 100
      End If


! end of reading input keywords
210   Continue
c      CLOSE(Input)
c---------------------------------------------------------------------------




C=====  Perform some PROCESSING of the data =======================================

      If(nosym .eqv. .false.) Then
         If(check_symm_presence .eqv. .false.) Then
         Write(6,'(A)') 'The SYMM keyword is mandatory if the system '
         Write(6,'(A)') 'contains more than one center of the same '//
     &                  'type!'
         Call quit(_RC_INPUT_ERROR_)
         End If !check_symm_presence
      End If !nosym





c preparing the info for computation of molar magnetization
c calculate the total number of directions for the average procedure




c--------  definintion g and D tensors ------------------------------
      If(compute_g_tensors) Then
         nst=0
         Do i=1, nMult
            nst=nst+ndim(i)
            If( nst > exch ) Then
               Write(6,'(A)') 'You have requested the computation of '//
     &                        'properties of more states'
               Write(6,'(A)') 'than the total number of exchange '//
     &                        'states.'
               Write(6,'(A)') 'The program doesn"t know how to compute'
               Write(6,'(A)') 'properties of inexistent states'
               Write(6,'(A)') 'New settings:'
               Write(6,'(A,i6)') 'NMULT = ',i-1

               nmult=i-1
               Go To 14
            End If
         End Do
14       continue
      End If
      If(DBG) Write(6,*) 'READIN:  after proc g and D'
      If(DBG) Call xFlush(6)







c--------  definintion of exchange ------------------------------
      If (npair > 0 ) Then
         ASUM=0
         Do lp=1,nPair
            ASUM = ASUM + i_pair(lp,1) + i_pair(lp,2)
         End Do
         If( Dipol .AND. (ASUM.eq.0) ) Then
            lp=0
            Do i=1,nCenter-1
               Do j=i+1,nCenter
                  If(i>=j) goto 17
                  lp=lp+1
                  i_pair(lp,1)=i
                  i_pair(lp,2)=j
                  If(DBG) Write(6,'(A,i3,A,2I3)') 'lp=',lp,
     &                            ' i_pair(lp,1:2)=',
     &                              i_pair(lp,1),i_pair(lp,2)
  17              Continue
               End Do
            End Do
         End If

         Duplicate_check(1:nPair)=0
         Do i=1,npair
            Duplicate_check(i)=1000*i_pair(i,1)+i_pair(i,2)
            If(DBG) Write(6,*) 'Duplicate_check: ',Duplicate_check(i)

            If ( i_pair(i,1) == i_pair(i,2) ) Then
               Write(6,'(A,i2,a)') 'The center ',i_pair(i,1),
     &                             ' interacts with itself.'
               Write(6,'(A)') 'This is not possible. '//
     &                        'The program has to stop.'
               Call quit(_RC_INPUT_ERROR_)

            Else If( i_pair(i,1) > i_pair(i,2) ) Then
               Write(6,'(A)') 'The convention of this program '//
     &                        'enforces the label of the first '//
     &                        'magnetic center of a given '
               Write(6,'(A)') 'interaction to be smaller than the '//
     &                        'label of the second magnetic center.'
               Write(6,'(A)') 'In order to avoid confusions, please '//
     &                        'respect this convention.'
               Write(6,'(A)') 'Avoid duplicate interactions!'
               Write(6,'(A)') 'The program will stop now.'
               Call quit(_RC_INPUT_ERROR_)
            End If
         End Do
         If(DBG) Call xFlush(6)

         ! check on the duplicate
         Do i=1,npair
            Do j=i+1,npair
               If (j==i) Go To 197
               If (Duplicate_check(i) == Duplicate_check(j) ) Then
                  Write(6,'(A)') 'Some interactions are declared '//
     &                           'twice in the input. Please declare '//
     &                           'all interactions only once!'
                  Write(6,'(A)') 'The program has to stop.'
                  Call quit(_RC_INPUT_ERROR_)
               End If
  197          Continue
            End Do
         End Do


         ! check on the indices of the exchange couplings:
         Do i=1, nPair
            If ( (i_pair(i,1) > nCenter) .OR.
     &           (i_pair(i,2) > nCenter) ) Then
               Write(6,'(A)') 'The numbering of the magnetic '//
     &                        'centers within NPAIR keyword is wrong.'
               Write(6,'(A,i2,a)')
     &                        'For the interaction Nr. = ',i,' you '//
     &                        'numerate individual centers with '//
     &                        'numbers larger than the total number '//
     &                        'of centers in this molecule.'
               Write(6,'(A)') 'NPAIR keyword is wrong.'
               Write(6,'(A,I4)') 'i_pair(i,1) =',i_pair(i,1)
               Write(6,'(A,I4)') 'i_pair(i,2) =',i_pair(i,2)
               Write(6,'(A)') 'The program has to stop.'
               Call quit(_RC_INPUT_ERROR_)
            End If
         End Do


!        check on the size of MxRank1 and MxRank2 wrt nexch(1) and nexch(2)
         If(JITO_exchange) Then
           l=0
           Do i=1,nneq
             Do j=1,neq(i)
               l=l+1
               nind(l,1)=i
               nind(l,2)=j
             End Do
           End Do

           Do lp=1,npair
             lb1=i_pair(lp,1)
             lb2=i_pair(lp,2)
             i1=nind(lb1,1) ! indices of non-equivalent sites
             i2=nind(lb2,1) ! indices of non-equivalent sites
             If(nExch(i1) < (imaxrank(lp,1)+1)) Then
               Write(6,'(100A)') ('#',i=1,100)
               Write(6,'(A)') 'interacting pair = ', lp
               Write(6,'(A)') 'type of site 1 = ', i1
               Write(6,'(A,i2,A,i2)') 'nExch(',i1,') = ',nExch(i1)
               Write(6,'(A,i2,A,i2)') 'Rank2(',lp,',1) = ',
     &                                              imaxrank(lp,1)
               Write(6,'(A)')'nExch < Rank+1 !!!'
               Write(6,'(A)') 'Rank of ITO operators for site 1 is '//
     &                 'larger than the number of defined exchange '//
     &                 'states for this site'
               Write(6,'(A)') 'The program will use only parameters '//
     &                 'which bring non-zero contribution to exchange.'
               Write(6,'(100A)') ('#',i=1,100)
             End If
             If(nExch(i2) < (imaxrank(lp,2)+1)) Then
               Write(6,'(100A)') ('#',i=1,100)
               Write(6,'(A)') 'interacting pair = ', lp
               Write(6,'(A)') 'type of site 2 = ', i2
               Write(6,'(A,i2,A,i2)') 'nExch(',i2,') = ',nExch(i2)
               Write(6,'(A,i2,A,i2)') 'Rank2(',lp,',2) = ',
     &                                              imaxrank(lp,2)
               Write(6,'(A)')'nExch < Rank+1 !!!'
               Write(6,'(A)') 'Rank of ITO operators for site 2 is '//
     &                 'larger than the number of defined exchange '//
     &                 'states for this site'
               Write(6,'(A)') 'The program will use only parameters '//
     &                 'which bring non-zero contribution to exchange.'
               Write(6,'(100A)') ('#',i=1,100)
             End If
           End Do
         End If ! JITO_exchange
       End If ! nPair


      If(DBG) Write(6,*) 'READIN:  before 200 '
      If(DBG) Call xFlush(6)
c--------  definintion of exchange ------------------------------









      Go To 200
C--------------------------------------------
      Write(6,*)' The following input line was not understood:'
      Write(6,'(A)') LINE
      Go To 999

997   continue
      Write(6,*)' READIN: Error reading "poly_aniso.input" '
      Write(6,*)' near line nr.',LINENR+1
      Go To 999
998   continue
      Write(6,*)' READIN: Unexpected End of input file.'
999   continue
      Call quit(_RC_INPUT_ERROR_)
590   continue
      Write(6,*) 'READIN: THE TINT command is incompatible with TEXP'
      Call quit(_RC_INPUT_ERROR_)
591   continue
      Write(6,*) 'READIN: THE HINT command is incompatible with HEXP'
      Call quit(_RC_INPUT_ERROR_)
595   continue
      Write(6,*) 'READIN: THE NCUT, ENCU, and ERAT are mutually'//
     & ' exclusive. You cannot use more than one keyword at the'//
     & ' same time.'
      Call quit(_RC_INPUT_ERROR_)
c ===============   NORMAL EndING  ===============================
200   continue

      If(IPRINT.gt.2) Then
      Write(6,'(5X,A)') 'NO ERORR WAS LOCATED WHILE READING INPUT'
      End If
      Call qExit('PA_readinp')

      Return
      End !Subroutine

