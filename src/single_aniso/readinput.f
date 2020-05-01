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
      Subroutine readin_single(iprint,nmult,ndim,ldim,ndimcf,ldimcf,
     & nlanth,axisoption,poly_file,Ifrestart,input_to_read, nk, mg,
     & zmagn,Do_structure_abc,cryst,coord,encut_definition,
     & compute_g_tensors,compute_CF,nDirTot,nss,nstate,
     & compute_magnetization,compute_torque,smagn,tinput,hinput,
     & compute_Mdir_vector, zeeman_energy, LUZee, doplot,
     & encut_rate,ncut,nTempMagn,TempMagn,m_paranoid,
     & compute_barrier,nBlock,AngPoints,input_file_name,
     & nT,nH,texp,chit_exp,zJ,hexp,magn_exp,hmin,hmax,
     & nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,
     & thrs,H_torq,T_torq)
C
C  THIS ROUTINE READS THE FILE "SINGLE_ANISO.INPUT".
C
C
      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "warnings.fh"
#include "mgrid.fh"

c----------------------------------------------------------------
c   magnetization vectors:
      Integer            :: nDir,nDirZee
      Real(kind=8)      :: dirX(nDir), dirY(nDir), dirZ(nDir)
      Real(kind=8)      :: dir_weight(nDirZee,3)
      Logical            :: compute_Mdir_vector, zeeman_energy
c      common/MVL/ compute_Mdir_vector
c      common/MZEL/ zeeman_energy
c----------------------------------------------------------------
      Integer :: nss, nstate
      Integer :: iprint,nt,nh,nk,mg,igsm,l,jEnd
      Integer :: nlanth,ndimcf,ldimcf,axisoption, i_OxStat
      Integer :: input_to_read,encut_definition,ncut,ntempmagn
      Integer :: ndirtot
      Integer :: nBlock
      Integer :: nmult,ndim(nMult), ldim
      Integer :: AngPoints
      Integer :: LUZee(nDirZee)

      Real(kind=8) :: tmin,tmax,hmin,hmax,t1,t2,zj,
     &                 tempmagn(nTempMagn), encut_rate
      Real(kind=8) :: texp(nT),chit_exp(nT)
      Real(kind=8) :: hexp(nH),magn_exp(nH,ntempmagn)
      Real(kind=8) :: zmagn(3,3),sum,tmp
      Real(kind=8) :: cryst(6),coord(3)
      Real(kind=8) :: column_check(3,3), row_check(3,3)
      Real(kind=8) :: check_dir_weight(nDirZee)
      Real(kind=8) :: zr(3,3),det_zmagn
      Real(kind=8) :: FindDetR
      Real(kind=8) :: Xfield
      Real(kind=8), intent(out) :: thrs
      Real(kind=8), intent(out) :: H_torq, T_torq

      Logical :: hcheck,tcheck,poly_file
      Logical :: Ifrestart,Do_structure_abc
      Logical :: compute_magnetization,encut_check,compute_cf
      Logical :: compute_g_tensors
      Logical :: checktmag
      Logical :: compute_barrier
      Logical :: compute_torque
      Logical :: smagn
      Logical :: tinput,hinput
      Logical, intent(out) :: m_paranoid
      Logical              :: doplot

      Character(2)  :: cME,clanth(37)
      Character(21) :: namefile_energy
      Character(180):: input_file_name,tmpline,err_msg

      External      :: FindDetR

      Integer  :: IsFreeUnit
      External :: IsFreeUnit
c=======================================================================
c      COMMON/CHISUBR/ TMIN,TMAX,T1,T2
c      COMMON/CHISUBL/ TINPUT
c      COMMON/MAGNSUBI/ NK,MG
c      COMMON/MAGNSUBR/ HMIN,HMAX
c      COMMON/MAGNSUBL/ HINPUT

      Integer        :: I,LINENR,j
      Character(280) :: LINE


      Logical :: DBG

      Call qEnter('SA_readin')
      DBG=.false.
C============ Some default settings=====================================
c  variables in "mgrid.fh"
      Do i=1,32
         Do j=1,3
           get_nP(j,i)=0
         End Do
      End Do
      nsymm  =1
      ngrid  =15
      get_nP(1, 1)=   5
      get_nP(1, 2)=   9
      get_nP(1, 3)=  17
      get_nP(1, 4)=  25
      get_nP(1, 5)=  29
      get_nP(1, 6)=  45
      get_nP(1, 7)=  49
      get_nP(1, 8)=  61
      get_nP(1, 9)=  77
      get_nP(1,10)=  93
      get_nP(1,11)= 105
      get_nP(1,12)= 125
      get_nP(1,13)= 141
      get_nP(1,14)= 161
      get_nP(1,15)= 185
      get_nP(1,16)= 229
      get_nP(1,17)= 309
      get_nP(1,18)= 401
      get_nP(1,19)= 505
      get_nP(1,20)= 621
      get_nP(1,21)= 749
      get_nP(1,22)= 889
      get_nP(1,23)=1041
      get_nP(1,24)=1205
      get_nP(1,25)=1381
      get_nP(1,26)=1569
      get_nP(1,27)=1769
      get_nP(1,28)=1981
      get_nP(1,29)=2205
      get_nP(1,30)=2441
      get_nP(1,31)=2689
      get_nP(1,32)=2949
      get_nP(2, 1)=   4
      get_nP(2, 2)=   6
      get_nP(2, 3)=  11
      get_nP(2, 4)=  16
      get_nP(2, 5)=  17
      get_nP(2, 6)=  27
      get_nP(2, 7)=  28
      get_nP(2, 8)=  34
      get_nP(2, 9)=  41
      get_nP(2,10)=  51
      get_nP(2,11)=  57
      get_nP(2,12)=  68
      get_nP(2,13)=  75
      get_nP(2,14)=  86
      get_nP(2,15)=  98
      get_nP(2,16)= 121
      get_nP(2,17)= 162
      get_nP(2,18)= 209
      get_nP(2,19)= 262
      get_nP(2,20)= 321
      get_nP(2,21)= 386
      get_nP(2,22)= 457
      get_nP(2,23)= 534
      get_nP(2,24)= 617
      get_nP(2,25)= 706
      get_nP(2,26)= 801
      get_nP(2,27)= 902
      get_nP(2,28)=1009
      get_nP(2,29)=1122
      get_nP(2,30)=1241
      get_nP(2,31)=1366
      get_nP(2,32)=1497
      get_nP(3, 1)=  3
      get_nP(3, 2)=  4
      get_nP(3, 3)=  7
      get_nP(3, 4)= 10
      get_nP(3, 5)= 10
      get_nP(3, 6)= 16
      get_nP(3, 7)= 16
      get_nP(3, 8)= 19
      get_nP(3, 9)= 22
      get_nP(3,10)= 28
      get_nP(3,11)= 31
      get_nP(3,12)= 37
      get_nP(3,13)= 40
      get_nP(3,14)= 46
      get_nP(3,15)= 52
      get_nP(3,16)= 64
      get_nP(3,17)= 85
      get_nP(3,18)=109
      get_nP(3,19)=136
      get_nP(3,20)=166
      get_nP(3,21)=199
      get_nP(3,22)=235
      get_nP(3,23)=274
      get_nP(3,24)=316
      get_nP(3,25)=361
      get_nP(3,26)=409
      get_nP(3,27)=460
      get_nP(3,28)=514
      get_nP(3,29)=571
      get_nP(3,30)=631
      get_nP(3,31)=694
      get_nP(3,32)=760
c variables in "mvect.fh"
      compute_Mdir_vector=.false.
      zeeman_energy=.false.
      Do i=1,nDir
        DirX(i)=0.0_wp
        DirY(i)=0.0_wp
        DirZ(i)=0.0_wp
      End Do
      Do i=1,nDirZee
        dir_weight(i,1)=0.0_wp
        dir_weight(i,2)=0.0_wp
        dir_weight(i,3)=0.0_wp
      End Do
C========== Initializations of arrays ==================================
      Do I=1,nTempMagn
        TempMagn(i)=0.0_wp
      End Do
C============ Initializations of constants =============================

      thrs                  = 1.0D-10
      ldim                  = 1
      ncut                  = 1
      TMIN                  = 0.0_wp
      TMAX                  = 300.0_wp
      XFIELD                = 0.0_wp
c      NT                    = 301
      HMIN                  =  0.0_wp
      HMAX                  = 10.0_wp
c      NH                    =  21
      NK                    = 200
      MG                    = 200
c      NDIR                  = 0
c      nDirZee               = 0
c      nTempMagn             = 1
      If(nTempMagn>0) TempMagn(1) = 2.0_wp
      T1                    = 5.0_wp
      T2                    = 6.0_wp
      ZJ                    = 0.0_wp
      IGSM                  = 1
      m_paranoid            =  .true.
      checkTMAG             =  .FALSE.
      compute_g_tensors     =  .FALSE.
      compute_magnetization =  .FALSE.
      compute_Mdir_vector   =  .FALSE.
      compute_CF            =  .FALSE.
      compute_barrier       =  .FALSE.
      compute_torque        =  .FALSE.
      smagn                 =  .FALSE.
      TINPUT                =  .FALSE.
      HINPUT                =  .FALSE.
      TCHECK                =  .FALSE.
      HCHECK                =  .FALSE.
      POLY_FILE             =  .FALSE.
      Do_structure_abc      =  .FALSE.
      zeeman_energy         =  .false.
      ENCUT_check           =  .false.
      doplot                =  .false.
      nlanth                =  0
      nDIMcf                =  0
      cME                   =  '  '
c -- lanthanides
      clanth( 1)            =  'CE'
      clanth( 2)            =  'PR'
      clanth( 3)            =  'ND'
      clanth( 4)            =  'PM'
      clanth( 5)            =  'SM'
      clanth( 6)            =  'EU'
      clanth( 7)            =  'GD'
      clanth( 8)            =  'TB'
      clanth( 9)            =  'DY'
      clanth(10)            =  'HO'
      clanth(11)            =  'ER'
      clanth(12)            =  'TM'
      clanth(13)            =  'YB'
      clanth(14)            =  'LU'
c -- actinides
      clanth(15)            =  'TH'
      clanth(16)            =  'PA'
      clanth(17)            =  'U'
      clanth(18)            =  'NP'
      clanth(19)            =  'PU'
      clanth(20)            =  'AM'
      clanth(21)            =  'CM'
      clanth(22)            =  'BK'
      clanth(23)            =  'CF'
      clanth(24)            =  'ES'
      clanth(25)            =  'FM'
      clanth(26)            =  'MD'
      clanth(27)            =  'NO'
      clanth(28)            =  'LR'
c -- transition metals
      clanth(29)            =  'SC'
      clanth(30)            =  'TI'
      clanth(31)            =  'V'
      clanth(32)            =  'CR'
      clanth(33)            =  'MN'
      clanth(34)            =  'FE'
      clanth(35)            =  'CO'
      clanth(36)            =  'NI'
      clanth(37)            =  'CU'

      Do i=1,6
      cryst(i)=0.0_wp
      End Do
      Do i=1,3
      coord(i)=0.0_wp
      End Do
      axisoption       =  1
      input_to_read    =  0
      encut_definition =  2
      encut_rate       =  1
      Do i=1,3
      coord(i)=0.0_wp
       Do j=1,3
       zmagn(i,j)=0.0_wp
       End Do
      End Do
      AngPoints = 46

C=========== End of default settings====================================
      REWIND (5)
50    READ(5,'(A280)',End=998) LINE
      IF(DBG) Write(6,'(A)') TRIM(LINE)
      Call NORMAL(LINE)
      If(LINE(1:7).ne.'&SINGLE') Go To 50
      LINENR=0
100   READ(5,'(A280)',End=998) LINE
      IF(DBG) Write(6,'(A)') TRIM(LINE)
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If(LINE(1:1).eq.'*') Go To 100
      If(LINE.eq.' ') Go To 100
      If ((LINE(1:4).eq.'End ').OR.(LINE(1:4).eq.'    '))  Go To 200

C ------------------------------------------
c      If (LINE(1:4).eq.'TYPE') Then
c        READ(5,*,ERR=997) ICALC
c         If     (icalc.eq.1) Then
c           compute_g_tensors     =  .true.
c         Else If (icalc.eq.2) Then
c           compute_chiT          =  .true.
c         Else If (icalc.eq.3) Then
c           compute_magnetization =  .true.
c         Else If (icalc.eq.4) Then
c           compute_g_tensors     =  .true.
c           compute_chiT          =  .true.
c         Else If (icalc.eq.5) Then
c           compute_g_tensors     =  .true.
c           compute_magnetization =  .true.
c         Else If (icalc.eq.6) Then
c           compute_chiT          =  .true.
c           compute_magnetization =  .true.
c         Else If (icalc.eq.7) Then
c           compute_g_tensors     =  .true.
c           compute_chiT          =  .true.
c           compute_magnetization =  .true.
c         Else
c         Write(6,'(A)') 'ICALC: the maximum value is 7. However, '//
c     &    'the calculation will continue by computing the magnetism.'
c           compute_g_tensors     =  .true.
c           compute_chiT          =  .true.
c           compute_magnetization =  .true.
c         End If
c        LINENR=LINENR+1
c        Go To 100
c      End If
C ------------------------------------------
      If (LINE(1:4).eq.'MLTP') Then
        READ(5,*,ERR=997) NMULT
        IF(DBG) Write(6,*) 'MLTP:  NMULT=',NMULT
        compute_g_tensors     =  .true.
        READ(5,*,Err=997) (NDIM(i),i=1,NMULT)
        IF(DBG) Write(6,*) 'MLTP: NDIM()=',(NDIM(i),i=1,NMULT)
        IGSM=NDIM(1)
        LINENR=LINENR+2
        Go To 100
      End If
C ------------------------------------------
      If (LINE(1:4).eq.'REST') Then
        Ifrestart=.true.
        READ(5,*,ERR=997) input_to_read
        IF(DBG) Write(6,*) 'REST: input_to_read=',input_to_read
        If ( (input_to_read==2) .OR. (input_to_read==3) .OR.
     &       (input_to_read==4)  ) Then
          BACKSPACE(5)
          READ(5,*) input_to_read, tmpline
          input_file_name=trim(tmpline)
        End If
        If (input_to_read == 1) Then
          Write(6,*)'RESTART: -- The SINGLE_ANISO will take all '//
     &              'ab initio information from the binary '//
     &              '$Project.aniso" file.'
        Else If (input_to_read == 2) Then
          Write(6,*)'RESTART: -- The SINGLE_ANISO will take all '//
     &              'ab initio information from the ASCII '//
     &              trim(input_file_name)//' file.'
        Else If ( input_to_read .eq. 3 ) Then
          Write(6,*)'RESTART: -- The SINGLE_ANISO will take all '//
     &              'ab initio information from the RASSI-HDF5 '//
     &              'binary file.'
        Else If (input_to_read == 4) Then
          Write(6,*)'RESTART: -- The SINGLE_ANISO will take all '//
     &              'ab initio information from the ASCII '//
     &              trim(input_file_name)//' file -- molcas-8.0 format.'
        Else
          Call WarningMessage(2,'SINGLE_ANISO:: RESTART  '//
     &                          'option is not known.')
          Call Quit_OnUserError()
        End If
       Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'TINT') Then
        If(TINPUT.EQV..FALSE.) Then
          TCHECK=.TRUE.

            t1=0.0_wp
            t2=0.0_wp

            READ(5,*,ERR=997) t1, t2, nT

            If( (t1<0).OR.(t2<0)) Then
               Call WarningMessage(2,
     &                  'TINT: negative temperature requested! ')
               Call Quit_OnUserError()
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
               Call Quit_OnUserError()
            End If

          IF(DBG) Write(6,*) 'TINT: Tmin, Tmax, nT=',Tmin, Tmax, nT
        Else
          goto 590
        End If
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'XFIE') Then
        READ(5,*,ERR=997) Xfield
        IF(DBG) Write(6,*) 'XFIE: Xfield=',Xfield
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'HINT') Then
        If(HINPUT.EQV..FALSE.) Then
           HCHECK=.TRUE.
           compute_magnetization=.true.

            t1=0.0_wp
            t2=0.0_wp

            READ(5,*,ERR=997) t1, t2, nH

            If ( (t1<0).OR.(t2<0) ) Then
               Call WarningMessage(2,
     &                  'HINT: negative field requested! ')
               Call Quit_OnUserError()
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
               Call Quit_OnUserError()
            End If

           IF(DBG) Write(6,*) 'HINT: Hmin, Hmax, nH=',Hmin, Hmax, nH
        Else
           Go To 591
        End If
         LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'NCUT') Then
        If(ENCUT_check) Then
           Go To 595
        Else
           ENCUT_check=.true.
           compute_magnetization=.true. !request for computation of M(H)
           encut_definition=1

           READ(5,*,ERR=997) NCUT  !E_cut=ESO(Ncut)

            If ( NCUT<0 ) Then
               Call WarningMessage(2,
     &                  'NCUT: negative NCUT requested! ')
               Call Quit_OnUserError()
            Else If ( NCUT==0 ) Then
               Call WarningMessage(2,
     &                  'NCUT: zero NCUT requested! ')
               Call Quit_OnUserError()
            End If

           IF(DBG) Write(6,*) 'NCUT: NCUT=',NCUT
           LINENR=LINENR+1
           Go To 100
        End If
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'ENCU') Then
        If(ENCUT_check) Then
           Go To 595
        Else
          ENCUT_check=.true.
          compute_magnetization=.true.
          encut_definition=2

          READ(5,*,ERR=997) NK, MG

          If ( (NK<=0).OR.(MG<=0) ) Then
             Call WarningMessage(2,
     &                'ENCU: zero or negative NK,MG requested! ')
             Call Quit_OnUserError()
          End If

          IF(DBG) Write(6,*) 'ENCU: NK, MG=',NK, MG
          LINENR=LINENR+1
          Go To 100
        End If
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'ERAT') Then
        If(ENCUT_check) Then
           Go To 595
        Else
          ENCUT_check=.true.
          compute_magnetization=.true.
          encut_definition=3  !Ncut = INT(nss*encut_rate); E_cut=E(Ncut)

          READ(5,*,ERR=997) encut_rate

          If ( encut_rate<=0.0_wp ) Then
             Call WarningMessage(2,
     &                'ERAT: zero or negative encut rate requested! ')
             Call Quit_OnUserError()
          End If

          IF(DBG) Write(6,*) 'ERAT: encut_rate=',encut_rate
          LINENR=LINENR+1
          Go To 100
        End If
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'MVEC') Then
        compute_magnetization=.true.   ! request for computation of M(H)
        compute_Mdir_vector=.true.
        READ(5,*,ERR=997) nDir
        IF(DBG) Write(6,*) 'MVEC: nDir=',nDir
        Do i=1,nDir
          READ(5,*,ERR=997) DirX(i), DirY(i), DirZ(i)
        IF(DBG) Write(6,*) 'MVEC: DirX,DirY,DirZ=',
     &                            DirX(i),DirY(i),DirZ(i)
        End Do
c  some processing:
        Do i=1,nDir
          sum=0.0_wp
          sum=DirX(i)*DirX(i)+DirY(i)*DirY(i)+DirZ(i)*DirZ(i)
          If ( sum .eq. 0.0_wp ) Then
             Write(err_msg,'(a,i3,a)') 'error: MVEC  vector ',i,
     &                          'has the modulus = 0.0_wp.'
             Call WarningMessage(2,err_msg)
             Call Quit_OnUserError()
          End If
          If ( sum .ne. 1.0_wp) Then
            Write(6,'(a,i3,a)') 'the vector',i,'was re-normalized.'
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
C-------------------------------------------
      If (LINE(1:4).eq.'MAVE') Then
        compute_magnetization=.true.

        READ(5,*,ERR=997) nsymm, ngrid

        IF(DBG) Write(6,*) 'MAVE: nsymm, ngrid=',nsymm, ngrid
        If((nsymm.lt.1).OR.(nsymm.gt.3)) Then
        Write(6,'(A)') '"nsymm" must take Integer values 1, 2 or 3.'
        Write(6,'(A,i5)') '"nsymm" = ',nsymm
        Call Quit_OnUserError()
        End If
        If((ngrid.lt.1).OR.(ngrid.gt.32)) Then
        Write(6,'(A)') '"ngrid" must take Integer values 1, 2, ... 32.'
        Write(6,'(A,i5)') '"ngrid" = ',ngrid
        Call Quit_OnUserError()
        End If
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
c      If (LINE(1:4).eq.'TLIN') Then
c        READ(5,*,ERR=997) T1, T2
c        LINENR=LINENR+1
c        Go To 100
c      End If
C-------------------------------------------
      If (LINE(1:4).eq.'SMAG') Then
        smagn=.true.
        IF(DBG) Write(6,*) 'SMAG: =',smagn
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'PLOT') Then
        doplot=.true.
        IF(DBG) Write(6,*) 'PLOT: =',doplot
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'TEXP') Then
        If(TCHECK.EQV..FALSE.) Then
          TINPUT=.TRUE.
          READ(5,*,ERR=997) NT
          IF(DBG) Write(6,*) 'TEXP: nT=',nT
          Do i=1,NT
            texp(i)=0.0_wp
            chit_exp(i)=0.0_wp
            READ(5,*,ERR=997) texp(i), chit_exp(i)
            IF(DBG) Write(6,*) 'TEXP: texp(i), chit_exp(i)=',
     &                                texp(i), chit_exp(i)
            ! check and clean negative values:
            if(    texp(i)<0.0_wp)     texp(i)=abs(    texp(i))
            if(chit_exp(i)<0.0_wp) chit_exp(i)=abs(chit_exp(i))
          End Do
          tmin=texp(1)
          tmax=texp(nT)
        Else
          Go To 590
        End If
        LINENR=LINENR+NT+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'HEXP') Then
        compute_magnetization=.true.
        If(checkTMAG) Then
          Write(6,'(A)') 'The data provided in TMAG will be ignored.'
        End If
        If(HCHECK.EQV..FALSE.) Then
          HINPUT=.TRUE.
          READ(5,*) nTempMagn, (TempMagn(i),i=1,nTempMagn)
          IF(DBG) Write(6,*) 'HEXP: nTempMagn =',nTempMagn
          IF(DBG) Write(6,*) 'HEXP: TempMagn()=',
     &                             (TempMagn(i),i=1,nTempMagn)
          READ(5,*) nH
          IF(DBG) Write(6,*) 'HEXP: nH =',nH
          if(nH<0) nH=abs(nH)
          if(nH==0) Call Quit_OnUserError()
          Do i=1,nH
            hexp(i)=0.0_wp
            Do j=1,nTempMagn
              magn_exp(i,j)=0.0_wp
            End Do
          End Do
          Do i=1,nH
            READ(5,*,ERR=997) Hexp(i), (magn_exp(i,j),j=1,nTempMagn)
            IF(DBG) Write(6,*) 'HEXP: Hexp(i),  magn_exp(i,j)=',
     &                           Hexp(i), (magn_exp(i,j),j=1,nTempMagn)
            ! check and clean negative values:
            If(hexp(i)<0.0_wp) hexp(i)=abs(hexp(i))
            Do j=1,nTempMagn
              If(magn_exp(i,j)<0.0_wp) magn_exp(i,j)=abs(magn_exp(i,j))
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
C-------------------------------------------
      If (LINE(1:4).eq.'ZJPR') Then
        READ(5,*,ERR=997) ZJ
        IF(DBG) Write(6,*) 'ZJPR: zJ =',zJ
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'TORQ') Then
        compute_torque=.true.
        READ(5,*,ERR=997) AngPoints, H_torq, T_torq
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'TMAG') Then
         If(HINPUT.EQV..FALSE.) Then
            compute_magnetization=.true.
            checkTMAG=.true.

            READ(5,*,ERR=997) nTempMagn, (TempMagn(i),i=1,nTempMagn)

            Do i=1,nTempMagn
               If ( TempMagn(i)<=0.0_wp ) Then
                  Call WarningMessage(2,
     &                'TMAG: zero or negative temperature requested! ')
                  If(TempMagn(i) <0.0_wp) TempMagn(i)=abs(TempMagn(i))
                  If(TempMagn(i)==0.0_wp) TempMagn(i)=0.0001_wp
               End If
            End Do

            IF(DBG) Write(6,*) 'TMAG: nTempMagn =',nTempMagn
            IF(DBG) Write(6,*) 'TMAG: TempMagn()=',
     &                             (TempMagn(i),i=1,nTempMagn)
          ! check and clean negative values:
         Else
            Write(6,'(A)') 'TMAG data is taken from HEXP.'
         End If
         LINENR=LINENR+1
         Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'PRLV') Then
        READ(5,*,ERR=997) IPRINT
        IF(DBG) Write(6,*) 'PRLV: IPRINT =',iPrint
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'POLY') Then
        IF(DBG) Write(6,*) 'POLY:'
        POLY_FILE = .TRUE.
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'CRYS') Then
        compute_CF = .TRUE.
        Read(5,*,ERR=997) cME
        IF(DBG) Write(6,*) 'CRYS: cME =',cME

      If( (cME.eq.'ce') .OR. (cME.eq.'Ce') .OR.
     &    (cME.eq.'cE') .OR. (cME.eq.'CE') ) Then
      nlanth=1
      nDIMcf=6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'pr') .OR. (cME.eq.'Pr') .OR.
     &         (cME.eq.'pR') .OR. (cME.eq.'PR') ) Then
      nlanth=2
      nDIMcf=9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'nd') .OR. (cME.eq.'Nd') .OR.
     &         (cME.eq.'nD') .OR. (cME.eq.'ND') ) Then
      nlanth=3
      nDIMcf=10  ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'pm') .OR. (cME.eq.'Pm') .OR.
     &         (cME.eq.'pM') .OR. (cME.eq.'PM') ) Then
      nlanth=4
      nDIMcf=9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'sm') .OR. (cME.eq.'Sm') .OR.
     &         (cME.eq.'sM') .OR. (cME.eq.'SM') ) Then
      nlanth=5
      nDIMcf=6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'eu') .OR. (cME.eq.'Eu') .OR.
     &         (cME.eq.'eU') .OR. (cME.eq.'EU') ) Then
      nlanth=6
      nDIMcf=1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'gd') .OR. (cME.eq.'Gd') .OR.
     &         (cME.eq.'gD') .OR. (cME.eq.'GD') ) Then
      nlanth=7
      nDIMcf=8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
      lDIMCF=1 ! (L=0)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'tb') .OR. (cME.eq.'Tb') .OR.
     &         (cME.eq.'tB') .OR. (cME.eq.'TB') ) Then
      nlanth=8
      nDIMcf=13  ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'dy') .OR. (cME.eq.'Dy') .OR.
     &         (cME.eq.'dY') .OR. (cME.eq.'DY') ) Then
      nlanth=9
      nDIMcf=16  ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'ho') .OR. (cME.eq.'Ho') .OR.
     &         (cME.eq.'hO') .OR. (cME.eq.'HO') ) Then
      nlanth=10
      nDIMcf=17  ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'er') .OR. (cME.eq.'Er') .OR.
     &         (cME.eq.'eR') .OR. (cME.eq.'ER') ) Then
      nlanth=11
      nDIMcf=16  ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'tm') .OR. (cME.eq.'Tm') .OR.
     &         (cME.eq.'tM') .OR. (cME.eq.'TM') ) Then
      nlanth=12
      nDIMcf=13  ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'yb') .OR. (cME.eq.'Yb') .OR.
     &         (cME.eq.'yB') .OR. (cME.eq.'YB') ) Then
      nlanth=13
      nDIMcf=8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'lu') .OR. (cME.eq.'Lu') .OR.
     &         (cME.eq.'lU') .OR. (cME.eq.'LU') ) Then
      nlanth=14
      nDIMcf=1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
      lDIMCF=1 ! (L=0)


      !- - - - - - - - - - - - - - - - - - - -
      ! ACTINIDES
      Else If( (cME.eq.'th') .OR. (cME.eq.'Th') .OR.
     &         (cME.eq.'tH') .OR. (cME.eq.'TH') ) Then
      nlanth=15
      nDIMcf=6  ! f1; multiplet J=L-S=3-1/2=5/2  =>  J = 2F_5/2
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'pa') .OR. (cME.eq.'Pa') .OR.
     &         (cME.eq.'pA') .OR. (cME.eq.'PA') ) Then
      nlanth=16
      nDIMcf=9  ! f2; multiplet J=L-S=5-1=4  => J = 3H_4
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'U' ) .OR. (cME.eq.'u' ) .OR.
     &         (cME.eq.'u ') .OR. (cME.eq.'U ') ) Then
      nlanth=17
      nDIMcf=10  ! f3; multiplet J=L-S=6-3/2=9/2  => J = 4I_9/2
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'np') .OR. (cME.eq.'Np') .OR.
     &         (cME.eq.'nP') .OR. (cME.eq.'NP') ) Then
      nlanth=18
      nDIMcf=9  ! f4; multiplet J=L-S=6-2=4  => J = 5I_4
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'pu') .OR. (cME.eq.'Pu') .OR.
     &         (cME.eq.'pU') .OR. (cME.eq.'PU') ) Then
      nlanth=19
      nDIMcf=6  ! f5; multiplet J=L-S=5-5/2=5/2  => J = 6H_5/2
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'am') .OR. (cME.eq.'Am') .OR.
     &         (cME.eq.'aM') .OR. (cME.eq.'AM') ) Then
      nlanth=20
      nDIMcf=1  ! f6; multiplet J=L-S=3-3=0  => J = 3F_0
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'cm') .OR. (cME.eq.'Cm') .OR.
     &         (cME.eq.'cM') .OR. (cME.eq.'CM') ) Then
      nlanth=21
      nDIMcf=8  ! f7; multiplet J=L+S=0+7/2=0  => J = 8S_7/2
      lDIMCF=1 ! (L=0)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'bk') .OR. (cME.eq.'Bk') .OR.
     &         (cME.eq.'bK') .OR. (cME.eq.'BK') ) Then
      nlanth=22
      nDIMcf=13  ! f8; multiplet J=L+S=3+3=0  => J = 7F_6
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'cf') .OR. (cME.eq.'Cf') .OR.
     &         (cME.eq.'cF') .OR. (cME.eq.'CF') ) Then
      nlanth=23
      nDIMcf=16  ! f9; multiplet J=L+S=5+5/2=15/2  => J = 6H_15/2
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'es') .OR. (cME.eq.'Es') .OR.
     &         (cME.eq.'eS') .OR. (cME.eq.'ES') ) Then
      nlanth=24
      nDIMcf=17  ! f10; multiplet J=L+S=6+2=8  => J = 5I_8
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'fm') .OR. (cME.eq.'Fm') .OR.
     &         (cME.eq.'fM') .OR. (cME.eq.'FM') ) Then
      nlanth=25
      nDIMcf=16  ! f11; multiplet J=L+S=6+3/2=15/2  => J = 4I_15/2
      lDIMCF=13 ! (L=6)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'md') .OR. (cME.eq.'Md') .OR.
     &         (cME.eq.'mD') .OR. (cME.eq.'MD') ) Then
      nlanth=26
      nDIMcf=13  ! f12; multiplet J=L+S=5+1=6  => J = 3H_6
      lDIMCF=11 ! (L=5)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'no') .OR. (cME.eq.'No') .OR.
     &         (cME.eq.'nO') .OR. (cME.eq.'NO') ) Then
      nlanth=27
      nDIMcf=8  ! f13; multiplet J=L+S=3+1/2=7/2  => J = 2F_7/2
      lDIMCF=7 ! (L=3)
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'lr') .OR. (cME.eq.'Lr') .OR.
     &         (cME.eq.'lR') .OR. (cME.eq.'LR') ) Then
      nlanth=28
      nDIMcf=1  ! f14; multiplet J=L+S=0+0=0  => J = 1S_0
      lDIMCF=1 ! (L=0)





!------------------------ transition metals --------------------------!

      Else If( (cME.eq.'Sc') .OR. (cME.eq.'Sc') .OR.
     &         (cME.eq.'sC') .OR. (cME.eq.'SC') ) Then

         nlanth=29
         ! Sc2+ -- d^1
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) Then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) Then
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=5 ! (L=2) d^1
           Else
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If

      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'Ti') .OR. (cME.eq.'Ti') .OR.
     &         (cME.eq.'tI') .OR. (cME.eq.'TI') ) Then

         nlanth=30
         ! Ti2+ -- d^2
         ! Ti3+ -- d^1
         ! Ti4+ -- d^0
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=7 ! (L=3) d^2  3F
           Else If (i_OxStat == 3) Then
             lDIMCF=5 ! (L=2) d^1  2D
           Else
             lDIMCF=1 ! (L=0) d^4
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If


      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'v' ) .OR. (cME.eq.'V' ) .OR.
     &         (cME.eq.'V ') .OR. (cME.eq.'v ') ) Then

         nlanth=31
         ! V2+ -- d^3
         ! V3+ -- d^2
         ! V4+ -- d^1
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=7 ! (L=3) d^3
           Else If (i_OxStat == 3) Then
             lDIMCF=7 ! (L=3) d^2
           Else If (i_OxStat == 4) Then
             lDIMCF=5 ! (L=2) d^1
           Else
             lDIMCF=1 ! (L=0) d^4
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If


      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'cr') .OR. (cME.eq.'cR') .OR.
     &         (cME.eq.'Cr') .OR. (cME.eq.'CR') ) Then

         nlanth=32
         Write(6,'(A)') 'Crystal field will not be computed'
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'mn') .OR. (cME.eq.'mN') .OR.
     &         (cME.eq.'Mn') .OR. (cME.eq.'MN') ) Then

         nlanth=33
         ! Mn3+ -- d^4
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<3) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 3) Then
             lDIMCF=5 ! (L=2) d^4
           Else
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'fe') .OR. (cME.eq.'fE') .OR.
     &         (cME.eq.'Fe') .OR. (cME.eq.'FE') ) Then

         nlanth=34
         ! Co2+ -- d^6 or d^4
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=5 ! (L=2)  d^6  or  d^4
           Else
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'co') .OR. (cME.eq.'cO') .OR.
     &         (cME.eq.'Co') .OR. (cME.eq.'CO') ) Then

         nlanth=35
         ! Co2+ -- d^7
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=7 ! (L=3) d^7
           Else
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'ni') .OR. (cME.eq.'nI') .OR.
     &         (cME.eq.'Ni') .OR. (cME.eq.'NI') ) Then

         nlanth=36
         ! Ni2+ -- d^8
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=7 ! (L=2) d^8
           Else
             lDIMCF=1 ! (L=0)
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If
      !- - - - - - - - - - - - - - - - - - - -
      Else If( (cME.eq.'cu') .OR. (cME.eq.'cU') .OR.
     &         (cME.eq.'Cu') .OR. (cME.eq.'CU') ) Then

         nlanth=37
         ! Cu2+ -- d^9
         READ(5,*,ERR=997) i_OxStat
         IF(DBG) Write(6,*) 'CRYS: i_OxStat =',i_OxStat,'nlanth=',nlanth

           If(i_OxStat<0) then
             Write(6,'(3A,i5)') 'Oxidation state of',cME,'is negative:',
     &                           i_OxStat
             Write(6,'(A)')  'It was re-set to positive.'
             i_OxStat=abs(i_OxStat)
           End If
           If (i_OxStat<2) then
             lDIMCF=1 ! (L=0)
             Write(6,'(3A,i5)') 'Oxidation state of ', cME,' is:',
     &                           i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           Else If (i_OxStat == 2) Then
             lDIMCF=5 ! (L=2) d^9
           Else
             lDIMCF=1 ! (L=0) d^4
             Write(6,'(A)') 'Oxidation state of ', cME,' is:', i_OxStat
             Write(6,'(A)') 'Crystal field will not be computed'
           End If
      !- - - - - - - - - - - - - - - - - - - -
      Else
         Write(6,'(A)') 'Label of the metal is not understood.'
         Write(6,'(A)') 'Crystal field will not be computed'
      End If



      If (IPRINT > 2) Then
      Write(6,'(5x,3A)') 'SINGLE_ANISO will calculate the parameters'//
     & ' of the crystal field for Ln = ',clanth(nlanth),','
      Write(6,'(5x,A,I2,a)') 'for the ground multiplet J.'//
     & ' Multiplicity of J = ', nDIMcf, ' and'
      Write(6,'(5x,A,I2)') 'for the ground LS term.'//
     & ' Multiplicity of L = ', lDIMcf
      End If
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'QUAX') Then
c        If ( check_CRYS ) Then
        READ(5,*,ERR=997) axisoption
        IF(DBG) Write(6,*) 'QUAX: axisoption =',axisoption
        LINENR=LINENR+1

        If( (axisoption.lt.1) .OR. (axisoption.gt.3) ) Then
           Call WarningMessage(2,'QUAX: axisoption out of range!'//
     &          ' Calculation will continue by employing the default'//
     &          ' option.')
        End If
        If( axisoption == 3 ) Then
          Do j=1,3
            READ(5,*,ERR=997) (zmagn(i,j),i=1,3)
            IF(DBG) Write(6,*) 'QUAX: zmagn(i,j) =',(zmagn(i,j),i=1,3)
          End Do
          LINENR=LINENR+3
        End If
c        Else
c        Write(6,'(A)') 'The CRYS keyword must be declared '//
c     &                 'above QUAX in the input!'
c        End If
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'PREX') Then
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'UBAR') Then
        compute_barrier=.TRUE.
        IF(DBG) Write(6,*) 'UBAR:'
        LINENR=LINENR+1
        Go To 100
      End If
C-------------------------------------------
      If (LINE(1:4).eq.'ABCC') Then
      Do_structure_abc = .TRUE.
      Read(5,*,ERR=997) (cryst(i),i=1,6)

      Do i=1,6
         If ( cryst(i)<=0 ) Then
            Call WarningMessage(2,
     &               'ABCC: zero or negative crystallographic '//
     &               'parameters requested! ')
            Call Quit_OnUserError()
         End If
      End Do

      IF(DBG) Write(6,*) 'ABCC: (cryst(i),i=1,6)=',(cryst(i),i=1,6)
      Read(5,*,ERR=997) (coord(i),i=1,3)
      IF(DBG) Write(6,*) 'ABCC: (coord(i),i=1,3)=',(coord(i),i=1,3)
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
      If (LINE(1:4).eq.'ZEEM') Then
        zeeman_energy=.true.
        compute_magnetization=.true.

        READ(5,*,ERR=997) nDirZee
        IF(DBG) Write(6,*) 'ZEEM: nDirZee=',nDirZee

        Do i=1,nDirZee
!         open the zeeman_energy_xxx.txt file where Zeeman eigenstates will
!         be further written in mangetization() subroutine
          Write(namefile_energy,'(5A)') 'zeeman_energy_',
     &                     CHAR(48+mod( int((i)/100),10)),
     &                     CHAR(48+mod( int((i)/10 ),10)),
     &                     CHAR(48+mod( int( i     ),10)),'.txt'
          !print *, 'namefile_energy: ', namefile_energy
          LUZee(i)=IsFreeUnit(30+i)
          Call molcas_open(LUZee(i),namefile_energy)
c          OPEN(30+i, FILE=namefile_energy)

          READ(5,*,ERR=997) (dir_weight(i,l),l=1,3)
          IF(DBG) Write(6,*) 'ZEEM: (dir_weight(i,l),l=1,3)=',
     &                              (dir_weight(i,l),l=1,3)

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
C-------------------------------------------

200   continue
      If(IPRINT.gt.2) Then
      Write(6,'(5X,A)') 'NO ERORR WAS LOCATED WHILE READING INPUT'
      End If



      If(compute_CF) Then
        If(axisoption.eq.3) Then
c check the determinant of the ZMAGN
        Det_zmagn=0.0_wp
        Do I=1,3
          Do J=1,3
          ZR(I,J)=0.0_wp
          ZR(I,J)=zmagn(I,J)
          End Do
        End Do
        Det_zmagn = FindDetR(ZR,3)
          If( Det_zmagn .lt. 0.0_wp ) Then
        Write(6,'(A)') 'QUAX: The determinant of the rotation matrix '//
     &                 'provided in the input is NEGATIVE.'
        Write(6,'(A,F22.14)') 'Determinant = ', Det_zmagn
        Write(6,'(A)') 'This means that the matrix you have provided '//
     &                 'can be decomposed in a product of two '
        Write(6,'(A)') 'matrices: Rotation*Inversion'
        Write(6,'(A)') 'The determinant of the Rotation matrix must '//
     &                 'be POSITIVE.'
        Write(6,'(A)') 'The program will stop.'
        Return
          End If

c check the orthogonality of the ZMAGN:
       Do i=1,3
         Do j=1,3
        column_check(i,j)=0.0_wp
        row_check(i,j)   =0.0_wp
           Do L=1,3
        column_check(i,j) = column_check(i,j) + zmagn(i,L) * zmagn(j,L)
           row_check(i,j) =    row_check(i,j) + zmagn(L,i) * zmagn(L,j)
           End Do
         End Do
       End Do

       Do i=1,3
         Do j=i+1,3
         If (i.eq.j) Go To 112
           If ( (ABS(column_check(1,2)).gt.0.0001_wp).OR.
     &          (ABS(column_check(1,3)).gt.0.0001_wp).OR.
     &          (ABS(column_check(2,3)).gt.0.0001_wp).OR.
     &          (ABS(   row_check(1,2)).gt.0.0001_wp).OR.
     &          (ABS(   row_check(1,3)).gt.0.0001_wp).OR.
     &          (ABS(   row_check(2,3)).gt.0.0001_wp) ) Then
          Write(6,'(A)') 'QUAX: The rotation matrix is not UNITARY.'
          Write(6,'(A,F19.14)') 'column_check(1,2) = ',column_check(1,2)
          Write(6,'(A,F19.14)') 'column_check(1,3) = ',column_check(1,3)
          Write(6,'(A,F19.14)') 'column_check(2,3) = ',column_check(2,3)
          Write(6,'(A,F19.14)') '   row_check(1,2) = ',   row_check(1,2)
          Write(6,'(A,F19.14)') '   row_check(1,3) = ',   row_check(1,3)
          Write(6,'(A,F19.14)') '   row_check(2,3) = ',   row_check(2,3)
          Write(6,'(A)') 'All above values must be exact 0.0.'
          Write(6,'(A)') 'Or at least less than than 0.0001.'
          Write(6,'(A)') 'Did you employ enough digits for '//
     &                   'the rotation matrix?'
          Write(6,'(A)') 'The program will stop.'
          Return
           End If
 112     continue
         End Do
       End Do
        End If ! axisoption
      End If ! compute_CF

c  preparing the info for computation of molar magnetization
      If(compute_magnetization) Then
c calculate the total number of directions for the average procedure
        nDirTot=0
        If(zeeman_energy) Then
        nDirTot=nDirZee
        End If
        If(compute_Mdir_vector) Then
        nDirTot=nDirTot+nDir
        End If
        nDirTot=nDirTot+get_nP(nsymm,ngrid)
      End If

C------ CHECK the data from INPUT ------------------------------
c      If(iprint.gt.10) Then

       If(dbg) Write(6,'(A,  F9.5)') 'ZJPR :         = ', zJ
       If(dbg) Write(6,'(A,  I3  )') 'PRLV :         = ',iprint

c      If (.not. compute_g_tensors) then
c         !generate an array of 10 low-lying groups of states
c         !
c         ndim(:)=0
c         nmult_try=10
c         j=0
c         ndim(i)=1
c         Do i=1, nss
c           etmp=eso(i)
c           do k=i+1,nss
c             if (abs(eso(k)-etmp) < 0.01_wp) ndim(i)=ndim(i)+1
c           enddo
c         End Do
c      End If

      If(compute_g_tensors) Then
        If(NMULT.gt.0) Then
        Write(6,'(A,I3)')   'MLTP :         = ', NMULT
           If(NMULT.le.20) Then
             Write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,NMULT)
           Else
             Write(6,'(A,20I3)') '               = ',(NDIM(i),i=1,20)
             Do j=21,NMULT,20
             jEnd=MIN(NMULT,J+19)
             Write(6,'(A,20I3)') '                 ',(NDIM(i),i=j,jEnd)
             End Do
           End If
        Else
        Write(6,'(A)') 'MLTP :         =  No pseudospin Hamiltonians'//
     &     ' will be computed. Is MLTP defined?'
        End If
      End If



      If((compute_CF).AND.(nDIMcf<=nss).AND.(lDIMcf<=nstate)) Then
        If (nlanth<15) Then
           Write(6,'(3A)') 'The Crystal-Field acting on the '//
     &                     'ground atomic multiplet of Ln = ',
     &                      clanth(nlanth),' is computed.'
        Else If (nlanth>=15 .and. nlanth<29) Then
           Write(6,'(3A)') 'The Crystal-Field acting on the '//
     &                     'ground atomic multiplet of Ac = ',
     &                      clanth(nlanth) ,' is computed.'
        Else If (nlanth>=29) Then
           Write(6,'(3A)') 'The Crystal-Field acting on the '//
     &                     'ground atomic |L,ML> multiplet of TM = ',
     &                      clanth(nlanth) ,' is computed.'
      End If


      Write(6,'(A,A )') 'CHIT :         = ',' molar magnetic '//
     & 'susceptibility is computed'
      If(TINPUT) Write(6,'(A)') 'TEXP :         = the experimental'//
     & ' temperature interval is read from the file "chitexp.input"'
      Write(6,'(A, I3)')  'TINT :      nT = ', nT
      Write(6,'(A,F7.3)') '          Tmin = ', Tmin
      Write(6,'(A,F7.3)') '          Tmax = ', Tmax


      End If
!--------------------------------------------------------------------
      If(compute_magnetization) Then

        Write(6,'(A,A )') 'MAGN :         = ',' molar magnetization'//
     &                    ' is computed'
        Write(6,'(A, I3)') 'NDIRTOT        = ', nDirTot
        Write(6,'(A, I3)') 'TMAG :         = ',nTempMagn
        Write(6,'(6x,A,20F7.3)') 'TempMagn = ',
     &                                 (TempMagn(i),i=1,nTempMagn)
        Write(6,'(A, I3)')  'HINT :      nH = ', nH
        Write(6,'(A,F7.3)') '          Hmin = ', Hmin
        Write(6,'(A,F7.3)') '          Hmax = ', Hmax
        Write(6,'(A, I3)') 'MAVE :   nDir = ', get_nP(nsymm,ngrid)

        If(HINPUT) Write(6,'(A)') 'HEXP :         = the experimental '//
     &                            'field interval is read from the '//
     &                            'file "mexp.input"'
        If(encut_definition.eq.1) Then
          Write(6,'(A, I3)')  'NCUT :         = ', ncut
        Else If(encut_definition.eq.2) Then
          Write(6,'(A,I4,a,i4)')  'ECUT :         = ', nk,', ',mg
        Else If(encut_definition.eq.3) Then
          Write(6,'(A,F7.3)') 'ERAT :         = ', encut_rate
        End If


        If(compute_Mdir_vector) Then
          Write(6,'(A,20I3)') 'MVEC :         = ', nDir
          If(nDir.gt.0) Then
            Do i=1,nDir
              Write(6,'(A,I2,A,3F11.6)') '   Dir :',i,' : ',
     &          dirX(i), dirY(i), dirZ(i)
            End Do
          End If
        End If
        If(zeeman_energy) Then
          If (nDirZee == 1) Then
             Write(6,'(2A,I2,1x,A)')
     &            'ZEEM :         = ',' Zeeman splitting '//
     &            'for the following direction of the '
             Write(6,'(18x,A)')
     &            'applied magnetic field is given in the'//
     &            ' "zeeman_energy_xxx.txt" file in $WorkDir/'
          Else If (nDirZee > 1) Then
             Write(6,'(2A,I2,1x,A)')
     &            'ZEEM :         = ',' Zeeman splitting '//
     &            'for the following',nDirZee,' directions of the '
             Write(6,'(18x,A)')
     &            'applied magnetic field are given in the'//
     &            ' "zeeman_energy_xxx.txt" files in $WorkDir/.'
          Else
             Write(6,'(A)') 'Error in input processing. nDirZee<0!'
             Call Quit_OnUserError()
          End If
          Do i=1, nDirZee
            Write(6,'(17x,3F11.6)')  (dir_weight(i,l),l=1,3)
          End Do
        End If
      End If !magnetization
!--------------------------------------------------------------------

      If(compute_torque) Then
        Write(6,'(A,A )') 'TORQ :         = ',' torque magnetization '//
     &                    'is computed'
      End If !torque

      If(doplot) Then
        Write(6,'(A,A )') 'PLOT :         = ',' GNUPLOT scripts and '//
     &                    'corresponding XT, M and UBAR plots will'//
     &                    'be generated'
      End If !




      If(Do_structure_abc) Then
        Write(6,'(2A)') 'ABCC :         = ','the main magnetic axes '//
     &                  'for the computed pseudospins are written '//
     &                  'also in the '
        Write(6,'( A)') 'crystallographic "abc" axes'
        Write(6,'(10x,A,F9.4)') 'a       = ', cryst(1)
        Write(6,'(10x,A,F9.4)') 'b       = ', cryst(2)
        Write(6,'(10x,A,F9.4)') 'c       = ', cryst(3)
        Write(6,'(10x,A,F9.4)') 'alpha   = ', cryst(4)
        Write(6,'(10x,A,F9.4)') 'beta    = ', cryst(5)
        Write(6,'(10x,A,F9.4)') 'gamma   = ', cryst(6)
        Write(6,'(10x,a,3F9.4)') 'coords: = ',(coord(i),i=1,3)
      End If

      If(compute_g_tensors) Then
        If(compute_barrier) Then
          Nblock=0
          Do i=1,nmult
            Nblock=Nblock+ndim(i)
          End Do
          Write(6,'(A,i4)') 'nBlock = ', nBlock
        End If
      End If




      Go To 190
C------ errors ------------------------------
      Write(6,*)' The following input line was not understood:'
      Write(6,'(A)') LINE
      Go To 999

997   continue
      Write(6,*)' READIN: Error reading standard input.'
      Write(6,*)' SINGLE_ANISO input near line nr.',LINENR+1
      Go To 999

998   continue
      Write(6,*)' READIN: Unexpected End of input file.'

999   continue
      Call XFLUSH(6)
      Call ABEnd()

590   continue
      Write(6,*) 'READIN: the TINT command is incompatible with TEXP'
      Call ABEnd()

591   continue
      Write(6,*) 'READIN: the HINT command is incompatible with HEXP'
      Call ABEnd()

595   continue
      Write(6,*) 'READIN: NCUT, ERAT and ENCU are mutually exclusive.'
      Call ABEnd()


 190  continue
      Call qExit('SA_readin')
      Return
      End
