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
      Subroutine fetch_init_const( nneq, neqv, nmax, exch, nLoc,
     &                             nCenter, nT, nH, nTempMagn, nDir,
     &                             nDirZee, nMult, nPair, MxRank1,
     &                             MxRank2, iReturn )
c  this routine looks into the file "single_aniso.input" for the "RESTart" keyword
c
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "warnings.fh"
      Integer, intent(out) :: nneq, neqv, nmax, exch, nLoc,
     &                        nCenter, nT, nH, nTempMagn, nDir,
     &                        nDirZee, nMult, nPair, MxRank1, MxRank2,
     &                        iReturn
c local variables:
      Integer :: NMAXC
      Parameter (NMAXC=99)
      Integer :: i, j, linenr, Input, nTempMagn_HEXP, nTempMagn_TMAG,
     &           nH_HEXP, nH_HINT, nT_TEXP, nT_TINT
      Integer :: neqA(NMAXC), nexchA(NMAXC)
      Integer :: sfs_check(NMAXC)
      Integer :: sos_check(NMAXC)
      Integer :: imaxrank(NMAXC,2)
      Integer :: idummy
      Integer :: irank1, irank2, iline
      Real(kind=8):: rdummy
      Real(kind=8):: TempMagn(NMAXC)
      Logical :: ab_initio_all
      Logical :: KeyCoor, KeyPair, KeyHEXP, KeyTEXP, KeyHINT, KeyTINT,
     &           KeyTMAG, KeyMLTP, KeyMVEC, KeyNNEQ, KeyZEEM, KeyITOJ
      Integer :: LUANISO, Isfreeunit
      Character(Len=1)   :: itype(NMAXC)
      Character(Len=280) :: line
      Character(Len=180)  :: namefile_aniso
      Logical :: ifHDF
      Logical :: DBG
      External Isfreeunit

      Call qEnter('PA_fetch')

      iReturn=0
      nH=0
      nT=0
      nTempMagn=1
      nneq=0
      neqv=0
      nmax=0
      exch=0
      nCenter=0
      nDirZee=0
      nDir=0
      nMult=0
      nLoc=0
      nPair=0
      MxRank1=0
      MxRank2=0
      luaniso=0
      rdummy=0.0_wp
      neqA(1:nmaxc)=0
      nexchA(1:nmaxc)=0
      sfs_check(1:nmaxc)=0
      sos_check(1:nmaxc)=0
      ab_initio_all=.false.
      itype(1:nmaxc)=' '
      imaxrank(1:nmaxc,1:2)=0

c      namefile_aniso='              '
      ifHDF=.false.
      Input=5

      DBG=.false.


      KeyNNEQ=.false.
      KeyPair=.false.
      KeyCoor=.false.
      KeyHEXP=.false.
      KeyTEXP=.false.
      KeyTMAG=.false.
      KeyTINT=.false.
      KeyHINT=.false.
      KeyMLTP=.false.
      KeyMVEC=.false.
      KeyZEEM=.false.
      KeyITOJ=.false.
      nH_HEXP=0
      nH_HINT=0
      nT_TEXP=0
      nT_TINT=0
      nTempMagn_HEXP=0
      nTempMagn_TMAG=0
C=========== End of default settings====================================
      REWIND(Input)
50    READ(Input,'(A280)',End=998) LINE
      Call NORMAL(LINE)
      If(LINE(1:11).ne.'&POLY_ANISO') Go To 50
      LINENR=0
100   READ(Input,'(A280)',End=998) line
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If (LINE(1:1).eq.'*') Go To 100
      If (LINE.eq.' ') Go To 100
      If((LINE(1:4).ne.'NNEQ').AND.(LINE(1:4).ne.'TEXP').AND.
     &   (LINE(1:4).ne.'HEXP').AND.(LINE(1:4).ne.'END ').AND.
     &   (LINE(1:4).ne.'    ').AND.(LINE(1:4).ne.'HINT').AND.
     &   (LINE(1:4).ne.'TINT').AND.(LINE(1:4).ne.'TMAG').AND.
     &   (LINE(1:4).ne.'MVEC').AND.(LINE(1:4).ne.'ZEEM').AND.
     &   (LINE(1:4).ne.'MLTP').AND.(LINE(1:4).ne.'LIN1').AND.
     &   (LINE(1:4).ne.'LIN3').AND.(LINE(1:4).ne.'LIN9').AND.
     &   (LINE(1:4).ne.'PAIR').AND.(LINE(1:4).ne.'ALIN').AND.
     &   (LINE(1:4).ne.'COOR').AND.(LINE(1:4).ne.'ITOJ')) Go To 100
      If((LINE(1:4).eq.'END ').OR.(LINE(1:4).eq.'    '))  Go To 200

      If (line(1:4).eq.'NNEQ') Then

         KeyNNEQ=.true.
         READ(Input,*) nneq, ab_initio_all, ifHDF

         If(DBG) WRITE(6,*) nneq, ab_initio_all, ifHDF

         If(nneq<0) Then
            Write(6,'(A)') 'nneq<0! Must be positive!'
            Call Quit_OnUserError()
         Else If (nneq==0) Then
            Write(6,'(A)') 'nneq=0! Must be larger than zero!'
            Call Quit_OnUserError()
         Else If (nneq>NMAXC) Then
            Write(6,'(A)') 'nneq>99! Must be smaller than this!'
            Call Quit_OnUserError()
         End If

         READ(Input,*) (neqA(i),i=1,nneq)
         If(DBG) WRITE(6,*) (neqA(i),i=1,nneq)

         Do i=1, nneq
           If(neqA(i)<0) Then
              Write(6,'(A,i2,A)') 'neq(',i,')<0! Must be positive!'
              Call Quit_OnUserError()
           Else If (neqA(i)==0) Then
              Write(6,'(A,i2,A)') 'neq(',i,')=0! Must be larger '//
     &                            'than zero!'
              Call Quit_OnUserError()
           End If
         End Do

         READ(Input,*) (nexchA(i),i=1,nneq)
         If(DBG) WRITE(6,*) (nexchA(i),i=1,nneq)

         Do i=1, nneq
           If(nexchA(i)<0) Then
              Write(6,'(A,i2,A)') 'nexch(',i,')<0! Must be positive!'
              Call Quit_OnUserError()
           Else If (nexchA(i)==0) Then
              Write(6,'(A,i2,A)') 'nexch(',i,')=0! Must be larger '//
     &                            'than zero!'
              Call Quit_OnUserError()
           End If
         End Do

         neqv=0
         neqv=MAXVAL(neqA(1:nneq))
         If(DBG) Write(6,*) 'neqv = ', neqv
         nmax=0
         nmax=MAXVAL(nexchA(1:nneq))
         If(DBG) Write(6,*) 'nmax = ', nmax

         ! compute "exch"
         exch=1
         Do i=1,nneq
           Do j=1,neqA(i)
             exch=exch*nexchA(i)
           End Do
         End Do
         If(DBG) Write(6,*) 'exch=',exch

         If(ab_initio_all .eqv. .false.) Then
           READ(Input,*,ERR=997) (itype(i),i=1,Nneq)
         Else
           Do i=1,nneq
             itype(i)='A'
           End Do
         End If !ab_initio_all

         ! check the maximal number of local spin-orbit states
         sos_check=0
         sfs_check=0
         nCenter=0
         Do i=1,NNEQ
           nCenter=nCenter+neqA(I)
         End Do


         Do i=1,NNEQ
           If (itype(i) .eq. 'A') Then
             If( ifHDF ) Then
                ! generating the name of the "aniso_input file for
                ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
                If(i<10) Then
                   Write(namefile_aniso,'(4A)') 'aniso_hdf_',
     &                       CHAR(48+mod( int( i     ),10)),'.input'
                Else If(i>=10 .and. i<=99) Then
                   Write(namefile_aniso,'(4A)') 'aniso_hdf_',
     &                       CHAR(48+mod( int((i)/10 ),10)),
     &                       CHAR(48+mod( int( i     ),10)),'.input'
                End If

#ifdef _HDF5_
                Call read_hdf5_init( NAMEFILE_ANISO,sfs_check(I),
     &                               sos_check(I))
                If (DBG) Write(6,*) ' sfs(I) ',sfs_check(I),
     &                              ' sos(I) ',sos_check(I)
#else
                Call WarningMessage(2,'File '//trim(NAMEFILE_ANISO)//
     &                               ' cannot be opened. Molcas was'//
     &                               ' compiled without HDF5 option.')
                Call Quit_OnUserError()
#endif
             Else
                ! generating the name of the "aniso_input file for
                ! each center. Maxmimum 10 centers. CHAR(48)=0 (zero)
                If(i<10) Then
                   Write(namefile_aniso,'(4A)') 'aniso_',
     &                       CHAR(48+mod( int( i     ),10)),'.input'
                Else If(i>=10 .and. i<=99) Then
                   Write(namefile_aniso,'(4A)') 'aniso_',
     &                       CHAR(48+mod( int((i)/10 ),10)),
     &                       CHAR(48+mod( int( i     ),10)),'.input'
                End If
                LUANISO = Isfreeunit(20)
                Call molcas_open( LUANISO, NAMEFILE_ANISO)
                READ( LUANISO,*) sfs_check(I), sos_check(I)
                If (DBG) Write(6,*) ' sfs(I) ',sfs_check(I),
     &                              ' sos(I) ',sos_check(I)
                CLOSE(LUANISO)
             End If ! ifHDF
           Else If((itype(i).eq.'B').OR.(itype(i).eq.'C')) Then
             sfs_check(I)=1
             sos_check(I)=NexchA(i)
           End If


         End Do ! NNEQ

         nLoc=MAXVAL(sos_check(1:nneq))

         LINENR=LINENR+3
         Go To 100
      End If



      If (line(1:4).eq.'TEXP') Then

          KeyTexp=.true.
          READ(Input,*) nT_TEXP

          If ( nT_TEXP<=0 ) Then
             Call WarningMessage(2,
     &                 'TEXP: Number of temperature points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'TEXP:: = nT_TEXP', nT_TEXP

          LINENR=LINENR+1
          Go To 100
      End If



      If (line(1:4).eq.'HEXP') Then

          KeyHexp=.true.
          READ(Input,*) nTempMagn_HEXP, (TempMagn(i),i=1,nTempMagn)
          READ(Input,*) nH_HEXP

          If ( nH_HEXP<=0 ) Then
             Call WarningMessage(2,
     &                 'HEXP: Number of field points <= 0! ')
             Call Quit_OnUserError()
          End If

          If ( nTempMagn_HEXP<=0 ) Then
             Call WarningMessage(2,
     &                 'HEXP: Number of temperature points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'HEXP:: = nH_HEXP ', nH_HEXP
          If(DBG) WRITE(6,*) 'HEXP:: = nTempMagn_HEXP ', nTempMagn_HEXP
          LINENR=LINENR+2
          Go To 100
      End If



      If (line(1:4).eq.'HINT') Then

          KeyHINT=.true.
          READ(Input,*) rdummy, rdummy, nH_HINT

          If ( nH_HINT<=0 ) Then
             Call WarningMessage(2,
     &                 'HINT: Number of field points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'HINT:: = nH_HINT ', nH_HINT
          LINENR=LINENR+1
          Go To 100
      End If



      If (line(1:4).eq.'TINT') Then

          KeyTINT=.true.
          READ(Input,*) rdummy, rdummy, nT_TINT

          If ( nT_TINT<=0 ) Then
             Call WarningMessage(2,
     &                 'TINT: Number of temperature points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'TINT:: = nT_TINT ', nT_TINT
          LINENR=LINENR+1
          Go To 100
      End If




      If (line(1:4).eq.'TMAG') Then

          KeyTMAG=.true.
          READ(Input,*) nTempMagn_TMAG

          If ( nTempMagn_TMAG<=0 ) Then
             Call WarningMessage(2,
     &                 'TMAG: Number of temperatureMAGN points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'TMAG:: = nTempMagn_TMAG ',nTempMagn_TMAG
          LINENR=LINENR+1
          Go To 100
      End If



      If (line(1:4).eq.'MVEC') Then

          KeyMVEC=.true.
          READ(Input,*) nDir

          If ( nDir<=0 ) Then
             Call WarningMessage(2,
     &                 'MVEC: Number of nDir points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'MVEC:: = nDir ',nDir
          LINENR=LINENR+1
          Go To 100
      End If



      If (line(1:4).eq.'ZEEM') Then
          KeyZEEM=.false.
          READ(Input,*) nDirZee

          If ( nDirZee<=0 ) Then
             Call WarningMessage(2,
     &                 'ZEEM: Number of nDirZee points <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'ZEEM:: = nDirZee ',nDirZee
          LINENR=LINENR+1
          Go To 100
      End If



      If (line(1:4).eq.'MLTP') Then

          KeyMLTP=.true.
          READ(Input,*) nMult

          If ( nMult<=0 ) Then
             Call WarningMessage(2,
     &                 'MLTP: Number of multiplets <= 0! ')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'MLTP:: =nMult ',nMult
          LINENR=LINENR+1
          Go To 100
      End If



      If ( (LINE(1:4).eq.'LIN9') .OR. (LINE(1:4).eq.'LIN3') .OR.
     &     (LINE(1:4).eq.'LIN1') .OR. (LINE(1:4).eq.'ALIN') .OR.
     &     (LINE(1:4).eq.'PAIR') .OR. (LINE(1:4).eq.'ITOJ') ) Then

          KeyPair=.true.
          READ(Input,*,ERR=997) nPair

          If ( nPair<=0 ) Then
             Call WarningMessage(2,
     &             'PAIR OR LINx:: Number of interacting pairs <= 0!')
             Call Quit_OnUserError()
          End If

          If(DBG) WRITE(6,*) 'PAIR:: =nPair ',nPair


          If(LINE(1:4).eq.'ITOJ') Then
             iline=0
             Do i=1,npair
                imaxrank(i,1)=0
                imaxrank(i,2)=0

                READ(Input,*,ERR=997) idummy, idummy,
     &                              imaxrank(i,1),imaxrank(i,2)
                Do irank1=1,2*imaxrank(i,1)+1
                  Do irank2=1,2*imaxrank(i,2)+1
                    READ(Input,*,ERR=997) idummy,idummy,idummy,idummy,
     &                                    rdummy,rdummy
                    iline=iline+1
                  End Do
                End Do
             End Do ! i
             MxRank1=MAXVAL(imaxrank(1:npair,1))
             MxRank2=MAXVAL(imaxrank(1:npair,2))
             LINENR=LINENR+iline
          End If
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'COOR') Then
          KeyCoor=.true.
          LINENR=LINENR+1
          Go To 100
      End If

200   Continue

      If((.not.KeyPair).AND.(KeyCoor)) nPair=nCenter*(nCenter-1)/2

      If(KeyHexp) Then
         nTempMagn=nTempMagn_HEXP
         nH       =nH_HEXP
      Else
         nTempMagn=nTempMagn_TMAG
         nH       =nH_HINT
      End If

      If(KeyTEXP) Then
         nT       =nT_TEXP
      Else
         nT       =nT_TINT
      End If

      ! in case the user did not set up some of the above keywords
      ! assume the following default ones:
      If (     nMult == 0 ) nMult=1
      If (        nT == 0 ) nT = 31
      If (        nH == 0 ) nH = 11
      If ( nTempMagn == 0 ) nTempMagn = 1

c preliminary check the values:
      If (DBG) then
         Write(6,*) 'nneq     =',nneq
         Write(6,*) 'neqv     =',neqv
         Write(6,*) 'exch     =',exch
         Write(6,*) 'nLoc     =',nLoc
         Write(6,*) 'nmax     =',nmax
         Write(6,*) 'nCenter  =',nCenter
         Write(6,*) 'nT       =',nT
         Write(6,*) 'nH       =',nH
         Write(6,*) 'nTempMagn=',ntempMagn
         Write(6,*) 'nDir     =',nDir
         Write(6,*) 'nDirZee  =',nDirZee
         Write(6,*) 'nMult    =',nMult
         Write(6,*) 'nPair    =',nPair
         Write(6,*) 'MxRank1  =',MxRank1
         Write(6,*) 'MxRank2  =',MxRank2
      End If

      Go To 190
C------ errors ------------------------------
997   continue
      Write(6,*)' READIN: Error reading "poly_aniso.input" '
      Write(6,*)' near line nr.',LINENR+1
      Go To 999
998   continue
      Write(6,*)' READIN: Unexpected End of input file.'
999   continue
      Call quit(_RC_INPUT_ERROR_)


190   Continue
      Call qExit('PA_fetch')
      Return
      End
