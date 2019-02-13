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
      Subroutine restart_check( Ifrestart, input_to_read,
     &                          input_file_name, nT, nH, nTempMagn,
     &                          nDir, nDirZee, nMult, GRAD)
c  this routine looks into the file "single_aniso.input" for the "RESTart" keyword
c
      Implicit None
      Integer, Parameter        :: wp=selected_real_kind(p=15,r=307)
      Integer ::  linenr, input_to_read, Input, nT, nH, nTempMagn
      Integer ::  nDir, nDirZee, nMult, i
      Logical ::  Ifrestart
      Logical ::  GRAD
      Real    ::  rdummy
      Character(280) :: line, tmp
      Character(180) :: input_file_name
      Integer :: ncut,nk,mg
      Real    :: encut_rate
      Logical :: KeyREST,KeyTEXP,KeyHEXP,KeyHINT,KeyTINT,KeyTMAG,
     &           KeyMVEC,KeyZEEM,KeyMLTP,KeyNCUT,KeyENCU,KeyERAT,KeyGRAD
      Logical :: DBG

      Call qEnter('SA_rest_chk')
      DBG=.false.

      nH=0
      nT=0
      nMult=0
      nDirZee=0
      nDir=0
      nk=0
      mg=0
      ncut=0
      encut_rate=0.0_wp
      nTempMagn=0
      Input=5

      KeyREST=.false.
      KeyTEXP=.false.
      KeyHEXP=.false.
      KeyHINT=.false.
      KeyTINT=.false.
      KeyTMAG=.false.
      KeyMVEC=.false.
      KeyZEEM=.false.
      KeyMLTP=.false.
      KeyNCUT=.false.
      KeyENCU=.false.
      KeyERAT=.false.
      KeyGRAD=.false.

C=========== End of default settings====================================
      REWIND(Input)
50    READ(Input,'(A180)',End=998) LINE
      Call NORMAL(LINE)
      If(LINE(1:7).ne.'&SINGLE') Go To 50
      LINENR=0
100   READ(Input,'(A280)',End=998) line
      LINENR=LINENR+1
      Call NORMAL(LINE)
      If (LINE(1:1).eq.'*') Go To 100
      If (LINE.eq.' ') Go To 100

      If((LINE(1:4).ne.'REST').AND.(LINE(1:4).ne.'TEXP').AND.
     &   (LINE(1:4).ne.'HEXP').AND.(LINE(1:4).ne.'END ').AND.
     &   (LINE(1:4).ne.'    ').AND.(LINE(1:4).ne.'HINT').AND.
     &   (LINE(1:4).ne.'TINT').AND.(LINE(1:4).ne.'TMAG').AND.
     &   (LINE(1:4).ne.'MVEC').AND.(LINE(1:4).ne.'ZEEM').AND.
     &   (LINE(1:4).ne.'MLTP').AND.(LINE(1:4).ne.'NCUT').AND.
     &   (LINE(1:4).ne.'ENCU').AND.(LINE(1:4).ne.'ERAT').AND.
     &   (LINE(1:4).ne.'GRAD')                          ) Go To 100
      If((LINE(1:4).eq.'END ').OR. (LINE(1:4).eq.'    ')) Go To 200

      If (line(1:4).eq.'REST') Then
         Ifrestart=.true.
         KeyREST=.true.
         READ(Input,*) input_to_read
         input_file_name='aniso.input'
         If(DBG) WRITE(6,*) input_to_read
         If ( (input_to_read==2) .OR. (input_to_read==3) .OR.
     &        (input_to_read==4) ) Then
           BACKSPACE(Input)
           READ(Input,*) input_to_read, tmp
           If(DBG) WRITE(6,*) tmp
           input_file_name=trim(tmp)
           If(DBG) WRITE(6,*) 'restart_check: REST, input_file_name='
           If(DBG) WRITE(6,*) input_file_name
         End If
         LINENR=LINENR+1
         Go To 100
      End If

      If (line(1:4).eq.'TEXP') Then
          READ(Input,*) nT
          IF(DBG) WRITE(6,*) 'restart_check: TEXP, nT=', nT
          KeyTEXP=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'GRAD') Then
          KeyGRAD=.true.
          GRAD=.true.
          IF(DBG) WRITE(6,*) 'restart_check:  GRAD = ', GRAD
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'HEXP') Then
          READ(Input,*) nTempMagn, (rdummy,i=1,nTempMagn)
          READ(Input,*) nH
          IF(DBG) WRITE(6,*) 'restart_check: HEXP, nH=', nH
          IF(DBG) WRITE(6,*) 'restart_check: HEXP, nTempMagn=',nTempMagn
          KeyHEXP=.true.
          LINENR=LINENR+2
          Go To 100
      End If

      If (line(1:4).eq.'HINT') Then
          READ(Input,*) rdummy, rdummy, nH
          IF(DBG) WRITE(6,*) 'restart_check: HINT, nH=', nH
          KeyHINT=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'TINT') Then
          READ(Input,*) rdummy, rdummy, nT
          IF(DBG) WRITE(6,*) 'restart_check: HINT, nT=', nT
          KeyTINT=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'TMAG') Then
          READ(Input,*) nTempMagn
          IF(DBG) WRITE(6,*) 'restart_check: TMAG, nTempMagn=',nTempMagn
          KeyTMAG=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'MVEC') Then
          READ(Input,*) nDir
          IF(DBG) WRITE(6,*) 'restart_check: MVEC, nDir=',nDir
          KeyMVEC=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'ZEEM') Then
          READ(Input,*) nDirZee
          IF(DBG) WRITE(6,*) 'restart_check: ZEEM, nDirZee=',nDirZee
          KeyZEEM=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (line(1:4).eq.'MLTP') Then
          READ(Input,*) nMult
          IF(DBG) WRITE(6,*) 'restart_check: MLTP, nMult=',nMult
          KeyMLTP=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (LINE(1:4).eq.'NCUT') Then
          READ(Input,*) NCUT
          IF(DBG) WRITE(6,*) 'restart_check: NCUT, NCUT=',NCUT
          KeyNCUT=.true.
          LINENR=LINENR+1
          Go To 100
      End If

      If (LINE(1:4).eq.'ENCU') Then
          READ(Input,*) NK, MG
          IF(DBG) WRITE(6,*) 'restart_check: ENCU, NK, MG=',NK,MG
          LINENR=LINENR+1
          KeyENCU=.true.
          Go To 100
      End If

      If (LINE(1:4).eq.'ERAT') Then
          READ(Input,*) encut_rate
          IF(DBG) WRITE(6,*) 'restart_check: ERAT, encut_rate=',
     &                        encut_rate
          KeyERAT=.true.
          LINENR=LINENR+1
          Go To 100
      End If

200   Continue
      Write(6,'(5X,A)') 'NO ERORR WAS LOCATED WHILE READING INPUT'

c      print *,'KeyREST=',KeyREST
c      print *,'KeyTEXP=',KeyTEXP
c      print *,'KeyHEXP=',KeyHEXP
c      print *,'KeyHINT=',KeyHINT
c      print *,'KeyTINT=',KeyTINT
c      print *,'KeyTMAG=',KeyTMAG
c      print *,'KeyMVEC=',KeyMVEC
c      print *,'KeyZEEM=',KeyZEEM
c      print *,'KeyMLTP=',KeyMLTP
c      print *,'KeyNCUT=',KeyNCUT
c      print *,'KeyENCU=',KeyENCU
c      print *,'KeyERAT=',KeyERAT
c      print *,'KeyGRAD=',KeyGRAD

c      print *,'LOGLINE=',KeyTMAG.OR.KeyZEEM.OR.KeyMVEC.OR.KeyHINT.OR.
c     &                   KeyHEXP.OR.KeyNCUT.OR.KeyENCU.OR.KeyERAT

      If( KeyTMAG.OR.KeyZEEM.OR.KeyMVEC.OR.KeyHINT.OR.KeyHEXP.OR.
     &    KeyNCUT.OR.KeyENCU.OR.KeyERAT ) Then
          If(nTempMagn==0) nTempMagn=1
          If(       nH==0) nH=21
      End If
          If(       nT==0) nT=301

c      print *, 'nTempMagn=',nTempMagn
c      print *, 'nH       =',nH
c      print *, 'nT       =',nT




      Go To 190
C------ errors ------------------------------
998   continue
      Write(6,*)' -- READIN: Unexpected End of input file.'


190   Continue
      Call qExit('SA_rest_chk')
      Return
      End
