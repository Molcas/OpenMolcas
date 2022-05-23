************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine InpCtl_GENANO()
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      Character*72 Key,KeyIn
      Logical RdHead,done
*----------------------------------------------------------------------*
      Title='Atom'
      RdHead=.false.
      nSets=1
      center='ANO '
      thr=1.0d-8
      kRfSet=1
      rowise=.false.
      lftdeg=.false.
      iProj=0
      Do 100 i=0,MxLqn
         nPrim(i)=0
         nCore(i)=0
100   Continue
      Do 110 i=1,MxSets
         wSet(i)=0.0d0
110   Continue
      wSet(1)=1.0d0
      kSet=0
      isUHF=0
*----------------------------------------------------------------------*
      Call RdNlst(5,'GENANO')
*----------------------------------------------------------------------*
      done=.false.
200   Read(5,'(a)',end=300,err=300) KeyIn
         Key=KeyIn
*        Write(*,*) 'echo> ',KeyIn
         Call zlcase(Key)
         If(Key(1:1).eq.'*') Then
            Continue
         Else If(Key.eq.' ') Then
            Continue
         Else If(Key.eq.'title') Then
            RdHead=.true.
         Else If(Key.eq.'sets') Then
            RdHead=.false.
            Read(5,*,end=900,err=900) nSets
            If(nSets.gt.MxSets) Then
               Write(6,'(a,i3)') ' *** Too many sets specified:',nSets
               nSets=MxSets
               Write(6,'(a,i3)') ' *** Reduced to:             ',nSets
            End If
            Do 210 i=1,nSets
               wSet(i)=1.0d0/nSets
210         Continue
         Else If(Key.eq.'center') Then
            RdHead=.false.
            Read(5,'(a)',end=900,err=900) center
            Call Upcase(center)
         Else If(Key.eq.'no threshold'.or.
     &                 Key.eq.'nothreshold') Then
            read(5,*,end=900,err=900) thr
         Else If(Key.eq.'row wise'.or.Key.eq.'rowwise') Then
            RdHead=.false.
            rowise=.true.
         Else If(Key.eq.'lift degeneracy'.or.
     &              Key.eq.'liftdegeneracy') Then
            RdHead=.false.
            lftdeg=.true.
         Else If(Key.eq.'rydberg') Then
            RdHead=.false.
            rydgen=.true.
         Else If(Key.eq.'weights') Then
            RdHead=.false.
            Read(5,*) (wSet(i),i=1,nSets)
         Else If(Key.eq.'orbitalweights'.or.
     &           Key.eq.'orbital weights') Then
            Read(5,*) wc0,wc1
*           Write(6,'(a,f12.6)') 'inpctl: wc0',wc0
*           Write(6,'(a,f12.6)') 'inpctl: wc1',wc1
         Else If(Key.eq.'project') Then
            RdHead=.false.
            iProj=1
         Else If(Key.eq.'project 1') Then
            RdHead=.false.
            iProj=1
         Else If(Key.eq.'project 2') Then
            RdHead=.false.
            iProj=2
         Else If(Key.eq.'end of input') Then
            RdHead=.false.
            done=.true.
         Else If(RdHead) Then
            If(Title.eq.'Atom') Title=KeyIn
            Write(6,*) KeyIn(:mylen(KeyIn))
         Else
            Write(6,*) '*** Illegal keyword: ',KeyIn
            Call Quit_OnUserError()
         End If
         If(.not.done) Go To 200
*----------------------------------------------------------------------*
300   Continue
*----------------------------------------------------------------------*
      Return
*----------------------------------------------------------------------*
900   Continue
      Write(6,*) 'Error while reading input, keyword: ',key
      Call Quit_OnUserError()
*----------------------------------------------------------------------*
      End
