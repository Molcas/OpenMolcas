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
* Copyright (C) 1993, Kurt Pfingst                                     *
************************************************************************
      SubRoutine Radlc(Zeta,nZeta,lsum,Rnr)
************************************************************************
*                                                                      *
* Object: to compute the radial part of the continuum  multipole       *
*         integrals within an R-matrix run                             *
*                                                                      *
* Called from: Mltint                                                  *
*                                                                      *
* Author: K.Pfingst 21/5/93                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "rmat.fh"
#include "real.fh"
#include "nrmf.fh"
      external fradf
      Character*80 Label
      Real*8 Zeta(nZeta), Rnr(nZeta,0:lsum)
      Parameter(limit=200,lenw=4*limit)
      Integer iScrt(limit)
      Real*8 Scrt(lenw)
*
*     Statement function for Cartesian index
*
      iRout = 122
      iPrint = nPrint(iRout)
ccccccccccccccccccccccccccccccccccccccc
c     iPrint = 99
ccccccccccccccccccccccccccccccccccccccc
      Call QEnter('Radlc')
      Result=Zero
*
************************************************************************
*
*
       Do 150 ir=0,lsum
        Do 250 iZeta=1,nZeta
         expsum=Zeta(iZeta)
         l=ir
*        If (quadpack) then
         If (quadpack.and..not.testint) then
            ier2=-1
            call dqag(fradf,0.0d0,Rmatr,Epsabs,Epsrel,Keyr,result2,
     &                abser,neval,ier2,
     &                limit,lenw,last,iScrt,Scrt)
            if (ier2.ne.0) then
               Call WarningMessage(1,
     &             ' WARNING in Radlc; Consult the output for details!')
               write(6,*)
               write(6,*) ' WARNING in Radlc'
               write(6,*)
               write(6,*)
     &                  ' ier=',ier2,' Error in Dqag called from Radlc.'
               write(6,*) ' result=',result2
               write(6,*) ' abser=',abser
               write(6,*) ' neval=',neval
               write(6,*)
            end if
            result=result2
*        Else If(Nagint) then
         Else If(Nagint.and..not.testint) then
            Call WarningMessage(2,
     &                          'Radlc: Nagint option not implemented!')
            Call Abend()
C           ier1=-1
C           call d01ajf(fradf,0.0d0,Rmatr,Epsabs,Epsrel,result1,abser,
C    *                  wrk1,4*INtparm,iwrk1,INtparm,ier1)
C           if (ier1.ne.0) then
C              write(*,*)
C              write(*,*) ' WARNING in Radlc'
C              write(*,*)
C              write(6,*)
C    &                ' ier=',ier1,' Error in D01ajf called from Radlc.'
C              write(6,*) ' result=',result1
C              write(6,*) ' abser=',abser
C              write(6,*) ' intparm=',intparm
C              write(6,*)
C           end if
C           result=result1
         EndIf
         If (testint) then
            Call WarningMessage(2,
     &                  'Radlc: testint option not implemented!')
            Call Abend()
C           ier2=-1
C           call dqag(fradf,0.0d0,Rmatr,Epsabs,Epsrel,Keyr,result2,
C    &                abser,neval,ier2,
C    &                limit,lenw,last,iScrt,Scrt)
C           if (ier2.ne.0) then
C              write(*,*)
C              write(*,*) ' WARNING in Radlc'
C              write(*,*)
C              write(6,*)
C    &                  ' ier=',ier2,' Error in Dqag called from Radlc.'
C              write(6,*) ' result=',result2
C              write(6,*) ' abser=',abser
C              write(6,*) ' neval=',neval
C              write(6,*)
C           end if
*
C           ier1=-1
C           call d01ajf(fradf,0.0d0,Rmatr,Epsabs,Epsrel,result1,abser,
C    *                  wrk1,4*INtparm,iwrk1,INtparm,ier1)
C           if (ier1.ne.0.or.ier2.ne.0) then
C              write(*,*)
C              write(*,*) ' WARNING in Radlc'
C              write(*,*)
C              write(6,*)
C    &                ' ier=',ier1,' Error in D01ajf called from Radlc.'
C              write(6,*) ' result=',result1
C              write(6,*) ' abser=',abser
C              write(6,*) ' intparm=',intparm
C              write(6,*)
C           end if
C           result=result1
*
C           diff=abs(result2-result1)
C           diff1=abs((result2-result1)/result2)
C           If (diff.gt.epsabs.or.diff1.gt.epsrel) then
C              write(*,*)
C              write(*,*) ' WARNING in Radlc'
C              write(*,*)
C              write(6,*) ' DIFFabs =',diff
C              write(*,*)
C              write(6,*) ' DIFFrel =',diff1
C              write(*,*)
C              write(6,*) ' NAG result=',result1
C              write(*,*)
C              write(6,*) ' NAG error =',ier1
C              write(*,*)
C              write(6,*) ' QUAD result=',result2
C              write(*,*)
C              write(6,*) ' QUAD error =',ier2
C              write(*,*)
C           End If
         End If
         Rnr(iZeta,ir)=result
250     continue
150    continue
*
************************************************************************
*
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in Radlc'
         Write (Label,'(A)') ' Rnr'
         Call RecPrt(Label,' ',Rnr(1,0),nZeta,lsum+1)
      End If
*
      Call QExit('Radlc')
      Return
      End
