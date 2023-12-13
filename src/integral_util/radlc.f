!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1993, Kurt Pfingst                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Radlc(Zeta,nZeta,lsum,Rnr)
!***********************************************************************
!                                                                      *
! Object: to compute the radial part of the continuum multipole        *
!         integrals within an R-matrix run                             *
!                                                                      *
! Called from: Mltint                                                  *
!                                                                      *
! Author: K.Pfingst 21/5/93                                            *
!***********************************************************************
      use Constants, only: Zero
      use rmat, only: l, ExpSum, QuadPack, TestInt, NagInt, EpsAbs,
     &                EpsRel, KeyR, RMatR
      Implicit None
      Integer nZeta, lSum
      Real*8 Zeta(nZeta), Rnr(nZeta,0:lsum)

      Real*8, external :: fradf
      Integer, Parameter :: limit=200, lenw=4*limit
#ifdef _DEBUGPRINT_
      Character(LEN=80) Label
#endif
      Integer iScrt(limit)
      Real*8 Scrt(lenw)
      Real*8 Result, Result2, AbsEr
      Integer ir, iZeta, ier2, nEval, Last
!
!     Statement function for Cartesian index
!
      Result=Zero
      Call Untested('Radlc')
!
!***********************************************************************
!
!
       Do 150 ir=0,lsum
        Do 250 iZeta=1,nZeta
         expsum=Zeta(iZeta)
         l=ir
!        If (quadpack) then
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
!        Else If(Nagint) then
         Else If(Nagint.and..not.testint) then
            Call WarningMessage(2,
     &                          'Radlc: Nagint option not implemented!')
            Call Abend()
!           ier1=-1
!           call d01ajf(fradf,0.0d0,Rmatr,Epsabs,Epsrel,result1,abser,
!    *                  wrk1,4*INtparm,iwrk1,INtparm,ier1)
!           if (ier1.ne.0) then
!              write(*,*)
!              write(*,*) ' WARNING in Radlc'
!              write(*,*)
!              write(6,*)
!    &                ' ier=',ier1,' Error in D01ajf called from Radlc.'
!              write(6,*) ' result=',result1
!              write(6,*) ' abser=',abser
!              write(6,*) ' intparm=',intparm
!              write(6,*)
!           end if
!           result=result1
         EndIf
         If (testint) then
            Call WarningMessage(2,
     &                  'Radlc: testint option not implemented!')
            Call Abend()
!           ier2=-1
!           call dqag(fradf,0.0d0,Rmatr,Epsabs,Epsrel,Keyr,result2,
!    &                abser,neval,ier2,
!    &                limit,lenw,last,iScrt,Scrt)
!           if (ier2.ne.0) then
!              write(*,*)
!              write(*,*) ' WARNING in Radlc'
!              write(*,*)
!              write(6,*)
!    &                  ' ier=',ier2,' Error in Dqag called from Radlc.'
!              write(6,*) ' result=',result2
!              write(6,*) ' abser=',abser
!              write(6,*) ' neval=',neval
!              write(6,*)
!           end if
!
!           ier1=-1
!           call d01ajf(fradf,0.0d0,Rmatr,Epsabs,Epsrel,result1,abser,
!    *                  wrk1,4*INtparm,iwrk1,INtparm,ier1)
!           if (ier1.ne.0.or.ier2.ne.0) then
!              write(*,*)
!              write(*,*) ' WARNING in Radlc'
!              write(*,*)
!              write(6,*)
!    &                ' ier=',ier1,' Error in D01ajf called from Radlc.'
!              write(6,*) ' result=',result1
!              write(6,*) ' abser=',abser
!              write(6,*) ' intparm=',intparm
!              write(6,*)
!           end if
!           result=result1
!
!           diff=abs(result2-result1)
!           diff1=abs((result2-result1)/result2)
!           If (diff.gt.epsabs.or.diff1.gt.epsrel) then
!              write(*,*)
!              write(*,*) ' WARNING in Radlc'
!              write(*,*)
!              write(6,*) ' DIFFabs =',diff
!              write(*,*)
!              write(6,*) ' DIFFrel =',diff1
!              write(*,*)
!              write(6,*) ' NAG result=',result1
!              write(*,*)
!              write(6,*) ' NAG error =',ier1
!              write(*,*)
!              write(6,*) ' QUAD result=',result2
!              write(*,*)
!              write(6,*) ' QUAD error =',ier2
!              write(*,*)
!           End If
         End If
         Rnr(iZeta,ir)=result
250     continue
150    continue
!
!***********************************************************************
!
!
#ifdef _DEBUGPRINT_
      Write (6,*) ' Result in Radlc'
      Write (Label,'(A)') ' Rnr'
      Call RecPrt(Label,' ',Rnr(1,0),nZeta,lsum+1)
#endif
!
      End SubRoutine Radlc
