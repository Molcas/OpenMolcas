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
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine will populate according to the aufbau principle.        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************
      Subroutine GoPop(Eps,Occ,Scr,n,PrtEor,PrThr,GapThr)
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
      Real*8  Eps(*)
      Real*8  Occ(*)
      Real*8  Scr(*)
      Integer n
      Logical PrtEor
      Real*8  PrThr
      Real*8  GapThr
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Real*8  TotNucChg
!     Real*8  eFermi
      Real*8  eGap
      Real*8  OccNo
      Real*8  eLo,eHi
      Real*8  tmp
      Integer nElec
      Integer nAlpha
      Integer nBeta
      Integer nOcc
      Integer nAct
      Integer kLo,kHi
      Integer i,j,k,m
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
!     GapThr=1.0d-2
!----------------------------------------------------------------------*
! Sort orbital energies                                                *
!----------------------------------------------------------------------*
      Do i=1,n
         Scr(i)=Eps(i)
      End Do
      Do i=1,n
         j=i
         Do k=i,n
            If(Scr(k).lt.Scr(j)) j=k
         End Do
         tmp=Scr(i)
         Scr(i)=Scr(j)
         Scr(j)=tmp
      End Do
!----------------------------------------------------------------------*
! How many alpha/beta                                                  *
!----------------------------------------------------------------------*
      Call Get_dScalar('Total nuclear Charge',TotNucChg)
      nElec  = Int(TotNucChg+0.5d0)
      nBeta  = Int(0.5d0*(TotNucChg+0.5d0))
      nAlpha = nElec-nBeta
!     Write(6,'(a,i5)') 'nElec . . . . . . . . . .',nElec
!     Write(6,'(a,i5)') 'nAlpha  . . . . . . . . .',nAlpha
!     Write(6,'(a,i5)') 'nBeta . . . . . . . . . .',nBeta
!----------------------------------------------------------------------*
! Optonally print sorted orbitals energies                             *
!----------------------------------------------------------------------*
      If(PrtEor) Then
         m=0
         Do i=1,n
            If(Scr(i).le.PrThr) m=i
         End Do
         Write(6,*)
         Write(6,'(a)') 'Sorted orbital energies'
         Write(6,'(a)') '-----------------------'
         Write(6,*)
         Write(6,'(a,i5,a,i5)') 'Printing',m,' out of',n
         Write(6,'(a,f6.1)') 'Filled orbitals:',0.5d0*TotNucChg
         Write(6,*)
         Write(6,'(i5,1h-,i5,2x,10f12.4)')                              &
     &      (i,Min(i+9,m),(Scr(j),j=i,Min(i+9,m)),i=1,m,10)
         Write(6,*)
      End If
!----------------------------------------------------------------------*
! Populate alpha                                                       *
!----------------------------------------------------------------------*
      If(nAlpha.ge.n) Then
!        Write(6,'(a)') 'Alpha is MB'
         eLo=Min(Eps(n)+1.0d-6,0.0d0)
         eHi=eLo
         OccNo=0.0d0
      Else If(nAlpha.le.0) Then
!        Write(6,'(a)') 'Alpha is empty'
         eLo=Eps(1)-1.0d0
         eHi=eLo
         OccNo=0.0d0
      Else
!        Write(6,'(a)') 'Alpha is not MB'
         If(nAlpha.gt.0) Then
!           eFermi = (Scr(nAlpha+1)+Scr(nAlpha))/2.0d0
            eGap   = (Scr(nAlpha+1)-Scr(nAlpha))
         Else
!           eFermi = 0.0D0
            eGap   = 0.0D0
         End If
!        Write(6,'(a,f12.6)') 'eFermi (alpha)  . . . . .',eFermi
!        Write(6,'(a,f12.6)') 'eGap (alpha)  . . . . . .',eGap
         If(eGap.gt.GapThr) Then
!           Write(6,'(a)') 'Alpha have large gap'
            OccNo=0.0d0
            eLo=0.25d0*Scr(nAlpha)+0.75d0*Scr(nAlpha+1)
            eHi=0.75d0*Scr(nAlpha)+0.25d0*Scr(nAlpha+1)
         Else
!           Write(6,'(a)') 'Alpha have small gap'
            kLo=1
            Do i=2,nAlpha
               If(Scr(i)-Scr(i-1).gt.GapThr) kLo=i
            End Do
            kHi=n
            Do i=n-1,Max(1,nAlpha),-1
               If(Scr(i+1)-Scr(i).gt.GapThr) kHi=i
            End Do
            If(kLo.gt.1) Then
               eLo=(Scr(kLo)+Scr(kLo-1))/2.0d0
            Else
               eLo=Scr(1)-1.0d0
            End If
            If(kHi.lt.n) Then
               eHi=(Scr(kHi)+Scr(kHi+1))/2.0d0
            Else
               eHi=Scr(n)+1.0d0
            End If
!           Write(6,'(a,i5)') 'kLo (alpha) . . . . . . .',kLo
!           Write(6,'(a,i5)') 'kHi (alpha) . . . . . . .',kHi
!           Write(6,'(a,f12.6)') 'eLo (alpha) . . . . . . .',eLo
!           Write(6,'(a,f12.6)') 'eHi (alpha) . . . . . . .',eHi
            nOcc=kLo-1
            nAct=kHi-kLo+1
            OccNo=1.0d0*(nAlpha-nOcc)/nAct
         End If
!        Write(6,'(a,f12.6)') 'OccNo (alpha) . . . . . .',OccNo
      End If
      Do i=1,n
         If(Eps(i).lt.eLo) Then
            Occ(i)=Occ(i)+1.0d0
         Else If(Eps(i).lt.eHi) Then
            Occ(i)=Occ(i)+OccNo
         End If
      End Do
!----------------------------------------------------------------------*
! Populate beta                                                        *
!----------------------------------------------------------------------*
      If(nBeta.ge.n) Then
!        Write(6,'(a)') 'Beta is MB'
         eLo=Min(Eps(n)+1.0d-6,0.0d0)
         eHi=eLo
         OccNo=0.0d0
      Else If(nBeta.le.0) Then
!        Write(6,'(a)') 'Beta is empty'
         eLo=Eps(1)-1.0d0
         eHi=eLo
         OccNo=0.0d0
      Else
!        Write(6,'(a)') 'Beta is not MB'
         If(nBeta.gt.0) Then
!           eFermi = (Scr(nBeta+1)+Scr(nBeta))/2.0d0
            eGap   = (Scr(nBeta+1)-Scr(nBeta))
         Else
!           eFermi = 0.0D0
            eGap   = 0.0D0
         End If
!        Write(6,'(a,f12.6)') 'eFermi (beta) . . . . . .',eFermi
!        Write(6,'(a,f12.6)') 'eGap (beta) . . . . . . .',eGap
         If(eGap.gt.GapThr) Then
!           Write(6,'(a)') 'Beta have large gap'
            OccNo=0.0d0
            eLo=0.25d0*Scr(nBeta)+0.75d0*Scr(nBeta+1)
            eHi=0.75d0*Scr(nBeta)+0.25d0*Scr(nBeta+1)
         Else
!           Write(6,'(a)') 'Beta have small gap'
            kLo=1
            Do i=2,nBeta
               If(Scr(i)-Scr(i-1).gt.GapThr) kLo=i
            End Do
            kHi=n
            Do i=n-1,Max(1,nBeta),-1
               If(Scr(i+1)-Scr(i).gt.GapThr) kHi=i
            End Do
            If(kLo.gt.1) Then
               eLo=(Scr(kLo)+Scr(kLo-1))/2.0d0
            Else
               eLo=Scr(1)-1.0d0
            End If
            If(kHi.lt.n) Then
               eHi=(Scr(kHi)+Scr(kHi+1))/2.0d0
            Else
               eHi=Scr(n)+1.0d0
            End If
!           Write(6,'(a,i5)') 'kLo (beta)  . . . . . . .',kLo
!           Write(6,'(a,i5)') 'kHi (beta)  . . . . . . .',kHi
!           Write(6,'(a,f12.6)') 'eLo (beta)  . . . . . . .',eLo
!           Write(6,'(a,f12.6)') 'eHi (beta)  . . . . . . .',eHi
            nOcc=kLo-1
            nAct=kHi-kLo+1
            OccNo=1.0d0*(nBeta-nOcc)/nAct
         End If
!        Write(6,'(a,f12.6)') 'OccNo (beta)  . . . . . .',OccNo
      End If
      Do i=1,n
         If(Eps(i).lt.eLo) Then
            Occ(i)=Occ(i)+1.0d0
         Else If(Eps(i).lt.eHi) Then
            Occ(i)=Occ(i)+OccNo
         End If
      End Do
!----------------------------------------------------------------------*
! Print population (debug)                                             *
!----------------------------------------------------------------------*
!     Write(6,*)
!     Write(6,'(a)') 'Occupation of orbitals'
!     Write(6,'(a)') '----------------------'
!     Write(6,*)
!     Write(6,'(a,i5,a,i5)') 'Printing',m,' out of',n
!     Write(6,*)
!     Write(6,'(i5,1h-,i5,2x,10f12.4)')                                 &
!    &   (i,Min(i+9,m),(Occ(j),j=i,Min(i+9,m)),i=1,m,10)
!----------------------------------------------------------------------*
! Done!                                                                *
!----------------------------------------------------------------------*
      Return
      End
