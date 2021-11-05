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
! Copyright (C) 2009, Giovanni Ghigo                                   *
!***********************************************************************
      Subroutine ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,           &
     &              lNMAT,lnTabDim, nTabDim,nMaxQ,nMat0,lVec)
!!
!!    Generate Logical Vector of useful States
!!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
      Dimension nMat0(nOsc)
      Integer nTabDim(0:lnTabDim), lVec(0:lnTabDim)
      Integer nMaxQ(nOsc)

      If (iPrint.GE.3) then
        Write(6,*) ' Original number of States=',max_nOrd+1
      EndIf
      Rewind (lNMAT)
      iIndex = 0
      Do iOrd = 0, max_nOrd
          iIndex = nTabDim(iOrd)
          Call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
          nSumQ = 0
          lVec(iOrd) = 1
          Do iOsc = 1, nOsc
            If (nMat0(iOsc).GT.nMaxQ(iOsc)) lVec(iOrd) = 0
            nSumQ = nSumQ + nMat0(iOsc)
          EndDo
          If (nSumQ.LT.minQ) lVec(iOrd) = 0
      EndDo
!!
      nYes = 0
      Do iOrd = 0,max_nOrd
        If (lVec(iOrd).EQ.1) nYes = nYes + 1
      EndDo
!!
      If (iPrint.GE.3) then
        Write(6,*) ' Selected number of States=',nYes
      EndIf
!!
      Return
      End


      Subroutine ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT,lnTabDim,     &
     & GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,                 &
     & dMinWind,dRho, nMat0,nTabDim,lVec,lTVec,EneMat)
!!
!!    Calculate Energy of Levels  GG 30-Dec-08 - 08-Jan-09
!!
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Constants_mula.fh"
#include "WrkSpc.fh"
      Real*8 GE1, GE2, harmfreq1(nOsc), harmfreq2(nOsc)
      Real*8 x_anharm1(nOsc,nOsc),x_anharm2(nOsc,nOsc)
      Real*8 dMinWind, dRho, dWlow, dWup
      Real*8 dEne, EneMat(0:max_nOrd)
      Integer nMat0(nOsc), nTabDim(0:lnTabDim)
      Integer lVec(0:lnTabDim), lTVec(0:lnTabDim)
      Logical lUpdate

      If (dMinWind.EQ.0.0d0) then
        lUpDate=.True.
        dMinWind = 1.0d0
      else
        lUpDate=.False.
      EndIf
      Do iOrd = 0, max_nOrd
        lTVec(iOrd) = lVec(iOrd)
      EndDo
!!
!!    Energy calculation
!!
      If (iPrint.GE.4) then
        Write(6,*)
        Write(6,*) ' States in the preliminar window :'
        If (nOsc.LE.24) then
          Write(6,'(a,108a)') '  ',('=',i=1,108)
          Write(6,*)'     jOrd    ene/au    ene/cm-1 '//                &
     &    'Vibrational quantum numbers'
          Write(6,'(a,108a)') '  ',('-',i=1,108)
        else
          Write(6,'(a,36a)') '  ',('=',i=1,36)
          Write(6,*)'        #    jOrd   ene/au      ene/cm-1 '
          Write(6,'(a,36a)') '  ',('-',i=1,36)
        EndIf
        Call XFlush(6)
      EndIf
!!
      Call GetMem('level1','Allo','Inte',iplevel1,nOsc)
      Call GetMem('level2','Allo','Inte',iplevel2,nOsc)
      Do iv=1,nOsc
        iWork(iplevel1+iv-1) = 0
      EndDo
      Rewind (lNMAT)
      Do iOrd = 0, max_nOrd
        If (lVec(iOrd).EQ.1) then
          iIndex = nTabDim(iOrd)
          Call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
          Do iv=1,nOsc
            iWork(iplevel2+iv-1) = nMat0(iv)
          EndDo
          l_harm=nOsc
          Call TransEnergy(                                             &
     &    GE1,x_anharm1,harmfreq1,iWork(iplevel1),                      &
     &    GE2,x_anharm2,harmfreq2,iWork(iplevel2),                      &
     &    dEne,l_harm)
          EneMat(iOrd) = dEne
          If (iPrint.GE.4) then
            If (nOsc.LE.24) then
              loc_n_max = 0
              Do j=1,nOsc
                loc_n_max = loc_n_max + nMat0(j)
              EndDo
              Write(6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,  &
     &            dEne*HarToRcm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
            else
               Write(6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne, &
     &            dEne*HarToRcm,loc_n_max
            EndIf
          EndIf
        EndIf
      End Do
!!
!!    Energy selection
!!
      If (iPrint.GE.3) then
        Write(6,*)
        Write(6,*) ' States in the window :'
        If (nOsc.LE.24) then
          Write(6,'(a,108a)') '  ',('=',i=1,108)
          Write(6,*)'     jOrd    ene/au    ene/cm-1 '//                &
     &    'Vibrational quantum numbers'
          Write(6,'(a,108a)') '  ',('-',i=1,108)
        else
          Write(6,'(a,36a)') '  ',('=',i=1,36)
          Write(6,*)'        #    jOrd   ene/au      ene/cm-1 '
          Write(6,'(a,36a)') '  ',('-',i=1,36)
        EndIf
        Call XFlush(6)
      EndIf
!!
      nYes_start = nYes
 100  Continue
      dWlow = 0.5d0*dMinWind/dRho
      dWup  = dWlow
      Do iOrd = 0,max_nOrd
        lVec(iOrd) = lTVec(iOrd)
        If (lVec(iOrd).EQ.1) then
          dEne = EneMat(iOrd)
          If (dEne.LT.-dWlow .or. dEne.GT.dWup) then
            lVec(iOrd) =  0
            nYes = nYes - 1
          else
            lVec(iOrd) = 1
            If (iPrint.GE.3) then
              If (nOsc.LE.24) then
                iIndex = nTabDim(iOrd)
                Call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
                loc_n_max = 0
                Do j=1,nOsc
                  loc_n_max = loc_n_max + nMat0(j)
                EndDo
                Write(6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,&
     &              dEne*HarToRcm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
              else
                Write(6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,&
     &              dEne*HarToRcm,loc_n_max
              EndIf
            EndIf
          EndIf
        EndIf
      EndDo
      If (nYes.LT.1 .and. lUpDate) then
        dMinWind = dMinWind + 1.0d0
        nYes = nYes_start
        GoTo 100
      EndIf
!!
      Call GetMem('level2','Free','Inte',iplevel2,nOsc)
      Call GetMem('level1','Free','Inte',iplevel1,nOsc)

      If (iPrint.GE.3) then
        If (nOsc.LE.30) Write(6,'(a,108a)') '  ',('-',i=1,108)
        If (nOsc.GT.30) Write(6,'(a,36a)') '  ',('-',i=1,36)
        Write(6,'(a,f12.9,a,f12.9,a)')'  Window: ',                     &
     &   -dWlow,         ' / ',dWup,         ' (au)'
        Write(6,'(a,f12.6,a,f12.6,a)')'  Window: ',                     &
     &   -dWlow*HarToRcm,' / ',dWup*HarToRcm,' (cm-1)'
      EndIf
      If (iPrint.GE.2) then
        Write(6,*) ' Final number of States=',nYes
      EndIf
      If (dMinWind.GT.1.0d0 .and. lUpDate .and. iPrint.GE.1) then
        Write(6,*)
        Write(6,*) ' *** Warning: Expansion factor has been set to ',   &
     &                                                     dMinWind
        Write(6,*)
      EndIf
      Call XFlush(6)
      Return
      End
