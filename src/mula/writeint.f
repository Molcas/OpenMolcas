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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine WriteInt(IntMat,TermMat,mMat,nMat,                     &
     &  OccNumMat1,OccNumMat2,MatEl,ForceField,E1,E2,                   &
     &  T0,harmfreq1,harmfreq2,x_anharm1,x_anharm2,                     &
     &  l_IntMat_1,l_IntMat_2,l_TermMat_1,l_TermMat_2,                  &
     &  nDimTot, nOsc)
!!
!!  Purpose:
!!    Write vibrational levels, intensities to log.
!!
!!  Input:
!!    FC       : Real*8 two dimensional array - Franck-Condon
!!               factors.
!!    Int_Mat  : Real*8 two dimensional array - Intensities.
!!    Term_Mat : Real*8 two dimensional array - Vibronic levels.
!!    mMat     : Integer two dimensional array - oscillator quanta for
!!               ground state.
!!    nMat     : Integer two dimensional array - oscillator quanta for
!!               excited state.
!!    m_plot,
!!    n_plot   : Integer array - transitions wanted in output.
!!
!!  Uses
!!    TabMod
!!    VibMod
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1996.
!!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
#include "indims.fh"
#include "inputdata.fh"

      Parameter ( nfreq = 2000 )
      Real*8   Intensity,max_Intensity
      Real*8 IntMat( 0:l_IntMat_1,0:l_IntMat_2 )
      Real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
      Real*8 OccNumMat1(0:nDimTot-1,nOsc),                              &
     &  OccNumMat2(0:nDimTot-1,nOsc)
      Integer mMat  (0:mdim1,mdim2)
      Integer nMat  (0:ndim1,ndim2)
      Character*23 mMatChar ( 0:mdim1+1 )
      Character*23 nMatChar( 0:ndim1+1 )
      Character*23   CharTemp
      Integer   TermMin,TermMax, ivee_cm, ivee_nm
      Real*8 E1(nDimTot),E2(nDimTot)
      Real*8 harmfreq1(nOsc),harmfreq2(nOsc)
      Real*8 x_anharm1(nOsc,nOsc),x_anharm2(nOsc,nOsc)
      Logical  MatEl,ForceField
      Character*12  Format
      Integer  nvTabDim

#include "inout.fh"
#include "WrkSpc.fh"
!!
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',                                           &
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',                                           &
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',                                           &
     &       '|       Results from intensity calculations       |'
      Write(6,'(a27,a)') ' ',                                           &
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',                                           &
     &       ' ================================================='
      Write(6,*)
      Write(6,*)
!!
!!---- Dimensions.
      max_mOrd = l_IntMat_1
      max_nOrd = l_IntMat_2
      nvar = nOsc
!!
      If(nvar.eq.1)  Format = '(a1,  i3,a2)'
      If(nvar.eq.2)  Format = '(a1, 2i3,a2)'
      If(nvar.eq.3)  Format = '(a1, 3i3,a2)'
      If(nvar.eq.4)  Format = '(a1, 4i3,a2)'
      If(nvar.eq.5)  Format = '(a1, 5i3,a2)'
      If(nvar.eq.6)  Format = '(a1, 6i2,a2)'
      If(nvar.eq.7)  Format = '(a1, 7i2,a2)'
      If(nvar.eq.8)  Format = '(a1, 8i1,a2)'
      If(nvar.eq.9)  Format = '(a1, 9i1,a2)'
      If(nvar.eq.10)  Format ='(a1,10i1,a2)'
      If(nvar.eq.11)  Format ='(a1,11i1,a2)'
      If(nvar.eq.12)  Format ='(a1,12i1,a2)'
      If(nvar.eq.13)  Format ='(a1,13i1,a2)'
      If(nvar.eq.14)  Format ='(a1,14i1,a2)'
      If(nvar.eq.15)  Format ='(a1,15i1,a2)'
      If(nvar.eq.16)  Format ='(a1,16i1,a2)'
      If(nvar.eq.17)  Format ='(a1,17i1,a2)'
      If(nvar.eq.18)  Format ='(a1,18i1,a2)'
      If(nvar.eq.19)  Format ='(a1,19i1,a2)'
      If(nvar.eq.20)  Format ='(a1,20i1,a2)'
      If ( nvar.lt.21 ) Then
      Do i = 0,mdim1
      Write(CharTemp,fmt=Format) '(',(mMat(i,j),j=1,nvar),')'
      mMatChar(i) = CharTemp
      End Do
      Do i = 0,ndim1
      Write(CharTemp,fmt=Format) '(',(nMat(i,j),j=1,nvar),')'
      nMatChar(i) = CharTemp
      End Do
      Else
      Do i = 0,max_mOrd
      Write(CharTemp,'(i5)') i
      mMatChar(i) = CharTemp
      End Do
      Do i = 0,max_nOrd
      Write(CharTemp,'(i5)') i
      nMatChar(i) = CharTemp
      End Do
!vv
      iprintLevel=0
      if(iprintLevel.eq.10) Then
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ',                                                   &
     &    'The meaning of n and nprime in FC-table below:'
      Write(6,fmt='(a1,a,40a3,a)')                                      &
     &    ' ','=============',('===',i=1,nvar),'='
!          Write(6,fmt='(40a3)',advance='no')  ('===',i=1,nvar)
!          Write(6,fmt='(   a)',advance='yes') '='
      Write(6,*)
      Write(6,*) ' ',                                                   &
     &    '     n                                  Oscillator quanta'
      Write(6,fmt='(a1,a,40a3,a)')                                      &
     &    ' ','-------------',('---',i=1,nvar),'-'
!          Write(6,fmt='(40a3)',advance='no')  ('---',i=1,nvar)
!          Write(6,fmt='(   a)',advance='yes') '-'
      If ( max_mOrd.gt.max_nOrd ) Then
      Do i = 0,max_mOrd
      Write(6,'(a4,i4,a5,40i3)') ' ',i,' ',                             &
     &          (mMat(i,j),j=1,nvar)
      End Do
      Else
      Do i = 0,max_nOrd
      Write(6,'(a4,i4,a5,40i3)') ' ',i,' ',                             &
     &          (nMat(i,j),j=1,nvar)
      End Do
      End If
      Write(6,fmt='(a1,a,40a3,a)')                                      &
     &      ' ','=============',('===',i=1,nvar),'='
!          Write(6,fmt='(40a3)',advance='no')  ('===',i=1,nvar)
!          Write(6,fmt='(   a)',advance='yes') '='
      EndIf
      End If
!!
!!---- Write vibrational levels.
      If ( WriteVibLevels ) Then
      Call GetMem('VibLevel1','Allo','Real',ipVibLevel1,max_mOrd+1)
      Call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
      If ( MatEl ) Then
      k = 0
      Do iOrd = 1,max_mOrd+1
      Work(ipVibLevel1+k) = (E1(iOrd)-E1(1))
      Work(ipVibLevel2+k) = T0+(E2(iOrd)-E2(1))
      k = k+1
      End Do
      Else
      Call GetMem('level1','Allo','Inte',iplevel1,nvar)
      Call GetMem('level2','Allo','Inte',iplevel2,nvar)

      G1 = 0.0d0
      G2 = G1+T0
      k = 0
      Work(ipVibLevel1+k) = 0.0d0
      do iv=1,nvar
      iWork(iplevel1+iv-1) = mMat(0,iv)
      enddo
      If ( max_mOrd.gt.0 ) Then
      Do iOrd = 1,max_mOrd
      k = k+1
      do iv=1,nvar
      iWork(iplevel2+iv-1)= mMat(iOrd,iv)
      enddo
      l_harm=nOsc
      Call TransEnergy(G1,x_anharm1,                                    &
     &             harmfreq1,iWork(iplevel1),                           &
     &             G1,x_anharm1,harmfreq1,iWork(iplevel2),              &
     &             Work(ipVibLevel1+k),l_harm)
      End Do
      End If
      k = 0
      Do iOrd = 0,max_nOrd
      do iv=1,nvar
      iWork(iplevel2+iv-1) = nMat(iOrd,iv)
      enddo
      l_harm=nOsc
      Call TransEnergy(G1,x_anharm2,                                    &
     &          harmfreq2,iWork(iplevel1),                              &
     &           G2,x_anharm2,harmfreq2,iWork(iplevel2),                &
     &          Work(ipVibLevel2+k),l_harm)
      k = k+1
      End Do
      Call GetMem('level1','Free','Inte',iplevel1,nvar)
      Call GetMem('level2','Free','Inte',iplevel2,nvar)
      End If
!!
      l_m=max_mOrd+1
      l_n=max_nOrd+1
      Call GetMem('VibSort1','Allo','Inte',ipVibSort1,l_m*2)
      Call GetMem('VibSort2','Allo','Inte',ipVibSort2,l_n*2)

      Do iOrd = 0,max_mOrd
      iWork(ipVibSort1+iOrd) =                                          &
     &       int(Work(ipVibLevel1+iOrd)*HarToRcm)
      iWork(ipVibSort1+iOrd+l_m) = iOrd
      End Do
      Do iOrd = 0,max_nOrd
      iwork(ipVibSort2+iOrd) =                                          &
     &       int(Work(ipVibLevel2+iOrd)*HarToRcm)
      iWork(ipVibSort2+iOrd+l_n) = iOrd
      End Do
      Call GetMem('VibLevel1','Free','Real',ipVibLevel1,max_mOrd+1)
      Call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)
!!
      If ( max_mOrd.gt.0 ) Then
      Do i = 0,max_mOrd
      Do j = i+1,max_mOrd
      If ( iWork(ipVibSort1+j).lt.iWork(ipVibSort1+i)) Then
      itemp = iWork(ipVibSort1+j)
      iwork(ipVibSort1+j) = iWork(ipVibSort1+i)
      iWork(ipVibSort1+i) = itemp
      itemp = iWork(ipVibSort1+j+l_m)
      iwork(ipVibSort1+j+l_m) =                                         &
     &                iWork(ipVibSort1+i+l_m)
      iwork(ipVibSort1+i+l_m) = itemp
      End If
      End Do
      End Do
      End If
      If ( max_nOrd.gt.0 ) Then
      Do i = 0,max_nOrd
      Do j = i+1,max_nOrd
      If ( iwork(ipVibSort2+j).lt.iWork(ipVibSort2+i)) Then
      itemp = iwork(ipVibSort2+j)
      iWork(ipVibSort2+j) = iwork(ipVibSort2+i)
      iWork(ipVibSort2+i) = itemp
      itemp = iWork(ipVibSort2+j+l_n)
      iWork(ipVibSort2+j+l_n) =                                         &
     &                iWork(ipVibSort2+i+l_n)
      iWork(ipVibSort2+i+l_n) = itemp
      End If
      End Do
      End Do
      End If
!!
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Vibrational levels for ground state'
      Write(6,*) ' ',                                                   &
     &    '=============================================='
      Write(6,*)
      Write(6,*) ' ',                                                   &
     &    '     Level                     Energy(cm-1)'
      Write(6,*) ' ',                                                   &
     &    '----------------------------------------------'
      Do i = 0,max_mOrd
      iOrd = iwork(ipVibSort1+i+l_m)
      Write(6,'(a4,a15,a15,i7)') ' ',mMatChar(iOrd),' ',                &
     &       iWork(ipVibSort1+i)
      End Do
      Write(6,*) ' ',                                                   &
     &    '=============================================='
!!
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Vibrational levels for excited state'
      Write(6,*) ' ',                                                   &
     &    '=============================================='
      Write(6,*)
      Write(6,*) ' ','     Level                     Energy(cm-1)'
      Write(6,*) ' ',                                                   &
     &    '----------------------------------------------'
      Do i = 0,max_nOrd
      iOrd = iWork(ipVibSort2+i+l_n)
      Write(6,'(a4,a15,a15,i7)') ' ',nMatChar(iOrd),' ',                &
     &       iWork(ipVibSort2+i)
      End Do
      Write(6,*) ' ',                                                   &
     &    '=============================================='
      Call GetMem('VibSort1','Free','Inte',ipVibSort1,l_m*2)
      Call GetMem('VibSort2','Free','Inte',ipVibSort2,l_n*2)

      End If
!!
!!---- Find start and end points for a given total quanta in mMat and nMat.
      m_plot_max = l_m_plot
      n_plot_max = l_n_plot
      max_mQuanta = iWork(ipm_plot)
      If ( m_plot_max.gt.1 ) Then
      Do i = 2,m_plot_max
      If ( iWork(ipm_plot+i-1).gt.max_mQuanta ) Then
      max_mQuanta = iWork(ipm_plot+i-1)
      End If
      End Do
      End If
      max_nQuanta = iWork(ipn_plot)
      If ( n_plot_max.gt.1 ) Then
      Do i = 2,n_plot_max
      If ( iWork(ipn_plot+i-1).gt.max_nQuanta ) Then
      max_nQuanta = iWork(ipn_plot+i-1)
      End If
      End Do
      End If
      maxQuanta = max(max_mQuanta,max_nQuanta)
      Call GetMem('mMatStart','Allo','Inte',ipmMatStart,maxQuanta+1)
      Call GetMem('mMatStop','Allo','Inte',ipmMatStop,maxQuanta+1)
      iWork(ipmMatStart) = 0
      iWork(ipmMatStop) = 0
      Do i = 1,maxQuanta
      iWork(ipmMatStart+i) = iWork(ipmMatStop+i-1)+1
      Call TabDim2_drv(i,nvar,nvTabDim)
      iWork(ipmMatStop+i) = nvTabDim-1
      End Do
!!
!!---- Intensities.
      n=(l_TermMat_1+2)*(l_TermMat_2+2)
      l_TermSort=n+1
      Call GetMem('TermSort','Allo','Inte',ipTermSort,l_TermSort*3)
      nval = 0
      max_Intensity = 0.0d0
!!
      If ( MatEl.and.( .not.ForceField )) Then
      If ( max_mQuanta.le.max_nQuanta ) Then
      Do i = 1,m_plot_max
      Do iOrd = iWork(ipmMatStart+iwork(ipm_plot+i-1)),                 &
     &  iWork(ipmMatStop+iWork(ipm_plot+i-1))
      Do jOrd = 0,3*max_mOrd/4
      iWork(ipTermSort+nval) =                                          &
     &                int(TermMat(iOrd,jOrd)*HarToRcm)
      iWork(ipTermSort+nval+l_TermSort) = iOrd
      iWork(ipTermSort+nval+l_TermSort*2) = jOrd
      nval = nval+1
      End Do
      End Do
      End Do
      Else
      Do iOrd = 0,3*max_mOrd/4
      Do j = 1,n_plot_max
      Do jOrd = iWork(ipmMatStart+iWork(ipn_plot+j-1)),                 &
     &  iWork(ipmMatStop+iWork(ipn_plot+j-1))
      If ( IntMat(iOrd,jOrd).gt.max_Intensity )                         &
     &                max_Intensity = IntMat(iOrd,jOrd)
      iWork(ipTermSort+nval) =                                          &
     &                int(Abs(TermMat(iOrd,jOrd))*                      &
     &                HarToRcm)
      iWork(ipTermSort+nval+l_TermSort) = iOrd
      iWork(ipTermSort+nval+l_TermSort*2) = jOrd
      nval = nval+1
      End Do
      End Do
      End Do
      End If
      Else
      Do i = 1,m_plot_max
      Do j = 1,n_plot_max
      Do iOrd = iWork(ipmMatStart+iWork(ipm_plot+i-1)),                 &
     &          iWork(ipmMatStop+iWork(ipm_plot+i-1))
      Do jOrd = iwork(ipmMatStart+iWork(ipn_plot+j-1)),                 &
     &             iWork(ipmMatStop+iwork(ipn_plot+j-1))
      If ( IntMat(iOrd,jOrd).gt.max_Intensity )                         &
     &                max_Intensity = IntMat(iOrd,jOrd)
      End Do
      End Do
      End Do
      End Do
      Do i = 1,m_plot_max
      Do j = 1,n_plot_max
      Do iOrd = iWork(ipmMatStart+iWork(ipm_plot+i-1)),                 &
     &          iWork(ipmMatStop+iWork(ipm_plot+i-1))
      Do jOrd = iWork(ipmMatStart+iWork(ipn_plot+j-1)),                 &
     &             iWork(ipmMatStop+iWork(ipn_plot+j-1))
!!                If (( IntMat(iOrd,jOrd).gt.int_thrs*max_Intensity ).or. &
!!                If (( IntMat(iOrd,jOrd).gt.1.0d-3*max_Intensity ).or. &
!!                    (( iOrd.eq.0 ).and.( jOrd.eq.0 ))) Then
      iWork(ipTermSort+nval) =                                          &
     &                   int(Abs(TermMat(iOrd,jOrd))*HarToRcm)
      iWork(ipTermSort+nval+l_TermSort) = iOrd
      iWork(ipTermSort+nval+l_TermSort*2) = jOrd
      nval = nval+1
!!                End If
      End Do
      End Do
      End Do
      End Do
      End If
!!
      nval = nval-1
      Do i = 0,nval
        Do j = i+1,nval
          If ( iWork(ipTermSort+j).lt.iWork(ipTermSort+i) ) Then
            itemp = iWork(ipTermSort+j)
            iWork(ipTermSort+j) = iWork(ipTermSort+i)
            iWork(ipTermSort+i) = itemp
            itemp = iWork(ipTermSort+j+l_TermSort)
            iWork(ipTermSort+j+l_TermSort) =                            &
     &              iWork(ipTermSort+i+l_TermSort)
            iWork(ipTermSort+i+l_TermSort) = itemp
            itemp = iWork(ipTermSort+j+l_TermSort*2)
            iWork(ipTermSort+j+l_TermSort*2) =                          &
     &              iWork(ipTermSort+i+l_TermSort*2)
            iWork(ipTermSort+i+l_TermSort*2) = itemp
          End If
        End Do
      End Do

! MatEl,ForceField,m_plot_max,n_plot_max== F T 1 2/5
!!
      Write(6,*)
      Write(6,'(90A)') ' ',('=',iv=1,89)
      If ( OscStr ) Then
        Write(6,'(A)') '  Ground                        Excited'//      &
     &  '                         Energy         Oscillator            '
        Write(6,'(A)') '  State                         State'//        &
     &  '                      cm-1 /  eV / nm      Strength          '
      Else
        Write(6,'(A)') '  Ground                        Excited'//      &
     &  '                         Energy        Intensities           '
        Write(6,'(A)') '  State                         State'//        &
     &  '                      cm-1 /  eV / nm                        '
      End If
      Write(6,'(90A)') ' ',('-',iv=1,89)
      const = 1.0d0

      Do i = 0,nval
        iOrd = iWork(ipTermSort+i+l_TermSort)
        jOrd = iWork(ipTermSort+i+l_TermSort*2)
        vee_cm=iWork(ipTermSort+i)
        vee_nm=1.0d7/vee_cm
        vee_eV=vee_cm/8065.6d0
        vee=vee_eV
        If (Use_cm) vee=vee_cm
        If (Use_nm) vee=vee_nm
        ivee_cm=INT(vee_cm+0.5d0)
        ivee_nm=INT(vee_nm+0.5d0)
        Intensity = IntMat(iOrd,jOrd)/const
        If (( MatEl ).and.( .not.ForceField )) Then
          If ( m_plot_max.lt.n_plot_max ) Then
            Write(6,'(a1,a15,a5,a1,f5.2,a1,f5.2,a1,f5.2,'//             &
     &      'a1,a5,i7,a10,e12.3)') ' ',                                 &
     &      mMatChar(iOrd),'---> ','(',OccNumMat2(jOrd,1),',',          &
     &      OccNumMat2(jOrd,2),',',OccNumMat2(jOrd,3),')',' ',          &
     &      ivee_cm,' ',Intensity
            if(Intensity.gt.1D-6) Then
              Call Add_Info('Energy',[vee_cm],1,5)
              Call Add_Info('Intensity',[Intensity],1,5)
            endif
          Else
            Write(6,'(a1,a15,a5,a1,f5.2,a1,f5.2,'//                     &
     &      'a1,f5.2,a1,a5,i7,a10,e12.3)') ' ',                         &
     &      mMatChar(iOrd),'.lt.--- ','(',OccNumMat2(jOrd,1),',',       &
     &      OccNumMat2(jOrd,2),',',OccNumMat2(jOrd,3),')',' ',          &
     &      ivee_cm,' ',Intensity
            if(Intensity.gt.1D-6) Then
              Call Add_Info('Energy',[vee_cm],1,5)
              Call Add_Info('Intensity',[Intensity],1,5)
            endif
          End If
        Else
          If ( m_plot_max.lt.n_plot_max ) Then
            Write(6,'(a1,a23,a7,a23,a4,i6,a1,f5.2,a1,i4,a2,e12.3)') ' ',&
     &      mMatChar(iOrd),' --->  ',                                   &
     &      nMatChar(jOrd),' ',ivee_cm,'/',vee_eV,'/',ivee_nm,'  ',     &
     &      Intensity
            if(Intensity.gt.1D-6) Then
              Call Add_Info('Energy',[vee],1,5)
              Call Add_Info('Intensity',[Intensity],1,5)
            endif
          Else
            Write(6,'(a1,a23,a7,a23,a4,i6,a1,f5.2,a1,i4,a2,e12.3)') ' ',&
     &      mMatChar(iOrd),' .lt.---  ',                                &
     &      nMatChar(jOrd),' ',ivee_cm,'/',vee_eV,'/',ivee_nm,'  ',     &
     &      Intensity
            if(Intensity.gt.1D-6) Then
              Call Add_Info('Energy',[vee],1,5)
              Call Add_Info('Intensity',[Intensity],1,5)
            endif
          End If
        End If
      End Do
      Write(6,'(90A)') ' ',('=',iv=1,89)
      Write(6,*)
      Write(6,*)
!
!***********************************************************************
!
!!---- Write to plot file.
!!
      conv = 0.0d0
      If (.NOT.Use_nm .AND. .NOT.Use_cm) conv = 1.239842d-4
      If (Use_cm) conv = 1.0d0
      If (Use_nm) conv = 1.0d7

! TermMin/Max must be in cm-1
      If(plotwindow) then
        If (Use_nm) then                       ! nm -> cm-1
          TermMin=int(conv/cmend  +0.999999D0) ! Note the inversion
          TermMax=int(conv/cmstart-0.999999D0) ! cmstart<=>cmend
        else If (Use_cm) then                  ! alreay cm-1
          TermMin=int(cmstart)
          TermMax=int(cmend+0.999999D0)
        else                                   ! eV -> cm-1
          TermMin=int(cmstart/conv+0.999999D0)
          TermMax=int(cmend/conv-0.999999D0)
        EndIf
      else
        TermMin = int(iWork(ipTermSort))-nfreq
        TermMax = int(iWork(ipTermSort+nval))+nfreq
      EndIf
      l_plotvec=TermMax-TermMin+1
      Call GetMem('plotvec','Allo','Real',ipplotvec,l_plotvec)
      call dcopy_(l_plotvec,[0.0d0],0,Work(ipplotvec),1)

      call molcas_open(plotunit,'plot.intensity')
!      Open (Unit=plotUnit,File='plot.intensity')
      If( broadplot ) Then
        FWHM = hbarcm/(2.0d0*LifeTime)
        G2 = FWHM/2.0d0
        Do iTrans = 0,nval
          iOrd = iWork(ipTermSort+iTrans+l_TermSort)
          jOrd = iWork(ipTermSort+iTrans+l_TermSort*2)
          ifreq = iWork(ipTermSort+iTrans)
          Do ipoint = TermMin,TermMax
          Work(ipplotVec+ipoint-TermMin) =                              &
     &          Work(ipplotVec+ipoint-TermMin)+                         &
     &                     IntMat(iOrd,jOrd)*(G2**2)/((G2**2)+          &
     &     (ipoint-ifreq)**2)
          End Do
        End Do
        Do ipoint = TermMin,TermMax
          If ( Work(ipplotVec+ipoint-TermMin).gt.0.0D0 ) Then
            If (.NOT.Use_nm) then
              Write(plotUnit,'(f12.6,e15.6)')                           &
     &          ipoint*conv,Work(ipplotVec+ipoint-TermMin)
            else
              Write(plotUnit,'(f12.6,e15.6)')                           &
     &          conv/ipoint,Work(ipplotVec+ipoint-TermMin)
            EndIf
          End If
        End Do
      Else
        Do iTrans = 0,nval
          iOrd = iWork(ipTermSort+iTrans+l_TermSort)
          jOrd = iWork(ipTermSort+iTrans+l_TermSort*2)
          ifreq = iWork(ipTermSort+iTrans)
          If (TermMin.le.ifreq .and. ifreq.le.TermMax) then
            Work(ipplotVec+ifreq-TermMin) =                             &
     &      Work(ipplotVec+ifreq-TermMin)+IntMat(iOrd,jOrd)
          EndIf
        End Do
        Do ipoint = TermMin,TermMax
          If ( Work(ipplotVec+ipoint-TermMin).gt.0.0d0 ) Then
            If (.NOT.Use_nm) then
              Write(plotUnit,'(f12.6,e15.6)') ipoint*conv,0.0d0
              Write(plotUnit,'(f12.6,e15.6)') ipoint*conv,              &
     &                      Work(ipplotVec+ipoint-TermMin)
              Write(plotUnit,'(f12.6,e15.6)') ipoint*conv,0.0d0
            else
              Write(plotUnit,'(f12.6,e15.6)') conv/ipoint,0.0d0
              Write(plotUnit,'(f12.6,e15.6)') conv/ipoint,              &
     &                      Work(ipplotVec+ipoint-TermMin)
              Write(plotUnit,'(f12.6,e15.6)') conv/ipoint,0.0d0
            EndIf
          End If
        End Do
!        Call GetMem('plotvec','Free','Real',ipplotvec,l_plotvec)
      End If
      Call GetMem('plotvec','Free','Real',ipplotvec,l_plotvec)
      Close ( plotUnit )
!!
      Call GetMem('TermSort','Free','Inte',ipTermSort,l_TermSort*3)
      Call GetMem('mMatStart','Free','Inte',ipmMatStart,maxQuanta+1)
      Call GetMem('mMatStop','Free','Inte',ipmMatStop,maxQuanta+1)
!!
! Avoid unused argument warnings
      If (.False.) Call Unused_real_array(OccNumMat1)
      End
