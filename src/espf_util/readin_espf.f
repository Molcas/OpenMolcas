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
      Subroutine ReadIn_ESPF(natom,ipCord,ipExt,MltOrd,iRMax,DeltaR,
     &                       Forces,Show_espf,ipIsMM,StandAlone,iGrdTyp,
     &                       DoTinker,DoGromacs,DynExtPot,ipMltp,natMM,
     &                       lMorok,DoDirect,ipGradCl,EnergyCl)
      Implicit Real*8 (a-h,o-z)
*
#include "espf.fh"
#include "opt_mmo.fh"
*
#include "print.fh"
      Character*180 Key,Line,PotFile,UpKey
      Character*10 ESPFKey
      Character*12 ExtPotFormat
      Logical Convert,DoTinker_old,DoTinker,DoGromacs_old,DoGromacs,
     &        Exist,Forces,Show_espf,StandAlone,DynExtPot,lMorok_old,
     &        lMorok,NoExt,DoDirect_old,DoDirect
      Save fift
      Data fift/1.5d1/
      Character*180 Get_Ln
      External Get_Ln
*
      Call qEnter('ReadIn')
*
* If some keywords are not given, what are the defauts ?
* 3 cases:
*   1) ESPF.DATA does not exist:
*      MULT: MltOrd = 0 (monopole)
*      GRID: Type = PNT ; iRMax = 4 shells ; DeltaR = 1 angstrom
*      EXTE: MANDATORY
*   2) ESPF.DATA exists:
*      Get back all values from ESPF.DATA
*      a) If the keyword "Forces" is read, only the old values
*      are retained, so it is clever not to include any other
*      ESPF keyword
*      b) If the keyword "Forces" is not read, all the new values
*      are compared to the old ones.
*
*
* Initialize values
*
      Write(ExtPotFormat,'(a4,i2,a6)') '(I4,',MxExtPotComp,'F10.5)'
      MltOrd = 1
      MltOrd_old = MltOrd
      nChg = -1
      iGrdTyp = 1
      iRMax = 4
*      iGrdTyp = 2
*      iRMax = 1
      DeltaR = One/Angstrom
      Convert = .False.
      DoTinker = .False.
      DoTinker_old = DoTinker
      DoGromacs = .False.
      DoGromacs_old = DoGromacs
      Forces = .Not.StandAlone
      Show_espf = .False.
      nMult = 0
      DynExtPot = .False.
      iQMChg = 1
      natMM = 0
      lMorok = .False.
      lMorok_old = lMorok
      NoExt = .False.
      DoDirect = .False.
      DoDirect_old = DoDirect
      nOrd_ext = 0
      MMIterMax = 0
      ConvF = 2.0D-4*AuToKjPerMolNm
*
* Print level
*
      iPL = iPL_espf()
*
* If the ESPF.DATA file exists, retrieve all data
* from it in "*_old" variables and arrays.
*
      IPotFl=15
      PotFile='***'
      Call F_Inquire('ESPF.DATA',Exist)
      If (Exist) Then
         IPotFl = IsFreeUnit(IPotFl)
         Call Molcas_Open(IPotFl,'ESPF.DATA')
10       Line = Get_Ln(IPotFl)
         ESPFKey = Line(1:10)
         If (ESPFKey.eq.'MLTORD    ') Then
            Call Get_I1(2,MltOrd_old)
            ibla = 0
            Do ii = 0, MltOrd_old
               ibla = ibla + (ii+2)*(ii+1)/2
            End Do
            MltOrd_old = ibla
         Else If (ESPFKey.eq.'IRMAX     ') Then
            Call Get_I1(2,iRMax_old)
         Else If (ESPFKey.eq.'DELTAR    ') Then
            Call Get_F1(2,DeltaR_old)
         Else If (ESPFKey.eq.'GRIDTYPE  ') Then
            Call Get_I1(2,iGrdTyp_old)
         Else If (ESPFKey.eq.'MULTIPOLE ') Then
            Call Get_I1(2,nMult)
            Call GetMem('ESPFMltp','ALLO','REAL',ipMltp,nMult)
            Do iMlt = 1, nMult, MltOrd_old
               Line = Get_Ln(IPotFl)
               Call Get_I1(1,iAt)
               Call Get_F(2,Work(ipMltp+iMlt-1),MltOrd_old)
            End Do
         Else If (ESPFKey.eq.'TINKER    ') Then
            DoTinker_old = .True.
         Else if (ESPFKey.eq.'GROMACS   ') Then
            DoGromacs_old = .True.
         Else If (ESPFKey.eq.'DIRECT    ') Then
            DoDirect_old = .True.
         Else If (ESPFKey.eq.'LA_MOROK  ') Then
            lMorok_old = .True.
         Else If (ESPFKey.eq.'ENDOFESPF ') Then
            Goto 11
         EndIf
         Goto 10
11       Close (IPotFl)
         iRMax = iRMax_old
         DeltaR = DeltaR_old
         MltOrd = MltOrd_old
         iGrdTyp = iGrdTyp_old
         DoTinker = DoTinker_old
         DoGromacs = DoGromacs_old
         lMorok = lMorok_old
      EndIf
      If (.Not.StandAlone) Go To 1999
*
* Copy input from standard input to a local scratch file
*
      LuSpool = isFreeUnit(IPotFl)
      Call SpoolInp(LuSpool)
*
*---- Locate "start of input"
      Rewind(LuSpool)
      Call RdNLst(LuSpool,'espf')
*
  999 Continue
      Key = Get_Ln(LuSpool)
      Line = Key
      Call UpCase(Line)
      If (Line(1:4).eq.'MULT') GoTo 1000
      If (Line(1:4).eq.'EXTE') GoTo 1001
      If (Line(1:4).eq.'GRID') GoTo 1002
      If (Line(1:4).eq.'FORC') Goto 1003
      If (Line(1:4).eq.'SHOW') Goto 1004
      If (Line(1:4).eq.'LAMO') Goto 1005
      If (Line(1:4).eq.'MMIT') Goto 1006
      If (Line(1:4).eq.'MMCO') Goto 1007
      If (Line(1:4).eq.'END ') GoTo 1999
      If (.not.Exist) Then
         Write (6,*)' Unidentified keyword:', Key
         Call FindErrorLine
         Call Quit_OnUserError()
      EndIf
*
*>>>>>>>>>>>>> MULT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1000 Continue
      Key = Get_Ln(LuSpool)
      Call Get_I1(1,MltOrd)
      If (MltOrd.lt.0) Then
         Write(6,'(A)')' Error in espf/readin: MltOrd < 0!'
         Call Quit_OnUserError()
      End If
      If (DoGromacs.and.(MltOrd.gt.0)) Then
         Write(6,'(A)')' Error in espf/readin: Gromacs calculation'//
     &                 ' requested with MltOrd > 0'
         Write(6,'(A)')' Only MltOrd = 0 is currently allowed'
         Call Quit_OnUserError()
      End If
      If (MltOrd.gt.1) Then
         Write(6,'(A)')' Error in espf/readin: MltOrd > 1 NYI!'
         Call Quit_OnUserError()
      End If
      ibla = 0
      Do ii = 0, MltOrd
         ibla = ibla + (ii+2)*(ii+1)/2
      End Do
      MltOrd = ibla
      GoTo 999
*
*>>>>>>>>>>>>> EXTE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1001 Continue
      If (Forces) Goto 999
      Key = Get_Ln(LuSpool)
      UpKey = Key
      Call Upcase(UpKey)
      Call Get_iNumber(Key(1:(Index(Key,' ')-1)),ibla,iErr)
      If (iErr .eq. 0) Then
         PotFile='* *'
         nChg = ibla
*
* nChg < 0: error
* nChg > 0: external potential is given as point charges and dipoles,
* like for the seward xfield keyword
* nchg = 0: external potential directly given on atom centers as:
* pot field_x field_y field_z dfield_xx dfield_yy dfield_zz
* dfield_xy dfield_xz dfield_yz (ONE LINE per CENTER)
*
*
         If (nChg .lt. 0) Then
            Write(6,*) 'Error in readin_espf: nChg < 0!'
            Call Quit_OnUserError()
         Else If (nChg.gt.0) Then
            Convert = (Index(UpKey,'ANGSTROM').ne.0)
            Call GetMem('PtChg','ALLO','REAL',ipPC,nChg*7)
            Do iChg = 1, nChg
               Key = Get_Ln(LuSpool)
               iCurr = ipPC+(iChg-1)*7
               Call Get_F(1,Work(iCurr),7)
               If (Convert) Then
                  Work(iCurr  ) = Work(iCurr  )/Angstrom
                  Work(iCurr+1) = Work(iCurr+1)/Angstrom
                  Work(iCurr+2) = Work(iCurr+2)/Angstrom
                  Work(iCurr+4) = Work(iCurr+4)*Angstrom
                  Work(iCurr+5) = Work(iCurr+5)*Angstrom
                  Work(iCurr+6) = Work(iCurr+6)*Angstrom
               End If
            End Do
            ESelf = SelfEn(nChg,ipPC)
            Convert = .False.
         Else
            Do iAt = 1, natom
               Key = Get_Ln(LuSpool)
               Call Get_I1(1,jAt)
               If ((jAt.lt.1).or.(jAt.gt.natom)) Then
                  Write(6,'(A)')
     &               ' Error in espf/readin: atom out of range.'
                  Call Quit_OnUserError()
               End If
               Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),
     &                    MxExtPotComp)
            End Do
         End If
      Else
         iAt = Index(UpKey,'NONE ')
         NoExt = (iAt .ne. 0)
*
* Is it a QM/MM computation ?
*
         iAt = Index(UpKey,'TINKER ')
         DoTinker = (iAt .ne. 0)
*
         iAt = Index(UpKey,'GROMACS ')
         If (iAt .ne. 0) Then
#ifdef _GROMACS_
            DoGromacs = .True.
#else
            Write(6,*) 'Interface to Gromacs not installed'
            Call Quit_OnUserError()
#endif
         End If
*
         DoDirect = (Index(UpKey(iAt+7:120),'DIRECT').ne.0)
ctmp
         If (DoDirect) Then
           Write(6,*) 'Direct not yet implemented, abort.'
           Call Quit_OnUserError()
         End If
ctmp
*
* What kind of charges Tinker or Gromacs will use in the
* microiterations
*
         If (NoExt) Then
            PotFile = '* *'
         Else If (DoTinker.Or.DoGromacs) Then
            PotFile = 'ESPF.EXTPOT'
            If (Index(UpKey(iAt+7:120),'MULL').ne.0) Then
               iQMChg = 2
            Else If (Index(UpKey(iAt+7:120),'LOPR').ne.0) Then
               iQMChg = 3
            End If
            If (DoDirect .and. iQMChg .eq. 1) iQMChg = 2
         Else
            PotFile = Key(1:Len(Key))
         End If
      End If
      GoTo 999
*
*>>>>>>>>>>>>> GRID <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1002 Continue
      Key = Get_Ln(LuSpool)
      Call Upcase(Key)
      If (Index(Key,'GEPOL').ne.0) Then
         iGrdTyp = 2
         Call Get_I1(2,iRMax)
         If(iRMax.le.0 .or. iRMax.gt.4) Then
           Write(6,'(A)')'Error in readin_espf: 1 <= iRMax <= 4 !!!'
           Call Quit_OnUserError()
         End If
      Else If(Index(Key,'PNT').ne.0) Then
         iGrdTyp = 1
         Call Get_I1(2,iRMax)
         If(iRMax.le.0) Then
           Write(6,'(A)')'Error in espf/readin: iRMax < 1 !!!'
           Call Quit_OnUserError()
         End If
         Call Get_F1(3,DeltaR)
         If(DeltaR.le.Zero) Then
           Write(6,'(A)')'Error in espf/readin: DeltaR < 0.0 !!!'
           Call Quit_OnUserError()
         End If
         If (Index(Key,'ANGSTROM').ne.0) Convert=.True.
         If(Convert) DeltaR = DeltaR/Angstrom
         Convert = .False.
      Else
         Write(6,'(A)') 'Unrecognized GRID: GEPOL or PNT(default)'
         Call Quit_OnUserError()
      End If
      GoTo 999
*
*>>>>>>>>>>>>> FORC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1003 Continue
      If (.not.Exist) Then
         Write(6,*)'Error! Forces: the ESPF data are missing'
         Call Quit_OnUserError()
      EndIf
      If (DoTinker) Then
         Write(6,*)'Please erase the @Tinker call together with Forces'
         Call Quit_OnUserError()
      EndIf
      Forces = .True.
      Write(6,'(/,A,/,A)')' This ESPF run will compute energy gradient',
     &  ' Any other keyword is ignored !'
*
* Here I assume all I need can be retrieved from the $Project.espf file
*
      GoTo 999
*
*>>>>>>>>>>>>> SHOW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1004 Continue
      Show_espf = .True.
      GoTo 999
*
*>>>>>>>>>>>>> LAMO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1005 Continue
      lMorok = .True.
      If (iPL.ge.2) Write(6,'(A)') ' Morokuma scheme on'
      GoTo 999
*
*>>>>>>>>>>>>> MMIT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1006 Continue
      If (.Not.DoGromacs) Then
         Write(6,'(A)')
     &      ' MM microiterations only available with Gromacs'
         Call Quit_OnUserError()
      End If
      Line = Get_Ln(LuSpool)
      Call Get_I1(1,MMIterMax)
      GoTo 999
*
*>>>>>>>>>>>>> MMCO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1007 Continue
      Line = Get_Ln(LuSpool)
      Call Get_F1(1,ConvF)
      ConvF = ConvF*AuToKjPerMolNm
      GoTo 999
*
*>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1999 Continue
*
* Remove local copy of standard input
*
      If (StandAlone) Close(LuSpool)
*
* "Forces" case: retrieve all the data and update the MM gradient
*
      If (Forces) Then
         MltOrd = MltOrd_old
         iRMax = iRMax_old
         DeltaR = DeltaR_old
         iGrdTyp =iGrdTyp_old
         DoTinker = DoTinker_old
         DoGromacs = DoGromacs_old
         iQMChg = 0
         If (DoTinker) Call RunTinker(natom,ipCord,ipMltp,ipIsMM,
     &               MltOrd,DynExtPot,iQMchg,natMM,StandAlone,DoDirect)
#ifdef _GROMACS_
         If (DoGromacs) Call RunGromacs(natom,Work(ipCord),ipMltp,
     &                                  MltOrd,Forces,ipGradCl,EnergyCl)
#endif
         If (nAtMM.ne.0) Write(6,*) 'MM gradients have been updated'
         lMorok = lMorok_old
         IPotFl = IsFreeUnit(IPotFl)
         Call Molcas_Open(IPotFl,'ESPF.EXTPOT')
         Line = Get_Ln(IPotFl)
         Call Get_I1(1,nChg)
         Do iAt = 1, natom
            Line = Get_Ln(IPotFl)
            Call Get_I1(1,jAt)
            Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
         End Do
         Close(IPotFl)
*
* No external potential
*
      Else If (NoExt) Then
         nChg = -1
         Write(6,'(/,A)')' No external electrostatic potential'
*
* External potential read from a file
*
      Else If (nChg .eq. -1) Then
         If (DoTinker) Call RunTinker(natom,ipCord,ipMltp,ipIsMM,
     &                  MltOrd,DynExtPot,iQMChg,natMM,StandAlone,
     &                  DoDirect)
#ifdef _GROMACS_
         If (DoGromacs) Call RunGromacs(natom,Work(ipCord),ipMltp,
     &                                  MltOrd,Forces,ipGradCl,EnergyCl)
#endif
         LuSpool = IsFreeUnit(1)
         Call Molcas_Open (LuSpool,PotFile(1:(Index(PotFile,' ')-1)))
         If (iPL.ge.3) Write(6,'(/,A,A)')' External potential read in ',
     &                           PotFile(1:(Index(PotFile,' ')-1))
         Key = Get_Ln(LuSpool)
         UpKey = Key
         Call Upcase(UpKey)
         Call Get_I1(1,nChg)
         If (nChg .lt. 0) Then
            Write(6,*) 'Error in readin_espf: nChg < 0!'
            Call Quit_OnUserError()
         Else If (nChg.gt.0) Then
            Call Get_I1(2,nOrd_ext)
            iShift = 4+3*nOrd_ext
            Convert = (Index(UpKey,'ANGSTROM').ne.0)
            Call GetMem('PtChg','ALLO','REAL',ipPC,nChg*iShift)
            Do iChg = 1, nChg
               Key = Get_Ln(LuSpool)
               iCurr = ipPC+(iChg-1)*iShift
               Call Get_F(1,Work(iCurr),iShift)
               If (Convert) Then
                  Work(iCurr  ) = Work(iCurr  )/Angstrom
                  Work(iCurr+1) = Work(iCurr+1)/Angstrom
                  Work(iCurr+2) = Work(iCurr+2)/Angstrom
                  If (nOrd_ext .ne. 0) Then
                     Work(iCurr+4) = Work(iCurr+4)*Angstrom
                     Work(iCurr+5) = Work(iCurr+5)*Angstrom
                     Work(iCurr+6) = Work(iCurr+6)*Angstrom
                  End If
               End If
            End Do
            If (.not.(DoTinker.Or.DoGromacs)) ESelf = SelfEn(nChg,ipPC)
            Convert = .False.
         Else
            Do iAt = 1, natom
               Key = Get_Ln(LuSpool)
               Call Get_I1(1,jAt)
               If ((jAt.lt.1).or.(jAt.gt.natom)) Then
                  Write(6,'(A)')
     &               ' Error in espf/readin: atom out of range.'
                  Call Quit_OnUserError()
               End If
               Call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),
     &                    MxExtPotComp)
            End Do
         End If
         Close(LuSpool)
      End If
*
* If nChg > 0, 2 possibilities:
* a) read external point charges (only, no dipoles) for a direct
*    QM/MM electrostatic coupling
* b) external potential calculated from point charges and dipoles
*
      If (nChg .gt. 0 .and. DoDirect) Then
         lXF = .True.
         nXF = nChg
         nOrd_XF = nOrd_ext
         iXPolType = 0
         nXMolnr = 0
         ipXF = ipPC
         Call GetMem('PtChg','FREE','REAL',ipPC,nChg*(4+3*nOrd_ext))
      Else If (nChg .gt. 0) Then
         Do iAt = 1, natom
            Do iChg = 1, nChg
               dx = Work(ipCord+(iAt-1)*3  )-Work(ipPC+(iChg-1)*7  )
               dy = Work(ipCord+(iAt-1)*3+1)-Work(ipPC+(iChg-1)*7+1)
               dz = Work(ipCord+(iAt-1)*3+2)-Work(ipPC+(iChg-1)*7+2)
               qChg = Work(ipPC+(iChg-1)*7+3)
               dpxChg = Work(ipPC+(iChg-1)*7+4)
               dpyChg = Work(ipPC+(iChg-1)*7+5)
               dpzChg = Work(ipPC+(iChg-1)*7+6)
               rAtChg = sqrt(dx*dx+dy*dy+dz*dz)

               rAC2 = rAtChg * rAtChg
               rAC3 = rAtChg * rAC2
               rAC5 = rAC2 * rAC3
               rAC7 = rAC2 * rAC5
               iStart = ipExt+(iAt-1)*MxExtPotComp
*      Potential E
               Work(iStart) = Work(iStart)
     &                 + qChg/rAtChg
     &                 - (dpxChg*dx + dpyChg*dy + dpzChg*dz)/rAC3
*      Field F / x
               Work(iStart+1) = Work(iStart+1)
     &                 - qChg*dx/rAC3
     &                 + (dpxChg*(three*dx*dx-rAC2)
     &                  + dpyChg*(three*dx*dy)
     &                  + dpzChg*(three*dx*dz))/rAC5
*      Field F / y
               Work(iStart+2) = Work(iStart+2)
     &                 - qChg*dy/rAC3
     &                 + (dpxChg*(three*dy*dx)
     &                  + dpyChg*(three*dy*dy-rAC2)
     &                  + dpzChg*(three*dy*dz))/rAC5
*      Field F / z
               Work(iStart+3) = Work(iStart+3)
     &                 - qChg*dz/rAC3
     &                 + (dpxChg*(three*dz*dx)
     &                  + dpyChg*(three*dz*dy)
     &                  + dpzChg*(three*dz*dz-rAC2))/rAC5
*      Gradient G / xx
               Work(iStart+4) = Work(iStart+4)
     &                 + qChg*(three*dx*dx-rAC2)/rAC5
     &                 - (dpxChg*(fift*dx*dx*dx-rnine*dx*rAC2)
     &                  + dpyChg*(fift*dx*dx*dy-three*dy*rAC2)
     &                  + dpzChg*(fift*dx*dx*dz-three*dz*rAC2))/rAC7
*      Gradient G / yy
               Work(iStart+5) = Work(iStart+5)
     &                 + qChg*(three*dy*dy-rAC2)/rAC5
     &                 - (dpxChg*(fift*dy*dy*dx-three*dx*rAC2)
     &                  + dpyChg*(fift*dy*dy*dy-rnine*dy*rAC2)
     &                  + dpzChg*(fift*dy*dy*dz-three*dz*rAC2))/rAC7
*      Gradient G / zz
               Work(iStart+6) = Work(iStart+6)
     &                 + qChg*(three*dz*dz-rAC2)/rAC5
     &                 - (dpxChg*(fift*dz*dz*dx-three*dx*rAC2)
     &                  + dpyChg*(fift*dz*dz*dy-three*dy*rAC2)
     &                  + dpzChg*(fift*dz*dz*dz-rnine*dz*rAC2))/rAC7
*      Gradient G / xy
               Work(iStart+7) = Work(iStart+7)
     &                 + qChg*(three*dx*dy)/rAC5
     &                 - (dpxChg*(fift*dx*dy*dx-three*dx*rAC2)
     &                  + dpyChg*(fift*dx*dy*dy-three*dy*rAC2)
     &                  + dpzChg*(fift*dx*dy*dz))/rAC7
*      Gradient G / xz
               Work(iStart+8) = Work(iStart+8)
     &                 + qChg*(three*dx*dz)/rAC5
     &                 - (dpxChg*(fift*dx*dz*dx-three*dx*rAC2)
     &                  + dpyChg*(fift*dx*dz*dy)
     &                  + dpzChg*(fift*dx*dz*dz-three*dz*rAC2))/rAC7
*      Gradient G / yz
               Work(iStart+9) = Work(iStart+9)
     &                 + qChg*(three*dy*dz)/rAC5
     &                 - (dpxChg*(fift*dy*dz*dx)
     &                  + dpyChg*(fift*dy*dz*dy-three*dy*rAC2)
     &                  + dpzChg*(fift*dy*dz*dz-three*dz*rAC2))/rAC7
            End Do
         End Do
         Call GetMem('PtChg','FREE','REAL',ipPC,nChg*7)
      End If
*
* Check the compatibility between old and new keywords
*
      If (Exist) Then
         If (.not.Forces.and.
     &         ((MltOrd.ne.MltOrd_old).or.
     &          (iRMax.ne.iRMax_old).or.
     &          (iGrdTyp.ne.iGrdTyp_old).or.
     &          (lMorok.neqv.lMorok_old).or.
     &          (DoDirect.neqv.DoDirect_old).or.
     &          (Abs(DeltaR-DeltaR_old).gt.1.0d-6))) Then
            Write(6,*)'Conficts between some old and new ESPF keywords'
            Write(6,*)'      ',
     &       'MltOrd     iRMax    DeltaR    iGrdTyp   lMorok   DoDirect'
            Write(6,*)' OLD: ',MltOrd_old,iRMax_old,DeltaR_old,
     &                            iGrdTyp_old,lMorok_old,DoDirect_Old
            Write(6,*)' NEW: ',MltOrd,iRMax,DeltaR,iGrdTyp,lMorok,
     &                            DoDirect
            Write(6,'(A)')' Check these values or erase ESPF.DATA'
         End If
      Else
         If (PotFile(1:(Index(PotFile,' ')-1)).eq.'***') Then
            Write(6,*)'Error! The EXTE data are missing'
            Call Quit_OnUserError()
         EndIf
      End If
*
* Some output
*
      If (iPL.ge.2) Then
         If (DoDirect) Then
            Write(6,'(A)') ' DIRECT keyword found',
     &                  ' The ESPF scheme is switched off'
            Write(6,'(A,I5,A)')' External potential due to',
     &                                     nChg,' point charges'
         Else
            If (nChg.eq.0) Then
               Write(6,'(A)')' External potential:'
            Else if (nChg.gt.0) Then
               Write(6,'(A,I5,A)')' External potential due to',
     &                                     nChg,' point charges:'
            End If
            If (nChg.ge.0)
     &   Write(6,'(A)')' Atom     E         Fx        Fy        Fz   '//
     &      '     Gxx       Gyy       Gzz       Gxy       Gxz       Gyz'
         End If
      End If
*
* Write the external potential in ESPF.EXTPOT for later use
*
      IPotFl = IsFreeUnit(IPotFl)
      Call Molcas_Open(IPotFl,'ESPF.EXTPOT')
      Write(IPotFl,'(I1)') 0
      Do iAt = 1, natom
         If (.not.DoDirect .and. nChg.ge.0 .and. iPL.ge.2)
     &      Write(6,ExtPotFormat) iAt,
     &         (Work(ipExt+(iAt-1)*MxExtPotComp+j),j=0,MxExtPotComp-1)
         Write(IPotFl,ExtPotFormat) iAt,
     &         (Work(ipExt+(iAt-1)*MxExtPotComp+j),j=0,MxExtPotComp-1)
      End Do
      Close(IPotFl)
      Write (6,*)
      Call qExit('ReadIn')
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
#ifndef _GROMACS_
* Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(ipGradCl)
         Call Unused_real(EnergyCl)
      End IF
#endif
      End
