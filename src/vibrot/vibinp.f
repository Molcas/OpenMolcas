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
* Copyright (C) 1991, Bjorn O. Roos                                    *
************************************************************************
      Subroutine Vibinp(ncase,ngrid,nvib,Umin,Umax,
     *                  Rout,PotR,E0,dE0,Redm,Teas,Req,scale,temp)
************************************************************************
*  Object: Read input and construct default input parameters           *
*  Called from main program                                            *
*  subroutine calls: Isotope (tables of isotope masses)                *
*  written by Bjoern Roos in March 1991                                *
************************************************************************
*     Use Isotopes
      Implicit real*8 (A-H,O-Z)
#include "dimensions.fh"
#include "intinp.fh"
#include "observ.fh"

#include "SysDef.fh"
#include "constants.fh"
*     Character*80 line,Title1(10),Title2(10)
      Character*80 Title1(10),Title2(10)
C For storing character data using gather/scatter DAFILE operations, it
C is imperative that the strings are aligned on integers.
C Assume an 8-byte character string is always dividible into integers.
C Else, IntCh and (maybe) Title1 and Title2 must be changed.
      Character*8 IntCh
      Parameter (ntab=20)
      Character*4 word,tabinp(ntab),Diatom,Diatomx
      Character*2 At1x,At2x
      Character*180 l84,l84x
      Logical exist
      Dimension Rout(npoint+4),potR(npoint+4)
      Dimension Rin(npin),Ein(npin)
*
      Character*180 Line, Get_Ln,Get_Ln_EOF
      External Get_Ln
C
      Data tabinp/'TITL','ATOM','GRID','RANG','VIBR',
     *            'ROTA','ORBI','NOSP','OBSE','STEP',
     *            'POTE','ROVI','TRAN','ASYM','PRWF',
     *            'SCAL','TEMP','ALLR','DUMM','END '/
*
      LuIn=5
C
C     Set default values to input variables
C
      Atom1='  '
      Atom2='  '
      ipot=0
C     ******        Indicator for potential input
      ngrid=199
C     ********      Maximum number of grid points
      Rmin=1.D0
      Rmax=5.D0
C     *********     Integration range
      Umin=log(Rmin)
      Umax=log(Rmax)
      nvib=3
      nvib1=nvib-1
C     ******        Number of vibrational states
      J1A=0
      J2A=5
C     *****         Range of rotational quantum numbers
      lambda=0
C     ********      Orbital angular momentum quantum number
      dE0=0.0004d0
C     **********    Start value for step size in eigenvalue search
      ispc=1
C     ******        Spectroscopic constants will be computed
      iobs=0
C     ******        No calculation of matrix elements for operators
      Teas=0.D0
C     *********     Asymtotic energy difference between two potentials
      iplotp=0
C     ********      No plot file for the potential
      iscale=0
C     ********      No scaling of input potential (BOR in June 2001)
      temp=300.0d0
C     ************  Temperature for Boltzmann weighting
      iallrot=0
C     *********     No calculation of transition between all rot. levels
CPAM97: New key word PRWF, flagged by IFPRWF .gt.0, default=0:
      IFPRWF=0
      do i=1,10
       Title1(i)=' '
      enddo
C
C
C     Position input file
C
      Call RdNLst(5,'VibRot')
C
C     Input format statements
C
1000  FORMAT(A)
C
C     Read input data from unit 5
C
      ntit1=0
 200  Read(5,'(A)') line
      Call Upcase(line)
      If(line(1:1).eq.'*') go to 200
      word=line(1:4)
      If (Word.eq.'   ') Word='END'
210   Continue
      Do 220 i=1,ntab
       ii=i
       if(word.ne.tabinp(i)) go to 220
       If (ii.eq.1) Then
         Goto 1
       Else If (ii.eq.2) Then
         Goto 2
       Else If (ii.eq.3) Then
         Goto 3
       Else If (ii.eq.4) Then
         Goto 4
       Else If (ii.eq.5) Then
         Goto 5
       Else If (ii.eq.6) Then
         Goto 6
       Else If (ii.eq.7) Then
         Goto 7
       Else If (ii.eq.8) Then
         Goto 8
       Else If (ii.eq.9) Then
         Goto 9
       Else If (ii.eq.10) Then
         Goto 10
       Else If (ii.eq.11) Then
         Goto 11
       Else If (ii.eq.12) Then
         Goto 12
       Else If (ii.eq.13) Then
         Goto 13
       Else If (ii.eq.14) Then
         Goto 14
       Else If (ii.eq.15) Then
         Goto 15
       Else If (ii.eq.16) Then
         Goto 16
       Else If (ii.eq.17) Then
         Goto 17
       Else If (ii.eq.18) Then
         Goto 18
       Else If (ii.eq.19) Then
         Goto 19
       Else If (ii.eq.20) Then
         Goto 20
       End If
220   Continue
      Write(6,*)
      Write(6,*)'******************************************'
      Write(6,*)' VIBINP Error: Input line not recognized. '
      Write(6,*)' Input line, in upper case:               '
      Write(6,'(a)') line
      Write(6,*)' Extracted keyword: ',word
      Write(6,*)'******************************************'
      Call Quit_OnUserError()
C
C Read title lines. Maximum 10 is allowed.
C
1     Read(5,'(a)') line
      word=line(1:4)
      Call Upcase(word)
      Do i=1,ntab
        If(word.eq.tabinp(i)) go to 210
      End Do
      If(ntit1.lt.10) then
        ntit1=ntit1+1
        Read(line,1000) Title1(ntit1)
      end if
      go to 1
C
C     Read isotope numbers
C     isn1 and isn2 are the isotope numbers for atoms 1 and 2,
C     Atom1 and Atom2 the corresponding chemical symbols
C     Masses are read from the next one/two lines if one/two
C     of the isotope numbers are negative.
C     Else the masses are obtained from the tables in ISOTOPE.
C     if isn=0 the mass for the most abundant isotope will be used
C
*2    Read(5,'(a)') line
*     If(line(1:1).eq.'*') go to 2
 2    Line=Get_Ln(LuIn)
      Call Upcase(line)
      ii=0
      ist=0
      Do 23 i=1,80
       If(line(i:i).eq.' ') go to 23
       If(i.lt.ist) go to 23
       ii=ii+1
       kk=0
       Do 21 k=i,80
        If(line(k:k).eq.' ') go to 22
        kk=kk+1
21     Continue
22     Continue
       ist=i+kk
       If (ii.eq.1) Then
*         Read(line(i:i+kk-1),*) isn1
          Call Get_I(1,isn1,1)
       ElseIf (ii.eq.2) Then
*         Read(line(i:i+kk-1),'(A)') Atom1(1:kk)
          Call Get_S(2,Atom1(1:kk),1)
       ElseIf (ii.eq.3) Then
*         Read(line(i:i+kk-1),*) isn2
          Call Get_I(3,isn2,1)
       ElseIf (ii.eq.4) Then
*         Read(line(i:i+kk-1),'(A)') Atom2(1:kk)
          Call Get_S(4,Atom2(1:kk),1)
          Go To 24
       End If
23    Continue
24    Continue
*     If(isn1.lt.0) Read(5,*) xMass1
*     If(isn2.lt.0) Read(5,*) xMass2
*     If(isn1.ge.0) Call Isotope(isn1,Atom1,xMass1)
*     If(isn2.ge.0) Call Isotope(isn2,Atom2,xMass2)
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      Conv=CONV_AMU_TO_AU_
      If(Atom1.eq.'D') Then
         ich1=1
      Else If(Atom1.eq.'T') Then
         ich1=1
      Else
         ich1=iNuclearChargeFromSymbol(Atom1)
      End If
      If(isn1.eq.0) Then
         if(Atom1.eq.'D') Then
            isn1=2
         Else If(Atom1.eq.'T') Then
            isn1=3
         Else
            isn1=iMostAbundantIsotope(ich1)
         End If
         xMass1=dNuclearMass(ich1,isn1)
      Else If(isn1.gt.0) Then
         xMass1=dNuclearMass(ich1,isn1)
      Else
         Read(5,*) xMass1
         xMass1=Conv*xMass1
      End If
*
      If(Atom2.eq.'D') Then
         ich2=1
      Else If(Atom2.eq.'T') Then
         ich2=1
      Else
         ich2=iNuclearChargeFromSymbol(Atom2)
      End If
      If(isn2.eq.0) Then
         if(Atom2.eq.'D') Then
            isn2=2
         Else If(Atom2.eq.'T') Then
            isn2=3
         Else
            isn2=iMostAbundantIsotope(ich2)
         End If
         xMass2=dNuclearMass(ich2,isn2)
      Else If(isn2.gt.0) Then
         xMass2=dNuclearMass(ich2,isn2)
      Else
         Read(5,*) xMass2
         xMass2=Conv*xMass2
      End If
*
*
      Redm=xMass1*xMass2/(xMass1+xMass2)
      go to 200
C
C     Read number of grid points for numerical integration (max: npoint)
C
*3    Read(5,'(a)') line
*     Call Upcase(line)
*     If(line(1:1).eq.'*') go to 3
*     Read(line,*) ngrid
 3    Line=Get_Ln(LuIn)
      Call Get_I(1,ngrid,1)
      If((ngrid/2)*2.eq.ngrid) ngrid=ngrid-1
C     ************************************** ngrid should be odd
      If(ngrid.ge.npoint) ngrid=npoint-1
      go to 200
C
C     Read upper and lower integration range in atomic units
C
*4    Read(5,'(a)') line
*     Call Upcase(line)
*     If(line(1:1).eq.'*') go to 4
*     Read(line,*) Rmin,Rmax
 4    Line=Get_Ln(LuIn)
      Call Get_F(1,Rmin,1)
      Call Get_F(2,RMax,1)
      Umin=log(Rmin)
      Umax=log(Rmax)
      go to 200
C
C     Read number of vibrational quantum numbers
C
*5    Read(5,'(a)') line
*     Call Upcase(line)
*     If(line(1:1).eq.'*') go to 5
*     Read(line,*) nvib
 5    Line=Get_Ln(LuIn)
      Call Get_I(1,nvib,1)
      n0=0
      nvib1=nvib-1
      go to 200
C
C     Read range for rotational quantum numbers
C
*6    Read(5,'(a)') line
*     Call Upcase(line)
*     If(line(1:1).eq.'*') go to 6
*     Read(line,*) J1A,J2A
 6    Line=Get_Ln(LuIn)
      Call Get_I(1,J1A,1)
      Call Get_I(2,J2A,1)
      if(J2A.ge.nRot_Max) Then
        Write(6,*)
        Write(6,*)'********************************'
        Write(6,*)' VIBINP Error: J2A.gt.nRot_Max. '
        Write(6,'(1x,a,2i6)')'J2A:',J2A
        Write(6,*)'********************************'
        Call Quit_OnUserError()
      end if

      go to 200
C
C     Read orbital angular momentum quantum number
C
*7    Read(5,'(a)') line
*     Call Upcase(line)
*     If(line(1:1).eq.'*') go to 7
*     Read(line,*) lambda
 7    Line=Get_Ln(LuIn)
      Call Get_I(1,lambda,1)
      go to 200
C
C     Read flag for spectroscopic constants
C
8     ispc=0
      go to 200
C
C     Read input for calculation of matrix elements of observables
C     like the dipole operator, etc.
C
9     Continue
      iobs=iobs+1
      iplot(iobs)=0
      If(iobs.gt.10) Then
        Write(6,*)
        Write(6,*)'***************************'
        Write(6,*)' VIBINP Error: IOBS.gt.10. '
        Write(6,'(1x,a,2i6)')'IOBS:',IOBS
        Write(6,*)'***************************'
        Call Quit_OnUserError()
      end if
cVV: make it really nasty. if it looks like a file - read from the file
c    else-  this is a title

      Read(LuIn,'(a)') Line
      Call f_inquire(Line(1:index(Line,' ')-1), exist)
      if(exist) then
        LuIn1=isFreeUnit(15)
        Call molcas_open(LuIn1,Line(1:index(Line,' ')-1))
        Line=Get_Ln(LuIn1)
      else
        LuIn1=LuIn
      endif
        Titobs(iobs)=Trim(Line)
*      write(6,*)' In VIBINP. IOBS=',IOBS
*      write(6,*)' TITOBS just read:'
*      write(6,'(a80)') TITOBS(IOBS)
      nobsi=0
 91   Line=Get_Ln_EOF(LuIn1)
      if(Line(1:3).eq.'EOF'.and.LuIn1.ne.LuIn) then
        close(LuIn1)
        Line=Get_Ln(LuIn)
      endif

      word=line(1:4)
      Call Upcase(word)
      Do 92 i=1,ntab
       If(word.eq.tabinp(i)) go to 210
92    Continue
      If(word.eq.'PLOT') go to 93
      nobsi=nobsi+1
      If(NOBSI.gt.NPIN) Then
        Write(6,*)
        Write(6,*)'**********************************'
        Write(6,*)' VIBINP Error: NOBSI.gt.NPIN      '
        Write(6,'(1x,a,2i6)')'NOBSI,NPIN:',NOBSI,NPIN
        Write(6,*)'**********************************'
        Call Quit_OnUserError()
      end if
      npobs(iobs)=nobsi
*     Read(line,*) RinO(nobsi,iobs),Obsin(nobsi,iobs)
      Call Get_F(1,RinO(nobsi,iobs),1)
      Call Get_F(2,Obsin(nobsi,iobs),1)
      go to 91
*93   Read(5,'(a)') line
 93   Line=Get_Ln(LuIn)
      If(line(1:1).eq.'*') go to  93
*     Read(line,*) R0o(iobs),R1o(iobs),dRo(iobs)
      Call Get_F(1,R0o(iobs),1)
      Call Get_F(2,R1o(iobs),1)
      Call Get_F(3,dRo(iobs),1)
      iplot(iobs)=1
      go to 200
C
C     Read starting value for step size in eigenvalue search
C
 10   Line=Get_Ln(LuIn)
      Call Get_F(1,dE0,1)
      go to 200
C
C     Read potential
C
11    Continue
      ipot=1
      nop=0
      Line=Get_Ln(LuIn)
      xxx=999.97
      Read(line,*,err=666) xxx
      If(xxx.eq.999.97) go to 666
      backspace(LuIn)
      Luin1=Luin
      go to 111
 666  call f_inquire(line,exist)
      write(6,*) line
      If(.not.exist) then
       write(6,*) 'File with potential non-existent'
       write(6,*) 'File=',line(:mylen(line))
       call Abend
      Endif
      LuIn1=isFreeUnit(15)
      Call molcas_open(LuIn1,line)
 111  Line=Get_Ln_EOF(LuIn1)
      if(Line(1:3).eq.'EOF'.and.LuIn1.ne.LuIn) then
        close(LuIn1)
        Line=Get_Ln(LuIn)
      endif
      word=line(1:4)
      Call Upcase(word)
      Do 112 i=1,ntab
112   If(word.eq.tabinp(i)) go to 210
      If(word.eq.'PLOT') go to 113
      nop=nop+1
      If(NOP.gt.NPIN) Then
        Write(6,*)
        Write(6,*)'**********************************'
        Write(6,*)' VIBINP Error: NOP.gt.NPIN'
        Write(6,'(1x,a,2i6)')'NOP,NPIN:',NOP,NPIN
        Write(6,*)'**********************************'
        Call Quit_OnUserError()
      end if
*     Read(line,*) Rin(nop),Ein(nop)
      Call Get_F(1,Rin(nop),1)
      Call Get_F(2,Ein(nop),1)
      go to 111
c      if(Rin(1).gt.Rmax.or.Rin(nop).lt.Rmin) go to 900
*113  Read(5,'(a)') line
*     If(line(1:1).eq.'*') go to 113
*     Read(line,*) R0p,R1p,dRp
 113  Line=Get_Ln(LuIn)
      Call Get_F(1,R0p,1)
      Call Get_F(2,R1p,1)
      Call Get_F(3,dRp,1)
      iplotp=1
      if(LuIn1.ne.LuIn) Close(Unit=LuIn1)
C
C
C     Calculation of ro-vibrational wave functions (ncase=1)
C
12    Continue
      ncase=1
      go to 200
C
C     Calculation of transition moments (ncase=2)
C
13    Continue
      ncase=2
      go to 200
C
C     Asymptotic energy difference between two potentials
C
14    Continue
*141  Read(5,'(a)') line
*     If(line(1:1).eq.'*') go to 141
*     Read(line,*) Teas
*141  Continue
      Line=Get_Ln(LuIn)
      Call Get_F(1,Teas,1)
      go to 200
C
15    Continue
C Flag for printing the wave function.
      IFPRWF=1
      go to 200
C
16    Continue
C     Scaling of input potential such that the binding energy is 0.1 au.
      iscale=1
      go to 200
C
17    Continue
C     TEMPerature
C     Temperature for vibrational averaging
      Line=Get_Ln(LuIn)
      Call Get_F(1,Temp,1)
      go to 200
C
18    Continue
C     ALLRotational

      iallrot=1
      go to 200
C
19    Continue
      go to 200
C
20    Continue
      if(J1A.lt.lambda) Then
        Write(6,*)
        Write(6,*)'********************************'
        Write(6,*)' VIBINP Warning: J1A.lt.Lambda. '
        Write(6,'(1x,a,2i6)')'J1A,Lambda:',J1A,Lambda
        Write(6,*)' J1A is now reset=Lambda.       '
        Write(6,*)'********************************'
      end if

      if(J2A.lt.J1A) Then
        Write(6,*)
        Write(6,*)'***************************'
        Write(6,*)' VIBINP Error: J2A.lt.J1A. '
        Write(6,'(1x,a,2i6)')'J2A,J1A:',J2A,J1A
        Write(6,*)'***************************'
        Call Quit_OnUserError()
      end if
C
C      Check for input error
C
       If(ncase.eq.1) Then
        If(Atom1.eq.'  '.or.Atom2.eq.'  ') Then
         Write(6,*)
         Write(6,*)'**********************************'
         Write(6,*)' VIBINP Error: No atoms in input.'
         Write(6,*)'**********************************'
         Call Quit_OnUserError()
        End If
       End If
C
C     Retrieve data from Vibwvs1 and Vibwvs2 (ncase=2)
C
      If(ncase.eq.2) then
       iadvi1=0
       Call iDafile(Vibwvs1,2,iad12,100,iadvi1)
       iOpt=2
       Call WR_VibRot_Info1(Vibwvs1,iOpt,iadvi1,
     &                      ntit1,J1A,J2A,lambda,n0,
     &                      nvib1,Redm,Umax,Umin,ngrid,
     &                      isn1,isn2,Req,xMass1,xMass2)
       Call cDaFile(Vibwvs1,iOpt,Title1,10*80,iadvi1)
       Call cDaFile(Vibwvs1,iOpt,IntCh,     8,iadvi1)
*
       Atom1=IntCh(1:2)
       Atom2=IntCh(5:6)
       Rmin=exp(Umin)
       Rmax=exp(Umax)
       iadvi2=0
       Call iDafile(Vibwvs2,2,iad13,100,iadvi2)
       iOpt=2
       Call WR_VibRot_Info1(Vibwvs2,iOpt,iadvi2,
     &                      ntit2,J1B,J2B,lambdx,n02,
     &                      nvib21,Redmx,Umaxx,Uminx,ngridx,
     &                      isn1x,isn2x,Reqx,xMass1,xMass2)
       Call cDaFile(Vibwvs2,iOpt,Title2,10*80,iadvi2)
       Call cDaFile(Vibwvs2,iOpt,IntCh,     8,iadvi2)
       At1x=IntCh(1:2)
       At2x=IntCh(5:6)
#if 0
       if(iadvi2.eq.-666) then
cvv
cvv   This if statement is a dummy!
cvv   to fake f90 under HP-UX64
cvv
       write(6,*)' isn', isn2
       write(6,*)' Mass', xMass1
       write(6,*)' Mass', xMass2
       write(6,*)' Title', Title2
       write(6,*)' IntCh', IntCh
       endif
#endif

C
C      Check for consistency of data on files
C
       IERR=0
       If(Atom1.ne.At1x.or.Atom2.ne.At2x) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: Inconsistent data.'
         Write(6,'(1x,a,1x,a,1x,a)')'ATOM1,AT1X:',ATOM1,AT1X
         IERR=1
       End If
       If(abs(Redm-Redmx).gt.1.d-06) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: REDM.ne.REDMX'
         Write(6,'(1x,a,2f16.6)')'REDM,REDMX:',REDM,REDMX
         IERR=1
       End If
       If(Umax.ne.Umaxx) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: UMAX.ne.UMAXX'
         Write(6,'(1x,a,2f16.6)')'UMAX,UMAXX:',UMAX,UMAXX
         IERR=1
       End If
       If(Umin.ne.Uminx) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: UMIN.ne.UMINX'
         Write(6,'(1x,a,2f16.6)')'UMIN,UMINX:',UMIN,UMINX
         IERR=1
       End If
       If(NGrid.ne.NGridx) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: NGrid.ne.NGridX'
         Write(6,'(1x,a,2i8)')'NGrid,NGridX:',NGrid,NGridX
         IERR=1
       End If
       If(Isn1.ne.Isn1x) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: Isn1.ne.Isn1X'
         Write(6,'(1x,a,2i8)')'Isn1,Isn1X:',Isn1,Isn1X
         IERR=1
       End If
       If(Isn2.ne.Isn2x) Then
         Write(6,*)
         Write(6,*)'***************************************'
         Write(6,*)' VIBINP Error: Isn2.ne.Isn2X'
         Write(6,'(1x,a,2i8)')'Isn2,Isn2X:',Isn2,Isn2X
         IERR=1
       End If
       If(IERR.NE.0) Then
         Write(6,*)'*****************************************'
         Write(6,*)' VIBINP: Irrecoverable errors.'
         Write(6,*)' Transition calculation, but data on the '
         Write(6,*)' two VibWvs files are not compatible.'
         Write(6,*)'*****************************************'
         Call Quit_OnUserError()
       End If
      Endif
C
C
C     Section for print output of input information
C
C
      ii=0
      Diatom=' '
      Diatomx(1:2)=Atom1
      Diatomx(3:4)=Atom2
      Do 110 i=1,4
       If(Diatomx(i:i).ne.' ') then
        ii=ii+1
        Diatom(ii:ii)=Diatomx(i:i)
       Endif
110   Continue
      If(ncase.eq.1) Write(6,1100) Diatom
1100  Format(/1x,'Vibration-Rotation spectrum for the ',A4,' molecule.')
      If(ncase.eq.2) Write(6,1101) Diatom
1101  Format(/1x,'Transition moments for the ',A4,' molecule.')
C
      If(ncase.eq.2) Write(6,1190)
1190  Format(/1x,'State number 1')
      Do i=1,ntit1
        Write(6,*) Title1(i)
      End Do
      Write(6,1200) J1A,J2A,lambda,n0,nvib1,
     *              isn1,Atom1,xMass1,isn2,Atom2,xMass2,Redm
1200  Format(/1x,'Rotational quantum number range ',2I3
     *       /1x,'Electronic angular momentum     ',I3
     *       /1x,'Vibrational quantum number range',2I3
     *       /1x,'Mass of atom ',I3,1x,A2,D15.6,' au'
     *       /1x,'Mass of atom ',I3,1x,A2,D15.6,' au'
     *       /1x,'Reduced mass       ',D15.6,' au')
      If(ncase.eq.2) then
       Write(6,1191)
1191   Format(/1x,'State number 2')
        Do i=1,ntit2
          Write(6,*) Title2(i)
        End Do
       Write(6,1201) J1B,J2B,lambdx,n02,nvib21,
     *              isn1,At1x,xMass1,isn2,At2x,xMass2,Redmx
1201   Format(/1x,'Rotational quantum number range ',2I3
     *        /1x,'Electronic angular momentum     ',I3
     *        /1x,'Vibrational quantum number range',2I3
     *        /1x,'Mass of atom ',I3,1x,A2,D15.6,' au'
     *        /1x,'Mass of atom ',I3,1x,A2,D15.6,' au'
     *        /1x,'Reduced mass       ',D15.6,' au')
      Endif
C
C     Numerical integration data
C
      del=(Umax-Umin)/(ngrid-1)
      Write(6,1300) ngrid,del,Rmin,Rmax,Umin,Umax
1300  Format(/1x,'Statistics for numerical integration'
     *       /1x,'Number of steps         ',I5
     *       /1x,'Step length             ',F10.6
     *       /1x,'Radial integration range',2F10.6,' au'
     *       /1x,'logarithmic range       ',2F10.6)
C
      If(ispc.eq.0.and.ncase.eq.1) Write(6,1400)
1400  Format(/1x,'Spectroscopic constants will not be computed')
      If(ispc.ne.0.and.ncase.eq.1) Write(6,1500)
1500  Format(/1x,'Spectroscopic constants will be computed')
      If(iobs.eq.0) Write(6,1600)
1600  Format(1x,'Matrix elements of operators will not be computed')
      If(iobs.ne.0) Write(6,1700) iobs
1700  Format(1x,'Matrix elements of',I3,' operators will be computed')
      If(ipot.eq.0.and.ncase.eq.1) Then
        Write(6,*)
        Write(6,*)'**********************************'
        Write(6,*)' VIBINP Error: IPOT=0 and NCASE=1 '
        Write(6,*)'**********************************'
        Call Quit_OnUserError()
      End If
C
C     Print input potential
C
      If(ipot.ne.0) then
       Write(6,1800)
1800   Format(/1x,'Potential read from input'/6x,' R(au)     E(au)')
       Do 30 i=1,nop
        Write(6,1810) Rin(i),Ein(i)
1810    Format(1x,2F12.6)
30     Continue
       if(Rin(1).gt.Rmax.or.Rin(nop).lt.Rmin) then
         Write(6,*)
         Write(6,*)'****************************'
         Write(6,*)' VIBINP Error: Potential is '
         Write(6,*)' not within defined range.  '
         Write(6,*)'****************************'
         Call Quit_OnUserError()
       endif
       If(iplotp.ne.0) write(6,1820)
1820   Format(1x,'A plot file of the potential will be generated')
      Endif
C
C     Compute radial coordinates
C
      U=Umin-del
      Rout0=exp(U)
      Do 300 i=1,ngrid
         U=U+del
         Rout(i)=exp(U)
300   Continue
      Rn1=exp(U+del)
      Rout(ngrid+1)=Rout0
      Rout(ngrid+2)=Rn1
      If(ipot.ne.0)
     *Call POT(Rin,Ein,Rout,PotR,ngrid+2,1,E0,Req,R0p,R1p,dRp,
     *          nop,Title1(1),iplotp,Redm,scale,0)
      iadvib=0
      If(ncase.eq.1) then
C Store data on Vibwvs (ncase=1)
         Call iDafile(Vibwvs,1,iad12,100,iadvib)
         iOpt=1
         Redm1=Redm*scale
         Call WR_VibRot_Info1(Vibwvs,iOpt,iadvib,
     &                        ntit1,J1A,J2A,lambda,n0,
     &                        nvib1,Redm1,Umax,Umin,ngrid,
     &                        Abs(isn1),Abs(isn2),Req,xMass1,xMass2)
         IntCh(1:4)=Atom1
         IntCh(5:8)=Atom2
         Call cDaFile(Vibwvs,iOpt,Title1,10*80,iadvib)
         Call cDaFile(Vibwvs,iOpt,IntCh,     8,iadvib)
      Endif

C Fit observable input
      Rout0=Rout(ngrid+1)
      Do 32 i=1,iobs
       Rout(ngrid+1)=Req
       Call POT(RinO(1,i),Obsin(1,i),Rout,EoutO(1,i),ngrid+1,2,O0,
     * Oeq,R0o(i),R1o(i),dRo(i),npobs(i),Titobs(i),iplot(i),Redm,scale,
     * iobs)
32    Continue
      Rout(ngrid+1)=Rout0
C
C     Print observable input data
C
      If(iobs.ne.0) then
       Write(6,1900)
1900   Format(/1x,'Input data for observables')
       Do 40 i=1,iobs,4
        Write(6,1910) (Titobs(k)(1:18),k=i,min(i+3,iobs))
1910    Format(/1x,4A24)
        l84=' '
        Do 35 k=i,min(i+3,iobs)
         l84(24*(k-i)+1:24*(k-i)+10)='     R(au)'
         l84(24*(k-i)+11:24*(k-i+1))='        Value '
35      Continue
        Write(6,'(a)') l84(:mylen(l84))
        Do 37 j=1,npin
         l84=' '
         Do 36 k=i,min(i+3,iobs)
          If(j.le.npobs(k))
     *    Write(l84(24*(k-i)+1:24*(k-i+1)),'(1x,F9.4,F14.6)')
     *          RinO(j,k),Obsin(j,k)
36       Continue
         If(l84.ne.' ') Write(6,'(a)') l84(:mylen(l84))
37      Continue
        l84=' '
        l84x=' '
        Do 38 k=i,min(i+3,iobs)
         If(iplot(i).ne.0)
     *   l84(24*(k-i)+1:24*(k-i)+14)='     Plot file'
         If(iplot(i).eq.0)
     *   l84(24*(k-i)+1:24*(k-i)+17)='     No plot file'
         Write(l84x(24*(k-i)+1:24*(k-i+1)),'(F10.4,F14.6)')
     *          Req,EoutO(ngrid+1,k)
38      Continue
         Write(6,'(A)') l84(:mylen(l84))
         Write(6,*)' Interpolated value at equilibrium for state 1:'
         Write(6,'(A)') l84x(:mylen(l84x))
40     Continue
      Endif
      Return
      End
