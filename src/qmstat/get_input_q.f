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
*
*-- Process input to QMSTAT. All input varibles are stored in
*   qminp.fh which in turn are initialized in qmstat_init.
*
      Subroutine Get_Qmstat_Input(iQ_Atoms)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "warnings.fh"

      Character*180 Key
      Character*20 Kword
      Character*180 Get_Ln
      Character VecsQue*3
      Dimension CoTEMP1(3),CoTEMP2(3),CoTEMP3(3),CoTEMP4(3),CoTEMP5(3)
      Dimension SlFacTemp(6)
      External Get_Ln,iClast
      Logical YesNo(20),Changed

      Call QEnter('Get_Input')
*
*-- Say what is done and set all YesNo to false; their purpose is to
*   keep track on compulsory keywords and certain keyword combinations.
*
*      Write(6,*)
*      Write(6,*)'Input processed...'
      Do 1, i=1,20
        YesNo(i)=.false.
1     Continue

*
*-- Use some nice routines to collect input.
*
      LuRd=79
      LuRd=IsFreeUnit(LuRd)
      Call SpoolInp(LuRd)
      Rewind(LuRd)
      Call RdNlst(LuRd,'QMSTAT')

*
*-- The turning-point in this do-while loop.
*
1000  Continue

*
*-- Use Get_Ln to read the lines; it takes care of commented lines.
*
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)

*
*-- The keywords and their labels.
*
      If(Kword(1:4).eq.'TITL') Go To 101
      If(Kword(1:4).eq.'SIMU') Go To 103
      If(Kword(1:4).eq.'THRE') Go To 104
      If(Kword(1:4).eq.'STEP') Go To 105
      If(Kword(1:4).eq.'RUN ') Go To 106
      If(Kword(1:4).eq.'PRIN') Go To 107
      If(Kword(1:4).eq.'EXTE') Go To 108
      If(Kword(1:4).eq.'EDIT') Go To 109
      If(Kword(1:4).eq.'CONF') Go To 110
      If(Kword(1:4).eq.'QMSU') Go To 111
      If(Kword(1:4).eq.'SOLV') Go To 112
      If(Kword(1:4).eq.'RASS') Go To 113
      If(Kword(1:4).eq.'SCFS') Go To 114
      If(Kword(1:4).eq.'SING') Go To 115
      If(Kword(1:4).eq.'ANAL') Go To 116
      If(Kword(1:4).eq.'EXTR') Go To 117
      If(Kword(1:4).eq.'END ') Go To 99999

*
*-- This code is only reached if an illegal keyword in the
*   first tier is encountered.
*
**
      iChrct=Len(Kword)
      Last=iCLast(Kword,iChrct)
      Write(6,*)
      Write(6,*)'ERROR!'
      Write(6,'(1X,A,A)')Kword(1:Last),' is not a valid keyword!'
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<TITLe>>>   Read title.
*
101   Continue
      Key=Get_Ln(LuRd)
      Joblab=Trim(Key)
      ATitle=.true.
      Go To 1000

*
* <<<SIMUlation parameters>>>   Read a variety of simulation
*                               parameters.
*
103   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<RADIe>>>
      If(Kword(1:4).eq.'RADI') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Rstart)
        Go To 103
*---<<<PERMitivity>>>
      Elseif(Kword(1:4).eq.'PERM') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Diel)
        Go To 103
*---<<<TEMPerature>>>
      Elseif(Kword(1:4).eq.'TEMP') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Temp)
        Go To 103
*---<<<PRESsure>>>
      Elseif(Kword(1:4).eq.'PRES') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Pres)
        Go To 103
*---<<<SURFace>>>
      Elseif(Kword(1:4).eq.'SURF') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Surf)
        Go To 103
*---<<<TRANslation>>>
      Elseif(Kword(1:4).eq.'TRAN') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,DelX)
        Go To 103
*---<<<ROTAtion>>>
      Elseif(Kword(1:4).eq.'ROTA') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,DelFi)
        Go To 103
*---<<<CAVIty>>>
      Elseif(Kword(1:4).eq.'CAVI') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,DelR)
        Go To 103
*---<<<FORCe>>>
      Elseif(Kword(1:4).eq.'FORC') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Forcek)
        Go To 103
*---<<<BREPulsion>>>
      Elseif(Kword(1:4).eq.'BREP') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,dLJRep)
        Go To 103
*---<<<SEED>>>
      Elseif(Kword(1:4).eq.'SEED') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,iSeed)
        Go To 103
*---<<<PARAlleltemp>>>
      Elseif(Kword(1:4).eq.'PARA') then
        ParallelT=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nTemp)
        Key=Get_Ln(LuRd)
        Call Get_I(1,nStFilT,nTemp)
        Key=Get_Ln(LuRd)
        Call Get_F(1,ParaTemps,nTemp)
        Go To 103
*---<<<END simulation parameters>>>
      Elseif(Kword(1:4).eq.'END ') then
        Go To 1000
      Endif
*---Here we come if something gets wrong above
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the SIMUlation parameter sect'
     &//'ion:',Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<THREshold>>>    Get the polarization thresholds.
*
104   Continue
      Key=Get_Ln(LuRd)
      Call Get_F1(1,Pollim)
      Call Get_F1(2,Enelim)
      Call Get_I1(3,itMax)
      Go To 1000

*
* <<<STEPs>>>   Specify how many macro- and microsteps.
*
105   Continue
      Key=Get_Ln(LuRd)
      Call Get_I1(1,NMacro)
      Call Get_I1(2,NMicro)
      Go To 1000

*
* <<<RUN >>>   What type of simulation are we to run?
*
106   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
      If(Kword(1:4).eq.'ANAL') Anal=.true.
      If(Kword(1:4).eq.'QMEQ') Qmeq=.true.
      If(Kword(1:4).eq.'QMPR') then
        QmProd=.true.
        Read(LuRd,*)Inter,iNrUt
        Write(SaFilUt(6:6),'(i1.1)')iNrUt
        iLuSaUt=32+iNrUt
      Endif
      If(Kword(1:2).eq.'SM') then
        Write(6,*)
        Write(6,*)'No classical simulations are available.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif
      YesNo(8)=.true.
      Go To 1000

*
* <<<PRINt>>>   Specify print-level.
*
107   Continue
      Key=Get_Ln(LuRd)
      Call Get_I1(1,iPrint)
      Go To 1000

*
* <<<EXTErnal>>>   External one-electron perturbation
*                  should be added on the hamiltonian.
*
108   Continue
      AddExt=.true.
      Key=Get_Ln(LuRd)
      Call Get_I1(1,nExtAddOns)
      If(nExtAddOns.gt.MxExtAddOn) then
        Write(6,*)
        Write(6,*)'Too many external perturbations asked for.'
        Call Quit(_RC_INPUT_ERROR_)
      Endif
      Do 10801, i=1,nExtAddOns
        Read(LuRd,*)ScalExt(i),ExtLabel(i),iCompExt(i)
10801 Continue
      Go To 1000

*
* <<<EDITstartfile>>> Section for editing and displaying stuff on
*                     given startfile.
*
109   Continue
      EdSt=.true.
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*----<<<DELEte>>>   Delete solvent molecules.
      If(Kword(1:4).eq.'DELE') then
        DelOrAdd(1)=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NrStarti)
        Call Get_I1(2,NrStartu)
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nDel)
        Go To 109
      Endif
*----<<<ADD >>>    Add solvent molecules.
      If(Kword(1:4).eq.'ADD ') then
        DelOrAdd(2)=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NrStarti)
        Call Get_I1(2,NrStartu)
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nAdd)
        Go To 109
      Endif
*----<<<QMDElete>>>   Substitute all slots with non-water coordinates
*                     with water coordinates.
      If(Kword(1:4).eq.'QMDE') then
        DelOrAdd(3)=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NrStarti)
        Call Get_I1(2,NrStartu)
        Go To 109
      Endif
*----<<<DUMP coordinates>>>  Dump coordinates in a way suitable for
*                            graphical display.
      If(Kword(1:4).eq.'DUMP') then
        DelOrAdd(4)=.true.
        Key=Get_Ln(LuRd)
        Call UpCase(Key)
        cDumpForm=Key(1:4)
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NrStarti)
        Go To 109
      Endif
*----<<<END editstartfile>>>
      If(Kword(1:4).eq.'END ') then
        Write(StFilIn(6:6),'(i1.1)')NrStarti
        Write(StFilUt(6:6),'(i1.1)')NrStartu
        Go To 1000
      Endif
*---Here we come if something gets wrong above
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the EDITstartfile sect'
     &//'ion:',Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<CONFiguration>>>    Where is the initial configuration to be
*                        obtained. Also, if we wish to edit the
*                        startfile.
*
110   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<ADD >>>  How many to add at random.
      If(Kword(1:4).eq.'ADD ') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,iExtra)
        If(iExtra.gt.MxPut) then
          Write(6,*)
          Write(6,*)'The present limit of explicit solvent molecules i'
     &//'s',MxPut,'.'
          Call Quit(_RC_INPUT_ERROR_)
        Endif
        Go To 110
*---<<<FILE>>>  Read configuration/s from a file and put them there.
      Elseif(Kword(1:4).eq.'FILE') then
        Key=Get_Ln(LuRd)
        Kword=Trim(Key)
        Call UpCase(Kword)
*------<<<STARtfile>>>  Read from startfile.
        If(Kword(1:4).eq.'STAR') then
          Key=Get_Ln(LuRd)
          Kword=Trim(Key)
          Call UpCase(Kword)
*---------<<<SCRAtch>>> Just put QM as given on RUNFILE.
          If(Kword(1:4).eq.'SCRA') then
            Key=Get_Ln(LuRd)
            Call Get_I1(1,iNrIn)
            Call Get_I1(2,iNrUt)
            iRead=8
*---------<<<COPY>>> Collect place of QM from startfile. WARNING!
*                    You must use consistent startfile and RUNFILE!
          Elseif(Kword(1:4).eq.'COPY') then
            Key=Get_Ln(LuRd)
            Call Get_I1(1,iNrIn)
            Call Get_I1(2,iNrUt)
            iRead=7
*---------<<<CM  >>> Put QM in CM of QM-place on startfile.
          Elseif(Kword(1:2).eq.'CM') then
            Key=Get_Ln(LuRd)
            Call Get_I1(1,iNrIn)
            Call Get_I1(2,iNrUt)
            iRead=6
          Else
            Write(6,*)
            Write(6,*)'Illegal StartFile option.'
            Call Quit(_RC_INPUT_ERROR_)
          Endif
*------<<<SAMPfile>>>  Read configurations from sampfile and collect
*                      the extracted information in iNrExtr.
        Elseif(Kword(1:4).eq.'SAMP') then
          Key=Get_Ln(LuRd)
          Call Get_I1(1,iNrIn)
          Call Get_I1(2,iNrExtr)
          iRead=9
          Write(SimEx(6:6),'(i1.1)')iNrExtr
          YesNo(9)=.true.
*-----CRASH-BOOM-BANG!
        Else
          Write(6,*)
          Write(6,*)' Error in CONFiguration section, FILE subsection.'
          Call Quit(_RC_INPUT_ERROR_)
        Endif
        Write(StFilIn(6:6),'(i1.1)')iNrIn
        Write(StFilUt(6:6),'(i1.1)')iNrUt
        Write(SaFilIn(6:6),'(i1.1)')iNrIn
        iLuStIn=8+iNrIn
        iLuStUt=16+iNrUt
        iLuSaIn=24+iNrIn
        Go To 110
*---<<<INPUt>>>  Signify that the first configuration will be given
*                explicitly in input. The coordinates are then given
*                in the solvent section.
      Elseif(Kword(1:4).eq.'INPU') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,iNrUt)
        iLuStUt=16+iNrUt
        YesNo(5)=.true.  !User will give coord:s in input.
        Go To 110
      Elseif(Kword(1:4).eq.'END ') then
        YesNo(2)=.true.  !Signify that source of starting conf. is
                       !specified.
        Go To 1000
      Endif
*---Here we come if something gets wrong above
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the CONFiguration section:'
     &,Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<QMSUrrounding>>>  Give parameters for the QM-Stat.Mech.
*                      interaction.
*
111   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<DPARameters>>>  Dispersion
      If(Kword(1:4).eq.'DPAR') then
        Do 1111, i=1,iQ_Atoms
          Key=Get_Ln(LuRd)
          Call Get_F(1,Udisp(i,1),1)
          Call Get_F(2,Udisp(i,2),1)
1111    Continue
        Go To 111
*---<<<ELECtrostatic>>> Electrostatic Slater Numbers
      Elseif(Kword(1:4).eq.'ELEC') then
1131    Continue
        Key=Get_Ln(LuRd)
        Kword=Trim(Key)
        Call UpCase(Kword)
*------<<<THREsholds>>>  First is the Cutoff (distance Quantum Site-
*                        Classical molecule) to evaluate Penetration
*                        effects. Second, difference between two Slater
*                        exponents to not be consider the same value.
        If(Kword(1:4).eq.'THRE') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Cut_Elc)
          Call Get_F1(2,DifSlExp)
          Go To 1131
        Endif
*-----<<<NOPEnetration>>>  Electrostatic Penetration Not Computed
        If(Kword(1:4).eq.'NOPE') then
          lSlater=.false.
          Go To 1131
        Endif
*-----<<<QUADrupoles>>>  Electrostatic Penetration Computed in quadrupoles.
        If(Kword(1:4).eq.'QUAD') then
          lQuad=.true.
          Go To 1131
        Endif
*------<<<END Electrostatic>>>
        If(Kword(1:4).eq.'END ') then
          Go To 111
        Endif
        Write(6,*)
        Write(6,*)' Error in QMSUrrounding section, ELECtrostatic sub'
     &//'section.'
        Call Quit(_RC_INPUT_ERROR_)
*---<<<XPARameters>>>  Exchange repulsion
      Elseif(Kword(1:4).eq.'XPAR') then
1112    Continue
        Key=Get_Ln(LuRd)
        Kword=Trim(Key)
        Call UpCase(Kword)
*------<<<S2  >>>  The S2 parameter
        If(Kword(1:4).eq.'S2  ') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Exrep2)
          Go To 1112
        Endif
*------<<<S4  >>>  The S4 parameter
        If(Kword(1:4).eq.'S4  ') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Exrep4)
          Go To 1112
        Endif
*------<<<S6  >>>  The S6 parameter
        If(Kword(1:4).eq.'S6  ') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Exrep6)
          Go To 1112
        Endif
*------<<<S10 >>>  The S10 parameter
        If(Kword(1:4).eq.'S10 ') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Exrep10)
          Go To 1112
        Endif
*------<<<CUTOff>>> The cut-off radii for repulsion. The first is
*                   outer radius that says EX=0 if R.gt.Cut_Ex1, while
*                   the second is a EX=infinity if R.lt.Cut_Ex2.
        If(Kword(1:4).eq.'CUTO') then
          Key=Get_Ln(LuRd)
          Call Get_F1(1,Cut_Ex1)
          Call Get_F1(2,Cut_Ex2)
          Go To 1112
        Endif
*------<<<END xparameters>>>
        If(Kword(1:4).eq.'END ') then
          Go To 111
        Endif
        Write(6,*)
        Write(6,*)' Error in QMSUrrounding section, XPARameters subsec'
     &//'tion.'
        Call Quit(_RC_INPUT_ERROR_)
*---<<<DAMPing>>>  Damping parameters for Qm-Surrounding interaction.
      Elseif(Kword(1:4).eq.'DAMP') then
1119    Continue
        Key=Get_Ln(LuRd)
        Kword=Trim(Key)
        Call UpCase(Kword)
*------<<<DISPersion>>>  Dispersion damping parameters. This part
*                        should have a AUTO keyword which collects
*                        default parameters from the MpProp file.
        If(Kword(1:4).eq.'DISP') then
          Dispdamp=.true.
          Key=Get_Ln(LuRd)
*-- Damping numbers for solvent
          Call Get_F(1,CharDi(1),1)
          Call Get_F(2,QuaDi(1,1),1)
          Call Get_F(3,QuaDi(2,1),1)
          Call Get_F(4,QuaDi(3,1),1)
          Key=Get_Ln(LuRd)
          Call Get_F(1,CharDi(2),1)
          Call Get_F(2,QuaDi(1,2),1)
          Call Get_F(3,QuaDi(2,2),1)
          Call Get_F(4,QuaDi(3,2),1)
*-- Damping numbers for solute
          Do 11711, i=1,iQ_Atoms
            Key=Get_Ln(LuRd)
            Call Get_F(1,CharDiQ(i),1)
            Call Get_F(2,QuaDiQ(1,i),1)
            Call Get_F(3,QuaDiQ(2,i),1)
            Call Get_F(4,QuaDiQ(3,i),1)
11711     Continue
          Go To 1119
*------<<<FIELd>>> Parameters for damping electric field.
        Elseif(Kword(1:4).eq.'FIEL') then
          Fielddamp=.true.
          Key=Get_Ln(LuRd)
          Call Get_F1(1,CAFieldG)
          Call Get_F1(2,CBFieldG)
          Call Get_F1(3,CFExp)
          Go To 1119
*------<<<END damping>>>
        Elseif(Kword(1:4).eq.'END ') then
          Go To 111
        Endif
        Write(6,*)
        Write(6,*)' Error in QMSUrrounding section, DAMPing subsection.'
        Call Quit(_RC_INPUT_ERROR_)
*---<<<END qmsurrounding>>>
      Elseif(Kword(1:4).eq.'END ') then
        Go To 1000
      Endif
*---And here we only go if unrecognized keyword is encountered.
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the QMSUrrounding section:'
     &,Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<SOLVent>>>  Specify stuff about the solvent. Usually, these
*                parameters should not be altered.
*
112   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<EXCHange>>>  Exchange repulsion parameters to solvent-solvent.
      If(Kword(1:4).eq.'EXCH') then
        Do 1081, i=1,nAtom
          Do 1082, j=1,i
            Key=Get_Ln(LuRd)
            Call Get_F(1,Sexrep(i,j),1)
            Call Get_F(2,Sexre1(i,j),1)
            Call Get_F(3,Sexre2(i,j),1)
            Sexrep(j,i)=Sexrep(i,j)
            Sexre1(j,i)=Sexre1(i,j)
            Sexre2(j,i)=Sexre2(i,j)
1082      Continue
1081    Continue
        Go To 112
      Endif
*---<<<DISPersion>>>  Dispersion parameters to solvent-solvent.
      If(Kword(1:4).eq.'DISP') then
        Do 1091, i=1,nPol
          Do 1092, j=1,i
            Key=Get_Ln(LuRd)
            Call Get_F(1,Disp(i,j),1)
            Disp(j,i)=Disp(i,j)
1092      Continue
1091    Continue
        Go To 112
      Endif
*---<<<COORdinates>>>  Explicitly given coordinates of solvent
*                      molecules. Also need number of particles.
      If(Kword(1:4).eq.'COOR') then
        YesNo(6)=.true.  !Signify that user gives coordinates.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nPart)
        kaunt=0
        Do 1101, i=1,nPart
          Do 1102, j=1,nAtom
            kaunt=kaunt+1
            Key=Get_Ln(LuRd)
            Call Get_F(1,Cordst(kaunt,1),1)
            Call Get_F(2,Cordst(kaunt,2),1)
            Call Get_F(3,Cordst(kaunt,3),1)
1102      Continue
          Do 1103, kk=1,3
            CoTEMP1(kk)=Cordst(kaunt-2,kk)
            CoTEMP2(kk)=Cordst(kaunt-1,kk)
            CoTEMP3(kk)=Cordst(kaunt-0,kk)
1103      Continue
          Call OffAtom(CoTEMP1,CoTEMP2,CoTEMP3,CoTEMP4,CoTEMP5)
          kaunt=kaunt+1
          Do 1104, kk=1,3
            Cordst(kaunt,kk)=CoTEMP4(kk)
1104      Continue
          kaunt=kaunt+1
          Do 1105, kk=1,3
            Cordst(kaunt,kk)=CoTEMP5(kk)
1105      Continue
1101    Continue
        Go To 112
      Endif
*---<<<CAVRepulsion>>>  Repulsion parameters between solvent and
*                       cavity boundary.
      If(Kword(1:4).eq.'CAVR') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,Exdtal)
        Call Get_F1(2,Exdt1)
        Go To 112
      Endif
*---<<<OCORbitals>>>  Occupied Orbitals for the solvent molecule
      If(Kword(1:4).eq.'OCOR') then
        Key=Get_Ln(LuRd)
        Call Get_I(1,iOrb(2),1)
        Go To 112
      Endif
*---<<<ATCEchpol>>>  Number of atoms, centers, charges and
*                    polarizabilities.
* Jose               Slater Sites
      If(Kword(1:4).eq.'ATCE') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nAtom)
        Call Get_I1(2,nCent)
        Call Get_I1(3,nCha)
        Call Get_I1(4,nPol)
        Call Get_I1(5,nSlSiteC)
        Go To 112
      Endif
*---<<<CHARge>>>  Magnitude of the charges.
      If(Kword(1:4).eq.'CHAR') then
        Key=Get_Ln(LuRd)
        Call Get_F(1,Qsta,nCha)
        Go To 112
      Endif
*---<<<POLArizability>>>  Magnitude of polarizabilities.
      If(Kword(1:4).eq.'POLA') then
        Key=Get_Ln(LuRd)
        Call Get_F(1,Pol,nPol)
        Go To 112
      Endif
*Jose+++++++++++++
*---<<<SLATer>>> Magnitude of Slater PreFactors and Exponents.
      If(Kword(1:4).eq.'SLAT') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,lMltSlC)
        If(lMltSlC.gt.1) then
          Write(6,*)
          Write(6,*)'Too high order of multipole in classical system'
          Write(6,*)'              Higher order is 1'
          Call Quit(_RC_INPUT_ERROR_)
        Endif
        Do 2221, i=1,nSlSiteC
          Do 2223, j=0,lMltSlC
            nS=j*(j+1)*(j+2)/6
            nT=(j+1)*(j+2)*(j+3)/6
            Key=Get_Ln(LuRd)
            Call Get_F1(1,SlExpTemp)
            SlExpC(j+1,i)=SlExpTemp
            njhr=nT-nS
            Key=Get_Ln(LuRd)
            Call Get_F(1,SlFacTemp,njhr)
            njhr=1
            Do 2224, k=nS+1,nT
              SlFactC(k,i)=SlFacTemp(njhr)
              njhr=njhr+1
2224        Continue
2223      Continue
          Key=Get_Ln(LuRd)
          Call Get_F(1,SlPC(i),1)
2221    Continue
        Go To 112
      Endif

*---<<<END solvent>>>
      If(Kword(1:4).eq.'END ') then
        Go To 1000
      Endif
*---And the bla bla bla if something gets wrong.
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the SOLVent section:'
     &,Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<RASSisection>>>  Give some numbers specific for the handling
*                     of the RASSI-construction of the wave-func.
*
113   Continue
      QmType='RASSI'
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<JOBFiles>>>  How many jobfiles and how many states in them.
      If(Kword(1:4).eq.'JOBF') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NrFiles)
        Key=Get_Ln(LuRd)
        Call Get_I(1,NrStates,NrFiles)
        Go To 113
      Endif
*---<<<EQSTate>>> Which state is to be equilibrated.
      If(Kword(1:4).eq.'EQST') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nEqState)
        Go To 113
      Endif
*---<<<MOREduce>>> Work in reduced MO-basis.
      If(Kword(1:4).eq.'MORE') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,ThrsRedOcc)
        MoAveRed=.true.
        Go To 113
      Endif
*---<<<CONTract>>> Contract the RASSI state basis.
      If(Kword(1:4).eq.'CONT') then
        Key=Get_Ln(LuRd)
        Call Get_F1(1,ThrsCont)
        ContrStateB=.true.
        Go To 113
      Endif
*---<<<LEVElshift>>> Introduce levelshift of RASSI states.
      If(Kword(1:4).eq.'LEVE') then
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nLvlShift)
        Key=Get_Ln(LuRd)
        Call Get_I(1,iLvlShift,nLvlShift)
        Key=Get_Ln(LuRd)
        Call Get_F(1,dLvlShift,nLvlShift)
*----- Just a little sorting.
7485    Continue
        Changed=.false.
        Do 7484, i=1,nLvlShift-1
          If(iLvlShift(i).gt.iLvlShift(i+1)) then
            iTemp=iLvlShift(i)
            iLvlShift(i)=iLvlShift(i+1)
            iLvlShift(i+1)=iTemp
            dTemp=dLvlShift(i)
            dLvlShift(i)=dLvlShift(i+1)
            dLvlShift(i+1)=dTemp
            Changed=.true.
          Endif
7484    Continue
        If(Changed) GoTo 7485
        Go To 113
      Endif
*---<<<CISElect>>> Use overlap criterion in choosing state.
      If(Kword(1:4).eq.'CISE') then
        lCiSelect=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,nCIRef)
        Key=Get_Ln(LuRd)
        Call Get_I(1,iCIInd,nCIRef)
        Key=Get_Ln(LuRd)
        Call Get_F(1,dCIRef,nCIRef)
        GoTo 113
      Endif
*---<<<END rassisection>>>
      If(Kword(1:4).eq.'END ') then
        YesNo(3)=.true.  !Rassi section has been visited.
        Go To 1000
      Endif
*---HOW COULD IT GET WRONG HERE?
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the RASSisection section:'
     &,Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<SCFSection>>>   Numbers for a SCF-QmStat run.
*
114   Continue
      QmType='SCF'
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<ORBItals>>>  Specifiy the reduced orbital space in which the
*                   problem is solved.
      If(Kword(1:4).eq.'ORBI') then
        Key=Get_Ln(LuRd)
        Call Get_I(1,iOrb(1),1)
        Call Get_I1(2,iOcc1)
        If(iOrb(1).gt.MxOrb) then
          Write(6,*)
          Write(6,*)'The parameter MxOrb is set too low, or your '
     &//'total number of orbitals too high.'
          Call Quit(_RC_INPUT_ERROR_)
        Endif
        Go To 114
      Endif
*---<<<MP2Denscorr>>>
      If(Kword(1:4).eq.'MP2D') then
        Mp2DensCorr=.true.
        Go To 114
      Endif
*---<<<END scfsection>>>
      If(Kword(1:4).eq.'END ') then
        YesNo(4)=.true.  !Scf section has been visited.
        Go To 1000
      Endif
*---ETWAS FALSCH!
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the SCFSection:',Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
* <<<SINGle-point>>>   Signify that a set of single point calculations
*                      are to be done.
*
115   Continue
      SingPoint=.true.
      YesNo(7)=.true.
      Go To 1000

*
* <<<ANALyze section>>> Give details what analysis of the sampfile
*                        coordinates that is to be done.
*
116   Continue
      Go To 1000

*
* <<<EXTRact section>>> Give details what QM and QM/MM analysis
*                        that is to be done from the sampfile
*                        coordinates.
*
117   Continue
      Key=Get_Ln(LuRd)
      Kword=Trim(Key)
      Call UpCase(Kword)
*---<<<TOTAl energy>>>
      If(Kword(1:4).eq.'TOTA') then
        lExtr(1)=.true.
        Go To 117
      Endif
*---<<<DIPOle>>>
      If(Kword(1:4).eq.'DIPO') then
        lExtr(2)=.true.
        Go To 117
      Endif
*---<<<QUADrupole>>>
      If(Kword(1:4).eq.'QUAD') then
        lExtr(3)=.true.
        Go To 117
      Endif
*---<<<EIGEn things>>>
      If(Kword(1:4).eq.'EIGE') then
        lExtr(4)=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,iExtr_Eig)
        Call Get_S(2,VecsQue,1)
        Call UpCase(VecsQue)
        If(VecsQue(1:3).eq.'YES') lExtr(5)=.true.
        Go To 117
      Endif
*---<<<EXPEctation values>>>
      If(Kword(1:4).eq.'EXPE') then
        lExtr(6)=.true.
        Go To 117
      Endif
*---<<<ELOCal>>>
      If(Kword(1:4).eq.'ELOC') then
        lExtr(7)=.true.
        Key=Get_Ln(LuRd)
        Call Get_I1(1,NExtr_Atm)
        Key=Get_Ln(LuRd)
        Call Get_I(1,iExtr_Atm,NExtr_Atm)
        Go To 117
      Endif
******JoseMEP****************
*---<<<MESP>>>
* The Main Electrostatic potential, field and field gradients will
* be obtained in order to produce perturbation integrals that will
* be used to optimize the intramolecular geometry of the QM system.
*
      If(Kword(1:4).eq.'MESP') then
        lExtr(8)=.true.
        Go To 117
      Endif
******************************
*---<<<END extract section>>>
      If(Kword(1:4).eq.'END ') then
        YesNo(10)=.true.
        GoTo 1000
      Endif
*---ETWAS FALSCH!
      Write(6,*)
      Write(6,*)' Unrecognized keyword in the EXTRact section:'
     &                                              ,Kword(1:4)
      Call Quit(_RC_INPUT_ERROR_)

*
*-- Exit
*
99999 Continue

*
*-- Check if mandatory input was included and that no blatant
*   inconcistencies exist. Not fool-proof, fools!
*
      Call MandatoryInp(YesNo)

*
*-- Good bye.
*
      Call Qexit('Get_Input')

      Return
      End
