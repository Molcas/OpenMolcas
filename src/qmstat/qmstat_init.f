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
*-- Subroutine with purpose to initialize and set defaults for the
*   input section.
*
      Subroutine Qmstat_Init
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"

*IO_stuff
      StFilIn='STFIL0'
      SAFilIn='SAFIL0'
      StFilUt='STFIL0'
      SAFilUt='SAFIL0'
      BlockIn='BLOCIN'
      BlockUt='BLOCUT'
      SimEx='EXTRA0'
*--Jose: File for the optimization procedure
      FieldNuc='AVENUC'
*--------
      Do 101, i=1,MxJobs
        Write(JbName(i),'(A,i3.3)')'JOB',i
101   Continue
      RassiM='RASSIM'
      GammaO='GAMORB'
      EigV='EIGV'
      AddOns(1)='ADDON1'
      AddOns(2)='ADDON2'
      AddOns(3)='ADDON3'
      iNrIn=-1
      iNrUt=0
      iLuStIn=8+iNrIn
      iLuStUt=16+iNrUt
      iLuSaIn=24+iNrIn
      iLuSaUt=32+iNrUt
      iLuBlockIn=3
      iLuBlockUt=4
      iRead=0
*Defaults
      nEqState=1
      Cut_Ex1=10.0d0
      Cut_Ex2=0.0d0
      DelX = 0.00d0
      DelFi = 0.0d0
      DelR = 0.00d0
      Temp = 300.0d0
      ISEED = 791204
      IPrint = 1
      NMACRO = 1
      NMICRO = 1
      RSTART = 80.0d0
      NPART = 0
      nAtom=3
      NCENT = 5
      NPOL = 3
      NCHA = 4
*Jose.Slater Penetration
      nSlSiteC=5
      lMltSlC=0
*****
      nLvlShift=0
      Do 119, i=1,MxAt
        iExtr_Atm(i)=-1
119   Continue
      QSTA(1) = 0.5836d0
      QSTA(2) = 0.5836d0
      QSTA(3) =-0.5836d0
      QSTA(4) =-0.5836d0
      POL(1) = 5.932d0
      POL(2) = 0.641d0
      POL(3) = 0.641d0
*Jose.Slater Penetration
      Cut_Elc=6.0d0
      DifSlExp=0.001d0

      SlFactC(1,1)=-0.50d0
      SlFactC(1,2)=-0.4164d0
      SlFactC(1,3)=-0.4164d0
      SlFactC(1,4)=-0.5836d0
      SlFactC(1,5)=-0.5836d0
      Do 217, i=1,5
        Do 218, j=2,4
          SlFactC(j,i)=0.0d0
218     Continue
217   Continue
      SlExpC(1,1)=2.5552d0
      SlExpC(1,2)=2.6085d0
      SlExpC(1,3)=2.6085d0
      SlExpC(1,4)=2.5552d0
      SlExpC(1,5)=2.5552d0
      Do 219, i=1,5
        SlExpC(2,i)=0.00d0
219   Continue
      SlPC(1)=0.5d0
      SlPC(2)=1.0d0
      SlPC(3)=1.0d0
      SlPC(4)=0.0d0
      SlPC(5)=0.0d0
*******************************
      sExRep(1,1) = 2.092338000000000000d0
      sExRe1(1,1) = 158.998000000000000d0
      sExRe2(1,1) = 4.660090000000000d10
      sExRep(2,1) = 2.112447000000000000d0
      sExRe1(2,1) = 8.31922000000000000d0
      sExRe2(2,1) = 97560.62000000000000d0
      sExRep(2,2) = 1.075803000000000000d0
      sExRe1(2,2) = 0.06521000000000000d0
      sExRe2(2,2) = 1121941276d0
      sExRep(3,1) = 2.112447000000000000d0
      sExRe1(3,1) = 8.31922000000000000d0
      sExRe2(3,1) = 97560.62000000000000d0
      sExRep(3,2) = 1.075803000000000000d0
      sExRe1(3,2) = 0.06521000000000000d0
      sExRe2(3,2) = 1121941276d0
      sExRep(3,3) = 1.075803000000000000d0
      sExRe1(3,3) = 0.06521000000000000d0
      sExRe2(3,3) = 1121941276d0
      Disp(1,1) = 11.3380000000000000d0
      Disp(2,1) = 3.38283000000000000d0
      Disp(2,2) = 0.627068000000000000d0
      Disp(3,1) = 3.38283000000000000d0
      Disp(3,2) = 0.627068000000000000d0
      Disp(3,3) = 0.627068000000000000d0
      DO 10, I=1,NPOL
        DO 20, J=1,I
          DISP(J,I)=DISP(I,J)
          SEXREP(J,I)=SEXREP(I,J)
          SEXRE1(J,I)=SEXRE1(I,J)
          SEXRE2(J,I)=SEXRE2(I,J)
20      Continue
10    Continue
      CORDST(1,1) = 0.0d0
      CORDST(2,1) = 0.0d0
      CORDST(3,1) = 0.0d0
      CORDST(4,1) = 0.3126d0
      CORDST(5,1) = -0.3126d0
      CORDST(1,2) = 0.0d0
      CORDST(2,2) = 1.43d0
      CORDST(3,2) = -1.43d0
      CORDST(4,2) = 0.0d0
      CORDST(5,2) = 0.0d0
      CORDST(1,3) = 0.3d0
      CORDST(2,3) = -0.807d0
      CORDST(3,3) = -0.807d0
      CORDST(4,3) = -0.1191d0
      CORDST(5,3) = -0.1191d0
      ForceK = 0.001d0
      dLJrep=0.0d0
      Pres = 1.0d0
      PolLim = 0.0001d0
      EneLim = 0.0000001d0
      itMax=30
      Exdtal = 30.0d0
      Exdt1 = 0.060d0
      Surf = 30.0d0
      iOrb(2)=5
      Diel = 80.0d0
      iExtra = 0
      Smeq=.false.
      Qmeq=.false.
      Fielddamp=.false.
      Dispdamp=.false.
      Smprod=.false.
      QmProd=.false.
      ChargedQM=.false.
      ATitle=.false.
      Anal=.false.
      ParallelT=.false.
      Mp2DensCorr=.false.
      MoAveRed=.false.
      lCiSelect=.false.
      EdSt=.false.
*JoseMEP***** The dimension was increased from 8 to 12
      Do 41, i=1,12
        DelOrAdd(i)=.false.
        lExtr(i)=.false.
        lAnal(i)=.false.
41    Continue
      lSlater=.true.
      lQuad=.false.


      Return
      End
