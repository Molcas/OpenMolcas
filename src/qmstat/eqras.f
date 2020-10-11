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
      Subroutine EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C
     &                ,iBigForDeAll,nSizeBig,ip_UNCLE_MOE,nB)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "integral.fh"
#include "WrkSpc.fh"
#include "warnings.fh"
#include "constants.fh"

      Parameter (Conver1=1.0d10*CONST_BOHR_RADIUS_IN_SI_)
      Parameter (Conver2=2.0d0*Pi/360.0d0)
*     Parameter (BoltzK=1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_)
      Dimension nBas(1),nBas_C(1),nCnC_C(MxBasC),Coord(MxAt*3)
      Dimension Eint(MxQCen,10),Poli(MxQCen,10)
      Dimension iFP(3),iGP(3),iDT(3),iFi(3),iFil(MxQCen,10)
      Dimension Smat(MxStOT),SmatPure(MxStOT),Vmat(MxStOT)
      Dimension BoMaH(MxAt),BoMaO(MxAt)
*****Jose** Interaction with Slater type to consider Penetration
      Dimension Eint_Nuc(MxAt)
*****JoseMEP**New variables for the MEP calculation
      Dimension SumElcPot(MxQCen,10)
      Dimension PertElcInt(MxBas*(MxBas+1)/2)
      Character*4 Labjhr
*****
      Character Memlabel*20,Memlaabe*20,Memlaaab*20,MemLaaaa*20
      Character MemQFal*20,ChCo*2,ChCo2*2
      Logical DidWeAccept,Haveri,Exist,CalledBefore,SampleThis,InCutOff
      Character Head*200
      Parameter (ExLim=10) !Over how long distance the exchange rep.
                           !is computed, the solv-solv.
      External Ranf
      Dimension iDum(1)

*
*-- Enter eqras.
*

*
*-- Numbers, initializations, conversions.
*
      BoltzK=1.0d-3*CONST_BOLTZMANN_/CONV_AU_TO_KJ_
      Ract=Rstart              !Initial radie
      delX=delX/Conver1        !angstrom-->Bohr
      delFi=delFi*Conver2      !degree-->radian
      delR=delR/Conver1        !angstrom-->Bohrn
      CalledBefore=.false.
      Samplethis=.true.
      iBigForDeAll=iBigT
      If(ParallelT) nMacro=nTemp*nMacro

*
*---- If we have vacuum, then no volume-pressure work is done, nor do we
*     have any cavitation free-energy.
*
      If(abs(Diel-1).le.0.0001d0) then
        Gamma=0
        Gam=0
      Else
       If(SURF.le.0.0d0) then
        Gamma=0
        Gam=0
       Else
        Gamma=Pres*4.188790205d0*.52917d0**3*6.023d0*.00001d0
     &       *1.01325d0/627.52d0/4.184d0
        Gam=.0005973455d0/74d0*SURF
       Endif
      Endif
*
*--- More numbers.
*
      If(Temp.le.0.0d0) then
        BetaBol=1.0d23
      Else
        BetaBol=1.0d0/(Temp*BoltzK)
      Endif
      DiFac=-(Diel-1.0d0)/(Diel+1.0d0)
      Expran=0.0d0
      Expre=0.0d0
      ggsum=0.0d0
      HighS=0.0d0
      iHowMSampIN=0
      iHowMSampUT=0
      iD=1
      Adisp=Disp(1,2)
      nBaseC=nBas_C(1)
      nBaseQ=nBas(1)
      iTri=(iOrb(1)*(iOrb(1)+1))/2
      iTriBasQ=nBaseQ*(nBaseQ+1)/2
      iTriState=nState*(nState+1)/2
      timeCLAS=0
      timeEX=0
      timeEL=0
      timeMC=0


*
*---- If some damping has been requested, prepare it here and print.
*
      If((Dispdamp.or.FieldDamp).and.iPrint.ge.8) then
        Write(6,*)
        Write(6,*)'-----Various damping data.'
        Write(6,*)
      Endif
      If(Dispdamp) then

*
*---- Construct the Born-Mayer parameters, a la Brdarski-Karlstrom
*
        Call BornMayerBK(iQ_Atoms,BoMaH,BoMaO)
      Endif

*
*---- Damping of field.
*
      If(Fielddamp) then
        if(iPrint.ge.8) then
          Write(6,*)'  Damping the field between Qm-region and solvent.'
          Write(6,*)'E_damp=E_0*(1-exp(alpha*distance))^N'
          Write(6,*)'  alpha(QM-oxygen)   =',CAFieldG
          Write(6,*)'  alpha(QM-hydrogen) =',CBFieldG
          Write(6,*)'  N                  =',CFexp
        endif
      Endif

*
*-- Check what type of simulation to run, and generate some output of
*   outmost beauty.
*
      If(Qmeq.and.iRead.ne.9) then
        Call NiceOutPut('EIQ',Gam,Gamma,BetaBol)
      Elseif(QmProd.and.iRead.ne.9) then
        Call NiceOutPut('PIQ',Gam,Gamma,BetaBol)
        Call DaName(iLuSaUt,SaFilUt)
        iDisk=0 !Put some dummy on the sampfile so we have space for
                !the real number later.
        Call iDaFile(iLuSaUt,1,[iHowMSampUT],1,iDisk)
        !Below we make a check for extreme cases. Our algorithm to
        !select sampling configurations sets this limit.
        iProdMax=(2**30-1)*2
        If((nMacro*nMicro).ge.iProdMax) then
          Write(6,*)
          Write(6,*)'WARNING! Too large numbers for nMacro and nMicro t'
     &//'o run production!'
          Write(6,*)'   Their product must not be greater than 2**31!'
          Write(6,*)'   If you wish to make such large samples, you ca'
     &//'n run several samplings and collect several sampfiles.'
          Call Quit(_RC_INTERNAL_ERROR_)
        Endif
      Elseif(iRead.eq.9) then  !If we read from sampfile: open the
                          !sampfile and read how many configurations
                          !and open extract file.
        Call NiceOutPut('SSS',Gam,Gamma,BetaBol)
        Call DaName(iLuSaIn,SaFilIn)
        iDiskSa=0
        Call iDaFile(iLuSaIn,2,iDum,1,iDiskSa)
        iHowMSampIN=iDum(1)
        iLuExtr=54
        iLuExtr=IsFreeUnit(iLuExtr)
        Call OpnFl(SimEx,iLuExtr,Exist)
        Write(iLuExtr,*)'Extract-File'
        Write(iLuExtr,*)
        !And put some words in the output
        Write(6,*)
        Write(6,*)'   Total number of sampled configurations:'
     &           ,iHowMSampIN
        Write(6,*)'   Reading from the file ',SaFilIn
        Write(6,*)'   Summarizing data put on ',SimEx
******JoseMEP
      ! If we perform MEP calculation, first we make some zeros
      ! and allocate some memory.
        If(lExtr(8)) then
          Do ijhr=1,MxQCen
            Do jjhr=1,10
               SumElcPot(ijhr,jjhr)=0.0d0
               AvElcPot(ijhr,jjhr)=0.0d0
            End do
          End do
          NCountField=0

          iTriMaxBasQ=MxBas*(MxBas+1)/2
          call dcopy_(iTriMaxBasQ,[ZERO],iZERO,PertNElcInt,iONE)
          Call GetMem('SumOvlAOQ','Allo','Real',ipAOSum,iTriBasQ)
          call dcopy_(iTriBasQ,[ZERO],iZERO,Work(ipAOSum),iONE)
        Endif
**********
      Else
        Write(6,*)
        Write(6,*)'An invalid number of iRead detected.'
        Call Quit(_RC_INTERNAL_ERROR_)
      Endif
*----------------------------------------------------------------------*
* If we have input file, then read from it.                            *
*----------------------------------------------------------------------*
      iCStart=(((iQ_Atoms-1)/nAtom)+1)*nCent+1
      iCNum=(iCStart-1)/nCent
      i9=0 !i9 is active if iRead.eq.9 and we are collecting
            !configurations from the sampfile.
58886 Continue
      i9=i9+1
      If(iRead.le.8.and.iRead.ge.6) then
        Call Get8(Ract,Dum)
      Elseif(iRead.eq.9) then
        Call Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)
      Else
        If(iExtra.gt.0) then
          Call NyPart(iExtra,nPart,Cordst,Rstart,nCent,iSeed)
        Endif
        If(iPrint.ge.10) then
          Write(Head,*)'Coordinates of the initial distribution.'
          Call Cooout(Head,Cordst,nPart,nCent)
        Endif
      Endif
*
*-- Give a startvalue for the Total energy. The effect is that we
*   always accept the first microstep.
*
      Etot=1D+10

*
*-- Some numbers.
*
      ncpart=Ncent*nPart
      ncParm=ncPart-(nCent*icNum)
      nClas=nPart-iCNum
      indma=npart*npol
      iCi=(iQ_Atoms*(iQ_Atoms+1))/2

*
*-- Put QM-molecule in its place.
*
      If(iRead.eq.8.or.iRead.eq.0) then
        Call PlaceIt(Coord,iQ_Atoms,iCNum)
      Elseif(iRead.eq.6) then
        Call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
        If(iPrint.ge.10) then
          Write(Head,*)'CM-centred coordinates after substitution.'
          Call Cooout(Head,Cordst,nPart,nCent)
        Endif
      Endif

*----------------------------------------------------------------------*
*                                                                      *
*------------------------- START SIMULATION ---------------------------*
*                                                                      *
*----------------------------------------------------------------------*

      iSnurr=0  !How many steps taken totally.
      nSiffiB=0  !A number for blocking.
*
*---- The Macrosteps.
*
      Do 2000, iMacro=1,nMacro
        Esav=0.0d0
        If(iRead.eq.9) then
          iAcc=1
        Else
          iAcc=0
        Endif

*------ If we are running parallel tempering, then...
        If(ParallelT) Call ParaRoot(Ract,BetaBol,Etot,CalledBefore
     &                             ,SampleThis)
*
*------ The Microsteps.
*
        Do 2001, iMicro=1,nMicro
          Call Timing(Cpu1,Tim1,Tim2,Tim3)
          Eold=Etot
          iSnurr=iSnurr+1
*
*-------- Generate new configuration, both solvent and QM-region.
*
          Call GeoGen(Ract,Rold,iCNum,iQ_Atoms)
*
*-------- Compute Solvent-solvent interaction.
*
          Call ClasClas(iCNum,iCStart,ncParm,Coord,iFP,iGP,iDT,iFI,iDist
     &                 ,iDistIm,Elene,Edisp,Exrep,E2Die,ExDie)
          Call QMPosition(EHam,Cordst,Coord,Forcek,dLJrep,Ract,iQ_Atoms)
          Call Timing(Cpu2,Tim1,Tim2,Tim3)
          timeCLAS=timeCLAS+(Cpu2-Cpu1)
*--------------------------------------------------------------------------*
* Work a bit with the quantum part.                                        *
*--------------------------------------------------------------------------*
          Do 4002, i=1,3
            xyzMyQ(i)=0  !Dipoles for the QM-part, see polink.f.
            xyzMyI(i)=0
4002      Continue
          nSize=3*nPol*nPart
          Do 400, i=1,iCi  !Allocate memory for the field on the QM-mol.
            Do 4000, j=1,10 !iCi: number of quantum molecule sites.
              Write(MemQFal,'(A,i2.2,i2.2)')'Falt',i,j
              Call GetMem(MemQFal,'Allo','Real',iFil(i,j),nSize)
4000        Continue
400       Continue
          Do 401, i=1,iCi
            Do 402, j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
              Do 403, k=1,nPart*nPol !Classical polarisation sites
                                     !including quantum molecule.
                Work(iFil(i,j)-1+k)=0.0d0
                Work(iFil(i,j)-1+k+nPart*nPol)=0.0d0
                Work(iFil(i,j)-1+k+2*nPart*nPol)=0.0d0
403           Continue
              Eint(i,j)=0.0d0
402         Continue
          If (i.le.MxAt) Eint_Nuc(i)=0.0d0
401       Continue
*
*-------- Compute the exchange operator.
*
          Call ExRas(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms
     &              ,nAtomsCC,Ax,Ay,Az,iTriState,Smat,SmatPure
     &              ,InCutOff,ipAOSum)
          Call Timing(Cpu3,Tim1,Tim2,Tim3)
          timeEX=timeEX+(Cpu3-Cpu2)
*
*------ Electrostatics commencing.
*
*
*------ Compute various gradients of 1/r.
*
          If(lSlater) then
            Call OneOverR_Sl(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp
     &                      ,iCNum,Eint,iQ_Atoms,outxyzRAS
     &                      ,Eint_Nuc)
          Else
            Call OneOverR(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp
     &                   ,iCNum,Eint,iQ_Atoms,outxyzRAS)
          Endif

*
*----- Couple the point-charges in the solvent to the QM-region.
*
          Call HelState(Eint,nState,iCi,RasCha,RasDip,RasQua,Vmat
     &                 ,iPrint)
*
*----- Let QM-region and solvent polarize.
*
          Call PolRas(iDist,iDistIM,iDT,iFI,iFP,iFil,iCStart
     &               ,iTriState,VMat,Smat,DiFac,Ract,iCNum,Energy
     &               ,nVarv,iSTC,Haveri,iQ_Atoms,ip_ExpVal,Poli)
*
*----- Energy from QM-nuclei interacting with solvent field.
*
          If(lSlater) then
            Do 702, i=1,iQ_Atoms
             Energy=Energy-Eint_Nuc(i)*ChaNuc(i)
702        Continue
          Else
            Do 703, i=1,iQ_Atoms
              Energy=Energy-Eint(i,1)*ChaNuc(i)
703         Continue
          Endif
*
*----- Some additional boost of short-range repulsion.
*
          Call BoostRep(AddRep,SmatPure,iSTC,nState,InCutOff)
*
*----- Sum-up what we will call QM-region energy.
*
          Energy=Energy-EEdisp+AddRep
*----------------------------------------------------------------------*
* Final induction and reaction field energies.                         *
*----------------------------------------------------------------------*
          Call ReaInd(iGP,iDT,iDistIm,iCNum,IndMa,NcParm,Sum1,s90um)
          Call Timing(Cpu4,Tim1,Tim2,Tim3)
          timeEL=timeEL+(Cpu4-Cpu3)
*----------------------------------------------------------------------*
* Construct the final energy.                                          *
*----------------------------------------------------------------------*
          EnCLAS=Elene+EHam-Edisp+Exrep+E2Die+ExDie
          Etot=EnCLAS-0.5*S90um-Sum1+Gamma*Ract**3+Energy+Gam*Ract**2
          Dele=Etot-Eold
*----------------------------------------------------------------------*
* Printing and various if requsted.                                    *
*----------------------------------------------------------------------*
          If(iPrint.ge.10) then
            If(Haveri) Etot=999999
            Write(6,*)
            Write(6,*)'   ----Microstep',iMicro
            Write(6,*)'         Number of iterations:',nVarv
            Write(6,*)'         Total energy:',Etot
            Write(6,*)'         Of which is'
            Write(6,*)'            Pairwise solvent-solvent interaction'
     &//':',EnCLAS
            Write(6,*)'              Solvent-solvent Electrostatic'
     &,Elene
            Write(6,*)'              Harmonic Spring:',EHam
            Write(6,*)'              Solvent-solvent Dispersion:',-Edisp
            Write(6,*)'              Solvent-solvent Exchange:',Exrep
            Write(6,*)'            Energy of induced dipoles in field f'
     &//'rom explicit solvent:',-Sum1
            Write(6,*)'            Energy of solvent charge distributio'
     &//'n in reaction field:',-0.5*S90um
            Write(6,*)'            Solvent E-interaction with image:'
     &,E2Die
            Write(6,*)'            Solvent Repulsion with boundary:'
     &,ExDie
            Write(6,*)'            Solvent-Solute dispersion:',EEdisp
            Write(6,*)'            Energy of QM-region:',Energy
            Write(6,*)'              Higher order overlap exchange pair'
     &//'-term:',AddRep
            Write(6,*)'            Surface tension term:',Gam*Ract**2
            Write(6,*)'            Volume-pressure term:',Gamma*Ract**3
            Write(6,*)'         Previous accepted energy:',Eold
            Write(6,*)'         Difference:',Dele
            Write(6,*)'         Total dipole in QM-region:(',-xyzMyQ(1),
     &',',-xyzMyQ(2),',',-xyzMyQ(3),')'
            Write(6,*)'         Radie:',Ract
            If(Haveri) then
              Write(6,*)'    WARNING! SOME OF THE NUMBERS ABOVE HAVE N'
     &//'O MEANING SINCE THE POLARIZATION DID NOT CONVERGE!!!'
              GoTo 8194
            Endif
          Endif

*
*-- If we are collecting stuff from a sampfile, now is the time to put
*   data on the extract file. If center-specific expectation values
*   are requested, call Allen.
*
          If(iRead.eq.9) then
            If(lExtr(6)) then
              E_Nuc_Rubbet=0.0d0
              If(lSlater) then
                Do 6347, iAt=1,iQ_Atoms
                  E_Nuc_Rubbet=E_Nuc_Rubbet
     &                  -(Eint_Nuc(iAt)+Poli(iAt,1))*ChaNuc(iAt)
6347            Continue
              Else
                Do 6348, iAt=1,iQ_Atoms
                  E_Nuc_Rubbet=E_Nuc_Rubbet
     &                  -(Eint(iAt,1)+Poli(iAt,1))*ChaNuc(iAt)
6348            Continue
              Endif
            Endif
            If(lExtr(7)) then
              Call AllenGinsberg('RASSI',Eint,Poli,ChaNuc,RasCha,RasDip
     &                          ,RasQua,MxStOT,iSTC,nState,iExtr_Atm
     &                          ,lExtr(4),iExtr_Eig,iQ_Atoms
     &                          ,ip_ExpCento,E_Nuc_Part,lSlater
     &                          ,Eint_Nuc)
            Endif
            Call Extract(iLuExtr,i9,Etot,xyzMyQ,HMatState,iSTC,iDt
     &                  ,nState,HMatSOld,xyzQuQ,ip_ExpVal,ip_ExpCento
     &                  ,E_Nuc_Rubbet,E_Nuc_Part)
****JoseMEP**********
      ! If MEP option. Add electr. potential, field, etc
      ! for all solvent config.
            If(lExtr(8)) then
              Labjhr='Add '
              Call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot
     &                     ,NCountField,PertElcInt,iONE
     &                     ,iONE,[iONE],[iONE],iONE)
              NCountField=NCountField+1
            Endif
*******
            GoTo 9090
          Endif

*
*-- Resume the MC-wrap up.
*
          Dele=Dele*BetaBol
          DidWeAccept=.true.
          If(Dele.lt.0) then
            iAcc=iAcc+1
          Else
            Expe=Exp(-Dele)
            Expran=ranf(iseed)
            If(iPrint.ge.10) then
              Write(6,*)'         Positive energy change!'
              Write(6,*)'         Boltzmann weight:',Expe
              Write(6,*)'         Random number:',ExpRan
            Endif
            iAcc=iAcc+1
            If(Expe.lt.ExpRan) then
              Call Oldge(iAcc,Etot,Eold,Ract,Rold)
              DidWeAccept=.false.
              If(iPrint.ge.10) then
                Write(6,*)'         Not accepted!'
              Endif
            Endif
          Endif
          If(DidWeAccept)Esav=Esav+Etot
*----------------------------------------------------------------------*
* If this is a production run, then put stuff on the sampfile.         *
*----------------------------------------------------------------------*
          If(QmProd.and.iRead.ne.9) then
            If(SampleThis) then
              If(Inter.ne.0) then
                Inte=(iSnurr/Inter)*Inter
                If(Inte.eq.((iMacro-1)*nMicro+iMicro)) then
                  Call Put9(Etot,Ract,iDT,iHowMSampUT,Gamma,Gam,Esav
     &                     ,iDisk)
                Endif
              Endif
            Endif
          Endif

*
*-- Free memory.
*
9090      Continue
          nSize=(nClas*(nClas-1)/2)*(nCent**2)
          nSizeIm=(nClas*nCent)**2
          Call GetMem('DistMat','Free','Real',iDist,nSize)
          Call GetMem('DistMatIm','Free','Real',iDistIm,nSizeIm)
          Do 90001,i=1,3
            Write(ChCo,'(I1.1)')i
            Write(MemLabel,*)'FP'//ChCo
            Write(MemLaabe,*)'GP'//ChCo
            Write(MemLaaab,*)'DT'//ChCo
            Write(MemLaaaa,*)'FI'//ChCo
            Call GetMem(MemLabel,'Free','Real',iFP(i),IndMa)
            Call GetMem(MemLaabe,'Free','Real',iGP(i),IndMa)
            Call GetMem(MemLaaab,'Free','Real',iDT(i),IndMa)
            Call GetMem(MemLaaaa,'Free','Real',iFI(i),IndMa)
90001     Continue
          nSize=3*nPol*nPart
          Do 90002, i=1,iCi
            Do 90003, j=1,10
              Write(ChCo,'(I2.2)')i
              Write(ChCo2,'(I2.2)')j
              Write(MemQFal,*)'Falt'//ChCo//ChCo2
              Call GetMem(MemQFal,'Free','Real',iFil(i,j),nSize)
90003       Continue
90002     Continue
          Call GetMem('Coeff','Free','Real',iSTC,nState**2)
          Call Timing(Cpu5,Tim1,Tim2,Tim3)
          timeMC=timeMC+(Cpu5-Cpu4)

*---------------------------------------------------------------------*
*-- End of Microstep.                                                 *
*                                                                     *
*---------------------------------------------------------------------*
2001    Continue
*---------------------------------------------------------------------*
*-- Have we collected all sampled configurations? If no, go up again. *
*   If yes, then close some files and take a little jump downwards to *
*   the END!!!                                                        *
*---------------------------------------------------------------------*
*Jose***************************************
* This point is also used to perform the Average of the Potential,
* Field and Field gradients to obtain and average Electrostatic
* perturbation. The Non-Electr. perturbation is also obtained here
* it will be added directly to the One-electron file.
*******************************************

        If(i9.lt.iHowMSampIN.and.iRead.eq.9) then
          GoTo 58886
        Elseif(i9.ge.iHowMSampIN.and.iRead.eq.9) then
******JoseMEP**
      ! If MEP option. Obtain the mean Potential, Field
      ! and Field Gradients.
      ! It si also obtained the average of the Non-Electrostatic
      ! perturbation
          If(lExtr(8)) then
              Labjhr='Aver'
              Call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot
     &                     ,NCountField,PertElcInt,iONE
     &                     ,iONE,[iONE],[iONE],iONE)
              AverFact=1.0d0/Dble(NCountField)
              Call DaxPy_(iTriBasQ,AverFact,Work(ipAOSum),iONE
     &                   ,PertNElcInt,iONE)
              Call GetMem('SumOvlAOQ','Free','Real',ipAOSum,iTriBasQ)
          Endif
*********

          Call DaClos(iLuSaIn)
          Close(iLuExtr)
          GoTo 58887
        Endif
*----------------------------------------------------------------------*
* Write to startfile.                                                  *
*----------------------------------------------------------------------*
8194    Continue
        ESav=Esav/Dble(iAcc)
        Call Put8(Ract,Etot,Gamma,Gam,ESav)
        If(Haveri) Call Quit(_RC_NOT_CONVERGED_)
*----------------------------------------------------------------------*
* Print some things here at the end of the macrostep.                  *
*----------------------------------------------------------------------*
        If(.not.ParallelT) then
          jMacro=iMacro
        Else
          jMacro=1+(iMacro-1)/nTemp
        Endif
        Write(6,'(A,i4)')'---Macrostep ',jMacro
        Write(6,'(A,i4)')'     Number of microsteps:',nMicro
        Pr=100.0d0*(Dble(iAcc)/Dble(nMicro))
        Write(6,'(A,i4,A,f5.1,A)')'     Number of acceptances:',iAcc,'('
     &                          ,Pr,'%)'
        Write(6,'(A,f12.4)')'     Radie (a.u.):',Ract
        Write(6,'(A,f16.8)')'     Average Energy (a.u.) in Macrostep:'
     &                     ,Esav
        Write(6,'(A,3(f12.4))')'     Total dipole in QM-region last mic'
     &//'rostep (a.u.):',-xyzMyQ(1),-xyzMyQ(2),-xyzMyQ(3)
        Write(6,*)
        Call xFlush(6)
*--------------------------------------------------------------------------*
* End of Macrostep.                                                        *
*--------------------------------------------------------------------------*
2000  Continue
*----------------------------------------------------------------------*
* Put some things on info-file. Used to make tests.                    *
*----------------------------------------------------------------------*
      Call Add_Info('Total Energy',[Etot],1,6)
      Call Add_Info('Induction of system',[Sum1],1,6)
      Call Add_Info('React. field int.',[s90um],1,6)
      Call Add_Info('Solv-Solu Disp.',[EEdisp],1,6)
      Call Add_Info('QM-region Energy',[Energy],1,6)
      Call Add_Info('QM-region dipole',xyzMyQ,3,5)
      RRRnVarv=dble(nVarv)
      Call Add_Info('Pol.Iterations',[RRRnVarv],1,8)
*----------------------------------------------------------------------*
* Close some external files.                                           *
*----------------------------------------------------------------------*
58887 Continue
      If(QmProd.and.iRead.ne.9) then
        iDisk=0
        Call iDaFile(iLuSaUt,1,[iHowMSampUT],1,iDisk)
        Call DaClos(iLuSaUt)
      Endif
*--------------------------------------------------------------------------*
* The End... be happy!                                                     *
*--------------------------------------------------------------------------*
      If(MoAveRed) then
        nB=nRedMO
        ip_UNCLE_MOE=ipAvRed
      Else
        nB=nBaseQ
      Endif
      nSizeBig=nState*(nState+1)*nB*(nB+1)/4

      Write(6,*)
      Write(6,*)
      Write(6,*)' Time statistics. (hour:minute:second)'
      Write(6,*)' -----------------------------------------------------'
      it1h=int(timeCLAS)/3600
      it1m=int(timeCLAS-it1h*3600)/60
      t1s=timeCLAS-it1h*3600-it1m*60
      Write(6,9)'   Time spent on pair-wise solvent-solvent interaction'
     &//'s: ',it1h,':',it1m,':',t1s
      it2h=int(timeEX)/3600
      it2m=int(timeEX-it2h*3600)/60
      t2s=timeEX-it2h*3600-it2m*60
      Write(6,9)'   Time spent on solvent-solute overlap calculations: '
     &,it2h,':',it2m,':',t2s
      it3h=int(timeEL)/3600
      it3m=int(timeEL-it3h*3600)/60
      t3s=timeEL-it3h*3600-it3m*60
      Write(6,9)'   Time spent on solvent and solute electrostatic inte'
     &//'raction: ',it3h,':',it3m,':',t3s
      it4h=int(timeMC)/3600
      it4m=int(timeMC-it4h*3600)/60
      t4s=timeMC-it4h*3600-it4m*60
      Write(6,9)'   Time spent on the Metropolis-Monte Carlo decision: '
     &,it4h,':',it4m,':',t4s
9     Format(A,I4,A,I3,A,F5.2)

      Return
      End
