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
* Copyright (C) 1992,2000, Roland Lindh                                *
************************************************************************
      SubRoutine InpRct(LuSpool)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*                                                                      *
*             Modified for Langevin polarizabilities, March 2000 (RL)  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "rctfld.fh"
#include "covradt_data.fh"
      Character*180 KWord, Key, Get_Ln
      External Get_Ln
*
      iRout=1
      iPrint=nPrint(iRout)
*     Call qEnter('InpRct')
*                                                                      *
************************************************************************
*                                                                      *
*---- Initiate data for dielectric medium
*
      Eps=One
      EpsInf=One
      Eps_User=-One
      EpsInf_User=Zero
      rds=Zero
      lMax=-1
      lRF = .False.
      lRFCav = .False.
      RF_Basis=.False.
*                                                                      *
************************************************************************
*                                                                      *
*---- Initiate data for PCM
*
      PCM=.False.
      i_sph_inp = 0
*     Default PCM parameters
      Call PCMDef(ISlPar,RSlPar,iPrint)
*                                                                      *
************************************************************************
*                                                                      *
*---- Setting up default values for parameters used in Langevin
*     polarizabilities
*
*     latato: Gitter type (1 or 4)
      latato = 1
      Cordsi(1,latato) = Half
      Cordsi(2,latato) = Half
      Cordsi(3,latato) = Half
*
*     polsi: site polarizability
      polsi = Zero
*
*     dipsi: site dipole moment
      dipsi = Zero
*
*     radlat: maximum extension of the lattic
      radlat = Zero
*
*     scala,scalb,scalc: lengths of the cell dimensions
      scala = Zero ! should be set by input
      scalb = Zero
      scalc = Zero
*
*     scaaa: overall scaling
      scaaa = One

*     rotation of grid
      rotAlpha = Zero
      rotBeta  = Zero
      rotGamma = Zero

*     rsca: scaling of radii
      rsca = One

*     Controls wheather a sparse grid should be used outside
*     a distance distSparce from ALL atoms and XF points
*     nSparse = scala(sparse) / scala(normal)
      LSparse = .False.
      nSparse = 1
      distSparse = Zero

*     Minimum distance for handling dipole-dipole interactions
      dipCutoff = Zero

*     Damping of dipole-dipole interactions
      lDamping = .True.

*     Use a dipole-dipole cutoff of Amber type, i.e. same exclusions
*     as for the static field
      lAmberPol = .False.

*     Controls scaling of contributions from negative entries in
*     the exclusion list, typically 1-4 interactions
      scal14 = One

*     Controls wheather an average of different grids is done
      LGridAverage = .False.

*
*     diedel: controlls deletion of gitter polarizabilities
      diedel = 0.01D0
*
*     tK: Boltzman factor
      tK = 0.001
*
*     clim: convergence threshold for solver
      clim = 1.0D-15
*
*     afac: controlls equation solver (valid values 0.3-0.97)
      afac = Half

*     iterDamp: controls update of dipoles in equation solver
*     newDip = oldDip*iterDamp + estimatedDip*(1-iterDamp)
      dampIter=0.4D0
      lDiprestart=.False.
*
*     nexpo: exponent of the potential with which the QM system kills
*            the gitter polatizabilities
      nexpo = 12
*
*     prefac: scaling of the polarization energies
      prefac = One
*
      lLangevin = .False.
*                                                                      *
************************************************************************
*                                                                      *
      iPrint=5
*                                                                      *
************************************************************************
*                                                                      *
*-----Process the input
*
 998  Read (LuSpool,'(A72)',END=977,ERR=988) Key
      KWord = Key
      Call UpCase(KWord)
      If (KWord(1:1).eq.'*')    Go To 998
      If (KWord.eq.'')       Go To 998
      If (KWord(1:4).eq.'REAC') Go To 900
      If (KWord(1:4).eq.'PRIN') Go To 910
      If (KWord(1:4).eq.'LANG') Go To 911
      If (KWord(1:4).eq.'GITT') Go To 912
      If (KWord(1:4).eq.'POLA') Go To 913
      If (KWord(1:4).eq.'DIPO') Go To 914
      If (KWord(1:4).eq.'OSCA') Go To 915
      If (KWord(1:4).eq.'DIED') Go To 916
      If (KWord(1:4).eq.'MXLX') Go To 917
      If (KWord(1:4).eq.'AFAC') Go To 918
      If (KWord(1:4).eq.'RFBA') Go To 919
      If (KWord(1:4).eq.'PCM-') Go To 920
      If (KWord(1:4).eq.'SOLV') Go To 921
      If (KWord(1:4).eq.'DIEL') Go To 922
      If (KWord(1:4).eq.'COND') Go To 923
      If (KWord(1:4).eq.'AARE') Go To 925
      If (KWord(1:4).eq.'R-MI') Go To 926
      If (KWord(1:4).eq.'PAUL') Go To 927
      If (KWord(1:4).eq.'SPHE') Go To 928
      If (KWord(1:4).eq.'TEMP') Go To 929
      If (KWord(1:4).eq.'RSCA') Go To 930
      If (KWord(1:4).eq.'ROTA') Go To 931
      If (KWord(1:4).eq.'SPAR') Go To 932
      If (KWord(1:4).eq.'IDAM') Go To 933
      If (KWord(1:4).eq.'AVER') Go To 934
      If (KWord(1:4).eq.'CONV') Go To 935
      If (KWord(1:4).eq.'SKIP') Go To 936
      If (KWord(1:4).eq.'NODA') Go To 937
      If (KWord(1:4).eq.'LRAD') Go To 938
      If (KWord(1:4).eq.'AMBE') Go To 939
      If (KWord(1:4).eq.'SC14') Go To 940
      If (KWord(1:4).eq.'DRES') Go To 941
      If (KWord(1:4).eq.'END ') Go To 997
      iChrct=Len(KWord)
      Last=iCLast(KWord,iChrct)
      Write (6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      Call WarningMessage(2,'InpRct: Error in keyword.')
      Call Quit_OnUserError()
*                                                                      *
****** REAC ************************************************************
*                                                                      *
*-----Read reaction field parameters
*
 900  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,Eps)
      Call Get_F1(2,rds)
      Call Get_I1(3,lmax)
      If (nToken(KWord).gt.3) Call Get_F1(4,EpsInf)
*
      lRF = .True.
      lRFCav = .True.
      iRc=-1
      iOpt=0
      Write(KWord,'(A,F10.5,A,F10.5,A,I4)')
     & 'eps=',Eps,' radius=',rds,' higest moment=',lMax
      Go To 998
*                                                                      *
****** PRIN ************************************************************
*                                                                      *
*-----Print level
*
 910  KWord=Get_Ln(LuSpool)
      Call Get_I1(1,n)
      Do i = 1, n
         KWord=Get_Ln(LuSpool)
         Call Get_I1(1,jRout)
         Call Get_I1(2,iPrint)
         nPrint(jRout)=iPrint
      End Do
      Go To 998
*                                                                      *
****** LANG ************************************************************
*                                                                      *
 911  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,scala)
      Call Get_F1(2,scalb)
      Call Get_F1(3,scalc)
      lLangevin = .True.
      lRF = .True.
      Go To 998
*                                                                      *
****** GITT ************************************************************
*                                                                      *
 912  KWord=Get_Ln(LuSpool)
      Call Get_I1(1,latato)
      Do i = 1, latato
         KWord=Get_Ln(LuSpool)
         Call Get_F(1,Cordsi(1,i),3)
      End Do
      Go To 998
*                                                                      *
****** POLA ************************************************************
*                                                                      *
 913  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,polsi)
      Go To 998
*                                                                      *
****** DIPO ************************************************************
*                                                                      *
 914  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,dipsi)
      Go To 998
*                                                                      *
****** OSCA ************************************************************
*                                                                      *
 915  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,scaaa)
      Go To 998
*                                                                      *
****** DIED ************************************************************
*                                                                      *
 916  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,diedel)
      Go To 998
*                                                                      *
****** MXLX ************************************************************
*                                                                      *
 917  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,radlat)
      Go To 998
*                                                                      *
****** AFAC ************************************************************
*                                                                      *
 918  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,afac)
      If (afac.ge.One) Then
         Call WarningMessage(2,'InpRct: afac invalid value!;'//
     &               '        afac >= 1.0 !')
         Call Quit_OnUserError()
      End If
      Go To 998
*                                                                      *
****** RFBA ************************************************************
*                                                                      *
C% 919  RF_Basis=.True.
C%      Go To 998
 919  Go To 998
*                                                                      *
****** PCM- ************************************************************
*                                                                      *
 920  PCM=.True.
      lRF = .True.
      lLangevin = .False.
      Go To 998
*                                                                      *
****** SOLV ************************************************************
*                                                                      *
 921  Solvent=Trim(Get_Ln(LuSpool))
      ISlPar(15) = NumSolv(Solvent)
      Go To 998
*                                                                      *
****** DIEL ************************************************************
*                                                                      *
*
*
 922  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,Eps_User)
      If (nToken(KWord).gt.1) Call Get_F1(2,EpsInf_User)
      Go To 998
*                                                                      *
****** COND ************************************************************
*                                                                      *
 923  Conductor=.True.
      ISlPar(16) = 2
      Go To 998
*                                                                      *
****** AARE ************************************************************
*                                                                      *
 925  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,aArea)
      RSlPar(7) = aArea
      Go To 998
*                                                                      *
****** RMIN ************************************************************
*                                                                      *
 926  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,r_min_Sphere)
      RSlPar(3) = r_min_Sphere
      Go To 998
*                                                                      *
****** PAUL ************************************************************
*                                                                      *
 927  ITypRad = 2
      ISlPar(9) = ITypRad
      Go To 998
*                                                                      *
****** SPHE ************************************************************
*                                                                      *
 928  KWord=Get_Ln(LuSpool)
      Call Get_I1(1,I_Sph)
      Call Get_F1(2,Radius)
      ITypRad = 3
      ISlPar(9) = ITypRad
      i_sph_inp = i_sph_inp + 1
      NOrdInp(i_sph_inp) = I_Sph
      RadInp(i_sph_inp) = Radius
      ISlPar(14) = i_sph_inp
      Go To 998
*                                                                      *
****** TEMP ************************************************************
*                                                                      *
*     Temperature for Langevin model
 929  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,Temp)
      tK=3.1668D-6*Temp
      Go To 998
*                                                                      *
****** RSCA ************************************************************
*                                                                      *
*     Simultaneous scaling of all radii
 930  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,rsca)
      Go To 998

*                                                                      *
****** ROTA ************************************************************
*                                                                      *
*     Rotation of the grid with Euler angles alpha, beta, gamma (degrees)
 931  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,rotAlpha)
      Call Get_F1(2,rotBeta)
      Call Get_F1(3,rotGamma)

      rotAlpha=rotAlpha/180.0D0*3.1415926535897932D0
      rotBeta=rotBeta/180.0D0*3.1415926535897932D0
      rotGamma=rotGamma/180.0D0*3.1415926535897932D0
      Go To 998

*                                                                      *
****** SPAR ************************************************************
*                                                                      *
*     Use sparse grid outside a distance distSparse from all atoms
 932  KWord=Get_Ln(LuSpool)
      LSparse = .True.
      Call Get_I1(1,nSparse)
      Call Get_F1(2,distSparse)
      Go To 998

*                                                                      *
****** IDAM ************************************************************
*                                                                      *
*     change dampIter

 933  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,dampIter)
      Go To 998

*                                                                      *
****** AVER ************************************************************
*                                                                      *
*     requests average among a number of random grids in the Langevin
*     routine. Note: Average is taken for each SCF iteration and only
*     the last grid is used for updating the density
*     nGridSeed is RNG seed, 0=no seed, -1=system timer

 934  KWord=Get_Ln(LuSpool)
      LGridAverage=.True.
      Call Get_I1(1,nGridAverage)
      Call Get_I1(2,nGridSeed)
      Go To 998

*                                                                      *
****** CONV ************************************************************
*                                                                      *
*     Change convergence threshold

 935  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,clim)
      Go To 998

*                                                                      *
****** SKIP ************************************************************
*                                                                      *
*     Set minimal distance for dipole-dipole interactions

 936  KWord=Get_Ln(LuSpool)
      Call Get_F1(1,dipCutoff)
      Go To 998
*                                                                      *
****** NODA ************************************************************
*                                                                      *
*     Turn off damping of dipole-dipole interactions

 937  lDamping=.False.
      Go To 998

*                                                                      *
****** LRAD ************************************************************
*                                                                      *
*     Input non-standard Langevin radii

 938  Continue
      Kword=Get_Ln(LuSpool)
      Call UpCase(KWord)
      If (KWord(1:1).eq.'*')    Go To 938
      If (KWord.eq.'')       Go To 938
      If (KWord(1:3).eq.'END')  Go To 998
      Read(KWord,*,err=988) ii,val
      CovRadT_(ii)=val
c      write(6,*)'COVR',ii,val
      Go To 938
*                                                                      *
****** AMBE ************************************************************
*                                                                      *
*     Use dipole-dipole cutoff of AMBER type

 939  Continue
      lAmberPol=.True.
      Go To 998
*                                                                      *
****** SC14 ************************************************************
*                                                                      *
*     Scaling of negative entries in exclusion list

 940  Continue
      KWord=Get_Ln(LuSpool)
      Call Get_F1(1,scal14)

      Go To 998

*                                                                      *
****** DRES ************************************************************
*                                                                      *
*     Restart dipoles from scratch in each QM iteration
*     This sometimes gives better convergence.

 941  Continue
      lDiprestart=.True.
      Go To 998

*                                                                      *
****** END  ************************************************************
*                                                                      *
*-----End of input
*
 997  Continue
*                                                                      *
************************************************************************
*                                                                      *
      If (lLangevin) Then
*
         polsi=polsi*scaaa**Three
         scala=scala*scaaa
         scalb=scalb*scaaa
         scalc=scalc*scaaa
         dipsi=dipsi*scaaa**(Three/Two)
*
*        Analyze the lattice
*
         gatom=Anal_Gitt(cordsi,latato)
*
*        Scale the grid.
*
         Do i=1,latato
            cordsi(1,i)=cordsi(1,i)*scala
            cordsi(2,i)=cordsi(2,i)*scalb
            cordsi(3,i)=cordsi(3,i)*scalc
         End Do
*
         If (lRFCav) Then
            If (radlat.eq.Zero) radlat=rds-One/Ten
            If (radlat.gt.rds) Then
               Call WarningMessage(2,'InpRct: radlat.gt.rds')
               Call Abend()
            End If
         Else
            Call WarningMessage(1,
     &           'Running Langevin without reaction field')
         EndIf
*
         tK=One/tK
         poltot=gatom*(polsi+dipsi**2*tK/Three)
         tal=poltot*Four*Pi/(Three*scala*scalb*scalc)
         epscm=(One+Two*tal)/(One-tal)
         If (Eps.lt.One) Eps=epscm
*
         tk5=Half*tk
         r2=radlat**2
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Call qExit('InpRct')
      Return
*                                                                      *
************************************************************************
*                                                                      *
*-----Error handling
*
 977  Call WarningMessage(2,'InpRct: Premature end of input file.')
      Call Quit_OnUserError()
 988  Call WarningMessage(2,'InpRct: Error while reading input file.')
      Call Quit_OnUserError()
*
      End
