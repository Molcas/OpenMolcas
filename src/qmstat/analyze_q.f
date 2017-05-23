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
      Subroutine Analyze_Q(iQ_Atoms)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Parameter(iHUltraMax=1000)
      Dimension iCo(3)
      Dimension gR(MxAt,3,iHUltraMax)
      Data Dum/0.0d0/

*----------------------------------------------------------------------*
* Some numbers and defaults.                                           *
*----------------------------------------------------------------------*
      iHMax=0
      iCStart=(((iQ_Atoms-1)/nAtom)+1)*nCent+1
      iCNum=iCStart/nCent
      dR=0.1d0
*----------------------------------------------------------------------*
* Just say what we are doing.                                          *
*----------------------------------------------------------------------*
      Call NiceOutPut('AAA',Dum,Dum,Dum)
*----------------------------------------------------------------------*
* Open sampfile. Get some numbers about how many sampled etc.          *
*----------------------------------------------------------------------*
      Call DaName(iLuSaIn,SaFilIn)
      iDiskSa=0
      Call iDaFile(iLuSaIn,2,iHowMSamp,1,iDiskSa)
      iDiskTemp=iDiskSa
      Call WrRdSim(iLuSaIn,2,iDiskSa,iTCSim,64,Etot,Ract,nPart,Dum,Dum
     &            ,Dum)
      iDiskSa=iDiskTemp
*----------------------------------------------------------------------*
* Say something about these numbers to the user.                       *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,*)'The sampfile ',SaFilIn,' contains ',iHowMSamp,' sample'
     &//'d configurations.'
      Write(6,*)'Total number of particles:', nPart
*----------------------------------------------------------------------*
* BEGIN ANALYZING!                                                     *
*----------------------------------------------------------------------*
      Do 101, iSamp=1,iHowMSamp
*----------------------------------------------------------------------*
* Begin by getting the coordinates for this configuration. They are    *
* stored in Work(iCo(i)) where i=1 means x-coordinate, i=2 y-coordinate*
* and i=3 z-coordinate.                                                *
*----------------------------------------------------------------------*
        Call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Dum
     &              ,Dum,Dum)
        iDiskSa=iTcSim(1)
        Do 1001, i=1,3
          Call GetMem('Coordinates','Allo','Real',iCo(i),nPart*nCent)
          Call dDafile(iLuSaIn,2,Work(iCo(i)),nPart*nCent,iDiskSa)
1001    Continue
*----------------------------------------------------------------------*
* Once we have coordinates, lets compute some distances and start      *
* building various distribution functions.                             *
*----------------------------------------------------------------------*
        Do 111, i=1,iQ_Atoms
          Do 112, j=1,nAtom
            Do 113, k=1,nPart-iCNum
              dist2=0.0d0
              Do 114, l=1,3
                ind=iCStart+(j-1)+(k-1)*nCent
                dist2=dist2+(Work(iCo(l)+i-1)-Work(iCo(l)+ind-1))**2
114           Continue
              dist=sqrt(dist2)
              iH=int((dist+dR*0.5d0)/dR)
              if(iH.gt.iHMax) then
                iHMax=iH
                if(iH.gt.iHUltraMax) then
                  Write(6,*)
                  Write(6,*)'Too fine sections for g(r). Increase secti'
     &//'on size or allocate more memory.'
                  Call Quit(_RC_INTERNAL_ERROR_)
                Endif
              Endif
              gR(i,j,iH)=gR(i,j,iH)+1/dist2
113         Continue
112       Continue
111     Continue
*----------------------------------------------------------------------*
* End loop over sampled configurations.                                *
*----------------------------------------------------------------------*
        Do 1002, i=1,3
          Call GetMem('Coordinates','Free','Real',iCo(i),nPart*nCent)
1002    Continue
101   Continue
*----------------------------------------------------------------------*
* Time to generate a nice output.                                      *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,*)'SUMMARY OF RESULTS FOR SAMPFILE ANALYSIS.'
      Write(6,*)
      Do 9991, i=1,iQ_Atoms
        Write(6,*)
        Write(6,*)'Quantum atom ',i
        Write(6,'(5X,A,5X,5(A,I2,1X))')'Separation'
     &                  ,('Solvent atom',k,k=1,nAtom)
        Do 9992, iH=1,iHMax
          Write(6,'(F15.7,5(F15.7))')dR*iH,(gR(i,j,iH),j=1,nAtom)
9992    Continue
9991  Continue

      Call DaClos(iLuSaIn)

      Return
      End
