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
      Subroutine GF_on_the_Fly(iDo_dDipM)
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 DipM(3)
      Integer mDisp(8)
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
      lUt=6
*
      nX=3*nsAtom
      nAtom=nsAtom
      nInter=mInt
      mTR=mTROld
      mInter=nInter+mTR
*
      Call Allocate_Work(ipEVec,2*nX**2)
      Call Allocate_Work(ipEVal,2*nX)
      Call Allocate_Work(ipRedMas,nX)
      Call Allocate_Work(ipdDipM,3*mInter)
      Call Allocate_Work(ipTmp1,nX**2)
      Call Allocate_Work(ipTmp2,nX**2)
      Call FZero(Work(ipdDipM),3*mInter)
*                                                                      *
************************************************************************
*                                                                      *
      DipM(1)=0.0D0
      DipM(2)=0.0D0
      DipM(3)=0.0D0
      Call GF(nX,mInter,nInter,Work(ipTmp1),Work(ipTmp2),
     &        Work(ipEVec),Work(ipEVal),Work(ipRedMas),
     &        iNeg,iDo_dDipM,Work(ipdDipM),mTR,Smmtrc,nAtom,DipM)
*
      Call Free_Work(ipTmp2)
      Call Free_Work(ipTmp1)
#ifdef _DEBUG_
      Call RecPrt('EVec',' ',Work(ipEVec),2*nX,nX)
      Call RecPrt('EVal',' ',Work(ipEVal),2,nX)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Write out frequencies and IR intensities
*
      Write (LUt,10)
      Write (LUt,10) ' Observe that the harmonic oscillator analysis'
     &             //' is only valid at stationary points!'
      Write (LUt,10)
      Write (LUt,10) ' Note that rotational and translational degrees'
     &             //' have been automatically removed.'
      Write (LUt,10)
      Write (LUt,10)
      Write (LUt,10) ' Harmonic frequencies in cm-1'
      Write (LUt,10)
      If (iDo_dDipM.eq.1) Then
         Write (LUt,10) ' IR Intensities in km/mol'
         Write (LUt,10)
      End If
      Write (LUt,10) ' Normal modes in gf_on_the_fly.f '
10    Format(a)
*
      iOff=0
      iCtl=iDo_dDipM
      iEl =3
      jSym=1
*
      Call Allocate_Work(ipTemp,3*mInter)
*
      Call DGeTMO(Work(ipdDipM),3,3,nInter,Work(ipTemp),nInter)
*
      Call Free_Work(ipdDipM)
      Call Allocate_Work(ipIRInt,mInter)
*
      Lu_10=10
      Lu_10=IsFreeUnit(Lu_10)
      Call Molcas_Open(lu_10,'UNSYM')
*
      Write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',jsym
      Call GF_Print(Work(ipEVal),Work(ipEVec),Work(ipTemp),iEl,
     &              mInter,nInter,iCtl,Work(ipIRInt),Work(ipRedMas),
     &              Lu_10,iOff)
*
      Close(Lu_10)
      Call Free_Work(ipTemp)
*
      Call Add_Info('Approx. Freq.',Work(ipEVal),nInter,1)
*                                                                      *
************************************************************************
*                                                                      *
*-----Save normal modes for later generation of Molden input.
*
      nDisp=mInter
      Call GetMem('NMod','Allo','Real',ipNMod,nDisp**2)
      ipNx=ipNMod
*
      lModes=0
      nModes=0
      nX=mInter
      call dcopy_(nX*nInter,Work(ipEVec),2,Work(ipNMod),1)
      lModes=lModes+nInter*nX
      nModes=nModes+nInter
#ifdef _DEBUG_
      Call RecPrt('NModes',' ',Work(ipNMod),nX,nInter)
#endif
*
      If (nSym.eq.1) Then
         Call Print_Mode_Components(Work(ipNMod),Work(ipEVal),
     &                              nModes,lModes,mDisp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write stuff on Molden input file
*
      mSym=1
      Call ICopy(8,0,0,mDisp,1)
      mDisp(1)=nInter
      Call Freq_Molden(Work(ipEVal),nModes,Work(ipNMod),lModes,mSym,
     &                 Work(ipIRInt),mDisp,Work(ipRedMas))
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Work(ipNMod)
      Call Free_Work(ipIRInt)
      Call Free_Work(ipEVal)
      Call Free_Work(ipEVec)
      Call Free_Work(ipRedMas)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
