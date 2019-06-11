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
      Subroutine HrmFrq(nAtom,nInter,iNeg,dDipM,mTR,Smmtrc,DipM,
     &                  IRInt, UserT, UserP, nUserPT, nsRot, lTherm,
     &                  lDoubleIso)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 dDipM(3,nInter+mTR), DipM(3), IRInt(nInter+mTR)
      Logical Smmtrc(3,nAtom), lTherm, lDoubleIso
      Integer mDisp(8), nUserPT, nsRot
      Real*8  UserT(64), UserP
*                                                                      *
************************************************************************
*                                                                      *
      Call QEnter('HrmFrq')
      LUt=6
*                                                                      *
************************************************************************
*                                                                      *
      mInter = nInter + mTR
      nX = 3*nAtom
*
      Call GetMem('EVec','Allo','Real',ipEVec,2*mInter**2)
      Call GetMem('EVal','Allo','Real',ipEVal,2*mInter)
      Call GetMem('RedMas','Allo','Real',ipRedMas,mInter)
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute harmonic frequencies and dipole derivatives
*
      Call GetMem('tmp1','Allo','Real',iptmp1,nX**2)
      Call GetMem('tmp2','Allo','Real',iptmp2,nX**2)
*
      iDo_dDipM=1
      Call GF(nX,mInter,nInter,Work(ipTmp1),Work(ipTmp2),Work(ipEVec),
     &        Work(ipEVal),Work(ipRedMas),iNeg,iDo_dDipM,dDipM,mTR,
     &        Smmtrc,nAtom,DipM)
*
      Call GetMem('tmp2','Free','Real',iptmp2,(3*nAtom)**2)
      Call GetMem('tmp1','Free','Real',iptmp1,(3*nAtom)**2)
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
      Write (LUt,10) ' IR Intensities in km/mol'
      Write (LUt,10)
10    Format(a)
*
      iOff=0
      iCtl=1
      iEl =3
      iSym=1
*
      Call Allocate_Work(ipTemp,3*mInter)
      Call DGeTMO(dDipM,3,3,nInter,Work(ipTemp),nInter)
*
      Lu_10=10
      Lu_10=IsFreeUnit(Lu_10)
      Call Molcas_Open(lu_10,'UNSYM')
*
      Write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym
      Call GF_Print(Work(ipEVal),Work(ipEVec),Work(ipTemp),iEl,
     &              mInter,nInter,iCtl,IRInt,Work(ipRedMas),Lu_10,iOff)
*
      Close(Lu_10)
      Call Free_Work(ipTemp)
*
      If (lTherm) Call Thermo_Driver(UserT,UserP,nUserPT,nsRot,
     &                                    ipEVal,nInter,lTherm)
*                                                                      *
************************************************************************
*                                                                      *
*     Do single and double isotope substitutions
*
      Call Get_iScalar('NSYM',nIrrep)
      If (nIrrep.eq.1) Call IsoLoop(lDoubleIso)
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
*                                                                      *
************************************************************************
*                                                                      *
*     Write stuff on Molden input file
*
      nSym=1
      Call ICopy(8,[0],0,mDisp,1)
      mDisp(1)=nInter
!      If (nIrrep.eq.1) Then
         Call Print_Mode_Components(Work(ipNMod),Work(ipEVal),
     &                              nModes,lModes,mDisp)
!      End If
      Call Freq_Molden(Work(ipEVal),nModes,Work(ipNMod),lModes,nSym,
     &                 IRInt,mDisp,Work(ipRedMas))
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate memory
*
      Call GetMem('NMod','Free','Real',ipNMod,nDisp**2)
      Call GetMem('EVal','Free','Real',ipEVal,2*nInter)
      Call GetMem('EVec','Free','Real',ipEVec,2*nInter**2)
      Call GetMem('RedMas','Free','Real',ipRedMas,mInter)
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('HrmFrq')
      Return
      End
