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
      Subroutine HrmFrq(nAtom,nInter,iNeg,dDipM,mTR,DipM,IRInt)
      use thermochem
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 dDipM(3,nInter+mTR), DipM(3), IRInt(nInter+mTR)
      Integer mDisp(8)
      Real*8, Allocatable:: EVec(:), EVal(:), RedMas(:), Tmp1(:),
     &                      Tmp2(:), Temp(:), NMod(:)
*                                                                      *
************************************************************************
*                                                                      *
      LUt=6
*                                                                      *
************************************************************************
*                                                                      *
      nDoF = nInter + mTR
      nX = 3*nAtom
*
      Call mma_allocate(EVec,2*nDoF**2,Label='EVec')
      Call mma_allocate(EVal,2*nDoF,Label='EVal')
      Call mma_allocate(RedMas,nDoF,Label='RedMas')
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute harmonic frequencies and dipole derivatives
*
      Call mma_allocate(tmp1,nX**2,Label='tmp1')
      Call mma_allocate(tmp2,nX**2,Label='tmp2')
*
      Call GF(nX,nDoF,nInter,Tmp1,Tmp2,EVec,EVal,RedMas,iNeg,dDipM,
     &        mTR,nAtom,DipM)
*
      Call mma_deallocate(tmp2)
      Call mma_deallocate(tmp1)
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
     &             //' have been automatically removed,'
      Write (LUt,10) ' if the energy is invariant to these degrees'
     &             //' of freedom.'
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
      Call mma_allocate(Temp,3*nDoF,Label='Temp')
      Call DGeTMO(dDipM,3,3,nInter,Temp,nInter)
*
      Lu_10=10
      Lu_10=IsFreeUnit(Lu_10)
      Call Molcas_Open(lu_10,'UNSYM')
*
      Write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym
      Call GF_Print(EVal,EVec,Temp,iEl,nDoF,nInter,iCtl,IRInt,RedMas,
     &             Lu_10,iOff)
*
      Close(Lu_10)
      Call mma_deallocate(Temp)
*
      If (lTherm) Call Thermo_Driver(UserT,UserP,nUserPT,nsRot,EVal,
     &                               nInter,lTherm)
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
      nDisp=nDoF
      Call mma_allocate(NMod,nDisp**2,Label='NMod')
*
      lModes=0
      nModes=0
      nX=nDoF
      call dcopy_(nX*nInter,EVec,2,NMod,1)
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
      Call Print_Mode_Components(NMod,EVal,nModes,lModes,mDisp)
      Call Freq_Molden(EVal,nModes,NMod,lModes,nSym,IRInt,mDisp,RedMas)
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate memory
*
      Call mma_deallocate(NMod)
      Call mma_deallocate(RedMas)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
