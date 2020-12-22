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
      use Symmetry_Info, only: nIrrep
      use Slapaf_Info, only:  Coor
      use Slapaf_Parameters, only: mTROld
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 DipM(3)
      Integer mDisp(8)
      Real*8, Allocatable:: EVec(:), EVal(:), RedMas(:), dDipM(:),
     &                      Tmp1(:), Tmp2(:), IRInt(:), Temp(:),
     &                      NMod(:)
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      lUt=6
*
      nX=3*SIZE(Coor,2)
      nAtom=SIZE(Coor,2)
      nInter=mInt
      mTR=mTROld
      nDoF=nInter+mTR
*
      Call mma_allocate(EVec,2*nX**2,Label='EVec')
      Call mma_allocate(EVal,2*nX,Label='EVal')
      Call mma_allocate(RedMas,nX,Label='RedMas')
      Call mma_allocate(dDipM,3*nDoF,Label='dDipM')
      dDipM(:)=Zero
      Call mma_allocate(Tmp1,nX**2,Label='Tmp1')
      Call mma_allocate(Tmp2,nX**2,Label='Tmp2')
*                                                                      *
************************************************************************
*                                                                      *
      DipM(1)=0.0D0
      DipM(2)=0.0D0
      DipM(3)=0.0D0
      Call GF(nX,nDoF,nInter,Tmp1,Tmp2,EVec,EVal,RedMas,iNeg,dDipM,
     &        mTR,nAtom,DipM)
*
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp1)
#ifdef _DEBUGPRINT_
      Call RecPrt('EVec',' ',EVec,2*nX,nX)
      Call RecPrt('EVal',' ',EVal,2,nX)
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
     &             //' have been automatically removed,'
      Write (LUt,10) ' if the energy is invariant to these degrees'
     &             //' of freedom.'
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
      Call mma_allocate(Temp,3*nDoF,Label='Temp')
*
      Call DGeTMO(dDipM,3,3,nInter,Temp,nInter)
*
      Call mma_deallocate(dDipM)

      Call mma_allocate(IRInt,nDoF,Label='IRInt')
*
      Lu_10=10
      Lu_10=IsFreeUnit(Lu_10)
      Call Molcas_Open(lu_10,'UNSYM')
*
      Write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',jsym
      Call GF_Print(EVal,EVec,Temp,iEl,nDoF,nInter,iCtl,IRInt,RedMas,
     &              Lu_10,iOff)
*
      Close(Lu_10)
      Call mma_deallocate(Temp)
*
      Call Add_Info('Approx. Freq.',EVal,nInter,1)
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
#ifdef _DEBUGPRINT_
      Call RecPrt('NModes',' ',NMod,nX,nInter)
#endif
*
      If (nIrrep.eq.1) Then
         Call Print_Mode_Components(NMod,EVal,nModes,lModes,mDisp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write stuff on Molden input file
*
      mSym=1
      Call ICopy(8,[0],0,mDisp,1)
      mDisp(1)=nInter
      Call Freq_Molden(EVal,nModes,NMod,lModes,mSym,IRInt,mDisp,RedMas)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(NMod)
      Call mma_deallocate(IRInt)
      Call mma_deallocate(RedMas)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
