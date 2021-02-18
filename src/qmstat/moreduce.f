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
      Subroutine MoReduce(nBas,MOsToKeep,ipAvRedMO)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
*Removed by Valera or maybe it wasn't in his version #include "MolcasQmStat.fh"
#include "lenin.fh"

      Parameter (ReduceWarning=0.5d0)

      Dimension nBas(MxSym),iTocBig(MxStOT)

      Character Header*50,BsLbl*1

      Dimension BsLbl(LENIN8*MxBas)

      Logical LindMOs(MxBas),First

      Data First /.true./

      Dimension Dummy(1)

*
*--- A word of welcome.
*
      Write(6,*)'     ----- Transform to average natural MO-reduced'
     &//' basis.'
*
*--- First we accumulate the different density matrices.
*
      nSize=nBas(1)*(nBas(1)+1)/2
      weight=ONE/dble(nState)
      Call GetMem('DenM','Allo','Real',iDin,nSize)
      Call GetMem('DenA','Allo','Real',iDav,nSize)
      call dcopy_(nSize,[ZERO],iZERO,Work(iDav),iONE)
      Do 201, iS1=1,nState
        Do 202, iS2=1,iS1
          index=(iS1*(iS1-1)/2+iS2-1)*nSize
          call dcopy_(nSize,Work(iBigT+index),iONE,Work(iDin),iONE)
          If(iS1.ne.iS2) GoTo 202
          kaunt=0
          Do 203, iB1=1,nBas(1)
            Do 204, iB2=1,iB1
              If(iB1.eq.iB2) then
                Fac=1.0d0*weight
              Else
                Fac=0.5d0*weight
              Endif
              Work(iDav+kaunt)=Work(iDav+kaunt)+Work(iDin+kaunt)*Fac
              kaunt=kaunt+1
204         Continue
203       Continue
202     Continue
201   Continue
*
*--- Then since we are working in the non-orthogonal AO-basis,
*    it is necessary to orthogonalize before we diagonalize
*    accumulated density.
*
      Call GetMem('Vecs','Allo','Real',iVecs,nBas(1)**2)
      Call GetMem('Vecs2','Allo','Real',iVecs2,nBas(1)**2)
      Call GetMem('AuxS','Allo','Real',iAUX,nBas(1)**2)
      Call GetMem('OvlSs','Allo','Real',iSs,nBas(1)**2)
      Call GetMem('OvlSs','Allo','Real',iSx,nSize)
      Call GetMem('OvlSi','Allo','Real',iSst,nBas(1)**2)
      Call GetMem('OvlSi','Allo','Real',iSt,nSize)
      Call GetMem('DavSq','Allo','Real',iDavS,nBas(1)**2)
      Call GetMem('Trans','Allo','Real',iTrans,nBas(1)**2)
      Call GetMem('Trans','Allo','Real',iTransB,nBas(1)**2)
      Call GetMem('OrtoAvDen','Allo','Real',iOtD,nBas(1)**2)
      Call GetMem('OrtoAcDeT','Allo','Real',iOtDt,nSize)
      Call GetMem('OvlS','Allo','Real',iS,nSize+4)
      Call GetMem('Occs','Allo','Real',iOcc,nBas(1))
      kaunter=0
      Do 211, iB1=1,nBas(1)
        Do 212, iB2=1,nBas(1)
          Work(iVecs+kaunter)=0
          If(iB1.eq.iB2)Work(iVecs+kaunter)=1
          kaunter=kaunter+1
212     Continue
211   Continue
*--- Symmetric orthogonalization, hence get overlap matrix, S.
      Lu_One=92
      Call OpnOne(irc,0,'ONEINT',Lu_One)
      irc=-1
      iopt=0
      iSmLbl=0
      icomp=1
      Call RdOne(irc,iopt,'Mltpl  0',icomp,Work(iS),iSmLbl)
      Call Jacob(Work(iS),Work(iVecs),nBas(1),nBas(1))
      call dcopy_(nSize,[ZERO],iZERO,Work(iSx),iONE)
      call dcopy_(nSize,[ZERO],iZERO,Work(iSt),iONE)
      Do 221, i=1,nBas(1)
        Sqroot=sqrt(Work(iS+i*(i+1)/2-1))
        Work(iSx+i*(i+1)/2-1)=ONE/Sqroot
        Work(iSt+i*(i+1)/2-1)=Sqroot
221   Continue
      Call Square(Work(iSx),Work(iSs),iONE,nBas(1),nBas(1))
      Call Square(Work(iSt),Work(iSst),iONE,nBas(1),nBas(1))
*------S^(-1/2)
      Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iVecs)
     &          ,nBas(1),Work(iSs),nBas(1),ZERO,Work(iAUX),nBas(1))
      Call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),ONE,Work(iAUX)
     &          ,nBas(1),Work(iVecs),nBas(1),ZERO,Work(iTrans),nBas(1))
*------S^(1/2)
      Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iVecs)
     &          ,nBas(1),Work(iSst),nBas(1),ZERO,Work(iAUX),nBas(1))
      Call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),ONE,Work(iAUX)
     &          ,nBas(1),Work(iVecs),nBas(1),ZERO,Work(iTransB),nBas(1))
*--- The density matrix transforms 'inversly' from the matrix-elements,
*    thus let S^(1/2) transform it, not S^(-1/2) which applies to the
*    matrix elements.
      Call Square(Work(iDav),Work(iDavS),iONE,nBas(1),nBas(1))
      Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iTransB)
     &          ,nBas(1),Work(iDavS),nBas(1),ZERO,Work(iAUX),nBas(1))
      Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iAUX)
     &          ,nBas(1),Work(iTransB),nBas(1),ZERO,Work(iOtD),nBas(1))
      kaunter=0
      Do 213, iB1=1,nBas(1)
        Do 214, iB2=1,nBas(1)
          Work(iVecs2+kaunter)=0
          If(iB1.eq.iB2)Work(iVecs2+kaunter)=1
          kaunter=kaunter+1
214     Continue
213   Continue
      Call SqToTri_Q(Work(iOtD),Work(iOtDt),nBas(1))
      Call Jacob(Work(iOtDt),Work(iVecs2),nBas(1),nBas(1))
*--- With diagonalized density matrix, collect occupation numbers and
*    natural orbital coefficients.
      Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iTrans)
     &          ,nBas(1),Work(iVecs2),nBas(1),ZERO,Work(iAUX),nBas(1))
      kaunt=0
      kaunter=0
      Do 231, i=1,nBas(1)
        Do 232, j=1,i
          If(i.eq.j) then
            Work(iOcc+kaunt)=Work(iOtDt+kaunter)
            kaunt=kaunt+1
          Endif
          kaunter=kaunter+1
232     Continue
231   Continue
      TraceFull=0
      Do 233, i=1,nBas(1)
        TraceFull=TraceFull+Work(iOcc+i-1)
233   Continue
      If(iPrint.ge.10) then
        Call Get_cArray('Unique Basis Names',BsLbl,LENIN8*nBas(1))
        Write(Header,'(A)')'All average transition density orbitals'
        ThrOcc=-1D-0
        Call Primo(Header,.true.,.false.,ThrOcc,Dum,iONE,nBas(1),nBas(1)
     &            ,BsLbl,Dummy,Work(iOcc),Work(iAUX),-1)
        Write(6,*)
        Write(6,*)'  Trace = ',TraceFull
      Endif
*--- Deallocations.
      Call GetMem('DenM','Free','Real',iDin,nSize)
      Call GetMem('DenA','Free','Real',iDav,nSize)
      Call GetMem('Vecs','Free','Real',iVecs,nBas(1)**2)
      Call GetMem('Vecs2','Free','Real',iVecs2,nBas(1)**2)
      Call GetMem('OvlSs','Free','Real',iSs,nBas(1)**2)
      Call GetMem('OvlSs','Free','Real',iSx,nSize)
      Call GetMem('OvlSi','Free','Real',iSst,nBas(1)**2)
      Call GetMem('OvlSi','Free','Real',iSt,nSize)
      Call GetMem('DavSq','Free','Real',iDavS,nBas(1)**2)
      Call GetMem('Trans','Free','Real',iTrans,nBas(1)**2)
      Call GetMem('Trans','Free','Real',iTransB,nBas(1)**2)
      Call GetMem('OrtoAvDen','Free','Real',iOtD,nBas(1)**2)
      Call GetMem('OrtoAcDeT','Free','Real',iOtDt,nSize)
      Call GetMem('OvlS','Free','Real',iS,nSize+4)

*
*--- Jetzt far wir mal wieder. Reduce MO-basis according to input
*    criterion.
*
      MOsToKeep=0
      Do 301, iB=1,nBas(1)
        If(Work(iOcc+iB-1).ge.ThrsRedOcc) then
          LindMOs(iB)=.true.
          MOsToKeep=MOsToKeep+1
        Else
          LindMOs(iB)=.false.
        Endif
301   Continue
      nSize=nBas(1)*MOsToKeep
      Call GetMem('UncleMoe','Allo','Real',ipAvRedMO,nSize)
      Call GetMem('NewOccs','Allo','Real',iNewOcc,MOsToKeep)
      ind2=0
      ind3=0
*--- Loop to suck-out the nice MOs.
      Do 302, iB=1,nBas(1)
        If(LindMOs(iB)) then
          ind1=nBas(1)*(iB-1)
          call dcopy_(nBas(1),Work(iAUX+ind1),iONE,Work(ipAvRedMO+ind2)
     &              ,iONE)
          Work(iNewOcc+ind3)=Work(iOcc+iB-1)
          ind2=ind2+nBas(1)
          ind3=ind3+1
        Endif
302   Continue
      TraceRed=0
      Do 303, i=1,MOsToKeep
        TraceRed=TraceRed+Work(iNewOcc+i-1)
303   Continue
*--- Make a trace check.
      If((TraceFull-TraceRed).ge.ReduceWarning) then
        Write(6,*)
        Write(6,*)'WARNING!  With your occupation threshold, the densit'
     &//'y matrix trace'
        Write(6,*)'differs by ',TraceFull-TraceRed,'.'
        Write(6,*)'You should consider lowering the threshold!'
      Endif
      If(iPrint.ge.5) then
        Call Get_cArray('Unique Basis Names',BsLbl,LENIN8*nBas(1))
        Write(Header,'(A)')'Reduced average orbitals'
        ThrOcc=-1D-0
        Call Primo(Header,.true.,.false.,ThrOcc,Dum,iONE,nBas(1)
     &            ,[MOsToKeep],BsLbl,Dummy,Work(iNewOcc)
     &            ,Work(ipAvRedMO),-1)
        Write(6,*)
        Write(6,*)'  Trace = ',TraceRed,MOsToKeep
      Endif

*
*--- Time to reduce all individual density matrices in the big TDM to
*    the reduced MO-basis. Once more, observe that the density
*    transforms contravariantly. But sadly, we need to invert the
*    full square MO-matrix before we make reductions.
*
      Call GetMem('InverseC','Allo','Real',ipInv,nBas(1)**2)
      Call MInv(Work(iAUX),Work(ipInv),Ising,Det,nBas(1))

*
*--- Now all those transformations and density reductions. To check for
*    density losses, the overlaps are read. These partially transformed
*    transition density matrix is stored in a scratch file.
*
      DiffMegaMax=0.0d0
      nSize=nBas(1)*(nBas(1)+1)/2
      nMtK=MOsToKeep*(MOsToKeep+1)/2
      Call GetMem('Temporary','Allo','Real',ipTEMP,nBas(1)**2)
      Call GetMem('MOtrDen','Allo','Real',ipTmoD,nBas(1)**2)
      Call GetMem('MOreDen','Allo','Real',ipTreD,MOsToKeep**2)
      Call GetMem('MOreDen','Allo','Real',ipTreT,nMtK)
      Call GetMem('DenM','Allo','Real',iDin,nSize)
      Call GetMem('DenMsq','Allo','Real',iDsq,nBas(1)**2)
      Call GetMem('OvlS','Allo','Real',iS,nSize+4)
      Call GetMem('Ssquare','Allo','Real',iSsq,nBas(1)**2)
      Call GetMem('Strans','Allo','Real',iStrans,nBas(1)**2)
      Call GetMem('Stri','Allo','Real',iStri,nMtK)
      Lu_Scratch=57
      Lu_Scratch=IsFreeUnit(Lu_Scratch)
      Call DaName(Lu_Scratch,'TDMSCR')
      irc=-1
      iopt=0
      iSmLbl=0
      icomp=1
      Call RdOne(irc,iopt,'Mltpl  0',icomp,Work(iS),iSmLbl)
      Call ClsOne(irc,iopt)
      iDiskUt=0
      Call iDaFile(Lu_Scratch,1,iTocBig,MxStOT,iDiskUt)
      Do 401, iS1=1,nState
        Do 402, iS2=1,iS1
*------- Collect this particular density matrix.
          index=(iS1*(iS1-1)/2+iS2-1)*nSize
          call dcopy_(nSize,Work(iBigT+index),iONE,Work(iDin),iONE)
*------- Square it and correct the non-diagonal (recall convention)
          Call Dsq(Work(iDin),Work(iDsq),iONE,nBas(1),nBas(1))
*------- Contravariant transformation of density matrix.
          Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(ipInv)
     &              ,nBas(1),Work(iDsq),nBas(1),ZERO,Work(ipTEMP)
     &              ,nBas(1))
          Call Dgemm_('N','T',nBas(1),nBas(1),nBas(1),ONE,Work(ipTEMP)
     &              ,nBas(1),Work(ipInv),nBas(1),ZERO,Work(ipTmoD)
     &              ,nBas(1))
*------- Covariant transformation of overlap matrix.
          Call Square(Work(iS),Work(iSsq),iONE,nBas(1),nBas(1))
          Call Dgemm_('T','N',nBas(1),nBas(1),nBas(1),ONE,Work(iAUX)
     &              ,nBas(1),Work(iSsq),nBas(1),ZERO,Work(ipTEMP)
     &              ,nBas(1))
          Call Dgemm_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(ipTEMP)
     &              ,nBas(1),Work(iAUX),nBas(1),ZERO,Work(iStrans)
     &              ,nBas(1))
*------- How much charge ('overlap') is there in this density element?
          ChargeNonReduced=Ddot_(nBas(1)**2,Work(ipTmoD),iONE
     &                         ,Work(iStrans),iONE)
*------- Reduction of density matrix and overlap matrix.
          kaunter=0
          ind1=0
          Do 405, iM1=1,nBas(1)
            Do 406, iM2=1,nBas(1)
              If(LindMOs(iM1).and.LindMOs(iM2)) then
                Work(ipTreD+ind1)=Work(ipTmoD+kaunter)
                Work(iSsq+ind1)=Work(iStrans+kaunter)
                ind1=ind1+1
              Endif
              kaunter=kaunter+1
406         Continue
405       Continue
*------- Triangualize, jetzt!
          kaunter=0
          Do 407, iM1=1,MOsToKeep
            Do 408, iM2=1,MOsToKeep
              If(iM1.ne.iM2)Work(ipTreD+kaunter)=2*Work(ipTreD+kaunter)
              kaunter=kaunter+1
408         Continue
407       Continue
          Call SqToTri_Q(Work(ipTreD),Work(ipTreT),MOsToKeep)
          Call SqToTri_Q(Work(iSsq),Work(iStri),MOsToKeep)
*------- Compute total electronic charge of this reduced density.
          ChargeReduced=Ddot_(nMtK,Work(ipTreT),iONE,Work(iStri),iONE)
*------- Renormalize to get right charge ('overlap'); to safeguard
*        against zero overlaps, make check.
          If(abs(ChargeNonReduced).le.1D-7.or.
     &          abs(ChargeReduced).le.1D-7) then
            Fac=1.0D0
          Else
            Fac=ChargeNonReduced/ChargeReduced
          Endif
          kaunter=0
          Do 409, iM1=1,MOsToKeep
            Do 410, iM2=1,iM1
              Work(ipTreT+kaunter)=Work(ipTreT+kaunter)*Fac
              kaunter=kaunter+1
410         Continue
409       Continue
*------- If sufficient printlevel, show moment modifications.
          If(iPrint.ge.10) then
            Call MomentMod(ipTreT,ipTmoD,iAUX,MOsToKeep,nBas(1),LindMOs
     &                    ,iS1,iS2,First,DiffMax)
            If(DiffMax.gt.DiffMegaMax) DiffMegaMax=DiffMax
          Endif
*------- Add previous disk address to T-o-C.
          ind=iS1*(iS1+1)/2-iS1+iS2
          iTocBig(ind)=iDiskUt
*------- 'Because I will take a gigant dump on you!'
          Call dDaFile(Lu_Scratch,1,Work(ipTreT),nMtk,iDiskUt)
402     Continue
401   Continue
*--- The real table-of-content
      iDiskUt=0
      Call iDaFile(Lu_Scratch,1,iTocBig,MxStOT,iDiskUt)
*--- Deallocations and closing.
      Call GetMem('NewOccs','Free','Real',iNewOcc,MOsToKeep)
      Call GetMem('InverseC','Free','Real',ipInv,nBas(1)**2)
      Call GetMem('Temporary','Free','Real',ipTEMP,nBas(1)**2)
      Call GetMem('MOtrDen','Free','Real',ipTmoD,nBas(1)**2)
      Call GetMem('MOreDen','Free','Real',ipTreD,MOsToKeep**2)
      Call GetMem('MOreDen','Free','Real',ipTreT,nMtK)
      Call GetMem('DenM','Free','Real',iDin,nSize)
      Call GetMem('DenMsq','Free','Real',iDsq,nBas(1)**2)
      Call GetMem('OvlS','Free','Real',iS,nSize+4)
      Call GetMem('Ssquare','Free','Real',iSsq,nBas(1)**2)
      Call GetMem('Strans','Free','Real',iStrans,nBas(1)**2)
      Call GetMem('Stri','Free','Real',iStri,nMtK)
      Call GetMem('AuxS','Free','Real',iAUX,nBas(1)**2)
      Call GetMem('Occs','Free','Real',iOcc,nBas(1))
      Call DaClos(Lu_Scratch)

*
*--- Report on the reduction.
*
      Write(6,*)
      Write(6,90)'AO-basis ---> MO-basis reduction complete.'
      Write(6,91)'From ',nBas(1),' functions to ',MosToKeep,'.'
      Write(6,90)'Reduced basis renormalized to have same overlap as no'
     &//'n-reduced.'
      If(iPrint.ge.10) then
        Write(6,92)'Largest dipole difference is ',DiffMegaMax
      Endif
90    Format('        ',A)
91    Format('        ',A,I3,A,I3,A)
92    Format('        ',A,F10.7)

      Return
      End


      Subroutine MomentMod(ipRe,ipNRe,iCmo,nBRe,nBNRe,LindMOs,iS1,iS2
     &                    ,First,DiffMax)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "files_qmstat.fh"
#include "qminp.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension DipRe(3),DipNRe(3)

      Logical LindMOs(MxBas),First

      If(First.and.iPrint.ge.5) then
        Write(6,*)
        Write(6,*)'     Modifications of dipoles by renormalization and'
     &//' basis reduction.'
        Write(6,*)
        Write(6,*)'     State pair    |  Difference '
        Write(6,*)'     --------------|---------------------'
        First=.false.
      Endif

      nSize1=nBNRe*(nBNRe+1)/2
      nSize2=nBRe*(nBRe+1)/2
      Call GetMem('DipX','Allo','Real',ipDx,nSize1+4)
      Call GetMem('DipY','Allo','Real',ipDy,nSize1+4)
      Call GetMem('DipZ','Allo','Real',ipDz,nSize1+4)
      Call GetMem('DipXre','Allo','Real',ipDxRe,nSize2)
      Call GetMem('DipYre','Allo','Real',ipDyRe,nSize2)
      Call GetMem('DipZre','Allo','Real',ipDzRe,nSize2)
      Call GetMem('DipXsq','Allo','Real',ipDxsq,nBNRe**2)
      Call GetMem('DipYsq','Allo','Real',ipDysq,nBNRe**2)
      Call GetMem('DipZsq','Allo','Real',ipDzsq,nBNRe**2)
      Call GetMem('DipXm','Allo','Real',ipDxM,nBNRe**2)
      Call GetMem('DipYm','Allo','Real',ipDyM,nBNRe**2)
      Call GetMem('DipZm','Allo','Real',ipDzM,nBNRe**2)
      Call GetMem('TEMP','Allo','Real',ipTEMP,nBNRe**2)
      irc=-1
      iopt=0
      iSmLbl=0
*--- X
      icomp=1
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDx),iSmLbl)
      Call Square(Work(ipDx),Work(ipDxsq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)
     &          ,nBNRe,Work(ipDxsq),nBNRe,ZERO,Work(ipTEMP)
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDxM)
     &          ,nBNRe)
*--- Y
      icomp=2
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDy),iSmLbl)
      Call Square(Work(ipDy),Work(ipDysq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)
     &          ,nBNRe,Work(ipDysq),nBNRe,ZERO,Work(ipTEMP)
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDyM)
     &          ,nBNRe)
*--- Z
      icomp=3
      Call RdOne(irc,iopt,'Mltpl  1',icomp,Work(ipDz),iSmLbl)
      Call Square(Work(ipDz),Work(ipDzsq),iONE,nBNRe,nBNRe)
      Call Dgemm_('T','N',nBNRe,nBNRe,nBNRe,ONE,Work(iCmo)
     &          ,nBNRe,Work(ipDzsq),nBNRe,ZERO,Work(ipTEMP)
     &          ,nBNRe)
      Call Dgemm_('N','N',nBNRe,nBNRe,nBNRe,ONE,Work(ipTEMP)
     &          ,nBNRe,Work(iCmo),nBNRe,ZERO,Work(ipDzM)
     &          ,nBNRe)
*--- Triangualize and reduce.
      kaunt1=0
      kaunt2=0
      Do 2001, i=1,nBNRe
        Do 2002, j=1,nBNRe
          If(j.le.i) then
            If(LindMOs(i).and.LindMOs(j)) then
              Work(ipDxRe+kaunt1)=Work(ipDxM+kaunt2)
              Work(ipDyRe+kaunt1)=Work(ipDyM+kaunt2)
              Work(ipDzRe+kaunt1)=Work(ipDzM+kaunt2)
              kaunt1=kaunt1+1
            Endif
          Endif
          kaunt2=kaunt2+1
2002    Continue
2001  Continue
*--- Density
      DipNRe(1)=Ddot_(nBNRe**2,Work(ipDxM),iONE,Work(ipNRe),iONE)
      DipNRe(2)=Ddot_(nBNRe**2,Work(ipDyM),iONE,Work(ipNRe),iONE)
      DipNRe(3)=Ddot_(nBNRe**2,Work(ipDzM),iONE,Work(ipNRe),iONE)
      DipRe(1)=Ddot_(nSize2,Work(ipDxRe),iONE,Work(ipRe),iONE)
      DipRe(2)=Ddot_(nSize2,Work(ipDyRe),iONE,Work(ipRe),iONE)
      DipRe(3)=Ddot_(nSize2,Work(ipDzRe),iONE,Work(ipRe),iONE)
      Diffx=abs(DipRe(1)-DipNRe(1))
      Diffy=abs(DipRe(2)-DipNRe(2))
      Diffz=abs(DipRe(3)-DipNRe(3))
      If(iPrint.ge.5) then
        Write(6,99)iS1,iS2,'(',Diffx,',',Diffy,',',Diffz,')'
      Endif
99    Format('     ',2I3,'          ',3(A,F10.7),A)
*--- Return number
      DiffMax=Diffy
      If(Diffx.ge.Diffy) DiffMax=Diffx
      If(Diffz.ge.Diffx.and.Diffz.ge.Diffy) DiffMax=Diffz
*--- Deallocate en masse.
      Call GetMem('DipX','Free','Real',ipDx,nSize1+4)
      Call GetMem('DipY','Free','Real',ipDy,nSize1+4)
      Call GetMem('DipZ','Free','Real',ipDz,nSize1+4)
      Call GetMem('DipXre','Free','Real',ipDxRe,nSize2)
      Call GetMem('DipYre','Free','Real',ipDyRe,nSize2)
      Call GetMem('DipZre','Free','Real',ipDzRe,nSize2)
      Call GetMem('DipXsq','Free','Real',ipDxsq,nBNRe**2)
      Call GetMem('DipYsq','Free','Real',ipDysq,nBNRe**2)
      Call GetMem('DipZsq','Free','Real',ipDzsq,nBNRe**2)
      Call GetMem('DipXm','Free','Real',ipDxM,nBNRe**2)
      Call GetMem('DipYm','Free','Real',ipDyM,nBNRe**2)
      Call GetMem('DipZm','Free','Real',ipDzM,nBNRe**2)
      Call GetMem('TEMP','Free','Real',ipTEMP,nBNRe**2)

      Return
      End
