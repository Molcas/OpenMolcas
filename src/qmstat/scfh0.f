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
*-- Read all MO-transformed integrals and construct zeroth Fock-matrix.
*
      Subroutine ScfH0(nBas)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "tratoc.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension nBasM(MxSym),nOrbM(MxSym),nDelM(MxSym),nFroM(MxSym)
      Dimension nBas(MxSym)
      Dimension iToc(64)
      Character NameM*40000,firstind*10

*
*-- Wilkommen.
*
      Write(6,*)
      Write(6,*)
      Write(6,*)'Reading MO-transformed integrals. Zero:th hamiltonian '
     &//'constructed.'

*
*-- Numbers and files.
*
      nSize=iOrb(1)*(iOrb(1)+1)/2
      Call GetMem('FockM','Allo','Real',iPointF,nSize)
      Call GetMem('SUPER','Allo','Real',iSupM,nSize**2)
      iLu1=56
      iLu2=58
      Call DaName(iLu1,'TRAONE')
      Call DaName(iLu2,'TRAINT')
      iDisk=0
*--- This is special utility to read header of TRAONE.
      Call Wr_Motra_Info(iLu1,2,iDisk,iToc,64,Ecor,nSymM,nBasM,nOrbM
     &,nFroM,nDelM,MxSym,NameM,4*2*5000) !Last argument depends on
                                         !mxorb in Molcas.fh.
      nOrbMFirst=nOrbM(1)

*
*-- One checks.
*
      If(nBasM(1).ne.nBas(1)) then
        Write(6,*)
        Write(6,*)'  ERROR! Conflict between one-electron file and MO-t'
     &//'ransformed one-electron file.'
        Write(6,*)'         nBas=',nBas(1),' MO-nBas=',nBasM(1)
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

*
*-- Read one-electron matrix elements.
*
      iDisk=iToc(2)
      Call dDaFile(iLu1,2,Work(iPointF),nSize,iDisk)
      call dcopy_(nSize,Work(iPointF),iONE,HHmat,iONE)
      Call GetMem('FockM','Free','Real',iPointF,nSize)
      Call DaClos(iLu1)

*
*-- Add external perturbation if requested.
*
      If(AddExt) then
        Write(6,*)'    -- Adding external perturbation.'
        nBTri=nBas(1)*(nBas(1)+1)/2
        Lu_One=49
        Lu_One=IsFreeUnit(Lu_One)
        Call OpnOne(irc,0,'ONEINT',Lu_One)
        Call GetMem('AOExt','Allo','Real',ipAOx,nBTri+4)
        Call GetMem('TEMP','Allo','Real',iTEMP,nBas(1)*iOrb(1))
        Call GetMem('Final','Allo','Real',iFine,iOrb(1)**2)
        Call GetMem('Squared','Allo','Real',iSqAO,nBas(1)**2)
        Call GetMem('MOExt','Allo','Real',ipMOx,nSize)
        Do 9901, iExt=1,nExtAddOns
          irc=-1
          iopt=0
          iSmLbl=0
          Call RdOne(irc,iopt,ExtLabel(iExt),iCompExt(iExt),Work(ipAOx)
     &              ,iSmLbl)
          Call DScal_(nBTri,ScalExt(iExt),Work(ipAOx),iONE)
          If(irc.ne.0) then
            Write(6,*)
            Write(6,*)'ERROR when reading ',ExtLabel(iExt),'.'
            Write(6,*)'Have Seward computed this integral?'
            Call Quit(_RC_IO_ERROR_READ_)
          Endif
          Call Square(Work(ipAOx),Work(iSqAO),iONE,nBas(1),nBas(1))
          Call Dgemm_('T','N',iOrb(1),nBas(1),nBas(1),ONE,Work(iV1)
     &              ,nBas(1),Work(iSqAO),nBas(1),ZERO,Work(iTEMP)
     &              ,iOrb(1))
          Call Dgemm_('N','N',iOrb(1),iOrb(1),nBas(1),ONE,Work(iTEMP)
     &              ,iOrb(1),Work(iV1),nBas(1),ZERO,Work(iFine),iOrb(1))
          Call SqToTri_Q(Work(iFine),Work(ipMOx),iOrb(1))
          Call DaxPy_(nSize,ONE,Work(ipMOx),iONE,HHmat,iONE)
9901    Continue
        Call GetMem('AOExt','Free','Real',ipAOx,nBTri+4)
        Call GetMem('TEMP','Free','Real',iTEMP,nBas(1)*iOrb(1))
        Call GetMem('Final','Free','Real',iFine,iOrb(1)**2)
        Call GetMem('Squared','Free','Real',iSqAO,nBas(1)**2)
        Call GetMem('MOExt','Free','Real',ipMOx,nSize)
        Call ClsOne(irc,Lu_One)
      Endif

*
*-- Now to the two-electron matrix elements.
*
      iDisk=0
      Call iDaFile(iLu2,2,iTraToc,nTraToc,iDisk)
      iDisk=iTraToc(1)
      iSup=0
*--- Ooohhhh, lets get crude! Read ALL (yes, you read right) integrals
*    and then order them.
      nBuf1=nOrbM(1)*(nOrbM(1)+1)/2
      nBuf2=nBuf1*(nBuf1+1)/2
*
*--- Lets check if this construct is possible. If not advise user what
*     to do.
*
      Call GetMem('MAX','Max','Real',iDum,nMAX)
      If(nMAX.lt.(nBuf2+nBuf1**2)) then
        Write(6,*)
        Write(6,*)'  Too many MO-transformed two-electron integrals'
     &          //' from Motra. Do you need all?'
        Write(6,*)'  If not, then use the DELEte keyword in Motra'
     &          //' to remove the superfluous ones.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif
*
*--- Proceed!
*
      Call GetMem('Buffer','Allo','Real',iBuff,nBuf2)
      Call GetMem('Temporary','Allo','Real',iTEMP,nBuf1**2)
      Call dDaFile(iLu2,2,Work(iBuff),nBuf2,iDisk)
      Do 311, i=1,nBuf1
        Do 312, j=i,nBuf1
          iSup=iSup+1
          If(i.le.nSize.and.j.le.nSize) then
            Work(iTEMP+(i-1)*nBuf1+j-1)=Work(iBuff+iSup-1)
            Work(iTEMP+(j-1)*nBuf1+i-1)=Work(iBuff+iSup-1)
          Endif
312     Continue
311   Continue

*
*   and sees to that right numbers get in the right place.
*
      Do 321, i=1,iOrb(1)
        Do 322, j=1,i
          Do 323, k=1,i
            llmax=k
            If(i.eq.k) llmax=j
            Do 324, l=1,llmax
              ij=ipair_qmstat(i,j)
              ik=ipair_qmstat(i,k)
              il=ipair_qmstat(i,l)
              jk=ipair_qmstat(j,k)
              jl=ipair_qmstat(j,l)
              kl=ipair_qmstat(k,l)
            Work(iSupM+nSize*(ij-1)+kl-1)=Work(iTEMP+(ij-1)*nBuf1+kl-1)
     &-(Work(iTEMP+(ik-1)*nBuf1+jl-1)+Work(iTEMP+(il-1)*nBuf1+jk-1))/4
            Work(iSupM+nSize*(kl-1)+ij-1)=Work(iSupM+nSize*(ij-1)+kl-1)
            Work(iSupM+nSize*(ik-1)+jl-1)=Work(iTEMP+(ik-1)*nBuf1+jl-1)
     &-(Work(iTEMP+(ij-1)*nBuf1+kl-1)+Work(iTEMP+(il-1)*nBuf1+jk-1))/4
            Work(iSupM+nSize*(jl-1)+ik-1)=Work(iSupM+nSize*(ik-1)+jl-1)
            Work(iSupM+nSize*(il-1)+jk-1)=Work(iTEMP+(il-1)*nBuf1+jk-1)
     &-(Work(iTEMP+(ik-1)*nBuf1+jl-1)+Work(iTEMP+(ij-1)*nBuf1+kl-1))/4
            Work(iSupM+nSize*(jk-1)+il-1)=Work(iSupM+nSize*(il-1)+jk-1)
324         Continue
323       Continue
322     Continue
321   Continue
      Call GetMem('Buffer','Free','Real',iBuff,nBuf2)
      Call GetMem('Temporary','Free','Real',iTEMP,nBuf1**2)
      Call DaClos(iLu2)

*
*-- Serious amount of printing!
*
      If(iPrint.ge.35) then
        Write(6,*)
        Write(6,*)'The Super Matrix in all its divine g(l)ory:'
        kaunter=0
        Do 331, i=1,iOrb(1)
          Do 332, j=1,i
            Write(firstind,'(I3,A,I3)')i,',',j
            Call TriPrt(firstind,' ',Work(iSupM+nSize*kaunter),iOrb(1))
            kaunter=kaunter+1
332       Continue
331     Continue
        Write(6,*)'Super Matrix End.'
      Endif
      Write(6,*)'...Done!'

*
*-- The end.
*
      Return
      End


*
*-- Little bastard.
*
      Integer Function iPair_qmstat(a,b)
      Implicit Integer (a-z)
      iPair_qmstat=(Max(a,b)*(Max(a,b)-1))/2+Min(a,b)
      Return
      End
