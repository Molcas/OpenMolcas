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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine PrepP()
************************************************************************
*                                                                      *
* Object: to set up the handling of the 2nd order density matrix.      *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              OpnOne                                                  *
*              GetMem                                                  *
*              RdOne                                                   *
*              PrMtrx                                                  *
*              ClsOne                                                  *
*              ErrOne                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '92                                              *
************************************************************************
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "pso.fh"
#include "etwas.fh"
#include "aces_gamma.fh"
#include "mp2alaska.fh"
#include "nsd.fh"
#include "dmrginfo_mclr.fh"
#include "nac.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
************ columbus interface ****************************************
#include "columbus_gamma.fh"
#include "setup.fh"
#include "shinf.fh"
************************************************************************
      Integer nFro(0:7)
      Integer Columbus
      Character*8 RlxLbl,Method, KSDFT*16
      Logical lPrint
      Logical DoCholesky
      Real*8 CoefX,CoefR
      Character*80 Fmt*60
*
*...  Prologue

      iRout = 250
      iPrint = nPrint(iRout)
      lPrint=iPrint.ge.6
      Call qEnter('PrepP')
#ifdef _CD_TIMING_
      Call CWTIME(PreppCPU1,PreppWall1)
#endif
*
      Call StatusLine(' Alaska:',' Prepare the 2-particle matrix')
*
      iD0Lbl=1
      iComp=1
*
      lsa=.False.
      Gamma_On=.False.
      Gamma_mrcisd=.FALSE.
      lPSO=.false.
      Case_2C=.False.
      Case_3C=.False.
      Case_mp2=.False.
*      ip_Z_p_k=ip_Dummy

      nDens = 0
      Do 1 iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
 1    Continue
*
*
*...  Get the method label
      Call Get_cArray('Relax Method',Method,8)
      Call Get_iScalar('Columbus',columbus)
      nCMo = n2Tot
      mCMo = n2Tot
      If (Method.eq. 'KS-DFT  ' .or.
     &    Method.eq. 'MCPDFT  ' .or.
     &    Method.eq. 'CASDFT  ' ) Then
         Call Get_iScalar('Multiplicity',iSpin)
         Call Get_cArray('DFT functional',KSDFT,16)
         Call Get_dScalar('DFT exch coeff',CoefX)
         Call Get_dScalar('DFT corr coeff',CoefR)
         ExFac=Get_ExFac(KSDFT)
         CoulFac=One
      Else
         iSpin=0
         ExFac=One
         CoulFac=One
      End If
*
*...  Check the wave function type
*
*                                                                      *
************************************************************************
*                                                                      *
      If ( Method.eq.'RHF-SCF ' .or.
     &     Method.eq.'UHF-SCF ' .or.
     &     Method.eq.'IVO-SCF ' .or.
     &     Method.eq.'MBPT2   ' .or.
     &     Method.eq.'KS-DFT  ' .or.
     &     Method.eq.'ROHF    ' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ',Method
            If (Method.eq.'KS-DFT  ') Then
               Write (6,'(2A)') ' Functional type:   ',KSDFT
               Fmt = '(1X,A26,20X,F18.6)'
               Write(6,Fmt)'Exchange scaling factor',CoefX
               Write(6,Fmt)'Correlation scaling factor',CoefR
            End If
            Write (6,*)
         End If
         If(Method.eq.'MBPT2   ') Then
            Case_mp2=.true.
            Call DecideOnCholesky(DoCholesky)
            If(.not.DoCholesky) Then
               iSeed = 10
               LuGam = IsFreeUnit(iSeed)
               Write(FnGam,'(A6)') 'LuGam'
               Call DaName_MF_WA(LuGam,FnGam)
               Gamma_on=.True.
               Call Aces_Gamma()
            End If
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'Corr. WF' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,*)
     &         ' Wavefunction type: an Aces 2 correlated wavefunction'
            Write (6,*)
         End If
         Gamma_On=.True.
         Call Aces_Gamma()
       Else if (Method(1:7).eq.'MR-CISD' .and. Columbus.eq.1) then
************ columbus interface ****************************************
*do not reconstruct the two-particle density from the one-particle
*density or partial two-particle densities but simply read them from
*file

        nShell=mSkal
        nPair=nShell*(nShell+1)/2
        nQuad=nPair*(nPair+1)/2

*---- Allocate Table Of Content for half sorted gammas.
*
      Call GetMem('G_Toc','Allo','Real',ipG_Toc,nQuad+2)

*  find free unit number
        lgtoc=61
        lgtoc=isFreeUnit(lgtoc)
        idisk=0
*  read table of contents of half-sorted gamma file
         Call DaName(lgtoc,'gtoc')
         Call ddafile(lgtoc,2,Work(ipg_toc),nQuad+2,idisk)
         Call Daclos(lgtoc)
         n=int(Work(ipg_toc+nQuad))
         lbin=int(Work(ipg_toc+nQuad+1))
         if (n.ne.nQuad) then
           Call WarningMessage(2,'n.ne.nQuad')
           Write (6,*) 'n,nQuad=',n,nQuad
           Call Abend()
         endif

c       open(unit=lgtoc,file='gtoc',form='unformatted')
c        read(lgtoc) n
c        if (n.ne.nquad) stop 'quad error'
c        call read_lgtoc(lgtoc,Work(ipg_toc),n)
c        read(lgtoc) lbin
c       close (lgtoc)
        Gamma_On=.True.
        Gamma_mrcisd=.TRUE.
*       open gamma file
         LuGamma=60
         LuGamma=isfreeunit(LuGamma)
*        closed in closep
         Call DaName_MF(LuGamma,'GAMMA')
*  allocate space for bins
         Call GetMem('Bin','Allo','Real',ipBin,2*lBin)
*  compute SO2cI array
         Call GetMem('SO2cI','Allo','Inte',ipSO2cI,2*nSOs)
         call so2ci(iWork(ipSO2cI),iWork(ipSOsh),nsos)
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'RASSCF  ' .or.
     &          Method.eq.'CASSCF  ' .or.
     &          Method.eq.'GASSCF  ' .or.
     &          Method.eq.'MCPDFT  ' .or.
     &          Method.eq.'DMRGSCF ' .or.
     &          Method.eq.'CASDFT  ') then
*
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
*
         nDSO = nDens
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ', Method
            If (Method.eq.'CASDFT  ' .or. Method.eq.'MCPDFT  ')
     &         Write (6,'(2A)') ' Functional type:   ',KSDFT
            Write (6,*)
         End If
         If (method.eq.'MCPDFT  ') lSA=.true.
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'CASSCFSA' .or.
     &          Method.eq.'DMRGSCFS' .or.
     &          Method.eq.'GASSCFSA' .or.
     &          Method.eq.'RASSCFSA' ) then
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
         nDSO = nDens
         Call Get_iScalar('SA ready',iGo)
         If (iGO.eq.1) lSA=.true.
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint) Then
            Write (6,*)
            If (lSA) Then
               Write (6,'(2A)') ' Wavefunction type: State average ',
     &                            Method(1:6)
            Else
               Write (6,'(2A)') ' Wavefunction type: ', Method
            End If
            Write (6,*)
         End If
         Method='RASSCF  '
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,'Alaska: Unknown wavefuntion type')
         Write (6,*) 'Wavefunction type:',Method
         Write (6,*) 'Illegal type of wave function!'
         Write (6,*) 'ALASKA can not continue.'
         Call Quit_OnUserError()
      End If
*
*...  Read the (non) variational 1st order density matrix
*...  density matrix in AO/SO basis
         nsa=1
         If (lsa) nsa=4
         If ( Method.eq.'MCPDFT  ') nsa=4
!AMS modification: add a fifth density slot
         Call GetMem('D0   ','Allo','Real',ipD0,nDens*nsa+nDens)
         call dcopy_(nDens*nsa+nDens,[0.0d0],0,Work(ipD0),1)
         Call GetMem('DVar ','Allo','Real',ipDVar,nDens*nsa)
         if (.not.gamma_mrcisd) then
         Call Get_D1ao(ipD1ao,length)
         If ( length.ne.nDens ) Then
            Call WarningMessage(2,'PrepP: length.ne.nDens')
            Write (6,*) 'length=',length
            Write (6,*) 'nDens=',nDens
            Call Abend()
         End If
         call dcopy_(nDens,Work(ipD1ao),1,Work(ipD0),1)
         Call Free_Work(ipD1ao)
         endif
*

         Call Get_D1ao_Var(ipD1ao_Var,length)
         call dcopy_(nDens,Work(ipD1ao_Var),1,Work(ipDVar),1)
*        if (gamma_mrcisd) then
*         call dcopy_(nDens,Work(ipD1ao_Var),1,Work(ipD0),1)
*        endif
         Call Free_Work(ipD1ao_Var)
*
         If (Method.eq.'UHF-SCF ' .or.
     &       Method.eq.'ROHF    ' .or.
     &    (Method.eq.'KS-DFT  '.and.iSpin.ne.1) .or.
     &       Method.eq.'Corr. WF' ) Then
            Call Get_D1sao(ipDS,Length)
            Call Get_D1sao_Var(ipDSVar,Length)
         Else
            Call GetMem('DS   ','Allo','Real',ipDS,nDens)
            Call GetMem('DSVar','Allo','Real',ipDSVar,nDens)
            Call FZero(Work(ipDS),nDens)
            Call FZero(Work(ipDSVar),nDens)
         End If
*
*   This is necessary for the ci-lag
*
*
*     Unfold density matrix
*
************ columbus interface ****************************************
*do not modify the effective density matrices and fock matrices
      if (.not.gamma_mrcisd) then
      ij = -1
      Do 10 iIrrep = 0, nIrrep-1
         Do 11 iBas = 1, nBas(iIrrep)
            Do 12 jBas = 1, iBas-1
               ij = ij + 1
               Work(ipDVar +ij) = Half*Work(ipDVar +ij)
               Work(ipD0   +ij) = Half*Work(ipD0   +ij)
               Work(ipDS   +ij) = Half*Work(ipDS   +ij)
               Work(ipDSVar+ij) = Half*Work(ipDSVar+ij)
 12         Continue
            ij = ij + 1
 11      Continue
 10   Continue
      endif

      If (iPrint.ge.99) Then
         RlxLbl='D1AO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[ipD0],Work)
         RlxLbl='D1AO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[ipDVar],Work)
         RlxLbl='DSAO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[ipDS],Work)
         RlxLbl='DSAO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[ipDSVar],Work)
      End If

*
*...  Get the MO-coefficients
************ columbus interface ****************************************
*     not that the columbus mcscf MO coefficients have been written
*     to the RUNFILE !

         If (Method.eq.'UHF-SCF ' .or.
     &       Method.eq.'ROHF    ' .or.
     &    (Method.eq.'KS-DFT  '.and.iSpin.ne.1) .or.
     &       Method.eq.'Corr. WF'      ) Then
            nsa=2
         Else
            nsa=1
            If (lsa) nsa=2
         End If
         Call GetMem('CMO','Allo','Real',ipCMO,nsa*mCMO)
         ipCMOa=ipCMO
         ipCMOb=ipCMO+(nsa-1)*mCMO
         Call Get_CMO(ipTemp,nTemp)
         call dcopy_(nTemp,Work(ipTemp),1,Work(ipCMO),1)
         Call Free_Work(ipTemp)
         If (iPrint.ge.99) Then
            ipTmp1 = ipCMO
            Do iIrrep = 0, nIrrep-1
               Call RecPrt(' CMO''s',' ',
     &                     Work(ipTmp1),nBas(iIrrep),
     &                     nBas(iIrrep))
               ipTmp1 = ipTmp1 + nBas(iIrrep)**2
            End Do
         End If
*
*
*...  Get additional information in the case of a RASSCF wave function
*...  Get the number of inactive, active and frozen orbitals
************ columbus interface ****************************************
*  no need for MRCI gradient
         If (.not.lpso .or. gamma_mrcisd ) Goto 1000
         Call Get_iScalar('nSym',i)
         Call Get_iArray('nIsh',nIsh,i)
         Call Get_iArray('nAsh',nAsh,i)
         Call Get_iArray('nFro',nFro,i)
         If (iPrint.ge.99) Then
            Write (6,*) ' nISh=',nISh
            Write (6,*) ' nASh=',nASh
            Write (6,*) ' nFro=',nFro
         End If
         nAct = 0
         nTst = 0
         Do iIrrep = 0, nIrrep-1
!            write(*,*)"nAsh(iIrrep)",nAsh(iIrrep)  ! yma
            nAct = nAct + nAsh(iIrrep)
            nTst = nTst + nFro(iIrrep)
         End Do
         If (nTst.ne.0) Then
            Call WarningMessage(2,
     &                  '; No frozen orbitals are allowed!'//
     &                  '; ALASKA can not continue;')
            Call Quit_OnUserError()
         End If
*
*...  Get the one body density for the active orbitals
*     (not needed for SA-CASSCF)
         nG1 = nAct*(nAct+1)/2
         nsa=1
         If (lsa) nsa=0
         Call GetMem(' G1 ','Allo','Real',ipG1,nG1*nsa)
         If (nsa.gt.0) Then
            Call Get_D1MO(ipTemp,nTemp)
            call dcopy_(nTemp,Work(ipTemp),1,Work(ipG1),1)
            Call Free_Work(ipTemp)
            If (iPrint.ge.99) Call TriPrt(' G1',' ',Work(ipG1),nAct)
         End If
*
*...  Get the two body density for the active orbitals
         nG2 = nG1*(nG1+1)/2
         nsa=1
         if (lsa) nsa=2
         Call GetMem(' G2 ','Allo','Real',ipG2,nG2*nsa)
!       write(*,*) 'got the 2rdm, Ithink.'
         if(Method.eq.'MCPDFT  ') then
           Call Get_P2MOt(ipTemp,nTemp)!PDFT-modified 2-RDM
         else
           Call Get_P2MO(ipTemp,nTemp)
         end if
         call dcopy_(nTemp,Work(ipTemp),1,Work(ipG2),1)
         Call Free_Work(ipTemp)
         If (iPrint.ge.99) Call TriPrt(' G2',' ',Work(ipG2),nG1)
         If (lsa) Then

*  CMO1 Ordinary CMOs
*
*  CMO2 CMO*Kappa
*
           Call Get_LCMO(ipLCMO,Length)
           call dcopy_(Length,Work(ipLCMO),1,Work(ipCMO+mcmo),1)
           Call Free_Work(ipLCMO)
           If (iPrint.ge.99) Then
            ipTmp1 = ipCMO+mcmo
            Do iIrrep = 0, nIrrep-1
               Call RecPrt('LCMO''s',' ',
     &                     Work(ipTmp1),nBas(iIrrep),
     &                     nBas(iIrrep))
               ipTmp1 = ipTmp1 + nBas(iIrrep)**2
            End Do
           End If
*
* P are stored as
*                            _                     _
*   P1=<i|e_pqrs|i> + sum_i <i|e_pqrs|i>+<i|e_pqrs|i>
*   P2=sum_i <i|e_pqrs|i>
*
           Call Get_PLMO(ipPLMO,Length)
           call dcopy_(Length,Work(ipPLMO),1,Work(ipG2+ng2),1)
           ndim1=0
           if(doDMRG)then
             ndim0=0  !yma
             do i=1,8
               ndim0=ndim0+LRras2(i)
             end do
             ndim1=(ndim0+1)*ndim0/2
             ndim2=(ndim1+1)*ndim1/2
             do i=1,ng2
               if(i.gt.ndim2)then
                 Work(ipG2+ng2+i-1)=0.0d0
               end if
             end do
           end if
           Call Free_Work(ipPLMO)
           Call Daxpy_(ng2,One,Work(ipG2+ng2),1,Work(ipG2),1)
           If(iPrint.ge.99)Call TriPrt(' G2L',' ',Work(ipG2+ng2),nG1)
           If(iPrint.ge.99)Call TriPrt(' G2T',' ',Work(ipG2),nG1)
*
           Call Get_D2AV(ipD2AV,Length)
           If (ng2.ne.Length) Then
              Call WarningMessage(2,'Prepp: D2AV, ng2.ne.Length!')
              Write (6,*) 'ng2,Length=',ng2,Length
! DMRG with the reduced AS
              if(doDMRG)then
                Length=ng2 ! yma
              else
                Call Abend()
              end if
           End If
           call dcopy_(Length,Work(ipD2AV),1,Work(ipG2+ng2),1)
           Call Free_Work(ipD2AV)
           If (iPrint.ge.99) Call TriPrt('G2A',' ',Work(ipG2+ng2),nG1)
*
*
*  Densities are stored as:
*
*       ipd0 AO:
*
*       D1 = inactive diagonal density matrix
*                                _                 _
*       D2 = <i|E_pq|i> + sum_i <i|E_pq|i>+<i|E_pq|i> + sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> - 1/2 D2
*
*       D3 = sum_i <i|E_pq|i> (active)
*
*       D4 = sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> (inactive)
*
*       G1 = <i|e_ab|i>
*       G2 = sum i <i|e_ab|i>
*
!************************
         RlxLbl='D1AO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,ipD0,Work)


           Call Getmem('TMP','ALLO','REAL',ipT,2*ndens)
           Call Get_D1I(Work(ipCMO),Work(ipD0+0*ndens),Work(ipT),
     &                  nish,nbas,nIrrep)
           Call Getmem('TMP','FREE','REAL',ipT,2*ndens)

!************************
         RlxLbl='D1AO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,ipD0,Work)
*
           Call dcopy_(ndens,Work(ipDVar),1,Work(ipD0+1*ndens),1)
           If (.not.isNAC) call daxpy_(ndens,-Half,Work(ipD0+0*ndens),1,
     &                                             Work(ipD0+1*ndens),1)
!         RlxLbl='D1COMBO  '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,ipD0+1*ndens,Work)
*
*   This is necessary for the kap-lag
*
           Call Get_D1AV(ipD1AV,Length)
           nG1 = nAct*(nAct+1)/2
           If (ng1.ne.Length) Then
              Call WarningMessage(2,'PrepP: D1AV, nG1.ne.Length!')
              Write (6,*) 'nG1,Length=',nG1,Length
              if(doDMRG)then ! yma
                If (length.ne.ng1) length=nG1
              else
                Call Abend()
              end if
              Call Abend()
           End If
           Call Get_D1A(Work(ipCMO),Work(ipD1AV),Work(ipD0+2*ndens),
     &                 nIrrep,nbas,nish,nash,ndens)
           Call Free_Work(ipD1AV)
!************************
         RlxLbl='D1AOA   '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,ipD0+2*ndens,Work)
*
           Call Get_DLAO(ipDLAO,Length)
           call dcopy_(Length,Work(ipDLAO),1,Work(ipD0+3*ndens),1)

!ANDREW - modify D2: should contain only the correction pieces

         If ( Method.eq.'MCPDFT  ') then
!Get the D_theta piece
         Call Get_D1ao(ipD1ao,ndens)
      ij = -1
      Do  iIrrep = 0, nIrrep-1
         Do iBas = 1, nBas(iIrrep)
            Do jBas = 1, iBas-1
               ij = ij + 1
               Work(ipD1ao+ij) = Half*Work(ipD1ao+ij)
            end do
            ij = ij + 1
          end do
      end do
          call daxpy_(ndens,-1d0,Work(ipD0+0*ndens),1,
     &                                             Work(ipD1ao),1)
!          write(*,*) 'do they match?'
!          do i=1,ndens
!            write(*,*) Work(ipd1ao-1+i),Work(ipD0+2*ndens+i-1)
!          end do

          call daxpy_(ndens,-Half,Work(ipD0+0*ndens),1,
     &                                             Work(ipD0+1*ndens),1)
          call daxpy_(ndens,-1.0d0,Work(ipD1ao),1,
     &                                             Work(ipD0+1*ndens),1)
!ANDREW - Generate new D5 piece:
          call dcopy_(ndens,[0.0d0],0,
     &                                             Work(ipD0+4*ndens),1)
          call daxpy_(ndens,0.5d0,Work(ipD0+0*ndens),1,
     &                                             Work(ipD0+4*ndens),1)
          call daxpy_(ndens,1.0d0,Work(ipD1ao),1,
     &                                             Work(ipD0+4*ndens),1)
         Call Free_Work(ipD1ao)
          end if


!          call dcopy_(ndens*5,0.0d0,0,
!     &                                             Work(ipD0),1)
!          call dcopy_(nG2,0.0d0,0,
!     &                                             Work(ipG2),1)


           Call Free_Work(ipDLAO)
!************************
           !Call dscal_(Length,0.5d0,Work(ipD0+3*ndens),1)
           !Call dscal_(Length,0.0d0,Work(ipD0+3*ndens),1)

         RlxLbl='DLAO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,ipD0+3*ndens,Work)
! DMRG with the reduced AS
           if(doDMRG)then
             length=ndim1  !yma
           end if
         End If
         If (iPrint.ge.99) Call TriPrt(' G2',' ',Work(ipG2),nG1)
*
*...  Close 'RELAX' file
1000     Continue
*
*...  Epilogue, end
#ifdef _CD_TIMING_
      Call CWTIME(PreppCPU2,PreppWall2)
      Prepp_CPU  = PreppCPU2 - PreppCPU1
      Prepp_Wall = PreppWall2 - PreppWall1
#endif
      Call qExit('PrepP')

      Return
      End

      Subroutine Get_D1I(CMO,D1It,D1I,nish,nbas,nsym)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Dimension CMO(*), D1I(*),D1It(*)
      Integer nbas(nsym),nish(nsym)

      iOff1 = 0
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nIsh(iSym)
        If ( iBas.ne.0 ) then
          iOff2 = iOff1
          Do i = 1,iBas
            Do j = 1,iBas
              Sum = Zero
              Do k = 0,iOrb-1
                Sum = Sum + Two * CMO(iOff1+k*iBas+i)
     &                          * CMO(iOff1+k*iBas+j)
              End Do
              D1I(iOff2 + j) = Sum
            End Do
            iOff2 = iOff2 + iBas
          End Do
          iOff1 = iOff1 + iBas*iBas
        End If
      End Do
      Call Fold2(nsym,nBas,D1I,D1It)
      Return
      End

      Subroutine Get_D1A(CMO,D1A_MO,D1A_AO,
     &                    nsym,nbas,nish,nash,ndens)


      Implicit Real*8 (A-H,O-Z)

      Dimension CMO(*) , D1A_MO(*) , D1A_AO(*)
      Integer nbas(nsym),nish(nsym),nash(nsym)

#include "WrkSpc.fh"
#include "real.fh"
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Call qEnter('Get_D1A')

      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      ii=0
      Call GetMem('Scr1','Allo','Real',ip1,2*ndens)
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Call dCopy_(iBas*iBas,[Zero],0,Work(ip1-1+iOff3),1)
        If ( iAsh.ne.0 ) then
          Call GetMem('Scr1','Allo','Real',iTmp1,iAsh*iAsh)
          Call GetMem('Scr2','Allo','Real',iTmp2,iAsh*iBas)
          ij=0
          Do i=1,iash
           Do j=1,iash
            Work(iTmp1+ij)=D1A_MO(itri(i+ii,j+ii))
            ij=ij+1
           End do
          End do
          ii=ii+iash
          Call DGEMM_('N','T',
     &                iBas,iAsh,iAsh,
     &                One,CMO(iOff2+iIsh*iBas),iBas,
     &                Work(iTmp1),iAsh,
     &                Zero,Work(iTmp2),iBas)
          Call DGEMM_('N','T',
     &                iBas,iBas,iAsh,
     &                One,Work(iTmp2),iBas,
     &                CMO(iOff2+iIsh*iBas),iBas,
     &                Zero,Work(ip1-1+iOff3),iBas)
          Call GetMem('Scr2','Free','Real',iTmp2,iAsh*iBas)
          Call GetMem('Scr1','Free','Real',iTmp1,iAsh*iAsh)
        End If
        iOff1 = iOff1 + (iAsh*iAsh+iAsh)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + iBas*iBas
      End Do
      Call Fold2(nsym,nBas,work(ip1),D1A_AO)
      Call GetMem('Scr1','FREE','Real',ip1,ndens)
      Call qExit('Get_D1A')
      Return
      End


      Subroutine Fold2(nSym,nBas,A,B)

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =   A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      end
************ columbus interface ****************************************
*read table of contents for gamma file

        subroutine read_lgtoc(lgtoc,gtoc,n)
        integer n,lgtoc
        real*8 gtoc(n)
          read(lgtoc) gtoc
        return
         end
