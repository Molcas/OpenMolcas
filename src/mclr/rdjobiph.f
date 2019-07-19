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
      Subroutine RdJobIph
************************************************************************
*                                                                      *
*     Read the contents of the JOBIPH file.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Files_mclr.fh"
#include "glbbas_mclr.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
      Character*72 Line
      Character*8 Method
      real*8 dv_ci2  ! yma added
      Logical Found
      Dimension rdum(1)
*                                                                      *
************************************************************************
*                                                                      *
*     itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
      debug=.FALSE.
#ifdef _DEBUG_
      debug=.TRUE.
#endif

*                                                                      *
************************************************************************
*                                                                      *
*     Save the ROOT input parameter                                    *
*----------------------------------------------------------------------*
      kRoots=lRoots
*----------------------------------------------------------------------*
*     Read the table of disk adresses                                  *
*----------------------------------------------------------------------*
      Call DaName(LuJob,FnJob)
      iDisk=0
      Call iDaFile(LuJob,2,iToc,iTOCIPH,iDisk)
*----------------------------------------------------------------------*
*     Read the the system description                                  *
*----------------------------------------------------------------------*
      ndum=lenin8*mxorb
      Call Getmem('TMP','ALLO','CHAR',ipdum,ndum)
      iDisk=iToc(1)

!      write(*,*)"if dmrg, it should be something else "
      Call WR_RASSCF_Info(LuJob,2,iDisk,
     &                    nActEl,iSpin,nSym,State_sym,nFro,
     &                    nIsh,nAsh,nDel,
     &                    nBas,MxSym,cwork(ipdum),LENIN8*mxorb,
     &                    nConf,HeaderJP,144,
     &                    TitleJP,4*18*mxTit,PotNuc0,lRoots,
     &                    nRoots,iRoot,mxRoot,
     &                    nRs1,nRs2,nRs3,
     &                    nHole1,nElec3,iPt2,Weight)

      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
        call dmrg_spc_change_mclr(LRras2(1:8),nrs2)
      end if
!      do i=1,8
!        write(*,*)i,"-irrep",nIsh(i),nAsh(i),nRs1(i),nRs2(i),nRs3(i)
!        call xflush(6)
!      end do

      Call Getmem('TMP','FREE','CHAR',ipdum,ndum)
*----------------------------------------------------------------------*
*     Overwrite the variable lroots if approriate, i.e if lroot        *
*     was set by input.                                                *
*----------------------------------------------------------------------*
      If ( kRoots.ne.-1 ) then
         If ( iPt2.ne.0 ) then
            Write (6,*) 'RdJobiph: kRoots.ne.-1 .and. iPt2.ne.0'
            Call Abend()
         Else if ( kRoots.gt.lRoots ) then
            Write (6,*) 'RdJobiph: kRoots.ne.-1 .and. kRoots.gt.lRoots'
            Call Abend()
         End If
         lRoots=kRoots
         nRoots=1
      End If
*----------------------------------------------------------------------*
*     Precompute the total sum of variables and size of matrices       *
*----------------------------------------------------------------------*
      ntIsh=0
      ntItri=0
      ntIsqr=0
      ntAsh=0
      ntAtri=0
      ntAsqr=0
      ntBas=0
      ntBtri=0
      ntBsqr=0
      nna=0
      Do 10 iSym=1,nSym
         norb(isym)=nbas(isym)-ndel(isym)
         ntIsh=ntIsh+nIsh(iSym)
         ntItri=ntItri+nIsh(iSym)*(nIsh(iSym)+1)/2
         ntIsqr=ntIsqr+nIsh(iSym)*nIsh(iSym)
         ntAsh=ntAsh+nAsh(iSym)
         ntAtri=ntAtri+nAsh(iSym)*(nAsh(iSym)+1)/2
         ntAsqr=ntAsqr+nAsh(iSym)*nAsh(iSym)
         ntBas=ntBas+nBas(iSym)
         ntBtri=ntBtri+nBas(iSym)*(nBas(iSym)+1)/2
         ntBsqr=ntBsqr+nBas(iSym)*nBas(iSym)
         nA(iSym)=nna
         nnA=nnA+nAsh(isym)
10    Continue


!>  Generate the nr. of csf in each sub-sym, used in geom-opt with SA DMRG-SCF
      if(doDMRG)then  ! yma
        Call GugaCtl_dmrg()   ! generate the Nr. of csfs in each sym
!        do isym=1,8
!          write(*,*)"isym_ncsf in rdjobiph ",ncsf(isym)
!        end do
      end if

*
*----------------------------------------------------------------------*
*     Load the orbitals used in the last macro iteration               *
*----------------------------------------------------------------------*
*
      Call Get_CMO(ipCMO,Length)
C
C     Read state for geo opt
C
      Call Get_iScalar('Relax CASSCF root',irlxroot)
      Call Get_cArray('Relax Method',Method,8)
      iMCPD=.False.
      if(Method.eq.'MCPDFT  ') then
        iMCPD=.True.
        Do i=1,lroots
          if(iroot(i).eq.irlxroot)istate=i
        end do
      end if
      If (Method.eq.'CASSCFSA') Then
         Call Get_iScalar('SA ready',iGo)
         If (iGO.eq.-1) Then
            Write (6,*) 'MCLR not implemented for SA-CASSCF'//
     &                  ' with non-equivalent weights!'
            Call Abend()
         Else
            If (iGo.ne.2) SA=.true.
            Found=.true.
            If (override) Then
               If (isNAC) Then
                  Do j=1,2
                    NSSA(j)=0
                    Do i=1,lroots
                       If (iroot(i).eq.NACStates(j)) NSSA(j)=i
                    End Do
                    If (NSSA(j).eq.0) Found=.false.
                  End Do
               Else
                  irlxroot=iroot(istate)
               End If
            Else
               istate=0
               Do i=1,lroots
                  If (iroot(i).eq.irlxroot) istate=i
               End Do
               If (istate.eq.0) Found=.false.
            End If
            If (.not.Found) Then
               Call WarningMessage(2,
     &              'Cannot relax a root not included in the SA')
            End If
         End If
      Else If (irlxroot.eq.1.and..Not.(McKinley.or.PT2.or.iMCPD)) Then
         Write (6,*)
         Write (6,*) 'W A R N I N G !'
         Write (6,*)
         Write (6,*) 'Redundant rlxroot input in RASSCF!'
         Write (6,*) 'I''ll sign off here without a clean termination!'
         Write (6,*) 'However, I have to fix the epilogue file.'
         Write (6,*)
         irc=-1
         iopt=1
         Call OPNMCK(irc,iopt,FNMCK,LUMCK)
         Call WrMck(iRC,iOpt,'nSym',1,nBas,iDummer)
         Call ClsFls_MCLR()
         Call Finish(0)
      End if
C
*     iDisk=iToc(9)
*     IF(IPT2.EQ.0) iDisk=iToc(2)
*     Call dDaFile(LuJob,2,Work(ipCMO),ntBsqr,iDisk)
      If( .false. ) then
         jpCMO=ipCMO
         Do 15 iSym=1,nSym
            call dcopy_(nbas(isym)*ndel(isym),[0d0],0,
     *                 Work(jpCMO+norb(isym)*nbas(isym)),1)
            Write(Line,'(A,i2.2)') 'MO coefficients, iSym = ',iSym
            Call RecPrt(Line,' ',Work(jpCMO),nBas(iSym),nBas(iSym))
            jpCMO=jpCMO+nBas(iSym)*nBas(iSym)
15       Continue
      End If
*----------------------------------------------------------------------*
*     Load the CI vectors for the SA roots                             *
*----------------------------------------------------------------------*

! If doDMRG, introducing CI coeffieients later :
!    1) only coefficients of importants DETs using MPS2CI
!    2) and together with DET numbers from GUGA generation part

      if(doDMRG)then  ! yma
        !ipCI=1
        ! need to be deleted the 1000 yma
        Call GetMem('OCIvec','Allo','Real',ipCI,nConf*nroots+1000)
      ! So, nothing need to be done here ...
      else
        Call GetMem('OCIvec','Allo','Real',ipCI,nConf*nroots)
        Do i=1,nroots
          j=iroot(i)
          iDisk=iToc(4)
          Do k=1,j-1
            Call dDaFile(LuJob,0,rdum,nConf,iDisk)
          End Do
          Call dDaFile(LuJob,2,Work(ipCI+(i-1)*nconf),nConf,iDisk)
        End Do
!#ifdef _DEBUG_           ! yma umcomment
        Do i=0,nroots-1            !yma
          inum=0
          dv_ci2=0.0d0
          do j=1,nconf
!yma        CI-threshold
            if(abs(Work(ipCI+nconf*i+j-1)).lt.0.0d0)then
              inum=inum+1
              Work(ipCI+nconf*i+j-1)=0.0d0
            else
              dv_ci2=dv_ci2+Work(ipCI+nconf*i+j-1)**2
            end if
          end do
!          Call DVcPrt('CI coefficients',' ',Work(ipCI+nconf*i),nConf)!yma
!          write(*,*)"dismissed dets num", inum
!          write(*,*)"absolutely CI^2",dv_ci2
        End DO
      End If

!#endif
*----------------------------------------------------------------------*
*     Load state energy                                                *
*----------------------------------------------------------------------*
      Call GetMem('Temp2','Allo','Real',ipTmp2,mxRoot*mxIter)
      iDisk=iToc(6)
#ifdef _DEBUG_
      If (debug) Then
         Write(6,*) 'NROOTS: ',nroots
         Write(6,*) 'iROOTS: ',(iroot(i),i=1,nroots)
         Write(6,*) 'lROOTS: ',lroots
      End If
#endif
      Call dDaFile(LuJob,2,Work(ipTmp2),mxRoot*mxIter,iDisk)

      Do  iter=0,mxIter-1
        Do i=1,nroots
          j=iroot(i)
          ! It should be 0.0d0 in DMRG case
          Temp=Work(ipTmp2+iter*mxRoot+j-1)
          If ( Temp.ne.0.0D0 ) ERASSCF(i)=Temp
*          If (debug) Write(*,*) ERASSCF(i),i
         End Do
      End Do

#ifdef _DEBUG_
      If (debug) Then
          Write(6,*) (Work(ipTmp2+i),i=0,lroots)
          Write(6,*)'RASSCF energies=',(ERASSCF(i),i=1,nroots)
      End If
#endif
      Call GetMem('Temp2','Free','Real',ipTmp2,mxRoot*100)
*
      nAct  = 0    ! 1/2

      if(doDMRG)then  ! yma
        call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if

      nAct2 = 0
      nAct4 = 0
      Do iSym = 1, nSym
         nAct = nAct + nAsh(iSym)
        nAct2=nAct2+nAsh(iSym)**2
      End Do
      Do iS = 1, nSym
         Do jS = 1, nSym
            Do kS = 1, nSym
              lS=iEOr(iEOr(is-1,js-1),ks-1)+1
              nAct4=nAct4+nAsh(iS)*nAsh(jS)*nAsh(kS)*nAsh(lS)
            End Do
         End Do
      End Do
*     Call GetMem(' G1sq','Allo','Real',ipG1,nAct2)

      nG1 = nAct*(nAct+1)/2
      Call GetMem(' G1 ','Allo','Real',ipG1t,nG1)
      nG2=nG1*(nG1+1)/2
      Call GetMem(' G2 ','Allo','Real',ipG2t,nG2)
      Call RDDENS(Work(ipG1t),ng1,Work(ipG2t),ng2)
      ipg1=ipg1t
      ipG2=ipG2t
      ipg2tmm=ipg2
      ipg2tpp=ipg2

#ifdef _DEBUG_
      Call Triprt('G1',' ',Work(ipG1t),ntash)
      Call Triprt('G2',' ',Work(ipG2),ng1)
#endif

      nG1 = nAct*(nAct+1)/2
      nG2=nG1*(nG1+1)/2

!      Call Triprt('G1',' ',Work(ipG1t),ntash)
!      Call Triprt('G2',' ',Work(ipG2),ng1)

      if(doDMRG)then ! yma
        call dmrg_dim_change_mclr(LRras2(1:8),nact,0)
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
        call dmrg_spc_change_mclr(LRras2(1:8),nrs2)
      end if

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
