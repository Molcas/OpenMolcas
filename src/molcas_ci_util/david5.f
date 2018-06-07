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
      Subroutine David5(nDet,mxItr,nItr,CI_Conv,ThrEne,
     &                 iSel,ExplE,ExplV,HTUTRI,GTUVXTRI)
      use citrans
      use faroald
      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasrc.fh"
#include "rasscf.fh"
#include "general.fh"
#include "csfbas.fh"
#include "gugx.fh"
#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "rasscf_lucia.fh"
#include "output_ras.fh"

* lroots, maxjt in rasscf.fh
      Dimension CI_conv(2,lRoots,MAXJT)
* nsel in general.fh
      Dimension iSel(nSel), ExplE(nSel), ExplV(nSel,nSel)

      PARAMETER (ROUTINE='DAVID5  ')
      Dimension Alpha(mxRoot), Beta(mxRoot)
      Dimension HTUTRI(*), GTUVXTRI(*)
      real*8, allocatable :: sgm(:,:), psi(:,:)
      real*8, allocatable :: htu(:,:), gtuvx(:,:,:,:)
      integer nnew(mxkeep)
*-------------------------------------------------------------------
*MGD dec 2017 : When optimizing many states, the lowest ones tend to
*converge much faster than the rest. Changed the code so that the converged states
*are not optimize further, saving potentially a lot of time.

      if (DoFaro) then
        ! fill in the integrals from their triangular storage
        allocate(htu(my_norb,my_norb))
        allocate(gtuvx(my_norb,my_norb,my_norb,my_norb))
        htu=0.0d0
        gtuvx=0.0d0
        itu=0
        ituvx=0
        do it=1,my_norb
          do iu=1,it
            itu=itu+1
*         write(6,'(1x,3I4,F21.14)') it,iu,itu,htutri(itu)
            htu(iu,it) = htutri(itu)
            htu(it,iu) = htutri(itu)
            do iv=1,it
              ixmax=iv
              if (it==iv) ixmax=iu
              do ix=1,ixmax
                ituvx=ituvx+1
*             write(6,'(1x,5I4,F21.14)') it,iu,iv,ix,ituvx,
*    &         gtuvxtri(ituvx)
                GTUVX(IT,IU,IV,IX)=GTUVXTRI(ITUVX)
                GTUVX(IU,IT,IV,IX)=GTUVXTRI(ITUVX)
                GTUVX(IT,IU,IX,IV)=GTUVXTRI(ITUVX)
                GTUVX(IU,IT,IX,IV)=GTUVXTRI(ITUVX)
                GTUVX(IV,IX,IT,IU)=GTUVXTRI(ITUVX)
                GTUVX(IX,IV,IT,IU)=GTUVXTRI(ITUVX)
                GTUVX(IV,IX,IU,IT)=GTUVXTRI(ITUVX)
                GTUVX(IX,IV,IU,IT)=GTUVXTRI(ITUVX)
              end do
            end do
          end do
        end do
        ! Euhm, stuff needed for awkward conversions from a
        ! non-specified SYG to GUGA format befor converting to
        ! determinants. This is because for Lucia, CSFs have been
        ! converted to SYG format somewhere up in cistart.
        CALL GetMem('CIVEC','Allo','Real',IVECSVC, nconf)
        call getmem('kcnf','allo','inte',ivkcnf,nactel)
      end if

      Call qEnter('David')
      Call Timing(Alfex_1,Swatch,Swatch,Swatch)
      Rc_CI = 0
      IPRLEV=IPRLOC(3)

* allocate space for CI-vectors
      l1 = lRoots*mxKeep
      l2 = l1*l1
      l3 = (l2+l1)/2
C Trying to avoid writing out of bound in CSDTVC :::: JESPER :::: CHEAT
      Call GetMem('Vector1','Allo','Real',iVec1,ndet)
      Call GetMem('Vector2','Allo','Real',iVec2,ndet)
      Call GetMem('Vector3','Allo','Real',iVec3,ndet)
      Call GetMem('Esmall','Allo','Real',iEs,l1)
      Call GetMem('Hsmall','Allo','Real',iHs,l3)
      Call GetMem('Ssmall','Allo','Real',iSs,l3)
      Call GetMem('Csmall','Allo','Real',iCs,l2)
      Call GetMem('Scr1','Allo','Real',iScr1,l2)
      Call GetMem('Scr2','Allo','Real',iScr2,l2)
      Call GetMem('Scr3','Allo','Real',iScr3,lRoots*nSel)
      Call GetMem('Scr4','Allo','Real',iScr4,lRoots*nSel)
      Call GetMem('Scr5','Allo','Real',iScr5,lRoots*nSel)
      CALL GetMem('CTEMP','ALLO','REAL',kctemp,ndet)
      CALL GetMem('SIGTEM','ALLO','REAL',ksigtemp,ndet)
*-------------------------------------------------------------------

* Print convergence thresholds in ITERFILE
      ThrRes = Max(0.2d-6,SQRT(ThrEne))
      Write(IterFile,'(19X,A,F18.10)') '- Threshold for energy   ...:',
     &         ThrEne
      Write(IterFile,'(19X,A,F18.10)') '- Threshold for Residual ...:',
     &         ThrRes
      Write(IterFile,'(20A4)') ('****',i=1,20)
      Write(IterFile,*)
      Write(IterFile,'(1X,A4,4X,A4,4X,A18,4X,A14,4X,A14)')
     &       'Iter','Root','Energy','dE','Residual'
      Write(IterFile,'(72A1)') ('=',i=1,72)
      Call xFlush(IterFile)
*===================================================================
* start long loop over iterations
      nconverged=0
      iskipconv=1
      nnew(1)=lroots
      Do iterci=1,mxItr
*MGD Dynamically evaluate the size of nkeep (useful for large lroots
*when we don't want too large hamiltonians with strong linear dependencies)
         nhamil=nsel
         if (nsel.lt.1) nhamil=200
         istart=0
         if (iterci.eq.1) istart=1
         nkeep=nhamil/(nnew(2-istart))+1
         nkeep=max(min(nkeep,mxkeep),3)
*-------------------------------------------------------------------
*MGD for stability purposes recompute sigma vec from time to time
         idelta=0
         if ((mod(iterci-1,24).eq.0)) idelta=1
         Do i=istart,idelta
* New CI vectors (iterci,mroot) are available.
* compute new sigma vectors
           Do mRoot = 1,nnew(2-i)
              Call Load_CI_vec(iterci-1+i,mRoot,lRoots,nConf,
     &                        Work(iVec1),LuDavid)
              If ( iprlev.ge.DEBUG ) then
                 lPrint = Min(nConf,200)
                 Write (6,'(1X,A,I2,A,I2)')
     &                'CI vector, iter =',iterci,' mRoot =',mRoot
                 Write (6,'(1X,A)')
     &                '(max. 200 elements)'
                 Write (6,'(1X,A)')
     &                '-----------------------------'
                 Call dVcPrt(' ',' ',Work(iVec1),lPrint)
              End If

              Call Timing(Rolex_1,Swatch,Swatch,Swatch)
              IF (DOFARO) THEN
                ! determinant wavefunctions
                allocate(sgm(ndeta,ndetb))
                allocate(psi(ndeta,ndetb))

                CALL DCOPY_(NCONF, 0.0D0, 0, WORK(IVECSVC), 1)
                CALL REORD2(MY_NORB,NACTEL,1,0,
     &                      IWORK(KICONF(1)),IWORK(KCFTP),
     &                      WORK(IVEC1),WORK(IVECSVC),IWORK(IVKCNF))
                CALL CITRANS_SORT('C',WORK(IVECSVC),WORK(IVEC2))
                PSI = 0.0D0
                CALL CITRANS_CSF2SD(WORK(IVEC2),PSI)
                SGM = 0.0D0
                CALL SIGMA_UPDATE(HTU,GTUVX,SGM,PSI)
                CALL CITRANS_SD2CSF(SGM,WORK(IVEC2))
                CALL CITRANS_SORT('O',WORK(IVEC2),WORK(IVECSVC))
                Call Reord2(my_norb,NACTEL,1,1,
     &                      iWork(KICONF(1)),iWork(KCFTP),
     &                      Work(IVECSVC),Work(IVEC2),iWork(ivkcnf))

                If ( iprlev.ge.DEBUG ) then
                  FP=DNRM2_(NCONF,WORK(IVEC2),1)
                  WRITE(6,'(1X,A,F21.14)') 'sigma dnrm2_(faroald): ', FP
                End If

                ! free the arrays
                deallocate(sgm,psi)
              Else
C     Convert the CI-vector from CSF to Det. basis
                call dcopy_(nconf, work(ivec1), 1, work(kctemp),1)
                call dcopy_(ndet, 0.0d0, 0, work(ksigtemp), 1)
                CALL csdtvc(work(kctemp), work(ksigtemp), 1, work(kdtoc)
     &             ,iwork(kicts(1)), LSym, 1)
                call dcopy_(ndet, 0.0d0, 0, work(ksigtemp), 1)
                c_pointer = kctemp
C     Calling Lucia to determine the sigma vector
                CALL Lucia_Util('Sigma',iDummy,iDummy,Dummy)
C     Set mark so densi_master knows that the Sigma-vector exists on disk.
                iSigma_on_disk = 1
                CALL CSDTVC(work(iVec2), work(kctemp), 2, work(kdtoc),
     &             iWork(kicts(1)), LSym, 1)

                If ( iprlev.ge.DEBUG ) then
                  FP=DNRM2_(NCONF,WORK(IVEC2),1)
                  WRITE(6,'(1X,A,F21.14)') 'sigma dnrm2_(lucia):   ', FP
                End If
              End If

C  Add ECORE_HEX (different from zero when particle-hole formalism used)
              ECORE_HEX = GET_ECORE()
              call daxpy_(nconf,ecore_hex,work(iVec1),1,work(iVec2),1)
C Timings on generation of the sigma vector
              Call Timing(Rolex_2,Swatch,Swatch,Swatch)
              Rolex_2 = Rolex_2 - Rolex_1
              Rolex_3 = Rolex_3 + Rolex_2

              If ( iprlev.ge.DEBUG ) then
                 lPrint = Min(nConf,200)
                 Write (6,*) ' '
                 Write (6,'(1X,A,I2,A,I2)')
     &                'sigma vector, iter =',iterci,' mRoot =',mRoot
                 Write (6,'(1X,A)')
     &                '(max. 200 elements)'
                 Write (6,'(1X,A)')
     &                '--------------------------------'
                 Call dVcPrt(' ',' ',Work(iVec2),lPrint)
              End If
              Call Save_Sig_vec(iterci-1+i,mRoot,lRoots,nConf,
     &                 Work(iVec2),LuDavid)
           End Do
         End Do
* Sigma vectors (iterci,mroot) have been computed, for mroot=1..lroots
*-------------------------------------------------------------------
* compute Hsmall and Ssmall
* These are Hsmall(jtrial,ktrial), where jtrial is (jter,jroot), and
* ktrial is (kter,kroot), and similar Ssmall.
* jtrial=1..mxKeep*lroots correspond to jter=iterci-mxKeep+1..iterci
* (Fewer, at the beginning)

         jtrial = 0
         jstart=Max(1,iterci-nKeep+1)
         Do jter = iterci,jstart,-1
            Do jRoot = 1,nnew(iterci-jter+1)
               jtrial = jtrial+1
               Call Load_CI_vec(jter,jRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
               Call Load_Sig_vec(jter,jRoot,lRoots,nConf,Work(iVec2),
     &                    LuDavid)
               ktrial = 0
               Do kter = iterci,jter,-1
                  max_kRoot = nnew(iterci-kter+1)
                  If ( kter.eq.jter ) max_kRoot = jRoot
                  Do kRoot = 1, max_kRoot
                     ktrial = ktrial+1
                     Call Load_CI_vec(kter,kRoot,lRoots,nConf,
     &                    Work(iVec3),LuDavid)
                     ij = ktrial+(jtrial*jtrial-jtrial)/2
                     Sji = dDot_(nConf,Work(iVec1),1,Work(iVec3),1)
                     Hji = dDot_(nConf,Work(iVec2),1,Work(iVec3),1)
                     Work(iSs+ij-1) = Sji
                     Work(iHs+ij-1) = Hji
                  End Do
               End Do
            End Do
         End Do
         ntrial=jtrial
         If ( iprlev.ge.DEBUG ) then
            Call TriPrt('Hsmall',' ',Work(iHs),ntrial)
            Call TriPrt('Ssmall',' ',Work(iSs),ntrial)
         End If
* Hsmall and Ssmall have been computed (ntrial x ntrial, in triangular
* storage.)
*-------------------------------------------------------------------

* solve secular equation HC=SCE.

* PAM2009 nBasVec on input=min(ntrial,nconf)
* nBasVec returned as nr of orthonormal solutions to HC=SCE
         nBasVec=nConf
         Call HCSCE(ntrial,
     &        Work(iHs),Work(iSs),Work(iCs),Work(iEs),nBasVec)
         If ( nBasVec.lt.lRoots ) then
            Write(6,*) 'David: nBasVec less than lRoots'
            Write(6,*) 'nBasvec, lRoots = ',nBasVec, lRoots
            If ( ICIRST.EQ.1 ) Write(6,*) 'CIREstart was used. ',
     & 'Check the number of roots in the previous calculation'
            Call QTrace
            Call Abend
         Endif
         If ( iprlev.ge.DEBUG ) then
            Call dVcPrt('Eigenvalues of Hsmall',' ',
     &           Work(iEs),ntrial)
            Call RecPrt('Eigenvectors of Hsmall',' ',
     &           Work(iCs),ntrial,ntrial)
         End If
*-------------------------------------------------------------------
* compute the current 'best' CI, sigma and residual vector

* CI vector is Work(iVec1)
* sigma vector is saved in Work(iVec2)
* residual vector is saved in Work(iVec3)
         Do mRoot=1,lRoots
*...      initialize 'best' CI and sigma vector
            Call dCopy_(nConf,0.0d0,0,Work(iVec1),1)
            Call dCopy_(nConf,0.0d0,0,Work(iVec2),1)
*...      accumulate contributions
            jtrial = 0
            Do jter=iterci,jstart,-1
               Do jRoot=1,nnew(iterci-jter+1)
                  jtrial = jtrial+1
                  Cik=Work(iCs-1+jtrial+(mRoot-1)*ntrial)
                  Call Load_CI_vec(jter,jRoot,lRoots,nConf,Work(iVec3),
     &                       LuDavid)
                  Call Daxpy_(nConf,Cik,Work(iVec3),1,Work(iVec1),1)
                  Call Load_Sig_vec(jter,jRoot,lRoots,nConf,Work(iVec3),
     &                       LuDavid)
                  Call Daxpy_(nConf,Cik,Work(iVec3),1,Work(iVec2),1)
               End Do
            End Do
            RR = dDot_(nConf,Work(iVec1),1,Work(iVec1),1)
            scl=1.0d0/sqrt(RR)
            Call DScal_(nConf,scl,Work(iVec1),1)
            Call DScal_(nConf,scl,Work(iVec2),1)
            Call Save_tmp_CI_vec(mRoot,lRoots,nConf,Work(iVec1),LuDavid)
            Call Save_tmp_Sig_vec(mRoot,lRoots,nConf,Work(iVec2),
     &                       LuDavid)
*...      compute residual vector
            E0 = Work(iEs+mRoot-1)
            Call dCopy_(nConf,Work(iVec2),1,Work(iVec3),1)
            call daxpy_(nConf,-E0,Work(iVec1),1,Work(iVec3),1)
*...      save current best energy and residual
            RR = dDot_(nConf,Work(iVec3),1,Work(iVec3),1)
            CI_conv(1,mroot,iterci) = E0
            CI_conv(2,mroot,iterci) = SQRT(RR)
*...  print vectors
            If ( iprlev.ge.DEBUG ) then
               lPrint = Min(nConf,200)
               Write (6,'(1X,A,I2,A,I2)')
     &              'new best CI vector, iter =',iterci,' mRoot =',mRoot
               Write (6,'(1X,A)')
     &              '(max. 200 elements)'
               Write (6,'(1X,A)')
     &              '--------------------------------------'
               Call dVcPrt(' ',' ',Work(iVec1),lPrint)
               Write (6,'(1X,A,I2,A,I2)')
     &           'new best sigma vector, iter =',iterci,' mRoot =',mRoot
               Write (6,'(1X,A)')
     &              '(max. 200 elements)'
               Write (6,'(1X,A)')
     &              '-----------------------------------------'
               Call dVcPrt(' ',' ',Work(iVec2),lPrint)
               Write (6,'(1X,A,I2,A,I2)')
     &           'new residual vector, iter =',iterci,' mRoot =',mRoot
               Write (6,'(1X,A)')
     &              '(max. 200 elements)'
               Write (6,'(1X,A)')
     &              '-----------------------------------------'
               Call dVcPrt(' ',' ',Work(iVec3),lPrint)
            End If
*...      to improve the preconditioner select all elements in the
*...      subspace of the explicit Hamiltonian
            If ( nSel.gt.1 ) then
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  iConf = iSel(i)
                  Work(iScr3+iOff+i-1) = Work(iVec3+iConf-1)
                  Work(iScr4+iOff+i-1) = Work(iVec1+iConf-1)
               End Do
            End If
         End Do
* Current best CI & Sigma vectors have been stored in a temporary place
* for mroot=1..lroots.
* Also, the selected elements of the CI and Sigma vectors have been
* saved at Work(iScr3)(Sigma)  and Work(iScr4)(CI)
*-------------------------------------------------------------------

* check for convergence
         nItr = iterci
         If ( iterci.gt.1 ) then
            dE = CI_conv(1,1,iterci-1) - CI_conv(1,1,iterci)
         Else
            dE = 0.0d0
         End If
         Write(IterFile,'(1X,I4,4X,I4,4X,F18.10,4X,F14.10,4X,F14.10)')
     &          IterCi,1,CI_conv(1,1,iterci),dE,
     &          CI_conv(2,1,iterci)
         Do jRoot = 2,lRoots
            If ( iterci.gt.1 ) then
               dE = CI_conv(1,jroot,iterci-1) - CI_conv(1,jroot,iterci)
            End If
            Write(IterFile,'(9X,I4,4X,F18.10,4X,F14.10,4X,F14.10)')
     &              jRoot,CI_conv(1,jroot,iterci),dE,
     &              CI_conv(2,jroot,iterci)
         End Do
         If (lRoots .gt. 1) Write(IterFile,*)
         Call xFlush(IterFile)
         If ( iprlev.gt.DEBUG ) then
            Write (6,*)
            Write (6,'(1X,120A1)')('*',i=1,120)
            Write (6,'(1X,A,I2)') 'CI iteration ',iterci
            ThrRes = Max(0.2d-6,SQRT(ThrEne))
            Write (6,'(1X,A,2F18.10)') 'ThrEne,ThrRes=',ThrEne,ThrRes
            Do jRoot = 1,lRoots
               If ( iterci.gt.1 ) then
                  dE =CI_conv(1,jroot,iterci-1)-CI_conv(1,jroot,iterci)
               Else
                  dE = 0.0d0
               End If
               Write (6,'(1X,A,I2,A,F18.10,2(A,F14.10))')
     &          ' root ',jRoot,
     &          ' energy =',CI_conv(1,jroot,iterci),
     &          ' dE =',dE,
     &          ' residual =',CI_conv(2,jroot,iterci)
            End Do
            Write (6,'(1X,120A1)')('*',i=1,120)
            Write (6,*)
         End If
         ThrRes = Max(0.2d-6,SQRT(ThrEne))
         iConv = 0
         nconverged=0
*Do not check for convergence of hidden roots
         Do jRoot=1,lRoots-hroots
            If ( iterci.gt.1 ) then
               dE = CI_conv(1,jroot,iterci-1) - CI_conv(1,jroot,iterci)
            Else
               dE = 0.0d0
            End If
            dE = abs(dE)
            R  = CI_conv(2,jroot,iterci)
            If ( (dE.lt.ThrEne) .and. (R.lt.ThrRes) ) Then
               iConv = iConv+1
               If (jRoot.eq.nconverged+1) nconverged=nconverged+1
            EndIf
         End Do
         if (iskipconv.eq.0) nconverged=0
         If ( iConv.ge.lRoots-hroots ) Goto 100
*-------------------------------------------------------------------
* compute correction vectors q1 = r/(E0-H) and q2 = c/(E0-H)

         If ( nSel.gt.1 ) then
            Call DGEMM_('T','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr3),nSel,
     &                  0.0d0,Work(iScr5),nSel)
            Do mRoot=1,lRoots
               E0 = Work(iEs+mRoot-1)
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  Z = E0-ExplE(i)
                  If ( Abs(Z).lt.0.001d0 ) Z = 0.001d0
                  Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
               End Do
            End Do
            Call DGEMM_('N','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr5),nSel,
     &                  0.0d0,Work(iScr3),nSel)
            Call DGEMM_('T','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr4),nSel,
     &                  0.0d0,Work(iScr5),nSel)
            Do mRoot=1,lRoots
               E0 = Work(iEs+mRoot-1)
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  Z = E0-ExplE(i)
                  If ( Abs(Z).lt.0.001d0 ) Z = 0.001d0
                  Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
               End Do
            End Do
            Call DGEMM_('N','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr5),nSel,
     &                  0.0d0,Work(iScr4),nSel)
         End If
*-------------------------------------------------------------------
         Do mRoot=nconverged+1,lRoots
            E0 = -Work(iEs+mRoot-1)
            Call Load_tmp_Sig_vec(mRoot,lRoots,nConf,Work(iVec1),
     &                LuDavid)
            Call Load_tmp_CI_vec(mRoot,lRoots,nConf,Work(iVec2),LuDavid)
            call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
            Call Load_H_diag(nConf,Work(iVec3),LuDavid)
            E0 = Work(iEs+mRoot-1)
            Do i = 0,nConf-1
               Z = E0-Work(iVec3+i)
               If ( ABS(Z).lt.1.0d-4 ) Z = 1.0d-4
               Work(iVec3+i) = Work(iVec1+i)/Z
            End Do
            If ( nSel.gt.1 ) then
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  iConf = iSel(i)
                  Work(iVec3+iConf-1) = Work(iScr3+iOff+i-1)
               End Do
            End If
            Alpha(mRoot) = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
            Call Load_H_diag(nConf,Work(iVec3),LuDavid)
            E0 = Work(iEs+mRoot-1)
            Do i = 0,nConf-1
               Z = E0-Work(iVec3+i)
               If ( ABS(Z).lt.1.0d-4 ) Z = 1.0d-4
               Work(iVec3+i) = Work(iVec2+i)/Z
            End Do
            If ( nSel.gt.1 ) then
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  iConf = iSel(i)
                  Work(iVec3+iConf-1) = Work(iScr4+iOff+i-1)
               End Do
            End If
            Beta(mRoot) = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
         End Do
*-------------------------------------------------------------------

* compute correction vectors q3 = (r-E1*q2)/(E0-H)
         If ( nSel.gt.1 ) then
            Do mRoot=1,lRoots
               Call Load_tmp_Sig_vec(mRoot,lRoots,nConf,Work(iVec1),
     &                     LuDavid)
               Call Load_tmp_CI_vec(mRoot,lRoots,nConf,Work(iVec2),
     &                     LuDavid)
               E0 = -Work(iEs+mRoot-1)
               call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
               E1 = -Alpha(mRoot)/Beta(mRoot)
               call daxpy_(nConf,E1,Work(iVec2),1,Work(iVec1),1)
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  iConf = iSel(i)
                  Work(iScr3+iOff+i-1) = Work(iVec1+iConf-1)
               End Do
            End Do
            Call DGEMM_('T','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr3),nSel,
     &                  0.0d0,Work(iScr5),nSel)
            Do mRoot=1,lRoots
               E0 = Work(iEs+mRoot-1)
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  Z = E0-ExplE(i)
                  If ( Abs(Z).lt.0.001d0 ) Z = 0.001d0
                  Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
               End Do
            End Do
            Call DGEMM_('N','N',
     &                  nSel,lRoots,nSel,
     &                  1.0d0,ExplV,nSel,
     &                  Work(iScr5),nSel,
     &                  0.0d0,Work(iScr3),nSel)
         End If
         Do i=min(iterci,mxkeep-1),2,-1
           nnew(i+1)=nnew(i)
         End Do
         nnew(2)=0
         Do mRoot=nconverged+1,lRoots
            Call Load_tmp_Sig_vec(mRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
            Call Load_tmp_CI_vec(mRoot,lRoots,nConf,Work(iVec2),LuDavid)
            E0 = -Work(iEs+mRoot-1)
            call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
            E1 = -Alpha(mRoot)/Beta(mRoot)
            call daxpy_(nConf,E1,Work(iVec2),1,Work(iVec1),1)
            Call Load_H_diag(nConf,Work(iVec3),LuDavid)
            E0 = Work(iEs+mRoot-1)
            Do i = 0,nConf-1
               Z = E0-Work(iVec3+i)
               If ( ABS(Z).lt.1.0d-4 ) Z = 1.0d-4
               Work(iVec3+i) = Work(iVec1+i)/Z
            End Do
            If ( nSel.gt.1 ) then
               iOff = (mRoot-1)*nSel
               Do i = 1,nSel
                  iConf = iSel(i)
                  Work(iVec3+iConf-1) = Work(iScr3+iOff+i-1)
               End Do
            End If
*Orthonormalize wrt previous vectors
            updsiz=dnrm2_(nconf,Work(iVec3),1)
            scl=1.0D0/updsiz
            Call DScal_(nConf,scl,Work(iVec3),1)
            Do jter=iterci,jstart+1,-1
               Do jRoot=1,nnew(iterci-jter+1)
                 Call Load_CI_vec(jter,jRoot,lRoots,nConf,Work(iVec2),
     &                    LuDavid)
                 ovl = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
                 call daxpy_(nConf,-ovl,Work(iVec2),1,Work(iVec3),1)
               End Do
            End Do
            updsiz=dnrm2_(nconf,Work(iVec3),1)
            If (updsiz.gt.1.0d-6) then
              scl=1.0D0/updsiz
              Call DScal_(nConf,scl,Work(iVec3),1)
              nnew(2)=nnew(2)+1
              Call Save_CI_vec(iterci,nnew(2),lRoots,nConf,Work(iVec3),
     &                    LuDavid)
            EndIf
         End Do
*-------------------------------------------------------------------
* move the current best CI and sigma vectors to the last place
* in the list of retained CI vectors
         Do mRoot=1,lRoots
           Call Load_tmp_CI_vec(mRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
           Call Save_CI_vec(iterci+1,mRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
           Call Load_tmp_Sig_vec(mRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
           Call Save_Sig_vec(iterci+1,mRoot,lRoots,nConf,Work(iVec1),
     &                    LuDavid)
         End Do

* end of the long loop over iterations
      End Do
*===================================================================

      mxItr = Min(mxCiIt,mxItr+12)
      If (IPRLEV.ge.USUAL) Then
        Write (6,*) '       ',
     &     'No convergence in the CI section: ',
     &     'MAXJT will be increased to ',mxItr
      End If
      Rc_CI = 16
      nItr=nItr-1
* deallocate local temporary vectors
 100  CALL GetMem('CTEMP','Free','REAL',kctemp,ndet)
      CALL GetMem('SIGTEM','Free','REAL',ksigtemp,ndet)
      Call GetMem('Vector1','Free','Real',iVec1,ndet)
      Call GetMem('Vector2','Free','Real',iVec2,ndet)
      Call GetMem('Vector3','Free','Real',iVec3,ndet)
      Call GetMem('Esmall','Free','Real',iEs,l1)
      Call GetMem('Hsmall','Free','Real',iHs,l3)
      Call GetMem('Ssmall','Free','Real',iSs,l3)
      Call GetMem('Csmall','Free','Real',iCs,l2)
      Call GetMem('Scr1','Free','Real',iScr1,l2)
      Call GetMem('Scr2','Free','Real',iScr2,l2)
      Call GetMem('Scr3','Free','Real',iScr3,lRoots*nSel)
      Call GetMem('Scr4','Free','Real',iScr4,lRoots*nSel)
      Call GetMem('Scr5','Free','Real',iScr5,lRoots*nSel)

      Call Timing(Alfex_2,Swatch,Swatch,Swatch)
      Alfex_2 = Alfex_2 - Alfex_1
      Alfex_3 = Alfex_3 + Alfex_2
      Call qExit('David')

      Return
      End
