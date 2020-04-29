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
      Subroutine InpCtl_MCLR(iPL)
************************************************************************
*                                                                      *
*     Read all relevant input data and display them                    *
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
#include "WrkSpc.fh"
#include "Pointers.fh"
#include "sa.fh"
#include "negpre.fh"

#include "detdim.fh"
#include "csfbas_mclr.fh"
#include "spinfo_mclr.fh"
#include "dmrginfo_mclr.fh"
      logical ldisk,ipopen

! ==========================================================
      integer,allocatable::index_SD(:) ! not final version
      real*8,allocatable::vector_cidmrg(:)
! ==========================================================

*                                                                      *
************************************************************************
*                                                                      *
      Call Rd1Int_MCLR !Read in interesting info from RUNFILE and ONEINT
      Call RdAB   ! Read in orbitals, perturbation type, etc.
*                                                                      *
************************************************************************
*                                                                      *
      Call Rd2Int(iPL) ! Read in 2el header
*                                                                      *
************************************************************************
*                                                                      *
      Call RdInp_MCLR()  ! Read in input
*                                                                      *
************************************************************************
*                                                                      *
*     Default activate ippage utility

*

      ldisk  =ipopen(0,.True.)
*
      write(*,*) "iMethod:",iMethod,iCASSCF
      If (iMethod.eq.iCASSCF) Then
         If (TimeDep) Then
            Call RdJobIph_td
         Else
            Call RdJobIph
         End If

*        Write(6,*) 'Setup of Determinant tables'
         Call DetCtl   ! set up determinant tables
*....... Read in tables from disk
         Call InCsfSD(State_sym,State_sym,.true.)

             !Call GetMem('CIvec','Allo','Real',ipNEW,NCONF)
             !Call GetMem('OCIvec','Free','Real',ipCI,nConf)
             !ipCI=ipNEW
*                                                                      *
************************************************************************
*                                                                      *
*        Write(6,*) 'Transformation of CI vector to symmetric '
*    &             ,'group from GUGA pepresentation'

         !> scratch  ! yma testing
!         if(doDMRG.and.doMCLR)then
!           call xflush(117)
!           close(117)
!         end if

         Do i=1,nroots
           ! yma
!          No need to copy,since there are no CI-vectors
           if(doDMRG.and.doMCLR)then
             Call Getmem('CIROOT','ALLO','REAL',ipT,ndets_RGLR)
           else
             Call Getmem('CIROOT','ALLO','REAL',ipT,nconf)
             call dcopy_(nconf,Work(ipCI+(i-1)*nconf),1,Work(ipT),1)
           end if

           !> If doDMRG
           if(doDMRG.and.doMCLR)then ! yma
           else
             Call GugaCtl_MCLR(ipT,1)   ! transform to sym. group
           end if

! Here should be the position for introducing the CI(SR) coefficients
!           iSSM=1     ! yma
!           write(*,*)"Set ISSM eq 1 ",ISSM

           if(doDMRG)then !yma
             ! mma_allocate and mma_deallocate
             allocate(index_SD(ndets_RGLR))
             allocate(vector_cidmrg(ndets_RGLR))
             call ci_reconstruct(i,ndets_RGLR,vector_cidmrg,
     &                           index_SD)
             do ii=1,ndets_RGLR
               if(abs(vector_cidmrg(ii)).lt.0.0d0)then
                 vector_cidmrg(ii)=0.0d0
               end if
             end do
             call CSDTVC_dmrg(work(ipT),vector_cidmrg,2,WORK(KDTOC),
     &                     index_SD,ISSM,1,IPRDIA)
             ! mma_allocate and mma_deallocate
             deallocate(index_SD)
             deallocate(vector_cidmrg)
           end if

           call dcopy_(nconf,Work(ipT),1,Work(ipCI+(i-1)*nconf),1)

          if(doDMRG.and.doMCLR)then !yma
            Call Getmem('CIROOT','FREE','REAL',ipT,ndets_RGLR)
          else
            Call Getmem('CIROOT','FREE','REAL',ipT,nconf)
          end if
        End Do

*                                                                      *
************************************************************************
*                                                                      *
        ldisk  =ipopen(nconf,page)
*
*        If we are computing Lagrangian multipliers we pick up all CI
*        vectors. For Hessian calculations we pick up just one vector.
*
C        Write (*,*) 'iState,SA,nroots=',iState,SA,nroots
         If (SA.or.iMCPD) Then
            ipcii=ipget(nconf*nroots)
            call dcopy_(nconf*nroots,Work(ipCI),1,Work(ipin(ipcii)),1)
            nDisp=1
         Else
            ipcii=ipget(nconf)
            ipCI_ = ipCI + (iState-1)*nConf
            call dcopy_(nConf,Work(ipCI_),1,Work(ipin(ipcii)),1)
            If (iRoot(iState).ne.1) Then
               Write (6,*) 'McKinley does not support computation of'
     &                   //' harmonic frequencies of excited states'
               Call Abend()
            End If
         End If
C        Call RecPrt('CI vector',' ',Work(ipin(ipcii)),1,nConf)
         Call Getmem('CIVEC','FREE','REAL',ipci,idum)
         ipci=ipcii
         irc=ipout(ipci)
*                                                                      *
************************************************************************
*                                                                      *
         If (ngp) Call rdciv
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call InpOne           ! read in oneham
      Call PrInp_MCLR(iPL)  ! Print all info
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
