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
      use Str_Info, only: DTOC
      use negpre
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Files_mclr.fh"
#include "stdalloc.fh"
#include "Pointers.fh"
#include "sa.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
#include "dmrginfo_mclr.fh"
      logical ldisk,ipopen
      Real*8, Allocatable:: CIVec(:,:), CITmp(:)

! ==========================================================
      integer,allocatable::index_SD(:) ! not final version
      real*8,allocatable::vector_cidmrg(:)
! ==========================================================

*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine RdJobIph_td(CIVec)
        Real*8, Allocatable:: CIVec(:,:)
        End Subroutine RdJobIph_td
        Subroutine RdJobIph(CIVec)
        Real*8, Allocatable:: CIVec(:,:)
        End Subroutine RdJobIph
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      !Read in interesting info from RUNFILE and ONEINT
      Call Rd1Int_MCLR()
      Call RdAB()   ! Read in orbitals, perturbation type, etc.
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
      If (iMethod.eq.iCASSCF) Then
         If (TimeDep) Then
            Call RdJobIph_td(CIVec)
         Else
            Call RdJobIph(CIVec)
         End If

*        Write(6,*) 'Setup of Determinant tables'
         Call DetCtl   ! set up determinant tables
*....... Read in tables from disk
         Call InCsfSD(State_sym,State_sym,.true.)
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
             Call mma_allocate(CITmp,ndets_RGLR,Label='CITmp')
           else
             Call mma_allocate(CITmp,nconf,Label='CITmp')
             call dcopy_(nconf,CIVec(:,i),1,CITmp,1)
           end if

           !> If doDMRG
           if(doDMRG.and.doMCLR)then ! yma
           else
             Call GugaCtl_MCLR(CITmp,1)   ! transform to sym. group
           end if

! Here should be the position for introducing the CI(SR) coefficients
!           iSSM=1     ! yma
!           write(6,*)"Set ISSM eq 1 ",ISSM

           if(doDMRG)then !yma
             ! mma_allocate and mma_deallocate
             allocate(index_SD(ndets_RGLR))
             allocate(vector_cidmrg(ndets_RGLR))
             call ci_reconstruct(i,ndets_RGLR,vector_cidmrg,index_SD)
             do ii=1,ndets_RGLR
               if(abs(vector_cidmrg(ii)).lt.0.0d0)then
                 vector_cidmrg(ii)=0.0d0
               end if
             end do
             call CSDTVC_dmrg(CITmp,vector_cidmrg,2,DTOC,
     &                     index_SD,ISSM,1,IPRDIA)
             ! mma_allocate and mma_deallocate
             deallocate(index_SD)
             deallocate(vector_cidmrg)
           end if

           call dcopy_(nconf,CITmp,1,CIVec(:,i),1)
           Call mma_deallocate(CITmp)
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
            irc=ipin(ipcii)
            call dcopy_(nconf*nroots,CIVec,1,W(ipcii)%Vec,1)
            nDisp=1
         Else
            ipcii=ipget(nconf)
            irc=ipin(ipcii)
            call dcopy_(nConf,CIVec(:,iState),1,W(ipcii)%Vec,1)
            If (iRoot(iState).ne.1) Then
               Write (6,*) 'McKinley does not support computation of'
     &                   //' harmonic frequencies of excited states'
               Call Abend()
            End If
         End If
C        irc=ipin(ipcii)
C        Call RecPrt('CI vector',' ',W(ipcii)%Vec,1,nConf)
         Call mma_deallocate(CIVec)
*
*        At this point we change to ipci being the index of the CI
*        vector in the ipage utility.
*
         ipci=ipcii
         irc=ipout(ipci)
*                                                                      *
************************************************************************
*                                                                      *
         If (ngp) Call rdciv()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call InpOne()         ! read in oneham
      Call PrInp_MCLR(iPL)  ! Print all info
*                                                                      *
************************************************************************
*                                                                      *
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Then
         Call Unused_integer(irc)
         Call Unused_logical(ldisk)
      End If
#endif
      End
