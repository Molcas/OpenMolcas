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
* Copyright (C) 2012, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      subroutine xdwinit(Heff,H0,U0)

      use definitions, only: wp, iwp, u6
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_grad
      use caspt2_global, only: CMO, CMO_Internal, CMOPT2, NCMO
      use caspt2_global, only: FIFA, DREF
      use caspt2_global, only: LUONEM
      use PrintLevel, only: DEBUG, INSANE, USUAL, VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nState, IfChol, iSCF, nAshT, nConf,
     &                         STSym, iAd1m, mState
      use pt2_guga, only: CIThr

      implicit none

      Real(kind=wp),intent(inout) :: Heff(Nstate,Nstate)
      Real(kind=wp),intent(inout) :: H0(Nstate,Nstate)
      Real(kind=wp),intent(inout) :: U0(Nstate,Nstate)

      Real(kind=wp) :: wgt,FIJ

      Integer(kind=iwp) :: iState,iDisk,I,J
      Real(kind=wp), allocatable:: CI(:), DAVE(:), CIRef(:,:),
     &                              CIXMS(:)


* Allocate memory for CI array state averaged 1-RDM
      call mma_allocATE(CI,NCONF,Label='CI')
      call mma_allocate(DAVE,SIZE(DREF),Label='DAVE')
      DAVE(:)=0.0_wp

* Set the weight for the density averaging
      wgt = 1.0_wp/real(Nstate,kind=wp)

* Loop over all states to compute the state-average density matrix
      do Istate=1,Nstate

        if (ISCF.NE.0) then
* Special case for a single Slater determinant
          CI(1)=1.0_wp
        else
* Get the CI array
          call loadCI(CI, Istate)
        end if

* Compute 1-particle active density matrix GAMMA1
        call POLY1(CI,nConf)

* Restructure GAMMA1 as DREF array
        call GETDREF(DREF,SIZE(DREF))

* Average the density
        call DAXPY_(SIZE(DREF),wgt,DREF,1,DAVE,1)

      end do

      if (IPRGLB.GE.INSANE) then
        write(u6,*)' State-average 1-RDM'
        do I=1,NASHT
          write(u6,'(1x,14f10.6)')(DAVE((I*(I-1))/2+J),J=1,I)
        end do
        write(u6,*)
      end if

* Copy the state-average 1-RDM into DREF and release memory
* for both the CI array and DAVE
      DREF(:)=DAVE(:)
      call mma_deallocate(CI)
      call mma_deallocate(DAVE)

* Load CASSCF MO coefficients
      call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
      CMO=>CMO_Internal
      iDisk=IAD1M(1)
      call ddafile(LUONEM,2,CMO,NCMO,iDisk)

* Build the state-average Fock matrix in MO basis
      if (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
        call INTCTL2(.false.)
      else
* INTCTL1 uses TRAONE and FOCK_RPT2 to get the matrices in MO basis.
        call INTCTL1(CMO,SIZE(CMO))
      end if

* Loop again over all states to compute H0 in the model space
* Loop over ket functions
      do J=1,Nstate
* Loop over bra functions
        do I=1,Nstate
* Compute matrix element <I|F|J> and store it into H0
          FIJ = 0.0_wp
          call FOPAB(FIFA,SIZE(FIFA),I,J,FIJ)
          H0(I,J) = FIJ
        end do
      end do
* End of loop over states

      if (IPRGLB.ge.USUAL) then
        write(u6,*)
        write(u6,*)' H0 in the original model space basis:'
        call prettyprint(H0,Nstate,Nstate)
      end if

* Diagonalize H0 in the model space
      call eigen(H0,U0,Nstate)

* Transform the Fock matrix in the new basis
      call transmat(H0,U0,Nstate)
        if (IPRGLB.ge.USUAL) then
          write(u6,*)' H0 eigenvectors:'
          call prettyprint(U0,Nstate,Nstate)
        end if
        if (IPRGLB.ge.DEBUG) then
          write(u6,*)' H0 in the rotated model space basis:'
          call prettyprint(H0,Nstate,Nstate)
        end if

* As well as Heff
      call transmat(Heff,U0,Nstate)
      if (IPRGLB.ge.VERBOSE) then
        write(u6,*)' Heff[1] in the rotated model space basis:'
        call prettyprint(Heff,Nstate,Nstate)
      end if

* Mix the CI arrays according to the H0 eigenvectors. Assume we can
* put all the original ones in memory, but put the resulting vectors
* one by one in a buffer.
        if (IPRGLB.ge.VERBOSE) then
          write(u6,'(A)')' The CASSCF states are now rotated'//
     &                  ' according to the H0 eigenvectors'
          write(u6,*)
        end if

      call mma_allocate(CIref,nConf,Nstate,Label='CIRef')
* Load the CI arrays into memory
      do I=1,Nstate
        call loadCI(CIref(:,I),I)
      end do

      call mma_allocate(CIXMS,Nconf,Label='CIXMS')
      do J=1,Nstate
* Transform the states
        call dgemm_('N','N',Nconf,1,Nstate,
     &              1.0_wp,CIREF,Nconf,U0(:,J),Nstate,
     &              0.0_wp,CIXMS,Nconf)

* Write the rotated CI coefficients back into LUCIEX and REPLACE the
* original unrotated CASSCF states. Note that the original states
* are still available in the JobIph file
        call writeCI(CIXMS,J)

        if (IPRGLB.ge.VERBOSE) then
          write(u6,'(1x,a,i3)')
     &    ' The CI coefficients of rotated model state nr. ',MSTATE(J)
          call PRWF_CP2(STSYM,NCONF,CIXMS,CITHR)
        end if
      end do
C
      If (do_grad) call dcopy_(NCMO,CMO,1,CMOPT2,1)
      If (do_grad) CMOPT2(:)=CMO(:)

* Release all memory
      call mma_deallocate(CIRef)
      call mma_deallocate(CIXMS)
      call mma_deallocate(CMO_Internal)
      nullify(CMO)

      end subroutine xdwinit

