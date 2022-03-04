!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996,1997, Bernd Schimmelpfennig                       *
!***********************************************************************
!hange start
!     program amfi
!     implicit real*8 (a-h,o-z)
!hange else
      Subroutine amfi(IN,LUPROP,iCenter)
      implicit Real*8 (a-h,o-z)
!hange end
!###########################################################################
!
!          A M F I
!
!    Atomic Mean-Field Spin-Orbit Integral Program
!
! Integral-code to generate the one- and two-electron spin-orbit integrals
! in the no-pair approximation for an atom.
!
! basis set is built by atomic functions of the form:
!
!     f(r,Omega)= r**l Y_(lm) (Omega)
!
! Allthough the code is created with a lot of care and love for
! the details, the author doesn't give any warranty for it's
! correctness.
!
! B.Schimmelpfennig  Fysikum/Stockholm Summer 1996
!
! If you use this code, please honour the authors work
! by citing this work properly.
!
! The author would like to thank the Deutsche Forschungsgemeinschaft
! for financing this project by a Forschungsstipendium.
!
!
!   The spatial integrals are expected to be used with a spin part
!   expressed in Paulis spin-matrices rather than with the Spin-operator
!   itself. So if a factor of two is somehow missing, check whether the
!   same form of the operator is used.
!
!
!   WARNING !!!   WARNING !!   WARNING !!  WARNING !!   WARNING !!
!
!   when writing spin-same-orbit and spin-other-orbit with sigma_i:
!
!   For the spin-other-orbit-integrals particle 1 and 2 are exchanged
!   on the arrays carteXOO,carteYOO,carteZOO!!!!!!!!!
!
!   The reason is to use most of the same-orbit part again and to
!   have the same symmetry for the integrals on the arrays.
!
!
!   if the spin-other-orbit-part is used in the formulation with
!   sigma_j, the particles are of cause not interchanged.
!
!
!
!   (i|HSO_mean|j) = (ij) + 1/2 * sum_M  occ(M) {
!                   2(ij|MM)_same - (iM|jM)_same -2(iM|jM)_other
!                   + (jM|iM)_same +2(jM|iM)_other }
!
!   in the subroutines some signs are changed  to reorder indices
!   in the integrals to (iM|jM) or (Mi|Mj) accoding to the way they
!   were calculated before.
!
!
!
!   one-particle integrals (really one-particle or mean-field)
!   are written to files in CONTANDMULT. Look there for information on
!   the format of files.
!
!
!  BUGS:  There is still a strange sign-error in the two-electron-integrals
!  if one applies straight-forward the formulae of the documentation.
!  This problem has been solved by the the cheater...
!
!  Everybody is welcome to find the problem in the formulas ........
!
!  First reasonable results on Thallium (SD with frozen 5D) 14.10.96
!
!
!
!
!
!  Connection to MOLCAS:
!  How wonderful, they normalize the functions exactly as I do, which
!  means they use the correct linear combinations.
!
!  Exponents and coefficients are expected in the MOLCAS-Format
!  first exponents
!  coefficients afterwards
!
!                                           8.5.97
!
!
!###########################################################################
#include "para.fh"
      logical keep    ! parameter to decide about keeping angular
!bs                     ! integrals in memory
      logical keepcart    ! parameter to decide about keeping cartesian
!bs                         ! integrals in memory
      logical makemean   ! 'true' =   generating a meanfield
      logical bonn       ! 'true' = Bonn-approach for spin-other orbit
      logical breit      ! if breit is set, BREIT-PAULI only
      logical SAMEORB    ! parameter for same-orbit only
      logical AIMP       ! parameter to delete CORE for AIMP
      logical oneonly    ! parameter to use only oneelectron integrals
      character*4  symmetry
      parameter (Lpowmax=6)
      dimension ixyzpow(3*(Lpowmax+1)*(Lpowmax+1)) !
      data ixyzpow /                                                    &
!bs   the ones and zeros stand four odd and even powers of x,y,z
!bs   if you want to go higher than l=6, you have to look up
!bs   the powers yourself, and add them to the table
     &0,0,0,                                                            & ! s-function
     &0,1,0, 0,0,1, 1,0,0,                                              & ! p-functions
     &1,1,0, 0,1,1, 0,0,0,  1,0,1, 0,0,0,                               & ! d-functions
     &0,1,0, 1,1,1, 0,1,0,  0,0,1, 1,0,0,                               & ! f-functions
     &0,0,1, 1,0,0,                                                     & ! f-functions
     &1,1,0, 0,1,1, 1,1,0,  0,1,1, 0,0,0,                               & ! g-functions
     &1,0,1, 0,0,0, 1,0,1,  0,0,0,                                      & ! g-functions
     &0,1,0, 1,1,1, 0,1,0,  1,1,1, 0,1,0,                               & ! h-functions
     &0,0,1, 1,0,0, 0,0,1,  1,0,0, 0,0,1,                               & ! h-functions
     &1,0,0,                                                            & ! h-functions
     &1,1,0, 0,1,1, 1,1,0, 0,1,1, 1,1,0,                                & ! i-functions
     &0,1,1, 0,0,0, 1,0,1, 0,0,0, 1,0,1,                                & ! i-functions
     &0,0,0, 1,0,1, 0,0,0                                               & ! i-functions
     &/
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: oneoverR3(:), CartOne(:,:), OneContr(:),    &
     &                      CoulOvlp(:), PowExp(:), preY(:), preXZ(:)
      Integer, Allocatable:: checkxy(:), checkz(:), interxyz(:,:),      &
     &                       SgnProd(:)
!
#include "ipowxyz.fh"
!##########################################################################
!bs  #####################################################################
!bs         version with all angular integrals in memory
!         keep=.true.
!bs  #####################################################################
!bs         version without  all angular integrals in memory
          keep=.false.
!bs  #####################################################################
!bs         version without  all cartesian integrals in memory
          keepcart=.false.
!bs  #####################################################################
!bs         version with all cartesian integrals in memory
!         keepcart=.true.
!bs  #####################################################################
      ifinite=0
!bs   initialize tables with double facultatives...
      call inidf
!bs   move some powers of x,y,z to the right place   BEGIN
!bs   check if Lpowmax is high enough..
      if (Lpowmax.lt.Lmax) then
      Call SysAbendMsg('amfi', 'increase lpowmax and edit ixyzpow',' ' )
      endif
      jrun=1
      do irun=0,Lmax
      do Mval=-irun,irun
      ipowxyz(1,Mval,irun)=ixyzpow(jrun)
      ipowxyz(2,Mval,irun)=ixyzpow(jrun+1)
      ipowxyz(3,Mval,irun)=ixyzpow(jrun+2)
      jrun=jrun+3
      enddo
      enddo
!bs   move some powers of x,y,z to the right place   END
!bs   read the input
      call readbas(Lhigh,makemean,bonn,breit,symmetry,sameorb,AIMP,     &
     &             oneonly,ncont4,numballcart,IN,ifinite)
!bs
      icartdim=mxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax
      ionecontrdim=mxcontL*MxcontL*(2*Lmax+1)*3*Lmax
      ioneoverR3dim=Lmax*(MxprimL*MxprimL+MxprimL)/2
      ipowexpdim=MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*(Lmax+Lmax+6)
      icoulovlpdim=MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*10
      Call mma_allocate(oneoverR3,ioneoverR3dim,Label='oneoverR3')
      Call mma_allocate(cartone,icartdim,3,Label='cartone')
      Call mma_allocate(OneContr,ionecontrdim,Label='OneContr')
      Call mma_allocate(CoulOvlp,icoulovlpdim,Label='coulovlp')
      Call mma_allocate(PowExp,iPowExpDim,Label='PowExp')
      oneoverR3(:)=0.0D0
      cartone(:,:)=0.0D0
      OneContr(:)=0.0D0
      CoulOvlp(:)=0.0D0
      PowExp(:)=0.0D0
!bs
!bs
 123  if (ifinite.eq.2) call finite
!bs
!bs
! Lhigh is the highest l-value in the basis set
      if (makemean.and.(.not.oneonly).and.ifinite.le.1)                 &
     &   call getAOs(Lhigh)
      call genpowers(Lhigh,PowExp,CoulOvlp)
!generate powers of exponents and overlaps
!bs   start generating modified contraction coefficients
!bs   generate starting adresses of contraction coefficients  on
!bs   contrarray
      call genstar(Lhigh)
!bs   generate ovlp of normalized primitives
      call genovlp(Lhigh,CoulOvlp)
      do lrun=0,Lhigh
!bs      cont(L) arranges all the contraction coefficients for a given
!bs      L-value and renormalizes them
         call cont(lrun,breit,ifinite)
      enddo
!bs
!bs        beginning the angular part
      if (.not.oneonly) then
!BS   write(6,*) '***************************************************'
!BS   write(6,*) '********   beginning the 2e-part ******************'
!BS   write(6,*) '***************************************************'
!bs
!bs  ###################################################################
!bs  ###################################################################
!bs  ###################################################################
!bs
!bs
      idim1=(2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)
      idim2=(Lmax+1)*(Lmax+1)*(Lmax+1)*(Lmax+1)
      Call mma_allocate(preY,idim1,Label='preY')
      Call mma_allocate(preXZ,idim1,Label='preXZ')
      Call mma_allocate(checkxy,idim2,Label='CheckXY')
      Call mma_allocate(checkz,idim2,Label='CheckZ')
      Call mma_allocate(interxyz,16,idim2,Label='InterXYZ')
      Call mma_allocate(SgnProd,idim1,Label='SgnProd')
!
!     subroutine for angular part
!
      call angular(Lhigh,keep,keepcart,makemean,bonn,breit,             &
     &             sameorb,ifinite,                                     &
     &             cartone(1,1),cartone(1,2),cartone(1,3),              &
     &             PowExp,CoulOvlp,preXZ,preY,                          &
     &             checkxy,checkz,InterXYZ,SgnProd)
!
      Call mma_deallocate(SgnProd)
      Call mma_deallocate(InterXYZ)
      Call mma_deallocate(CheckZ)
      Call mma_deallocate(CheckXY)
      Call mma_deallocate(preXZ)
      Call mma_deallocate(preY)
      endif
      if (ifinite.eq.1) then ! redo everything for finite core
!BS      write(6,*) 'once more the two-electron integrals'
         ifinite=2
         goto 123
      endif
!bs ####################################################################
!bs ####################################################################
!bs ####################################################################
!BS   write(6,*) '***************************************************'
!BS   write(6,*) '*******   beginning the 1-electron-part  **********'
!BS   write(6,*) '***************************************************'
!
!bs   The one-electron spin-orbit integrals
!
      call gen1overR3(Lhigh,oneoverR3)
!
!     1/r**3  for normalized functions
!
      call contandmult(Lhigh,makemean,AIMP,oneonly,numballcart,LUPROP,  &
     &                 ifinite,CartOne,OneContr,oneoverR3,iCenter)
!
!bs   multiplies radial integrals with l,m-dependent
!bs   factors and contraction coefficients
      Call mma_deallocate(CoulOvlp)
      Call mma_deallocate(PowExp)
      Call mma_deallocate(OneContr)
      Call mma_deallocate(CartOne)
      Call mma_deallocate(oneoverR3)
!BS   write(6,*) '***************************************************'
!BS   write(6,*) '*******   end of  the 1-electron-part    **********'
!BS   write(6,*) '***************************************************'
!bs ####################################################################
!bs ####################################################################
!bs ####################################################################
      Return
      End Subroutine Amfi
