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

subroutine amfi(LUIN,LUPROP,iCenter)
!#######################################################################
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
!  Connection to MOLCAS:
!  How wonderful, they normalize the functions exactly as I do, which
!  means they use the correct linear combinations.
!
!  Exponents and coefficients are expected in the MOLCAS-Format
!  first exponents
!  coefficients afterwards
!
!                                           8.5.97
!#######################################################################

use AMFI_global, only: ipowxyz, MxcontL, MxprimL, Lmax
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LUIN, LUPROP, iCenter
integer(kind=iwp) :: icartdim, icoulovlpdim, idim1, idim2, ifinite, ionecontrdim, ioneoverR3dim, ipowexpdim, irun, jrun, lhigh, &
                     lrun, Mval, ncont4, numballcart
logical(kind=iwp) :: AIMP, bonn, breit, keep, makemean, oneonly, SAMEORB
character(len=4) :: symmetry
integer(kind=iwp), allocatable :: checkxy(:), checkz(:), interxyz(:,:), SgnProd(:)
real(kind=wp), allocatable :: CartOne(:,:), CoulOvlp(:), Energy(:), evec(:,:), eval(:), OneContr(:), oneoverR3(:), PowExp(:), &
                              preXZ(:), preY(:), scratch(:,:,:), TKIN(:,:), type1(:), type2(:)
!bs the ones and zeros stand four odd and even powers of x,y,z
!bs if you want to go higher than l=6, you have to look up
!bs the powers yourself, and add them to the table
integer(kind=iwp), parameter :: Lpowmax = 6, &
                                ixyzpow(3*(Lpowmax+1)**2) = [                                                   &
                                  0,0,0,                                                                        & ! s-function
                                  0,1,0,0,0,1,1,0,0,                                                            & ! p-functions
                                  1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,                                                & ! d-functions
                                  0,1,0,1,1,1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,                                    & ! f-functions
                                  1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,                        & ! g-functions
                                  0,1,0,1,1,1,0,1,0,1,1,1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,            & ! h-functions
                                  1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0 & ! i-functions
                                ]
!keep     : parameter to decide about keeping angular integrals in memory
!makemean : 'true' = generating a mean field
!bonn     : 'true' = Bonn-approach for spin-other orbit
!breit    : if breit is set, BREIT-PAULI only
!SAMEORB  : parameter for same-orbit only
!AIMP     : parameter to delete CORE for AIMP
!oneonly  : parameter to use only one-electron integrals

!#######################################################################
!bs ####################################################################
!bs        version with all angular integrals in memory
!keep = .true.
!bs ####################################################################
!bs        version without all angular integrals in memory
keep = .false.
!bs ####################################################################
!bs initialize tables with double factorials...
call inidf()
!bs move some powers of x,y,z to the right place   BEGIN
!bs check if Lpowmax is high enough..
if (Lpowmax < Lmax) call SysAbendMsg('amfi','increase lpowmax and edit ixyzpow',' ')
jrun = 1
do irun=0,Lmax
  do Mval=-irun,irun
    ipowxyz(:,Mval,irun) = ixyzpow(jrun:jrun+2)
    jrun = jrun+3
  end do
end do
!bs move some powers of x,y,z to the right place   END
!bs read the input
call readbas(Lhigh,makemean,bonn,breit,symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,LUIN,ifinite)

icartdim = MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax
ionecontrdim = MxcontL*MxcontL*(2*Lmax+1)*3*Lmax
ioneoverR3dim = Lmax*(MxprimL*MxprimL+MxprimL)/2
ipowexpdim = MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*(Lmax+Lmax+6)
icoulovlpdim = MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*10
call mma_allocate(oneoverR3,ioneoverR3dim,label='oneoverR3')
call mma_allocate(cartone,icartdim,3,label='cartone')
call mma_allocate(OneContr,ionecontrdim,label='OneContr')
call mma_allocate(CoulOvlp,icoulovlpdim,label='coulovlp')
call mma_allocate(PowExp,iPowExpDim,label='PowExp')
call mma_allocate(TKIN,MxprimL,MxprimL,label='TKIN')
call mma_allocate(evec,MxprimL,MxprimL,label='evec')
call mma_allocate(eval,MxprimL,label='eval')
call mma_allocate(Energy,MxprimL,label='Energy')
call mma_allocate(type1,MxprimL,label='type1')
call mma_allocate(type2,MxprimL,label='type2')
call mma_allocate(scratch,MxprimL,MxprimL,3,label='scratch')
oneoverR3(:) = Zero
cartone(:,:) = Zero
OneContr(:) = Zero
CoulOvlp(:) = Zero
PowExp(:) = Zero

do
  if (ifinite == 2) call finite()

  ! Lhigh is the highest l-value in the basis set
  if (makemean .and. (.not. oneonly) .and. (ifinite <= 1)) call getAOs(Lhigh)
  call genpowers(Lhigh,PowExp,CoulOvlp)
  ! generate powers of exponents and overlaps
  !bs generate ovlp of normalized primitives
  call genovlp(Lhigh,CoulOvlp,eval)
  do lrun=0,Lhigh
    !bs cont(L) arranges all the contraction coefficients for a given
    !bs L-value and renormalizes them
    call cont(lrun,breit,ifinite,TKIN,evec,eval,Energy,type1,type2,scratch)
  end do

  !bs beginning the angular part
  if (.not. oneonly) then
    !BS write(u6,*) '***************************************************'
    !BS write(u6,*) '********   beginning the 2e-part ******************'
    !BS write(u6,*) '***************************************************'

    !bs ################################################################
    !bs ################################################################
    !bs ################################################################

    idim1 = (2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)
    idim2 = (Lmax+1)*(Lmax+1)*(Lmax+1)*(Lmax+1)
    call mma_allocate(preY,idim1,label='preY')
    call mma_allocate(preXZ,idim1,label='preXZ')
    call mma_allocate(checkxy,idim2,label='CheckXY')
    call mma_allocate(checkz,idim2,label='CheckZ')
    call mma_allocate(interxyz,16,idim2,label='InterXYZ')
    call mma_allocate(SgnProd,idim1,label='SgnProd')

    ! subroutine for angular part

    call angular(Lhigh,keep,makemean,bonn,breit,sameorb,ifinite,cartone(1,1),cartone(1,2),cartone(1,3),PowExp,CoulOvlp,preXZ,preY, &
                 checkxy,checkz,InterXYZ,SgnProd)

    call mma_deallocate(SgnProd)
    call mma_deallocate(InterXYZ)
    call mma_deallocate(CheckZ)
    call mma_deallocate(CheckXY)
    call mma_deallocate(preXZ)
    call mma_deallocate(preY)
  end if
  if (ifinite /= 1) exit
  ! redo everything for finite core
  !BS write(u6,*) 'once more the two-electron integrals'
  ifinite = 2
end do
!bs ####################################################################
!bs ####################################################################
!bs ####################################################################
!BS write(u6,*) '***************************************************'
!BS write(u6,*) '*******   beginning the 1-electron-part  **********'
!BS write(u6,*) '***************************************************'

!bs The one-electron spin-orbit integrals

call gen1overR3(Lhigh,oneoverR3)

! 1/r**3 for normalized functions

call contandmult(Lhigh,AIMP,oneonly,numballcart,LUPROP,ifinite,CartOne,OneContr,oneoverR3,iCenter)

!bs multiplies radial integrals with l,m-dependent
!bs factors and contraction coefficients
call mma_deallocate(CoulOvlp)
call mma_deallocate(PowExp)
call mma_deallocate(OneContr)
call mma_deallocate(CartOne)
call mma_deallocate(oneoverR3)
call mma_deallocate(TKIN)
call mma_deallocate(evec)
call mma_deallocate(eval)
call mma_deallocate(Energy)
call mma_deallocate(type1)
call mma_deallocate(type2)
call mma_deallocate(scratch)
!BS write(u6,*) '***************************************************'
!BS write(u6,*) '*******   end of  the 1-electron-part    **********'
!BS write(u6,*) '***************************************************'
!bs ####################################################################
!bs ####################################################################
!bs ####################################################################

return

end subroutine amfi
