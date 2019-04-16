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
* Copyright (C) 1996,1997, Bernd Schimmelpfennig                       *
************************************************************************
Change start
c     program amfi
c     implicit real*8 (a-h,o-z)
Change else
      Subroutine amfi(IN,LUPROP,iCenter)
      implicit Real*8 (a-h,o-z)
Change end
c###########################################################################
c
c          A M F I
c
c    Atomic Mean-Field Spin-Orbit Integral Program
c
c Integral-code to generate the one- and two-electron spin-orbit integrals
c in the no-pair approximation for an atom.
c
c basis set is built by atomic functions of the form:
c
c     f(r,Omega)= r**l Y_(lm) (Omega)
c
c Allthough the code is created with a lot of care and love for
c the details, the author doesn't give any warranty for it's
c correctness.
c
c B.Schimmelpfennig  Fysikum/Stockholm Summer 1996
c
c If you use this code, please honour the authors work
c by citing this work properly.
c
c The author would like to thank the Deutsche Forschungsgemeinschaft
c for financing this project by a Forschungsstipendium.
c
c
c   The spatial integrals are expected to be used with a spin part
c   expressed in Paulis spin-matrices rather than with the Spin-operator
c   itself. So if a factor of two is somehow missing, check whether the
c   same form of the operator is used.
c
c
c   WARNING !!!   WARNING !!   WARNING !!  WARNING !!   WARNING !!
c
c   when writing spin-same-orbit and spin-other-orbit with sigma_i:
c
c   For the spin-other-orbit-integrals particle 1 and 2 are exchanged
c   on the arrays carteXOO,carteYOO,carteZOO!!!!!!!!!
c
c   The reason is to use most of the same-orbit part again and to
c   have the same symmetry for the integrals on the arrays.
c
c
c   if the spin-other-orbit-part is used in the formulation with
c   sigma_j, the particles are of cause not interchanged.
c
c
c
c   (i|HSO_mean|j) = (ij) + 1/2 * sum_M  occ(M) {
c                   2(ij|MM)_same - (iM|jM)_same -2(iM|jM)_other
c                   + (jM|iM)_same +2(jM|iM)_other }
c
c   in the subroutines some signs are changed  to reorder indices
c   in the integrals to (iM|jM) or (Mi|Mj) accoding to the way they
c   were calculated before.
c
c
c
c   one-particle integrals (really one-particle or mean-field)
c   are written to files in CONTANDMULT. Look there for information on
c   the format of files.
c
c
c  BUGS:  There is still a strange sign-error in the two-electron-integrals
c  if one applies straight-forward the formulae of the documentation.
c  This problem has been solved by the the cheater...
c
c  Everybody is welcome to find the problem in the formulas ........
c
c  First reasonable results on Thallium (SD with frozen 5D) 14.10.96
c
c
c
c
c
c  Connection to MOLCAS:
c  How wonderful, they normalize the functions exactly as I do, which
c  means they use the correct linear combinations.
c
c  Exponents and coefficients are expected in the MOLCAS-Format
c  first exponents
c  coefficients afterwards
c
c                                           8.5.97
c
c
c###########################################################################
#include "para.fh"
      logical keep    ! parameter to decide about keeping angular
cbs                     ! integrals in memory
      logical keepcart    ! parameter to decide about keeping cartesian
cbs                         ! integrals in memory
      logical makemean   ! 'true' =   generating a meanfield
      logical bonn       ! 'true' = Bonn-approach for spin-other orbit
      logical breit      ! if breit is set, BREIT-PAULI only
      logical SAMEORB    ! parameter for same-orbit only
      logical AIMP       ! parameter to delete CORE for AIMP
      logical oneonly    ! parameter to use only oneelectron integrals
      character*4  symmetry
#include "datapow.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      common /ipowxyz/ ipowxyz(3,-Lmax:Lmax,0:Lmax)
c##########################################################################
cbs  #####################################################################
cbs         version with all angular integrals in memory
c         keep=.true.
cbs  #####################################################################
cbs         version without  all angular integrals in memory
          keep=.false.
cbs  #####################################################################
cbs         version without  all cartesian integrals in memory
          keepcart=.false.
cbs  #####################################################################
cbs         version with all cartesian integrals in memory
c         keepcart=.true.
cbs  #####################################################################
      ifinite=0
cbs   initialize tables with double facultatives...
      call inidf
cbs   move some powers of x,y,z to the right place   BEGIN
cbs   check if Lpowmax is high enough..
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
cbs   move some powers of x,y,z to the right place   END
cbs   read the input
      call readbas(Lhigh,makemean,bonn,breit,
     *symmetry,sameorb,AIMP,oneonly,ncont4,numballcart,IN,
     *ifinite)
cbs
      icartdim=mxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)*Lmax
      ionecontrdim=mxcontL*MxcontL*(2*Lmax+1)*3*Lmax
      ioneoverR3dim=Lmax*(MxprimL*MxprimL+MxprimL)/2
      ipowexpdim=MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*(Lmax+Lmax+6)
      icoulovlpdim=MxprimL*MxprimL*(Lmax+1)*(Lmax+1)*10
      Call GetMem('oneoverR3','Allo','Real',ioneoverR3,ioneoverR3dim)
      Call GetMem('cartone','Allo','Real',icartx,3*icartdim)
      Call GetMem('onecontr','Allo','Real',ionecontr,ionecontrdim)
      Call GetMem('powexp','Allo','Real',ipowexp,ipowexpdim)
      Call GetMem('coulovlp','Allo','Real',icoulovlp,icoulovlpdim)
      call dzero(work(ioneoverR3),ioneoverR3dim)
      call dzero(work(icartx),3*icartdim)
      call dzero(work(ionecontr),ionecontrdim)
      call dzero(work(ipowexp),ipowexpdim)
      call dzero(work(icoulovlp),icoulovlpdim)
      icarty=icartx+icartdim
      icartz=icarty+icartdim
cbs
cbs
 123  if (ifinite.eq.2) call finite
cbs
cbs
! Lhigh is the highest l-value in the basis set
      if (makemean.and.(.not.oneonly).and.ifinite.le.1)
     *call getAOs(Lhigh)
      call genpowers(Lhigh,work(ipowexp),work(icoulovlp))
!generate powers of exponents and overlaps
cbs   start generating modified contraction coefficients
cbs   generate starting adresses of contraction coefficients  on
cbs   contrarray
      call genstar(Lhigh)
cbs   generate ovlp of normalized primitives
      call genovlp(Lhigh,work(icoulovlp))
      do lrun=0,Lhigh
cbs   cont(L) arranges all the contraction coefficients for a given L-value
cbs   and renormalizes them
      call cont(lrun,breit,ifinite)
      enddo
cbs
cbs        beginning the angular part
      if (.not.oneonly) then
CBS   write(6,*) '***************************************************'
CBS   write(6,*) '********   beginning the 2e-part ******************'
CBS   write(6,*) '***************************************************'
cbs
cbs  #####################################################################################
cbs  #####################################################################################
cbs  #####################################################################################
cbs
cbs
      idim1=(2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)*(2*Lmax+1)
      idim2=(Lmax+1)*(Lmax+1)*(Lmax+1)*(Lmax+1)
      Call GetMem('preY','Allo','Real',ipreY,idim1)
      Call GetMem('preXZ','Allo','Real',ipreXZ,idim1)
      Call GetMem('icheckxy','Allo','Inte',iicheckxy,idim2)
      Call GetMem('icheckz','Allo','Inte',iicheckz,idim2)
      Call GetMem('interxyz','Allo','Inte',iinterxyz,16*idim2)
      Call GetMem('isgnprod','Allo','Inte',iisgnprod,idim1)
      call angular(Lhigh,keep,keepcart,makemean,bonn,breit,
     *sameorb,ifinite,work(icartx),work(icarty),work(icartz),
     *work(ipowexp),work(icoulovlp),work(ipreXZ),work(ipreY),
     *iwork(iicheckxy),iwork(iicheckz),iwork(iinterxyz),
     *iwork(iisgnprod)) ! subroutine for angular part
      Call GetMem('isgnprod','Free','Inte',iisgnprod,idim1)
      Call GetMem('interxyz','Free','Inte',iinterxyz,16*idim2)
      Call GetMem('icheckz','Free','Inte',iicheckz,idim2)
      Call GetMem('icheckxy','Free','Inte',iicheckxy,idim2)
      Call GetMem('preXZ','Free','Real',ipreXZ,idim1)
      Call GetMem('preY','Free','Real',ipreY,idim1)
      endif
      if (ifinite.eq.1) then ! redo everything for finite core
CBS   write(6,*) 'once more the two-electron integrals'
      ifinite=2
      goto 123
      endif
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################
CBS   write(6,*) '***************************************************'
CBS   write(6,*) '*******   beginning the 1-electron-part  **********'
CBS   write(6,*) '***************************************************'
cbs    the one-electron spin-orbit integrals
      call gen1overR3(Lhigh,work(ioneoverR3))
! 1/r**3  for normalized functions
      call contandmult(Lhigh,makemean,AIMP,oneonly,numballcart,LUPROP,
     *ifinite,work(icartx),work(icarty),work(icartz),work(ionecontr),
     *work(ioneoverR3),iCenter)
cbs   multiplies radial integrals with l,m-dependent
cbs   factors and contraction coefficients
      Call GetMem('coulovlp','free','Real',icoulovlp,icoulovlpdim)
      Call GetMem('powexp','free','Real',ipowexp,ipowexpdim)
      Call GetMem('onecontr','free','Real',ionecontr,ionecontrdim)
      Call GetMem('cartone','free','Real',icartx,3*icartdim)
      Call GetMem('oneoverR3','free','Real',ioneoverR3,ioneoverR3dim)
CBS   write(6,*) '***************************************************'
CBS   write(6,*) '*******   end of  the 1-electron-part    **********'
CBS   write(6,*) '***************************************************'
cbs ########################################################################################
cbs ########################################################################################
cbs ########################################################################################
      Return
      end
      subroutine finite
cbs
cbs   subroutine to set up parameters for finite nucleus. The s-functions are replaced
cbs   by just one exponent which models the nucleus.
cbs
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
      common /nucleus/ charge,Exp_finite
      noccorb(0)=1
      do l=1,lmax_occ
      noccorb(l)=0
      enddo
      occup(1,0)=-charge
      nprimit_keep=nprimit(0)
      ncontrac_keep=ncontrac(0)
      nprimit(0)=1
      ncontrac(0)=1
      exponents(1,0)=0.5d0*Exp_finite
      return
      end
