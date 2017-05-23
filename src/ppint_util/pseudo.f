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
      Subroutine Pseudo(ai,xi,yi,zi,lit,
     &                  aj,xj,yj,zj,ljt,gout,intmax,lmn1u,
     &                  ccr,zcr,nkcrl,nkcru,lcr,ncr,xc,yc,zc,npot)
************************************************************************
c This collection of subroutines, adapted from the ARGOS suite of programs,
c calculates integrals of an atomic semi-local one-component (i.e., non-
c relativistic or scalar-relativistic) pseudopotential (PP) between two shells
c of primitive cartesian GTOs. The integrals are contained in field gout(*).
c lmn1u is the maximum l value + 1 of the basis functions
c lproju is the maximum l value + 1 of the pseudopotential
c GTO parameters
c
c The GTOs are given by Ni (x-xi)^li (y-yi)^mi (z-zi)^ni exp(-ai(r-ri)^2)
c and Nj (x-xj)^lj (y-yj)^mj (z-zj)^nj exp(-aj(r-rj)^2). The Ni, Nj are
c chosen such that the integrals are normalized to (2l-1)!!(2m-1)!!(2n-1)!!.
c Input data are
c ai, aj : exponents
c xi,yi,zi, xj,yj,zj : centers
c lit = li+mi+ni+1, ljt = lj+mj+nj+1 : angular-momentum values + 1
c PP parameters
c
c The pseudopotential is of the form V_loc + sum_l V_l P_l, where V_loc and the V_l
c are radial potentials and P_l is the projector onto angular-momentum l.
c Both V_loc and the V_l are given by expansions of the form
c V = sum_i ccr(i) ( r^(ncr(i)-2) exp(-zcr(i)*r^2) ),
c where r is the distance from the PP center (xc, yc, zc).
c The PP parameters are assumed to be taken from a list of PPs for several atoms
c containing the following input data:
c npot : number of radial potentials in the list
c ccr(i) : coefficients
c ncr(i), zcr(i) : exponential parameters
c nkcrl(1,k) : lower index i for V_loc at nucleus k
c nkcru(1,k) : upper index i for V_loc at nucleus k
c nkcrl(l+2,k) : lower index i for V_l at nucleus k
c nkcru(l+2,k) : upper index i for V_l at nucleus k
c lcr(k) : maximum l value + 1 of the PP at nucleus k
c k=kcrs : specifies the PP for which integrals are calculated
c xc,yc,zc : position of nucleus kcrs
c
************************************************************************
      implicit real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      parameter(lproju=9)
      parameter (imax=100,kcrs=1)
cAOM      parameter(iptmax=100000)
      parameter(facij=5.56832799683170784522d0)
      real*8 gout(2*intmax)
      Real*8 ccr(npot),zcr(npot), crda(lproju,3), crdb(lproju,3)
      Integer nkcrl(lproju+1,kcrs),nkcru(lproju+1,kcrs),lcr(kcrs),
     &        ncr(imax)
*
cAOM      dimension ipt(55),a(iptmax)
      dimension ipt(55)
*
C     Write (*,*) 'ncr',(ncr(i),i=1,npot)
C     Write (*,*) 'zcr',(zcr(i),i=1,npot)
C     Write (*,*) 'ccr',(ccr(i),i=1,npot)
C     Write (*,*) 'nkcrl',(nkcrl(i,1),i=1,lcr(1)+1)
C     Write (*,*) 'nkcru',(nkcru(i,1),i=1,lcr(1)+1)
C     Write (*,*)
*                                                                      *
************************************************************************
*                                                                      *
* precalculation of tables for PP integrals
* iptmax is the dimension of the work space
*
      ipt(1) = 1
      lmnvmx = (lmn1u*(lmn1u + 1)*(lmn1u + 2))/6
      ipt(9) = ipt(1) +  3*lmnvmx
      ipt(10) = ipt(9) + 35
      ipt(11) = ipt(10) + 35
      ncru = 0
      do  1 i=1,npot
        ncru=max(ncru,ncr(i))
  1   continue
      ndfac = max(4*lmn1u+2*lproju-3,6*lproju+3,4*lmn1u-1,
     &    2*lmn1u+2*lproju+1,4,ncru+4*lmn1u+2*lproju-1)
      ipt(12) = ipt(11) +  ndfac
      ipt(13) = ipt(12) +  (lmn1u*(lmn1u+1))/2
      lmax =  max(1,lmn1u-1 + max(lmn1u-1,lproju))
      l1max2 = (lmax+1)**2
      ipt(14) = ipt(13) +  l1max2
      ipt(15) = ipt(14) +  l1max2
      lmnpwr = (((lmax*(lmax+2)*(lmax+4))/3)*(lmax+3)+
     &    (lmax+2)**2*(lmax+4))/16
      ipt(16) = ipt(15) +  lmnpwr
      ipt(17) = ipt(16) +  lmnpwr
      ipt(18) = ipt(17) +  lmnpwr
      ipt(19) = ipt(18) +  lmnpwr
      ipt(20) = ipt(19) +  3*lproju**2
      l2m1 = 2*lproju-1
      ipt(21) = ipt(20) +  3*l2m1
      ipt(31) = ipt(21) +  3*l2m1
cAOM      if ( ipt(31) .gt. iptmax ) stop 'work space too small'
      ltot1 = lit + ljt - 1
      ltot12 = ltot1**2
      ipt(33) = ipt(31) +  ltot12
      ipt(34) = ipt(33) +  ltot12
      ipt(35) = ipt(34) +  ltot1
      ipt(36) = ipt(35) +  ltot1
      ipt(37) = ipt(36) +  ltot1
cAOM      if ( ipt(37) .gt. iptmax ) stop 'work space too small'
      mproju = 2*lproju+1
      lamau = lit+lproju
      ipt(38) = ipt(37) +  lit*mproju*lamau
      lambu = ljt+lproju
      ipt(39) = ipt(38) +  ljt*mproju*lambu
      ipt(40) = ipt(39) +  ltot1*lambu*lamau
      ipt(41) = ipt(40) +  ljt
      ipt(42) = ipt(41) +  ljt
      ipt(43) = ipt(42) +  ljt
      ipt(44) = ipt(43) +  ljt
      ipt(45) = ipt(44) +  lit
      ipt(46) = ipt(45) +  lit
      ipt(47) = ipt(46) +  ltot1+lproju
      ipt(48) = ipt(47) +  lamau
      ipt(49) = ipt(48) +  lambu
      ipt(50) = ipt(49) +  ltot1
cAOM      ipt(51) = ipt(50) +  lambu*lamau
      lWork = ipt(50) +  lambu*lamau

cAOM      if ( ipt(51) .gt. iptmax ) stop 'work space too small'
      Call GetMem('Wrk','Allo','Real',mWrk,lWork)
      Do I=1,50
        ipt(i)=ipt(i)+mWrk-1
      Enddo
      call lmnvgn_molcas(lmn1u,Work(ipt(1)))
      eps = 1.0d-12
      call cortab_molcas(Work(ipt(12)),Work(ipt(11)),eps,Work(ipt(19)),
     &      Work(ipt(9)),Work(ipt(10)),Work(ipt(13)),Work(ipt(14)),
     &      Work(ipt(15)),Work(ipt(16)),Work(ipt(17)),lmax,lmn1u,lproju,
     &      Work(ipt(20)),Work(ipt(21)),ndfac,Work(ipt(18)))
*                                                                      *
************************************************************************
*                                                                      *
*     integrals for V_loc part of PP
*
      lproju1=lproju+1
      call pseud1_molcas(Work,Work(ipt(31)),ccr,gout,ipt,
     &        Work(ipt(1)),ltot1,ncr,nkcrl,nkcru,Work(ipt(33)),
     &        Work(ipt(34)),Work(ipt(35)),Work(ipt(36)),zcr,lit,ljt,
     &        ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lproju1,crda,crdb)
*                                                                      *
************************************************************************
*                                                                      *
*     integrals for V_l P_l parts of PP
*
      lcru=lcr(kcrs)
      call pseud2_molcas(Work,Work(ipt(37)),Work(ipt(38)),ccr,
     &        gout,ipt,lambu,ltot1,mproju,ncr,nkcrl,nkcru,Work(ipt(39)),
     &        zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,
     &        kcrs,lcru,lproju1,crda,crdb)
      Call GetMem('Wrk','Free','Real',mWrk,lWork)
*                                                                      *
************************************************************************
*                                                                      *
*     final results with normalization factors
*
*     faci=sqrt( (4.D0*ai)**lit/2.D0*sqrt(2.D0*ai) )
*    &   /  ( (4.D0*ai)**((lit-1)/2.D0) * (Two*ai/Pi)**(Three/Four) )
*     facj=sqrt( (4.D0*aj)**ljt/2.D0*sqrt(2.D0*aj) )
*    &   /  ( (4.D0*aj)**((ljt-1)/2.D0) * (Two*aj/Pi)**(Three/Four) )
*     facij=faci*facj
      litot=lit*(lit+1)/2
      ljtot=ljt*(ljt+1)/2
*
      Call DScal_(litot*ljtot,facij,gout,1)
*
*     Do  i=1,litot
*        Do  j=1,ljtot
*          gout(i+(j-1)*litot) = gout(i+(j-1)*litot)*facij
*        End Do
*     End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Reorder integrals from Argos to Molcas order
*
      Call Molcas_Order(GOut,litot,ljtot)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
