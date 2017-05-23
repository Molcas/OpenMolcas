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
      subroutine gencoul(l1,l2,l3,l4,makemean,
     *bonn,breit,sameorb,cont4SO,cont4OO,icont4,powexp
     *,coulovlp)
      implicit real*8(a-h,o-z)
cbs   SUBROUTINE to generate all required radial
cbs   integrals for the four angular momenta l1-l4
#include "para.fh"
#include "param.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      logical makemean,bonn,breit,sameorb
      dimension cont4SO(*),cont4OO(*),
     *powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5))
     *,coulovlp(*)
      max1=1  !starting values for limits of precalculated
c             ! powers of function Cfunct(X)
      max2=1
cbs   first of all, this routine determines, for which L
cbs   values the radial integrals have to be solved
cbs   initialize the number of blocks for the different
cbs   l-combinations
cbs   no (ss|ss) contributions
      if (l1.eq.0.and.l2.eq.0.and.l3.eq.0.and.l4.eq.0) return
c      ! no integrals for <ss|ss>
      if (makemean) then
            nblock=1  ! sp sp are the first, so the first block
            Lstarter(1)=1
      else
      Call SysAbendMsg('gencoul',
     & 'only mean-field with this version',' ')
      endif
cbs   keep track of L-values for later purposes
      Lvalues(1)=l1
      Lvalues(2)=l2
      Lvalues(3)=l3
      Lvalues(4)=l4
cbs   now nanz is given the new value
      nanz=ncontrac(l1)*ncontrac(l2)*ncontrac(l3)*ncontrac(l4)
      nprimprod=nprimit(l1)*nprimit(l2)*nprimit(l3)*nprimit(l4)
      Call GetMem('QUOT1','Allo','Real',iquot1,7*nprimprod)
      iquot2=iquot1+nprimprod
      iquotp1=iquot2+nprimprod
      iquotp2=iquotp1+nprimprod
      iprim=iquotp2+nprimprod
      iscr1=iprim+nprimprod
      iscr2=iscr1+nprimprod
c
      call initfrac(nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4),
     *Work(iquot1),Work(iquot2),exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4))
cbs   prepare the powers needed for cfunctx
c
c
c     There are seven different CASES of integrals following
c       (   A  --  C)
c
c     The structure is the same for all cases, therefore comments can be found only on case A
c
c
c
cbs   ###########################################################################################################
cbs   the (+2) cases          CASE A
cbs   ##########################################################################################################
      incl1=1  !  Those increments define the case
      incl3=1
cbs   determine the possible L-values for the integrals by checking for triangular equation
c
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
c
cbs   returns first and last L-values (Lanf,Lend), for which
cbs   radial integrals have to be calculated
      if(Lend-Lanf.ge.0) then
cbs   if there are blocks
              Lblocks(1)=(Lend-Lanf)/2+1 ! L increases in steps of 2,
cbs                                       due to parity conservation
              Lfirst(1)=Lanf
              Llast(1)=Lend
      else
              Lblocks(1)=0
      endif
      if (Lblocks(1).gt.0) then    ! integrals have to be calculated
cbs### check, whether integrals fit on array ################
      if  (Lstarter(1)+nanz*Lblocks(1).gt.icont4) then
      write(6,*) 'end at: ',Lstarter(1)+nanz*Lblocks(1)
       Call SysAbendMsg('gencoul','increase icont4 in amfi.f',' ')
      endif
cbs### check, whether integrals fit on array ################
      istart=Lstarter(1)
c      ! gives the address, where to write the contracted integrals
cbs   ipow1 and ipow2 are the the numbers of powers in the prefactor
cbs   of the function Cfunct
cbs   now loop over possible L-values
      do Lrun= Lfirst(1),Llast(1),2
                      ipow1=2+(l2+l4+Lrun)/2
                      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
cbs   those powers have to be generated...
      call getpow(ipow1,Work(iquot1),Work(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
cbs   those powers have to be generated...
      call getpow(ipow2,Work(iquot2),Work(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
c     in buildcoul the radial integrals are calculated
                      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *         Lrun,Work(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *                nprimit(l4),
     *                exponents(1,l1),exponents(1,l2),
     *                exponents(1,l3),exponents(1,l4),
     *                powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *                Work(iquotp1),Work(iquotp2),coulovlp)
cbs   in the contcas_ routines the integrals are contracted, including exponents as prefactors...
                      if (bonn.or.breit.or.sameorb) then
                      call contcasASO(l1,l2,l3,l4,istart,Work(iprim),
     *           Work(iscr1),Work(iscr2),cont4SO)
                else
                      call contcasASO(l1,l2,l3,l4,istart,Work(iprim),
     *           Work(iscr1),Work(iscr2),cont4SO)
                call contcasAOO(l1,l2,l3,l4,istart,Work(iprim),
     *           Work(iscr1),Work(iscr2),cont4OO)
                endif
                istart=istart+nanz! start for next block  contr integr.
      enddo
      endif
cbs   ##########################################################################################################
cbs   the (0) cases         CASE  B
cbs   ##########################################################################################################
      incl1=0
      incl3=0
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then
      Lblocks(2)=(Lend-Lanf)/2+1
      Lfirst(2)=Lanf
      Llast(2)=Lend
      Lblocks(3)=(Lend-Lanf)/2+1
      Lfirst(3)=Lanf
      Llast(3)=Lend
      else
      Lblocks(2)=0
      Lblocks(3)=0
      endif
      Lstarter(2)=Lstarter(1)+
     *nanz*Lblocks(1)
      Lstarter(3)=Lstarter(2)+
     *nanz*Lblocks(2)
cbs   primitive integrals are the same for type 2 and 3  !!!!!
      if (Lblocks(2).gt.0) then
cbs### check, whether integrals fit on array ################
      if  (Lstarter(2)+2*nanz*Lblocks(2).gt.icont4) then
      write(6,*) 'end at: ',Lstarter(2)+2*nanz*Lblocks(2)
       Call SysAbendMsg ('gencoul','increase icont4 in amfi.f',' ' )
      endif
cbs### check, whether integrals fit on array ################
      istart=Lstarter(2)
      istart2=Lstarter(3)
      do Lrun= Lfirst(2),Llast(2),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,Work(iquot1),Work(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,Work(iquot2),Work(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,Work(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4),
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *Work(iquotp1),Work(iquotp2),coulovlp)
      if (bonn.or.breit.or.sameorb) then
      call contcasB1SO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      else
      call contcasB1SO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      call contcasB2SO(l1,l2,l3,l4,istart2,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      Call contcasB1OO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4OO)
      Call contcasB2OO(l1,l2,l3,l4,istart2,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4OO)
      endif
      istart=istart+nanz
      istart2=istart2+nanz
      enddo
      endif
cbs   ##########################################################################################################
cbs   the (-2) cases      CASE C
cbs   ##########################################################################################################
      if (l1.eq.0.or.l3.eq.0) then
      Lblocks(4)=0
      else
      incl1=-1
      incl3=-1
      call getlimit(l1+incl1,l2,l3+incl3,l4,Lanf,Lend)
      if(Lend-Lanf.ge.0) then
      Lblocks(4)=(Lend-Lanf)/2+1
      Lfirst(4)=Lanf
      Llast(4)=Lend
      else
      Lblocks(4)=0
      endif
      endif
      Lstarter(4)=Lstarter(3)+
     *nanz*Lblocks(3)
      if (Lblocks(4).gt.0) then
cbs### check, whether integrals fit on array ################
      if  (Lstarter(4)+nanz*Lblocks(4).gt.icont4) then
      write(6,*) 'end at: ',Lstarter(4)+nanz*Lblocks(4)
      Call SysAbendMsg('gencoul', 'increase icont4 in amfi.f', ' ' )
      endif
cbs### check, whether integrals fit on array ################
      istart=Lstarter(4)
      do Lrun= Lfirst(4),Llast(4),2
      ipow1=2+(l2+l4+Lrun)/2
      ipow2=2+(l1+l3+incl1+incl3+Lrun)/2
      call getpow(ipow1,Work(iquot1),Work(iquotp1),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call getpow(ipow2,Work(iquot2),Work(iquotp2),
     *nprimit(l1),nprimit(l2),nprimit(l3),nprimit(l4))
      call buildcoul(l1,l2,l3,l4,incl1,incl3,
     *Lrun,Work(iprim),nprimit(l1),nprimit(l2),nprimit(l3),
     *nprimit(l4),
     *exponents(1,l1),exponents(1,l2),
     *exponents(1,l3),exponents(1,l4),
     *powexp(1,1,l3,l1,lrun),powexp(1,1,l4,l2,lrun),
     *Work(iquotp1),Work(iquotp2),coulovlp)
      if (bonn.or.breit.or.sameorb) then
      call contcasCSO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      else
      call contcasCSO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4SO)
      call contcasCOO(l1,l2,l3,l4,istart,Work(iprim),
     *Work(iscr1),Work(iscr2),cont4OO)
      endif
      istart=istart+nanz
      enddo
      endif
      Call GetMem('QUOT1','Free','Real',iquot1,7*nprimprod)
      return
      end
