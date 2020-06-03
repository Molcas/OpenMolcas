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
      subroutine angular(Lhigh,keep,keepcart,makemean,bonn,
     *                   breit,sameorb,ifinite,
     *                   onecartx,onecarty,onecartz,powexp,coulovlp,
     *                   preXZ,preY,icheckxy,icheckz,interxyz,isgnprod)

c
cbs   COMBINES THE RADIAL INTEGRALS WITH THE ANGULAR FACTORS
c
cbs   if keep=.true. then
cbs   all the integrals will be kept in memory.
cbs   Perhaps, there will be the option to make the
cbs   transformation to the cartesian basis-sets
cbs   everytime, they are required.
cbs   Therefore, the integrals are kept in memory and
cbs   can be further transformed, whenever required.
cbs   in order not to waste to much memory, the atomic
cbs   integrals are thrown away after each l,l,l,l-block
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: ConOO(:), ConSO(:), CartOO(:), CartSO(:)
      logical keep,keepcart,makemean,bonn,
     *        breiT,sameorb,cleaner,NFINI
cbs   NFINI means not finite nucleus
      dimension l2block(0:Lmax,0:Lmax,0:Lmax,0:Lmax),
     *          onecartX(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *          onecartY(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *          onecartZ(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *powexp(MxprimL,MxprimL,0:Lmax,0:Lmax,0:(Lmax+Lmax+5)),coulovlp(*),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *icheckxy(0:Lmax,0:Lmax,0:Lmax,0:Lmax),
     *icheckz(0:Lmax,0:Lmax,0:Lmax,0:Lmax),
     *interxyz(16,0:Lmax,0:Lmax,0:Lmax,0:Lmax),
     *isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
cbs ####################################################################
cbs   some preparation of factors needed later on..                    #
cbs ####################################################################
      ipnt(i,j)=(max(i,j)*max(i,j)-max(i,j))/2+min(i,j)
      roottwo=sqrt(2d0)
cbs   calculate some prefactors that will be needed quite often
      call prefac(Lmax,preroots,clebsch)
      if (ifinite.ne.2) then
cbs      clean array for one electron integrals
         onecartX(:,:,:,:)=0.0D0
         onecartY(:,:,:,:)=0.0D0
         onecartZ(:,:,:,:)=0.0D0
         NFINI=.true.
      else
         NFINI=.false.
      endif
*
cbs   generate an array with sign for (even/odd) m-values
      isignM(0)=1
      do I=2,Lmax,2
      isignM(I)=1
      isignM(-I)=1
      enddo
      do I=1,Lmax,2
      isignM(I)=-1
      isignM(-I)=-1
      enddo
      call genprexyz(preXZ)
      call genprexyz2(preXZ)
      call genprexyz3(preXZ)
      call genprexyz4(preXZ)
      call genprexyz5(preXZ)
      call genprexyz6(preY,preXZ)
      call genprexyz7(preXZ)
      call genprexyz8(preXZ)
      call genprexyz9(preXZ)
      call genprexyz10(preXZ)
      call genprexyz11(preY)
      call genprexyz12(preY)
      call genprexyz13(icheckxy)
      call genprexyz14(icheckz,interxyz)
      call genprexyz15a(icheckxy,icheckz,interxyz)
cbs #####################################################################
cbs   isgnprod gives the sign due to powers (-1)**M  this are again
cbs   angular m-values
cbs #####################################################################
      do M4=-Lmax,Lmax
      if (M4.gt.0) then
      inter4=isignM(M4)
      else
      inter4=1
      endif
      do M3=-Lmax,Lmax
      if (M3.gt.0) then
      inter3=inter4*isignM(M3)
      else
      inter3=inter4
      endif
      do M2=-Lmax,Lmax
      if (M2.gt.0) then
      inter2=inter3*isignM(M2)
      else
      inter2=inter3
      endif
      do M1=-Lmax,Lmax
      if (M1.gt.0) then
      isgnprod(m1,m2,m3,m4)=inter2*isignM(M1)
      else
      isgnprod(m1,m2,m3,m4)=inter2
      endif
      enddo
      enddo
      enddo
      enddo
cbs #####################################################################
cbs   some preparation of factors needed later on..  finished           #
cbs #####################################################################
c
c
c
cbs   counter for total number of cartesian integrals                   !  set some counters
      numbcart=0
cbs   same orbit integrals integrals  on carteXSO carteYSO and carteSO
cbs   other orbit integrals  on carteXOO carteYOO and carteOO
      iangfirst=0 ! first block of angular integrals
cbs #####################################################################
cbs   loop over all (l,l,l,l) blocks generated in the radial part       #
cbs #####################################################################
      do lrun4=0,Lmax
      do lrun3=0,Lmax
      do lrun2=0,Lmax
      do lrun1=0,Lmax
      l2block(lrun1,lrun2,lrun3,lrun4)=0
      enddo
      enddo
      enddo
      enddo
cbs   loop over all possible < l1 l2, l3 l4 > blocks
CBS   write(6,'(A)') '   L1   L2   L3   L4'
      do l1=0,Lhigh   ! improving is probably possible...
      do l2=0,Lhigh
      do l3=0,l1
      do l4=0,l2
cbs   check parity
      if (mod(l1+l2+l3+l4,2).eq.0) then
cbs   check that Lleft and Lright do not always differ by more than one
cbs   a difference of two means two spin flips and is therefore not allowed
      Lleftmax=l1+l2
      Lrightmax=l3+l4
      Lleftmin=iabs(l1-l2)
      Lrightmin=iabs(l3-l4)
      if ((Lrightmin-Lleftmax.le.1.and.Lrightmax-Lleftmin.gt.-1).or.
     *(Lleftmin-Lrightmax.le.1.and.Lleftmax-Lrightmin.gt.-1)) then
cbs   additional check for mean-field
      if ((l1.eq.l3.and.l2.eq.l4).or.(l1.eq.l2.and.l3.eq.l4)) then
      if (l1+l3.ne.0) then
CBS   write(6,'(4I5)') l1,l2,l3,l4
CBS   now I determine the size of the angular integral arrays
        jblock=0
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then
        jblock=jblock+1
        endif
        endif
        enddo
        enddo
        enddo
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3
        m4=m1+m2-m3
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !  all m=0 make no sense
        if (iabs(m4).le.l4)  then
        jblock=jblock+1
        endif
        endif
        endif
        enddo
        enddo
        enddo
CBS   done !!
cbs     number of contracted integrals for each block
        ncont=ncontrac(l1)*ncontrac(l2)*
     *  ncontrac(l3)*ncontrac(l4)
      mxangint=jblock*ncont
cbs   determine the size icont4 for the radial integrals
      call gencoulDIM(l1,l2,l3,l4,makemean,bonn,breit,
     *                sameorb,icont4)
      Call GetMem('ANGSO','Allo','Real',iangSO,2*mxangint)
      iangOO=iangSO+mxangint
      Call mma_allocate(CartSO,nCont,Label='CartSO')
      Call mma_allocate(CartOO,nCont,Label='CartOO')
      Call mma_allocate(ConSO,iCont4,Label='ConSO')
      Call mma_allocate(ConOO,iCont4,Label='ConOO')
*
      call gencoul(l1,l2,l3,l4,makemean,bonn,breit,
     *             sameorb,conSO,conOO,icont4,powexp,coulovlp)
!gen and trans integrals
        l2block(l1,l2,l3,l4)=1  ! can be used for getting the
cbs   local counter for integral adresses
        mblock=0 ! counter of (m,m,m,m)-blocks for (l1,l2,l3,l4)
cbs     if keep is set to false, the angular integrals are
cbs     thrown away after each block of l-values
cbs     which means integrals start at address 0
        if (.not.keep) iangfirst=0
        locstar=iangfirst ! local starting adress counter
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
        do m4=-l4,l4
        mcombina(1,m1,m2,m3,m4)=0  ! will hold type of integrals (1,2,3)
        mcombina(2,m1,m2,m3,m4)=0  ! will hold number of block
        enddo
        enddo
        enddo
        enddo
        do m1=-l1,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3
c####################################################################################
c####################################################################################
c########## the L- -type block to be combined with sigma+ ###########################
c####################################################################################
c####################################################################################
        m4=m1+m2-m3+1
        if (iabs(m4).le.l4) then !the L- -block to  combine with sigma+
cbs     not all m-combinations are needed for the mean-field
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then
        mcombina(1,m1,m2,m3,m4)=1
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then
        write(6,*) 'not enough space allocated for angular integrals'
        write(6,*) 'increase mxangint to at least ',
     *  locstar+ncont
        Call Abend()
        endif
cbs mkangLmin = make_angular_integrals_for_L- type operator
cbs really generates  the angular prefactors and combines them with
cbs the radial integrals
        call mkangLmin(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       work(iangSO+locstar),
     *       work(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       ConSO(Lstarter(1)),
     *       ConSO(Lstarter(2)),
     *       ConSO(Lstarter(3)),
     *       ConSO(Lstarter(4)),
     *       ConOO(Lstarter(1)),
     *       ConOO(Lstarter(2)),
     *       ConOO(Lstarter(3)),
     *       ConOO(Lstarter(4)),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb)
        locstar=locstar+ncont ! increase starting address
        mcombina(2,m1,m2,m3,m4)=mblock  ! set the block number
c####################################################################################
c####################################################################################
c########## the L+ -type block to be combined with sigma- ###########################
c####################################################################################
c####################################################################################
c
c   these integrals are obtained by changing the signs of the m-values.
c   As the integrals are the same, the pointer points to the same integrals...
c
c
        mcombina(1,-m1,-m2,-m3,-m4)=3
        mcombina(2,-m1,-m2,-m3,-m4)=mblock
        endif
        Endif
        enddo
        enddo
        enddo
c####################################################################################
c####################################################################################
c########## the L0 -type block to be combined with sigma0 ###########################
c####################################################################################
c####################################################################################
        do m1=  0,l1
        do m2=-l2,l2
        do m3=-l3,l3
cbs     m4 is more or less fixed by m1-3
        m4=m1+m2-m3 ! the L0-block to be combined with sigma0
cbs     not all m-combinations are needed for the mean-field
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.iabs(m2).eq.iabs(m4)).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.
     *  (iabs(m1).eq.iabs(m2).or.iabs(m3).eq.iabs(m4)))) then
c
        if (m1.ne.0.or.m2.ne.0.or.m3.ne.0) then !all m=0 make no sense
        if (iabs(m4).le.l4)  then
        mcombina(1,m1,m2,m3,m4)=2
        mblock=mblock+1
        if (locstar+ncont.gt.mxangint) then
        write(6,*) 'not enough space allocated for angular integrals'
        write(6,*) 'increase mxangint to at least ',
     *  locstar+ncont
        Call Abend()
        endif
        call mkangL0(Lmax,l1,l2,l3,l4,m1,m2,m3,m4,
     *       work(iangSO+locstar),
     *       work(iangOO+locstar),
     *       Lfirst(1),Llast(1),Lblocks(1),
     *       ncontrac(l1),ncontrac(l2),ncontrac(l3),ncontrac(l4),
     *       ConSO(Lstarter(1)),
     *       ConSO(Lstarter(2)),
     *       ConSO(Lstarter(3)),
     *       ConSO(Lstarter(4)),
     *       ConOO(Lstarter(1)),
     *       ConOO(Lstarter(2)),
     *       ConOO(Lstarter(3)),
     *       ConOO(Lstarter(4)),
     *       preroots,clebsch,scratch4,bonn,breit,
     *       sameorb)
        locstar=locstar+ncont
        mcombina(2,m1,m2,m3,m4)=mblock
        endif
        endif
        endif
        enddo
        enddo
        enddo
cbs  ##################################################################################
cbs  ##################################################################################
cbs     transformation to l,m dependent integrals is finished
cbs  ##################################################################################
c
c
c
c
cbs  ##################################################################################
cbs     begin transformation to cartesian integrals
cbs  ##################################################################################
cbs  ##################################################################################
cbs     check out, which combinations of m-values will
cbs     contribute to cartesian integrals
        do m1=-l1,l1       !
        do m2=-l2,l2  ! these indices now run over the real harmonics
        do m3=-l3,l3  !
        do m4=-l4,l4  !
        mcombcart(1,m1,m2,m3,m4)=0     ! will hold the type  x=1 y=2 z=3
        mcombcart(2,m1,m2,m3,m4)=0     ! will hold the block number
        enddo
        enddo
        enddo
        enddo
        mblockx=0
        mblocky=0
        mblockz=0
        do m3=-l3,l3
        do m4=-l4,l4
cbs     if the l-values are the same : triangular matrix over m-values is sufficient
        if (l1.eq.l3) then
        m1upper=m3
        else
        m1upper=l1
        endif
        if (makemean) m1upper=l1
cbs     if the l-values are the same : triangular matrix over m-values is sufficient
        if (l2.eq.l4) then
        m2upper=m4
        else
        m2upper=l2
        endif
        if (makemean) m2upper=l2
        do m1=-l1,m1upper
        If (l1.eq.l3.and.m1.eq.m3) then ! clean real zeros by symmetry
cbs     this a problem of the spin-other-orbit integrals, as they are by formula
cbs     not antisymmetric in the indices for particle 1.
        cleaner=.true.
        else
        cleaner=.false.
        endif
        do m2=-l2,m2upper
cbs     not all m-combinations are needed for the mean-field
        if ((.not.makemean).or.
     *  (l1.eq.l3.and.l2.eq.l4.and.m2.eq.m4).or.
     *  (l1.eq.l2.and.l3.eq.l4.and.(m1.eq.m2.or.m3.eq.m4))) then
C
        indx=ipowxyz(1,m1,l1)+ipowxyz(1,m2,l2)+
     *  ipowxyz(1,m3,l3)+ipowxyz(1,m4,l4)
        indy=ipowxyz(2,m1,l1)+ipowxyz(2,m2,l2)+
     *  ipowxyz(2,m3,l3)+ipowxyz(2,m4,l4)
        indz=ipowxyz(3,m1,l1)+ipowxyz(3,m2,l2)+
     *  ipowxyz(3,m3,l3)+ipowxyz(3,m4,l4)
        indx=mod(indx,2)
        indy=mod(indy,2)
        indz=mod(indz,2)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA X      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (indx.eq.0.and.indy.eq.1.and.indz.eq.1.and.
     *  icheckxy(iabs(m1),iabs(m2),iabs(m3),iabs(m4)).gt.0) then
! Y*Z ->  transforms like  L_x (B1)
cbs     integrals for sigma_x
        mblockx=mblockx+1
        mcombcart(1,m1,m2,m3,m4)=1
        mcombcart(2,m1,m2,m3,m4)=mblockx
        call tosigX(m1,m2,m3,m4,work(iangSO+iangfirst),
     *              mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *              ncontrac(l4),CartSO,preXZ,
     *              interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),
     *              isgnprod,cleaner)
c
        if (.not.bonn.and.(.not.breiT))
     *  call tosigX(m1,m2,m3,m4,work(iangOO+iangfirst),
     *  mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *  ncontrac(l4),cartOO,preXZ,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then
             if (m2.eq.m4.and.m1.lt.m3.and.
     *       iabs(m1+m3).eq.1.and.l1.ne.0) then
             call two2mean13(CartSO,occup(1,l2),AOcoeffs(1,1,l2),
     *                       onecartx(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *                       ncontrac(l1),ncontrac(l2),noccorb(l2))
             endif
             endif
             if (l1.eq.l2.and.l3.eq.l4) then
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then
             if (m3.lt.m4.and.iabs(m4+m3).eq.1) then
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34a(cartSO,cartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if(NFINI) call two2mean34a(cartSO,cartOO,
     *                                  occup(1,l1),AOcoeffs(1,1,l1),
     *                           onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *                                  ncontrac(l3),ncontrac(l1),
     *                                  noccorb(l2),sameorb)
             endif
             endif
             if (m3.gt.m4.and.iabs(m4+m3).eq.1) then
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartx(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             endif
             endif
             elseif(m3.eq.m4.and.l1.ne.0) then
             if (m1.lt.m2.and.iabs(m1+m2).eq.1) then
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12a(CartSO,cartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             if (m1.gt.m2.and.iabs(m1+m2).eq.1) then
cbs   for the "Bonn-approach"   exchange cartexOO by cartexSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12b(cartSO,CartSO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12b(CartSO,CartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartx(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             endif
             endif
c#######################################################################
c############  mean-field-part #########################################
c#######################################################################
        endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Y      ++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.0.and.indz.eq.1.and.
     *  icheckxy(iabs(m1),iabs(m2),iabs(m3),iabs(m4)).gt.0) then
! X*Z transforms like L_y  (B2)
cbs     integrals for sigma_y
        mblocky=mblocky+1
        mcombcart(1,m1,m2,m3,m4)=2
        mcombcart(2,m1,m2,m3,m4)=mblocky
        call tosigY(m1,m2,m3,m4,work(iangSO+iangfirst),
     *              mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *              ncontrac(l4),CartSO,preY,
     *              interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),
     *              isgnprod,cleaner)
c
        if (.not.bonn.and.(.not.breit))
     *  call tosigY(m1,m2,m3,m4,work(iangOO+iangfirst),
     *              mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *              ncontrac(l4),cartOO,preY,
     *  interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),isgnprod,
     *  cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.iabs(m3-m1).eq.1.and.l1.ne.0) then
             call two2mean13(CartSO,occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartY(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))
             endif
             endif
             if (l1.eq.l2.and.l3.eq.l4) then
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then
             if (m3.lt.m4.and.iabs(m3-m4).eq.1) then
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34a(CartSO,CartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if (NFINI) call two2mean34a(CartSO,CartOO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             endif
             endif
             if (m3.gt.m4.and.iabs(m3-m4).eq.1) then
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartY(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             endif
             endif
             elseif(m3.eq.m4.and.l1.ne.0) then
             if (m1.lt.m2.and.iabs(m1-m2).eq.1) then
cbs   for the "Bonn-approach"   exchange carteOO by carteSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12a(CartSO,CartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             if (m1.gt.m2.anD.Iabs(m1-m2).eq.1) then
cbs   for the "Bonn-approach"   exchange carteYOO by carteYSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12b(CartSO,CartSO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12b(CartSO,CartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartY(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             endif
             endif
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++      SIGMA Z      ++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        elseif (indx.eq.1.and.indy.eq.1.and.indz.eq.0.and.
     *  icheckz(iabs(m1),iabs(m2),iabs(m3),iabs(m4)).gt.0) then
! X*Y transforms like L_z  (A2)
cbs     integrals for sigma_z
        mblockz=mblockz+1
        mcombcart(1,m1,m2,m3,m4)=3
        mcombcart(2,m1,m2,m3,m4)=mblockz
        call tosigZ(m1,m2,m3,m4,work(iangSO+iangfirst),
     *              mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *              ncontrac(l4),CartSO,preXZ,
     *              interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),
     *              isgnprod,cleaner)
c
        if (.not.bonn.and.(.not.breit))
     *  call tosigZ(m1,m2,m3,m4,work(iangOO+iangfirst),
     *              mcombina,ncontrac(l1),ncontrac(l2),ncontrac(l3),
     *              ncontrac(l4),CartOO,preXZ,
     *              interxyz(1,iabs(m1),iabs(m2),iabs(m3),iabs(m4)),
     *              isgnprod,cleaner)
        if (makemean) then ! generate mean-field-contributions
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
             if (l1.eq.l3.and.l2.eq.l4) then
             if (m2.eq.m4.and.m1.lt.m3.
     *       and.m1.eq.-m3.and.l1.ne.0) then
             call two2mean13(CartSO,occup(1,l2),
     *       AOcoeffs(1,1,l2),onecartz(1,1,ipnt(m1+l1+1,m3+l3+1),l1),
     *       ncontrac(l1),ncontrac(l2),noccorb(l2))
             endif
             endif
             if (l1.eq.l2.and.l3.eq.l4) then
             if (m1.eq.m2.and.l3.ne.0.and.l3.ne.l1) then
             if (m3.lt.m4.and.m3.eq.-m4) then
cbs   for the "Bonn-approach"   exchange carteOO by carteSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34a(CartSO,CartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if (NFINI) call two2mean34a(CartSO,CartOO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             endif
             endif
             if (m3.gt.m4.and.m3.eq.-m4) then
cbs   for the "Bonn-approach"   exchange carteOO by carteSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean34b(CartSO,CartSO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             else
             if (NFINI) call two2mean34b(CartSO,CartOO,occup(1,l1),
     *       AOcoeffs(1,1,l1),onecartz(1,1,ipnt(m3+l3+1,m4+l4+1),l3),
     *       ncontrac(l3),ncontrac(l1),noccorb(l2),sameorb)
             endif
             endif
             elseif(m3.eq.m4.and.l1.ne.0) then
             if (m1.lt.m2.and.m1.eq.-m2) then
cbs   for the "Bonn-approach"   exchange carteOO by carteSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12a(CartSO,CartSO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12a(CartSO,CartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             if (m1.gt.m2.and.m1.eq.-m2) then
cbs   for the "Bonn-approach"   exchange carteOO by carteSO
             if (bonn.or.breiT) then
             if (NFINI) call two2mean12b(cartSO,CartSO,
     *       occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             else
             if (NFINI) call two2mean12b(cartSO,cartOO,occup(1,l3),
     *       AOcoeffs(1,1,l3),onecartz(1,1,ipnt(m1+l1+1,m2+l2+1),l1),
     *       ncontrac(l1),ncontrac(l3),noccorb(l3),sameorb)
             endif
             endif
             endif
             endif
c##########################################################################
c############  mean-field-part ############################################
c##########################################################################
        endif
        endif
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        endif ! for check of significance for meanfield.
        enddo
        enddo
        enddo
        enddo
        numbcart=numbcart+(mblockx+mblocky+mblockz)*ncont
cbs   just controlling if x and y integrals have the same number of blocks
      if (mblockx.ne.mblocky) then
      write(6,*)
     *'numbers of integrals for sigma_x and sigma_y not equal!'
      write(6,'(A12,4I3,2(A3,I5))')
     *'l1,l2,l3,l4 ',l1,l2,l3,l4,' X:',mblockx,' Y:',mblocky
      write(6,*) ' check the ipowxyz-array'
      Call Abend()
      endif
cbs   start adresses for the next <ll|ll> block of integrals
      Call GetMem('ANGSO','free','real',iangSO,2*mxangint)
      Call mma_deallocate(CartSO)
      Call mma_deallocate(CartOO)
      Call mma_deallocate(ConSO)
      Call mma_deallocate(ConOO)
      endif
      endif
      endif
      endif
      enddo
      enddo
      enddo
      enddo
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(keepcart)
      end
