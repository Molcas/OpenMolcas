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
      Subroutine contandmult(Lhigh,makemean,AIMP,oneonly,numballcart,
     &                       LUPROP,ifinite,onecart,
     &                       onecontr,oneoverR3,iCenter)
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "ired.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Real*8, Allocatable:: Dummy(:), OCA(:,:), OCA2(:,:), OCA3(:,:)
      logical makemean,AIMP,oneonly
      character*8 xa,ya,za
      dimension xa(4),ya(4),za(4),
     *          onecart(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,3),
     *          onecontr(mxcontL,MxcontL,-Lmax:Lmax,3,Lmax),
     *          oneoverR3((MxprimL*MxprimL+MxprimL)/2,Lmax)
      common /nucleus/ charge,Exp_Finite
*
      IPNT(I,J)=(J*J-J)/2+I
*
cbs   get back the real number of functions for the finite nucleus
      if (ifinite.eq.2) ncontrac(0)=ncontrac_keep
c###############################################################################
cbs   subroutine to contract radial one-electron integrals
cbs   and multiply them with angular factors
c###############################################################################
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1SPNORB'
      ya(4)='Y1SPNORB'
      ZA(4)='Z1SPNORB'
c
cbs   clean the arrays for cartesian integrals
C
      length3=(numbalLcart*numbalLcart+numbalLcart)/2
      Call mma_allocate(OCA,Length3,3,Label='OCA')
      Call mma_allocate(OCA2,Length3,3,Label='OCA2')
      Call mma_allocate(Dummy,MxContL**2,Label='Dummy')
      Dummy(:)=0.0D0
      OCA(:,:)=0.0D0
      OCA2(:,:)=0.0D0
c
c
c
c
cbs   one-electron-integrals:
cbs   1. index: number of first contracted function
cbs   2. index: number of second contracted function
cbs   3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
cbs   4. index: L-value
cbs    onecart(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,1),
cbs    onecart(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,2),
cbs    onecart(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax,3)
c
c
c
cbs   generate one-electron integrals for all L greater/equal 1
      if (ifinite.eq.2) charge=0d0 ! nuclear integrals
cbs                                  are modelled for finite nucleus somewhere else
      do L=1,Lhigh
              call contone(L,oneoverr3(1,L),onecontr(1,1,-Lmax,1,L),
     *                     Lmax,contrarray(iaddtyp3(L)),nprimit(L),
     *                     ncontrac(L),MxcontL,Dummy,
     *                     onecart(1,1,1,L,1),
     *                     onecart(1,1,1,L,2),
     *                     onecart(1,1,1,L,3),
     *                     charge,oneonly)
      Enddo
c
cbs   ***********************************************************************
cbs   now move all integrals to one big arrays for X,Y,Z
cbs   ***********************************************************************
      do Lrun=1,Lhigh  !loop over L-values (integrals are diagonal in L)
      mrun=0
      do Msec=-Lrun,Lrun    ! cartesian M-values  (Mfirst,Msec) with
      do Mfirst=-Lrun,Msec  ! Mfirst <= Msec (actually '=' does never
c                             appear as there is no L-component  in Ag
C
c
cbs   determine  if L_X L_Y or L_Z
        ipowx=ipowxyz(1,mfirst,Lrun)+ipowxyz(1,msec,Lrun)
        ipowy=ipowxyz(2,mfirst,Lrun)+ipowxyz(2,msec,Lrun)
        ipowz=ipowxyz(3,mfirst,Lrun)+ipowxyz(3,msec,Lrun)
c
        mrun=mrun+1
cbs     now determine the irreducable representations
        iredfirst=iredLM(Mfirst,Lrun)
        iredsec=iredLM(Msec,Lrun)
cbs     check out which IR is the lower one.
        if (iredfirst.le.iredsec) then
*
cbs     calculate shift to get to the beginning of the block
           iredired= shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)
     *             + incrlm(Mfirst,Lrun)*itotalperIR(iredsec)
     *             + incrLM(Msec,Lrun)
           if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     *         mod(ipowz,2).eq.1) then
              do icartfirst=1,ncontrac(Lrun) ! loop first index
                 do icartsec=1,ncontrac(Lrun)   ! loop second index
                    oca(iredired+icartsec,1)=oca(iredired+icartsec,1)
     *                +onecart(icartfirst,icartsec,mrun,Lrun,1)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredsec)
             enddo
          endif
          if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     *        mod(ipowz,2).eq.1) then
              do icartfirst=1,ncontrac(Lrun) ! loop first index
                 do icartsec=1,ncontrac(Lrun)   ! loop second index
                    oca(iredired+icartsec,2)=oca(iredired+icartsec,2)
     *                +onecart(icartfirst,icartsec,mrun,Lrun,2)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredsec)
              enddo
           endif
           if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     *         mod(ipowz,2).eq.0) then
              do icartfirst=1,ncontrac(Lrun) ! loop first index
                 do icartsec=1,ncontrac(Lrun)   ! loop second index
                    oca(iredired+icartsec,3)=oca(iredired+icartsec,3)
     *                +onecart(icartfirst,icartsec,mrun,Lrun,3)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredsec)
              enddo
           endif
        elseif (iredfirst.gt.iredsec) then
cbs        In this case, indices are exchanged with respect to former
cbs        symmetry of blocks. Therefore, there will be a minus sign
c
cbs        calculate shift to get to the beginning of the block
           iredired=shiftIRIR((iredfirst*iredfirst-iredfirst)/2+iredsec)
     *             + incrLM(Msec,Lrun)*itotalperIR(iredfirst)
     *             + incrLM(Mfirst,Lrun)
           if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     *         mod(ipowz,2).eq.1) then
              do icartsec=1,ncontrac(Lrun) !loopsecond index
                 do icartfirst=1,ncontrac(Lrun) !loop first index
                    oca(iredired+icartfirst,1)=
     *                oca(iredired+icartfirst,1)
     *               -onecart(icartsec,icartfirst,mrun,Lrun,1)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredfirst)
              enddo
           endif
           if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     *         mod(ipowz,2).eq.1) then
              do icartsec=1,ncontrac(Lrun) !loop second index
                 do icartfirst=1,ncontrac(Lrun) !loop first index
                    oca(iredired+icartfirst,2)=
     *                oca(iredired+icartfirst,2)
     *               -onecart(icartsec,icartfirst,mrun,Lrun,2)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredfirst)
              enddo
           endif
           if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     *         mod(ipowz,2).eq.0) then
              do icartsec=1,ncontrac(Lrun) !loop  second index
                 do icartfirst=1,ncontrac(Lrun) !loop first index
                    oca(iredired+icartfirst,3)=
     *                oca(iredired+icartfirst,3)
     *               -onecart(icartsec,icartfirst,mrun,Lrun,3)
                 enddo
cbs              shift pointer by number of functions in IR
                 iredired=iredired+itotalperIR(iredfirst)
              enddo
           endif
        endif
      enddo
      enddo
      enddo
C
cbs   copy integrals on arrays with no symmetry blocking at all
cbs   which means huge triangular matrices
      irun=0
      do norb2=1,numballcarT
         ired2=iredoffunctnew(norb2)
         norbsh2=norb2-shiftIRED(ired2)
         do norb1=1,norb2
            ired1=iredoffunctnew(norb1)
            norbsh1=noRb1-shiftIRED(ired1)
            irun=irun+1
            iredired=shiftIRIR((ired2*ired2-ired2)/2+ired1)
            if (ired1.ne.ired2) then
               oca2(irun,1)=
     &            oca(iredired+norbsh2+(norbsH1-1)*itotalperIR(IREd2),1)
               oca2(irun,2)=
     &            oca(iredired+norbsh2+(norbsH1-1)*itotalperIR(IREd2),2)
               oca2(irun,3)=
     &            oca(iredired+norbsh2+(norbsH1-1)*itotalperIR(IREd2),3)
            else
               oca2(irun,1)=
     &            oca(iredired+norbsh2*(norbsH2-1)/2+norbsh1,1)
               oca2(irun,2)=
     &            oca(iredired+norbsh2*(norbsH2-1)/2+norbsh1,2)
               oca2(irun,3)=
     &            oca(iredired+norbsh2*(norbsH2-1)/2+norbsh1,3)
            endif
         Enddo
      enddo
      if (.not.AIMP) then
c     write a hermit-like file   b.s. 4.10.96
CBS   write(6,*) 'number of orbitals ',numbalLcarT
CBS   write(6,*) 'length of triangular matrix ', length3
CBS    This was removed and will be done in SEWARD
CBS           OPEN(LUPROP,STATUS='UNKNOWN',FORM='UNFORMATTED',
CBS  *        FILE='AOPROPER_MF')
CBS    rewind LUPROP
              write(LUPROP)  iCenter
              write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *        i=1,numbofsym),
     *        numballcart,(Loffunction(I),I=1,numballcart),
     *        (Moffunction(I),I=1,numballcart),
     *        Lhigh,(ncontrac(I),I=0,Lhigh)
              write(LUPROP) (oca2(irun,1),irun=1,length3)
              write(LUPROP)  Ya
              write(LUPROP) (oca2(irun,2),irun=1,length3)
              write(LUPROP)  Za
              write(LUPROP) (oca2(irun,3),irun=1,length3)
CBS   close(luprop)
      else
cbs   reorder for AIMP
cbs   write(6,*) 'reorder integrals for AIMP'
      length3=ikeeporb*(ikeeporb+1)/2
      Call mma_allocate(OCA3,length3,3,Label='OCA3')
      OCA3(:,:)=0.0D0
cbs   write(6,*) 'number of orbitals ',ikeeporb
cbs   write(6,*) 'length of triangular matrix ', length3
      do irun2=1,ikeeporb
         do irun1=1,irun2
            ind2=ikeeplist(irun2)
            ind1=ikeeplist(irun1)
            ipntold=ipnt(ind1,ind2)
            ipntnew=ipnt(irun1,irun2)
*
            oca3(ipntnew,1)=oca2(ipntold,1)
            oca3(ipntnew,2)=oca2(ipntold,2)
            oca3(ipntnew,3)=oca2(ipntold,3)
         enddo
      enddo
CBS   write(6,*) 'transfered to new blocks'
CBS   Luprop=19
CBS   OPEN(LUPROP,STATUS='UNKNOWN',FORM='UNFORMATTED',
CBS  *     FILE='AOPROPER_MF')
CBS   rewind LUPROP
*
      write(LUPROP)  iCenter
      write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *               i=1,numbofsym),
     *               ikeeporb,(Loffunction(ikeeplist(i)),i=1,ikeeporb),
     *               (Moffunction(ikeeplist(i)),I=1,ikeeporb),
     *               Lhigh,((NContrac(I)-icore(I)),I=0,Lhigh)
      write(LUPROP) (oca3(irun,1),irun=1,length3)
      write(LUPROP)  Ya
      write(LUPROP) (oca3(irun,2),irun=1,length3)
      write(LUPROP)  Za
      write(LUPROP) (oca3(irun,3),irun=1,length3)
*
      Call mma_deallocate(OCA3)
CBS   close(luprop)
      endif
cbs
cbs   that is it!!
cbs
      Call mma_deallocate(OCA2)
      Call mma_deallocate(OCA)
      Call mma_deallocate(Dummy)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(makemean)
      end
