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
     &                       LUPROP,ifinite,onecartX,onecartY,onecartZ,
     &                       onecontr,oneoverR3,iCenter)
      implicit real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "ired.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      logical makemean,AIMP,oneonly
      character*8 xa,ya,za
      dimension xa(4),ya(4),za(4),
     *onecartX(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecartY(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecartZ(mxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1),Lmax),
     *onecontr(mxcontL,MxcontL,-Lmax:Lmax,3,Lmax),
     *oneoverR3((MxprimL*MxprimL+MxprimL)/2,Lmax)
      common /nucleus/ charge,Exp_Finite
      IPNT(I,J)=(J*J-J)/2+I
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
      iloca=length3
      Call GetMem('OCAX','Allo','Real',iocax,6*iloca+MxContL*MxContL)
      iocay=iocax+iloca
      iocaz=iocay+iloca
      iocax2=iocaz+iloca
      iocay2=iocax2+iloca
      iocaz2=iocay2+iloca
      idummy=iocaz2+iloca
      call dzero(work(iocax),6*length3+MxContL*MxContL)
c
c
c
c
cbs   one-electron-integrals:
cbs   1. index: number of first contracted function
cbs   2. index: number of second contracted function
cbs   3. index: pointer(m1,m2)    m1< m2 otherwise change sign of integral
cbs   4. index: L-value
cbs    onecartX(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartY(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax),
cbs    onecartZ(mxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1),Lmax)
c
c
c
cbs   generate one-electron integrals for all L greater/equal 1
      if (ifinite.eq.2) charge=0d0 ! nuclear integrals
cbs                                  are modelled for finite nucleus somewhere else
      do L=1,Lhigh
              call contone(L,oneoverr3(1,L),onecontr(1,1,-Lmax,1,L),
     *  Lmax,contrarray(iaddtyp3(L)),nprimit(L),ncontrac(L),
     *        MxcontL,work(idummy),
     *        onecartx(1,1,1,L),onecartY(1,1,1,L),onecartZ(1,1,1,L),
     *        charge,oneonly)
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
cbs     calculate shift to get to the beginning of the block
       iredired=shiftIRIR((iredsec*iredsec-iredsec)/2+iredfirst)
     *                +incrlm(Mfirst,Lrun)*itotalperIR(iredsec)+
     *                incrLM(Msec,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then
                      do icartfirst=1,ncontrac(Lrun) ! loop first index
                      do icartsec=1,ncontrac(Lrun)   ! loop second index
                      work(iocax+iredired+(icartsec-1))=
     *                work(iocax+iredired+(icartsec-1))
     *                +onecartx(icartfirst,icartsec,mrun,Lrun)
                      enddo
cbs             shift pointer by number of functions in IR
                      iredired=iredired+itotalperIR(iredsec)
                      enddo
        endif
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then
                      do icartfirst=1,ncontrac(Lrun) ! loop first index
                      do icartsec=1,ncontrac(Lrun)   ! loop second index
                      work(iocay+iredired+(icartsec-1))=
     *                work(iocay+iredired+(icartsec-1))
     *                +onecarty(icartfirst,icartsec,mrun,Lrun)
                      enddo
cbs             shift pointer by number of functions in IR
                      iredired=iredired+itotalperIR(iredsec)
                      enddo
        endif
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then
                      do icartfirst=1,ncontrac(Lrun) ! loop first index
                      do icartsec=1,ncontrac(Lrun)   ! loop second index
                      work(iocaz+iredired+(icartsec-1))=
     *                work(iocaz+iredired+(icartsec-1))
     *                +onecartz(icartfirst,icartsec,mrun,Lrun)
                      enddo
cbs             shift pointer by number of functions in IR
                      iredired=iredired+itotalperIR(iredsec)
                      enddo
        endif
              elseif (iredfirst.gt.iredsec) then
cbs     In this case, indices are exchanged with respect to former
cbs     symmetry of blocks. Therefore, there will be a minus sign
c
cbs     calculate shift to get to the beginning of the block
            iredired=shiftIRIR((iredfirst*iredfirst-iredfirst)/2+
     *                iredsec)+
     *                incrLM(Msec,Lrun)*itotalperIR(iredfirst)+
     *                incrLM(Mfirst,Lrun)
       if (mod(ipowx,2).eq.0.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.1) then
                      do icartsec=1,ncontrac(Lrun) !loopsecond index
                      do icartfirst=1,ncontrac(Lrun) !loop first index
                      work(iocax+iredired+(icartfirst-1))=
     *                work(iocax+iredired+(icartfirst-1))
     *         -onecartx(icartsec,icartfirst,mrun,Lrun)
                      enddo
cbs             shift pointer by number of functions in IR
                      iredired=iredired+itotalperIR(iredfirst)
                      enddo
        endif
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.0.and.
     * mod(ipowz,2).eq.1) then
                      do icartsec=1,ncontrac(Lrun) !loop second index
                      do icartfirst=1,ncontrac(Lrun) !loop first index
                      work(iocay+iredired+(icartfirst-1))=
     *                work(iocay+iredired+(icartfirst-1))
     *               -onecarty(icartsec,icartfirst,mrun,Lrun)
                      enddo
cbs             shift pointer by number of functions in IR
                      iredired=iredired+itotalperIR(iredfirst)
                      enddo
        endif
       if (mod(ipowx,2).eq.1.and.mod(ipowy,2).eq.1.and.
     * mod(ipowz,2).eq.0) then
                      do icartsec=1,ncontrac(Lrun) !loop  second index
                      do icartfirst=1,ncontrac(Lrun) !loop first index
                      work(iocaz+iredired+(icartfirst-1))=
     *                work(iocaz+iredired+(icartfirst-1))
     *               -onecartz(icartsec,icartfirst,mrun,Lrun)
                      enddo
                      iredired=iredired+itotalperIR(iredfirst)
                      enddo
        endif
      endif
      enddo
      enddo
      enddo
C
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
      iredirEd=shiftIRIR((ired2*ired2-ired2)/2+
     *                ired1)
      if (ired1.ne.ired2) then
        work(iocax2+irun-1)=work(iocax-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2))
        work(iocay2+irun-1)=work(iocay-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2))
        work(iocaz2+irun-1)=work(iocaz-1+iredired+norbsh2+
     * (norbsH1-1)*itotalperIR(IREd2))
      else
       work(iocax2+irun-1)=work(iocax-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)
       work(iocay2+irun-1)=work(iocay-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)
       work(iocaz2+irun-1)=work(iocaz-1+iredired+norbsh2*
     * (norbsH2-1)/2+norbsh1)
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
              write(LUPROP) (work(iocax2+irun),irun=0,length3-1)
              write(LUPROP)  Ya
              write(LUPROP) (work(iocay2+irun),irun=0,length3-1)
              write(LUPROP)  Za
              write(LUPROP) (work(iocaz2+irun),irun=0,length3-1)
CBS   close(luprop)
      else
cbs   reorder for AIMP
cbs   write(6,*) 'reorder integrals for AIMP'
      length3=ikeeporb*(ikeeporb+1)/2
      Call GetMem('OCAX3','Allo','Real',iocax3,3*length3)
      call dzero(work(iocax3),3*length3)
      iocay3=iocax3+length3
      iocaz3=iocay3+length3
cbs   write(6,*) 'number of orbitals ',ikeeporb
cbs   write(6,*) 'length of triangular matrix ', length3
      do irun2=1,ikeeporb
      do irun1=1,irun2
      ind2=ikeeplist(irun2)
      ind1=ikeeplist(irun1)
      ipntold=ipnt(ind1,ind2)-1
      ipntnew=ipnt(irun1,irun2)-1
      work(iocax3+ipntnew)=work(iocax2+ipntold)
      work(iocay3+ipntnew)=work(iocay2+ipntold)
      work(iocaz3+ipntnew)=work(iocaz2+ipntold)
      enddo
      enddo
CBS   write(6,*) 'transfered to new blocks'
CBS    Luprop=19
CBS           OPEN(LUPROP,STATUS='UNKNOWN',FORM='UNFORMATTED',
CBS  *        FILE='AOPROPER_MF')
CBS           rewind LUPROP
              write(LUPROP)  iCenter
              write(LUPROP)  xa,numbofsym,(nrtofiperIR(I),
     *        i=1,numbofsym),
     *        ikeeporb,(Loffunction(ikeeplist(i)),i=1,ikeeporb),
     *        (Moffunction(ikeeplist(i)),I=1,ikeeporb),
     *        Lhigh,((NContrac(I)-icore(I)),I=0,Lhigh)
              write(LUPROP) (work(iocax3+irun),irun=0,length3-1)
              write(LUPROP)  Ya
              write(LUPROP) (work(iocay3+irun),irun=0,length3-1)
              write(LUPROP)  Za
              write(LUPROP) (work(iocaz3+irun),irun=0,length3-1)
      Call GetMem('OCAX3','free','Real',iocax3,3*length3)
CBS   close(luprop)
      endif
cbs
cbs   that is it!!
cbs
      Call GetMem('OCAX','free','Real',iocax,6*iloca+MxContL*MxContL)
      return
c Avoid unused argument warnings
      if (.false.) call Unused_logical(makemean)
      end
