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
      subroutine contone(L,oneoverR3,onecontr,Lmax,
     *contcoeff,nprim,ncont,MxcontL,dummy,
     *onecartx,onecartY,onecartZ,charge,oneonly)
cbs   contracts one-electron integrals and multiplies with l,m-dependent
cbs   factors for L-,L0,L+
      implicit real*8 (a-h,o-z)
      dimension oneoverR3(*),
     *onecontr(MxcontL,MxcontL,-Lmax:Lmax,3),
     *contcoeff(nprim,ncont),dummy(ncont,ncont),
     *onecartx(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),
     *onecarty(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1)),
     *onecartz(MxcontL,MxcontL,
     *(Lmax+Lmax+1)*(Lmax+1))
      logical oneonly
      ipnt(I,J)=(max(i,j)*(max(i,j)-1))/2+min(i,j)
cbs   first of all cleaning dummy and onecontr
      do jrun=1,ncont
      do irun=1,ncont
      dummy(irun,jrun)=0d0
      enddo
      enddo
      if (oneonly) then
      iprod=MxcontL*MxcontL*(Lmax+Lmax+1)*(Lmax+1)
      call dzero(onecartx,iprod)
      call dzero(onecarty,iprod)
      call dzero(onecartz,iprod)
      endif
      iprod=3*(Lmax+lmax+1)*MxcontL*MxcontL
      call dzero(onecontr,iprod)
cbs   contract onto dummy
      do icont2=1,ncont
      do icont1=1,ncont
      do iprim2=1,nprim
      do iprim1=1,nprim
      dummy(icont1,icont2)=dummy(icont1,icont2)+
     *contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*
     *oneoverR3(ipnt(iprim1,iprim2))
      enddo
      enddo
      enddo
      enddo
      do icont2=1,ncont
      do icont1=1,ncont
      dummy(icont1,icont2)=dummy(icont1,icont2)*charge
      enddo
      enddo
cbs   start to add l,m dependent factors
      do M=-L,L
      factormin=sqrt(DBLE(L*L-M*M+L+M))
      factor0=DBLE(M)
      factorplus=sqrt(DBLE(L*L-M*M+L-M))
      do irun=1,ncont
      do jrun=1,ncont
      onecontr(irun,jrun,M,1)=dummy(jrun,irun)*factormin  ! L-minus
      enddo
      enddo
      do irun=1,ncont
      do jrun=1,ncont
      onecontr(irun,jrun,M,2)=dummy(jrun,irun)*factor0    ! L-0
      enddo
      enddo
      do irun=1,ncont
      do jrun=1,ncont
      onecontr(irun,jrun,M,3)=dummy(jrun,irun)*factorplus ! L-plus
      enddo
      enddo
      enddo
cbs   make the final cartesian integrals
      call cartoneX(L,Lmax,onecontr,ncont,
     *MxcontL,onecartX(1,1,1))
      call cartoneY(L,Lmax,onecontr,ncont,
     *MxcontL,onecartY(1,1,1))
      call cartoneZ(L,Lmax,onecontr,ncont,
     *MxcontL,onecartZ(1,1,1))
      return
      end
