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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine str2vb2_cvb(bikcof,ikcoff,cvb,cvbdet,iway,
     >  idetvb,
     >  i2s,nS,nalf1,nMs,
     >  ifnss,ndetvbs,
     >  absym,mnion,mxion,nconf,ndetvb,nvb,kbasis,nel,nalf,neltot,
     >  work,nconfion)
      implicit real*8 (a-h,o-z)
#include "malloc_cvb.fh"
      logical absym
      dimension bikcof(*),ikcoff(0:neltot,0:neltot,0:neltot)
      dimension cvb(nvb),cvbdet(ndetvb)
      dimension idetvb(ndetvb)
      dimension i2S(nS),nalf1(nMs)
      dimension ifnss(0:neltot,0:neltot),ndetvbs(0:neltot,0:neltot)
      dimension work(ndetvb),nconfion(0:*)
      save one,sqp5,sq2
      data one/1d0/,sqp5/.70710678118654752440d0/
      data sq2/1.41421356237309504880d0/

c
      i2s_keep      = 0 ! dummy initialize
      nalfsing_keep = 0 ! dummy initialize
c  Determinant to structure transformation
      if(iway.eq.1)then
        call fzero(cvb,nvb)
        do 100 idet=1,ndetvb
        work(idet)=cvbdet(idetvb(idet))
100     continue
      elseif(iway.eq.2)then
        call fzero(work,ndetvb)
      endif
      idadd=0
      isadd=0
      iconfadd=0
      do 200 ion=0,nel/2
      if(nconfion(ion).eq.0)goto 201
      nelsing=nel-2*ion
c  Investigate different S and Ms possibilities and
c  prepare to collect different BIKCOF matrices if
c  necessary :
      n_spin=0
      n_spin_values=0
      do 300 iS=1,nS
      if(i2s(iS).le.nelsing)then
        n_spin=n_spin+ifnss(nelsing,i2s(iS))
        n_spin_values=n_spin_values+1
        i2s_keep=i2s(iS)
      endif
300   continue
      n_det=0
      n_det_values=0
      do 400 iMs=1,nMs
      nalfsing=nalf1(iMs)-ion
      if(nalfsing.ge.0)then
        n_det=n_det+ndetvbs(nelsing,nalfsing)
        n_det_values=n_det_values+1
        nalfsing_keep=nalfsing
      endif
400   continue
      if(kbasis.eq.6)then
        do 500 iS=1,nS
        nalfsing_det=(nelsing+i2s(iS))/2
        if(i2s(iS).le.nelsing)then
          do 600 iMs=1,nMs
          nalfsing=nalf1(iMs)-ion
          if(nalfsing.ge.0)then
            if(nalfsing.ne.nalfsing_det)goto 600
            if(iway.eq.1)then
              do 700 idet=1,ifnss(nelsing,i2s(iS))
              if(i2s(iS).eq.0.and.absym.and.
     >          ndetvbs(nelsing,nalfsing).ne.1)then
                call daxpy_(nconfion(ion),sq2,work(idet+idadd),n_det,
     >            cvb(idet+isadd),n_spin)
              else
                call daxpy_(nconfion(ion),one,work(idet+idadd),n_det,
     >            cvb(idet+isadd),n_spin)
              endif
700           continue
            elseif(iway.eq.2)then
              do 800 idet=1,ifnss(nelsing,i2s(iS))
              if(i2s(iS).eq.0.and.absym.and.
     >          ndetvbs(nelsing,nalfsing).ne.1)then
                call daxpy_(nconfion(ion),sqp5,cvb(idet+isadd),n_spin,
     >            work(idet+idadd),n_det)
                call daxpy_(nconfion(ion),sqp5,cvb(idet+isadd),n_spin,
     >            work(ndetvbs(nelsing,nalfsing)-idet+1+idadd),n_det)
              else
                call daxpy_(nconfion(ion),one,cvb(idet+isadd),n_spin,
     >            work(idet+idadd),n_det)
              endif
800           continue
            endif
          endif
600       continue
        endif
500     continue
c  Skip collection if not necessary ...
      elseif(n_spin_values.eq.1.and.n_det_values.eq.1)then
        if(iway.eq.1)then
          call mxattbp_cvb(
     >      bikcof(1+ikcoff(nelsing,nalfsing_keep,i2s_keep)),
     >      work(1+idadd),n_spin,n_det,nconfion(ion),cvb(1+isadd))
        elseif(iway.eq.2)then
          call mxatbp_cvb(
     >      bikcof(1+ikcoff(nelsing,nalfsing_keep,i2s_keep)),
     >      cvb(1+isadd),n_det,n_spin,nconfion(ion),work(1+idadd))
        endif
      else
        i1 = mstackrz_cvb(n_det*n_spin)
        i_spin=0
        do 900 iS=1,nS
        if(i2s(iS).le.nelsing)then
          i_det=0
          do 1000 iMs=1,nMs
          nalfsing=nalf1(iMs)-ion
          if(nalfsing.ge.0)then
            if(ikcoff(nelsing,nalfsing,i2s(iS)).ne.-1)then
              ioff_bikcof=1+ikcoff(nelsing,nalfsing,i2s(iS))
              ioff_i1=i1+i_spin*n_det+i_det
              do 1100 j_spin=1,ifnss(nelsing,i2s(iS))
              call fmove_cvb(bikcof(ioff_bikcof),w(ioff_i1),
     >          ndetvbs(nelsing,nalfsing))
              ioff_bikcof=ioff_bikcof+ndetvbs(nelsing,nalfsing)
              ioff_i1=ioff_i1+n_det
1100          continue
            endif
            i_det=i_det+ndetvbs(nelsing,nalfsing)
          endif
1000    continue
        i_spin=i_spin+ifnss(nelsing,i2s(iS))
        endif
900     continue

        if(iway.eq.1)then
          call mxattbp_cvb(w(i1),work(1+idadd),
     >     n_spin,n_det,nconfion(ion),cvb(1+isadd))
        elseif(iway.eq.2)then
          call mxatbp_cvb(w(i1),cvb(1+isadd),
     >      n_det,n_spin,nconfion(ion),work(1+idadd))
        endif
        call mfreer_cvb(i1)
      endif
      isadd=isadd+nconfion(ion)*n_spin
      idadd=idadd+nconfion(ion)*n_det
201   iconfadd=iconfadd+nconfion(ion)
200   continue
      if(iway.eq.2)then
        do 1200 idet=1,ndetvb
        cvbdet(idetvb(idet))=work(idet)
1200    continue
      endif
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(mnion)
        call Unused_integer(mxion)
        call Unused_integer(nconf)
        call Unused_integer(nalf)
      end if
      end
