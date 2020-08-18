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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine dpci2vb2_cvb(civec,cvbdet,dvbdet,evbdet,ic1,ret,ic,
     >  nda,ndb,ndetvb,
     >  nfrag,nda_fr,ndb_fr,
     >  ia12ind,ib12ind,nc_facalf,nc_facbet,ncindalf,ncindbet,
     >  istack,mxstack,coeff,idetind,
     >  iapr,ixapr,
     >  ipr_off,ixapr_off,ixbpr_off,
     >  ndetvb_fr,ndavb)
      implicit real*8 (a-h,o-z)
      dimension civec(nda,ndb)
      dimension cvbdet(ndetvb),dvbdet(ndetvb),evbdet(ndetvb)
      dimension istack(mxstack)
      dimension nda_fr(nfrag),ndb_fr(nfrag)
      dimension ncindalf(0:nfrag),ncindbet(0:nfrag)
      dimension nc_facalf(nfrag),nc_facbet(nfrag)
      dimension ia12ind(*),ib12ind(*)
      dimension coeff(0:nfrag),idetind(nfrag)
      dimension iapr(ndetvb),ixapr(ndavb)
      dimension ipr_off(nfrag),ixapr_off(nfrag),ixbpr_off(nfrag)
      dimension ndetvb_fr(nfrag)

      if(ic.eq.0)then
        call fzero(cvbdet,ndetvb)
      elseif(ic.eq.1.or.ic.eq.4)then
        call fzero(civec,nda*ndb)
      elseif(ic.eq.3)then
        ret=0d0
      elseif(ic.eq.5)then
        call fzero(evbdet,ndetvb)
      endif

      do 100 ifr=1,nfrag
      if(ifr.eq.1)then
        ipr_off(ifr)=0
        ixapr_off(ifr)=0
        ixbpr_off(ifr)=0
      else
        ipr_off(ifr)=ipr_off(ifr-1)+ndetvb_fr(ifr-1)
        ixapr_off(ifr)=ixapr_off(ifr-1)+nda_fr(ifr-1)+1
        ixbpr_off(ifr)=ixbpr_off(ifr-1)+ndb_fr(ifr-1)+1
      endif
100   continue

      cinrm=0d0
      coeff(0)=1d0
      ncindalf(0)=1
      ncindbet(0)=1
      do i=1,nfrag
      if(i.eq.1)then
        nc_facalf(i)=1
        nc_facbet(i)=1
      else
        nc_facalf(i)=nc_facalf(i-1)*nda_fr(i-1)
        nc_facbet(i)=nc_facbet(i-1)*ndb_fr(i-1)
      endif
      enddo

      nloop=nfrag
c  MXITERS -> NDA_FR

c  Following is code for a set of nested loops. To deal with the
c  complication that the number of nested loops is not known at
c  compile time, a simple integer stack is used.
c  NESTLEVEL=1 signifies we are doing outermost loop and so on.

      nestlevel=0
      call istkinit_cvb(istack,mxstack)

c  Here we go to the beginning of the next loop in the sequence :
1000  continue
      if(nestlevel.lt.nloop)then
        nestlevel=nestlevel+1
        iter=0
        mxiter=ndetvb_fr(nestlevel)
        call istkpush_cvb(istack,iter)
        call istkpush_cvb(istack,mxiter)
      endif

c  Here we do the next loop iteration of the current loop :
2000  if(nestlevel.eq.0)goto 3000
      call istkpop_cvb(istack,mxiter)
      call istkpop_cvb(istack,iter)
      iter=iter+1
      if(iter.gt.mxiter)then
        nestlevel=nestlevel-1
        goto 2000
      else
        call istkpush_cvb(istack,iter)
        call istkpush_cvb(istack,mxiter)
      endif

c  Here goes the code specific to this loop level.
      idetvb=0
      ixa = 0 ! dummy initialize
      do 2100 ia=1,nda_fr(nestlevel)
      do 2101 ixa=ixapr(ia+ixapr_off(nestlevel)),
     >  ixapr(ia+1+ixapr_off(nestlevel))-1
      idetvb=idetvb+1
      if(idetvb.eq.iter)goto 2200
2101  continue
2100  continue
      write(6,*)' Error in DPCI2VB '
      call abend_cvb()
2200  ib=iapr(ixa+ipr_off(nestlevel))

      idetind(nestlevel)=iter+ipr_off(nestlevel)
      if((ic.eq.1.and.ic1.eq.0).or.ic.eq.3)then
        coeff(nestlevel)=coeff(nestlevel-1)
     >    *cvbdet(iter+ipr_off(nestlevel))
      elseif(ic.ne.4)then
        coeff(nestlevel)=dvbdet(iter+ipr_off(nestlevel))
      endif

      ncindalf(nestlevel)=ncindalf(nestlevel-1)
     >  +nc_facalf(nestlevel)*(ia-1)
      ncindbet(nestlevel)=ncindbet(nestlevel-1)
     >  +nc_facbet(nestlevel)*(ib-1)

      if(nestlevel.eq.nfrag)then
        ia_ci=ia12ind(ncindalf(nestlevel))
        if(ia_ci.ne.0)then
          ib_ci=ib12ind(ncindbet(nestlevel))
          if(ib_ci.ne.0)then
            if(ic.eq.0)then
              if(ic1.eq.0)then
c  --  CI2VBC  --
                cinrm=cinrm+civec(abs(ia_ci),abs(ib_ci))
     >            *civec(abs(ia_ci),abs(ib_ci))
                if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                  do ifr=1,nfrag
                  cvbdet(idetind(ifr))=cvbdet(idetind(ifr))
     >              +civec(abs(ia_ci),abs(ib_ci))
                  enddo
                else
                  do ifr=1,nfrag
                  cvbdet(idetind(ifr))=cvbdet(idetind(ifr))
     >              -civec(abs(ia_ci),abs(ib_ci))
                  enddo
                endif
              elseif(ic1.eq.2)then
c  --  CI2VBG  --
                do ifr=1,nfrag
                cf=1d0
                do jfr=1,nfrag
                if(jfr.ne.ifr)cf=cf*coeff(jfr)
                enddo
                if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                  cvbdet(idetind(ifr))=cvbdet(idetind(ifr))
     >              +cf*civec(abs(ia_ci),abs(ib_ci))
                else
                  cvbdet(idetind(ifr))=cvbdet(idetind(ifr))
     >              -cf*civec(abs(ia_ci),abs(ib_ci))
                endif
                enddo
              endif
            elseif(ic.eq.1)then
              if(ic1.eq.0)then
c  --  VB2CIC  --
                if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                  civec(abs(ia_ci),abs(ib_ci))=coeff(nestlevel)
                else
                  civec(abs(ia_ci),abs(ib_ci))=-coeff(nestlevel)
                endif
              elseif(ic1.eq.1)then
c  --  VB2CIF  --
                do ifr=1,nfrag
                cf=1d0
                do jfr=1,nfrag
                if(jfr.ne.ifr)cf=cf*coeff(jfr)
                enddo
                if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                  civec(abs(ia_ci),abs(ib_ci))=
     >              civec(abs(ia_ci),abs(ib_ci))
     >              +cf*cvbdet(idetind(ifr))
                else
                  civec(abs(ia_ci),abs(ib_ci))=
     >              civec(abs(ia_ci),abs(ib_ci))
     >              -cf*cvbdet(idetind(ifr))
                endif
                enddo
              endif
            elseif(ic.eq.2)then
c  --  VB2CIAF  --
              do ifr=1,nfrag
              cf=1d0
              do jfr=1,nfrag
              if(jfr.ne.ifr)cf=cf*coeff(jfr)
              enddo
              if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                civec(abs(ia_ci),abs(ib_ci))=
     >            civec(abs(ia_ci),abs(ib_ci))
     >            +cf*cvbdet(idetind(ifr))
              else
                civec(abs(ia_ci),abs(ib_ci))=
     >            civec(abs(ia_ci),abs(ib_ci))
     >            -cf*cvbdet(idetind(ifr))
              endif
              enddo
            elseif(ic.eq.3)then
c  --   VB2CIDOT  --
              if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                ret=ret+civec(abs(ia_ci),abs(ib_ci))*coeff(nestlevel)
              else
                ret=ret-civec(abs(ia_ci),abs(ib_ci))*coeff(nestlevel)
              endif
            elseif(ic.eq.4)then
c  --  SETIAPR  --
              civec(abs(ia_ci),abs(ib_ci))=1d0
            elseif(ic.eq.5)then
c  --  CI2ORDR  --
              do ifr=1,nfrag
              cf=0d0
              do jfr=1,nfrag
              if(jfr.ne.ifr)then
                cf2=cvbdet(idetind(jfr))
                do kfr=1,nfrag
                if(kfr.ne.ifr.and.kfr.ne.jfr)
     >            cf2=cf2*dvbdet(idetind(kfr))
                enddo
                cf=cf+cf2
              endif
              enddo
              if((ia_ci.gt.0).eqv.(ib_ci.gt.0))then
                evbdet(idetind(ifr))=evbdet(idetind(ifr))
     >            +cf*civec(abs(ia_ci),abs(ib_ci))
              else
                evbdet(idetind(ifr))=evbdet(idetind(ifr))
     >            -cf*civec(abs(ia_ci),abs(ib_ci))
              endif
              enddo
            endif
          endif
        endif
      endif

      goto 1000

c  This is the end ...
3000  continue

      if(ic.eq.0.and.ic1.eq.0)then
c  "Normalize" the coefficients for each fragment :
        fac=1d0/sqrt(cinrm**(1d0/DBLE(nfrag)))
        ndetvb_add=1
        do ifr=1,nfrag
        cnrm=ddot_(ndetvb_fr(ifr),cvbdet(ndetvb_add),1,
     >    cvbdet(ndetvb_add),1)
        fac1=fac*cnrm
        call dscal_(ndetvb_fr(ifr),fac1*sqrt(cnrm),
     >    cvbdet(ndetvb_add),1)
        ndetvb_add=ndetvb_add+ndetvb_fr(ifr)
        enddo
      endif
      return
      end
