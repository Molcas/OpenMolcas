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
c  **************************************************************
c  ** Transformation between VB structures and VB determinants **
c  ** Also generation of spin-function to determinant transf.  **
c  ** and of perfect-pairing guess.                            **
c  **************************************************************
c  *********************************************************************
c  *                                                                   *
c  *  Transformation between VB structures (VB CSFs) and determinants. *
c  *                                                                   *
c  *  The VB structures are ordered according to spatial               *
c  *  configurations, and are normally defined in terms of spin        *
c  *  eigenfunctions.                                                  *
c  *  VB determinants are ordered with alpha & beta indices in         *
c  *  increasing order, with alpha being the slower index.             *
c  *                                                                   *
c  *  STR2VB[CFG] : Structure-to-determinant transformation.           *
c  *  VB2STR[CFG] : Determinant-to-structure transformation.           *
c  *                                                                   *
c  *  [C] : Transforms coefficients.                                   *
c  *  [F] : Transforms a first-order change of the coefficients.       *
c  *  [G] : Transforms a gradient-type quantity.                       *
c  *                                                                   *
c  *  At present there is no actual difference between [C] and [F].    *
c  *                                                                   *
c  *  The difference between [CF] and [G] is mainly important when     *
c  *  non-orthogonal spin functions are used (Rumer/Project.)          *
c  *                                                                   *
c **********************************************************************
      subroutine str2vbc_cvb(cvb,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)

      kab=2

      kbs=nint(w(lb(kab)))
      if(kbs.ne.kbasiscvb)then
        call mkbiks_cvb()
        kbs=kbasiscvb
      endif
      ioffs_cvb=1
      ioffs_cvbdet=1
      idetvb_add=ll(17)
      ifnss_add=lb(4)
      if(kbasiscvb.eq.6)ifnss_add=lb(5)
      ndetvbs_add=lb(6)
      do 200 ifrag=1,nfrag
      iwrk=mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
      call str2vb2_cvb(w(lb(kab)+1),iw(lb(3)),cvb(ioffs_cvb),
     >  cvbdet(ioffs_cvbdet),2,
     >  iw(idetvb_add),
     >  i2s_fr(1,ifrag),nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),
     >  iw(ifnss_add),iw(ndetvbs_add),
     >  absym(1),
     >  mnion_fr(ifrag),mxion_fr(ifrag),nconf_fr(ifrag),
     >  ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,
     >  nel_fr(ifrag),nalf_fr(1,ifrag),nel,
     >  w(iwrk),nconfion_fr(0,ifrag))
      call mfreer_cvb(iwrk)
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
200   ioffs_cvbdet=ioffs_cvbdet+ndetvb_fr(ifrag)
      return

      end
c
      subroutine vb2strg_cvb(cvbdet,cvb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)
      kab=2
      kbs=nint(w(lb(kab)))
      if(kbs.ne.kbasiscvb)then
        call mkbiks_cvb()
        kbs=kbasiscvb
      endif
      ioffs_cvb=1
      ioffs_cvbdet=1
      idetvb_add=ll(17)
      ifnss_add=lb(4)
      if(kbasiscvb.eq.6)ifnss_add=lb(5)
      ndetvbs_add=lb(6)
      do 400 ifrag=1,nfrag
      iwrk=mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
      call str2vb2_cvb(w(lb(kab)+1),iw(lb(3)),cvb(ioffs_cvb),
     >  cvbdet(ioffs_cvbdet),1,
     >  iw(idetvb_add),
     >  i2s_fr(1,ifrag),nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),
     >  iw(ifnss_add),iw(ndetvbs_add),
     >  absym(1),
     >  mnion_fr(ifrag),mxion_fr(ifrag),nconf_fr(ifrag),
     >  ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,
     >  nel_fr(ifrag),nalf_fr(1,ifrag),nel,
     >  w(iwrk),nconfion_fr(0,ifrag))
      call mfreer_cvb(iwrk)
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
400   ioffs_cvbdet=ioffs_cvbdet+ndetvb_fr(ifrag)
      return
      end


      subroutine vb2strc_cvb(cvbdet,cvb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)

      kab=1

      kbs=nint(w(lb(kab)))
      if(kbs.ne.kbasiscvb)then
        call mkbiks_cvb()
        kbs=kbasiscvb
      endif
      ioffs_cvb=1
      ioffs_cvbdet=1
      idetvb_add=ll(17)
      ifnss_add=lb(4)
      if(kbasiscvb.eq.6)ifnss_add=lb(5)
      ndetvbs_add=lb(6)
      do 400 ifrag=1,nfrag
      iwrk=mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
      call str2vb2_cvb(w(lb(kab)+1),iw(lb(3)),cvb(ioffs_cvb),
     >  cvbdet(ioffs_cvbdet),1,
     >  iw(idetvb_add),
     >  i2s_fr(1,ifrag),nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),
     >  iw(ifnss_add),iw(ndetvbs_add),
     >  absym(1),
     >  mnion_fr(ifrag),mxion_fr(ifrag),nconf_fr(ifrag),
     >  ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,
     >  nel_fr(ifrag),nalf_fr(1,ifrag),nel,
     >  w(iwrk),nconfion_fr(0,ifrag))
      call mfreer_cvb(iwrk)
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
400   ioffs_cvbdet=ioffs_cvbdet+ndetvb_fr(ifrag)
      return
      end

      subroutine str2vbg_cvb(cvb,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)

      kab=1

      kbs=nint(w(lb(kab)))
      if(kbs.ne.kbasiscvb)then
        call mkbiks_cvb()
        kbs=kbasiscvb
      endif
      ioffs_cvb=1
      ioffs_cvbdet=1
      idetvb_add=ll(17)
      ifnss_add=lb(4)
      if(kbasiscvb.eq.6)ifnss_add=lb(5)
      ndetvbs_add=lb(6)
      do 200 ifrag=1,nfrag
      iwrk=mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
      call str2vb2_cvb(w(lb(kab)+1),iw(lb(3)),cvb(ioffs_cvb),
     >  cvbdet(ioffs_cvbdet),2,
     >  iw(idetvb_add),
     >  i2s_fr(1,ifrag),nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),
     >  iw(ifnss_add),iw(ndetvbs_add),
     >  absym(1),
     >  mnion_fr(ifrag),mxion_fr(ifrag),nconf_fr(ifrag),
     >  ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,
     >  nel_fr(ifrag),nalf_fr(1,ifrag),nel,
     >  w(iwrk),nconfion_fr(0,ifrag))
      call mfreer_cvb(iwrk)
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
200   ioffs_cvbdet=ioffs_cvbdet+ndetvb_fr(ifrag)
      return

      end

      subroutine str2vbf_cvb(cvb,cvbdet)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


#include "frag_cvb.fh"
#include "malloc_cvb.fh"
      dimension cvb(nvb),cvbdet(ndetvb)

      kab=2

      kbs=nint(w(lb(kab)))
      if(kbs.ne.kbasiscvb)then
        call mkbiks_cvb()
        kbs=kbasiscvb
      endif
      ioffs_cvb=1
      ioffs_cvbdet=1
      idetvb_add=ll(17)
      ifnss_add=lb(4)
      if(kbasiscvb.eq.6)ifnss_add=lb(5)
      ndetvbs_add=lb(6)
      do 200 ifrag=1,nfrag
      iwrk=mstackr_cvb(max(ndetvb_fr(ifrag),nvb_fr(ifrag)))
      call str2vb2_cvb(w(lb(kab)+1),iw(lb(3)),cvb(ioffs_cvb),
     >  cvbdet(ioffs_cvbdet),2,
     >  iw(idetvb_add),
     >  i2s_fr(1,ifrag),nS_fr(ifrag),nalf_fr(1,ifrag),nMs_fr(ifrag),
     >  iw(ifnss_add),iw(ndetvbs_add),
     >  absym(1),
     >  mnion_fr(ifrag),mxion_fr(ifrag),nconf_fr(ifrag),
     >  ndetvb_fr(ifrag),nvb_fr(ifrag),kbs,
     >  nel_fr(ifrag),nalf_fr(1,ifrag),nel,
     >  w(iwrk),nconfion_fr(0,ifrag))
      call mfreer_cvb(iwrk)
      idetvb_add=idetvb_add+ndetvb_fr(ifrag)
      ioffs_cvb=ioffs_cvb+nvb_fr(ifrag)
200   ioffs_cvbdet=ioffs_cvbdet+ndetvb_fr(ifrag)
      return

      end
