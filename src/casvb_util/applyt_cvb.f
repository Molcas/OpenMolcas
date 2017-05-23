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
c  ************************************
c  ** Routines involving CI and ORBS **
c  ************************************
c  *********************************************************************
c  *                                                                   *
c  *  APPLYT := transform CI vector IVEC according to orbital          *
c  *         := transformation in IGJORB (from GAUSSJ)                 *
c  *                                                                   *
c  *  :> IVEC:      starting CI vector                                 *
c  *  :> IGJORB:    simple orbital updates                             *
c  *  :> N_APPLYT   stats                                              *
c  *  :> I1ALF:     string handling                                    *
c  *  :> I1BET:           do.                                          *
c  *  :> IATO:            do.                                          *
c  *  :> IBTO:            do.                                          *
c  *  :> PHATO:           do.                                          *
c  *  :> PHBTO:           do.                                          *
c  *                                                                   *
c  *  :< IVEC:      transformed CI vector                              *
c  *  :< N_APPLYT   stats                                              *
c  *                                                                   *
c  *********************************************************************
      subroutine applyt_cvb(cvec,igjorb)
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension igjorb(*),cvec(*)

      ivec=nint(cvec(1))
      n_applyt=n_applyt+1
      ioff=idbl_cvb(norb*norb)
      if(iform_ci(ivec).eq.0)then
        call permci_cvb(cvec,igjorb(1+ioff))
        call applyt2_cvb(w(iaddr_ci(ivec)),
     >    igjorb,igjorb(1+norb+ioff),
     >    iw(ll(1)),iw(ll(2)),iw(ll(5)),iw(ll(6)),w(ll(9)),w(ll(10)))
      else
        write(6,*)' Unsupported format in APPLYT :',iform_ci(ivec)
        call abend_cvb()
      endif

      call setcnt2_cvb(ivec,0)
      return
      end
