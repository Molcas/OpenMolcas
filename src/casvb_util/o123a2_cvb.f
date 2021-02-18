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
      subroutine o123a2_cvb(nparm,grad,eigvec,eigval,gradp)
      implicit real*8 (a-h,o-z)
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

      dimension grad(nparm)
      dimension eigvec(nparm,nparm),eigval(nparm)
      dimension gradp(nparm)

      call gethess_cvb(eigvec)
      call mxdiag_cvb(eigvec,eigval,nparm)
      call mxatb_cvb(grad,eigvec,1,nparm,nparm,gradp)
      if(ip.ge.2)then
        write(6,'(a)')' Gradient in basis of Hessian eigenvectors :'
        call vecprint_cvb(gradp,nparm)
      endif
      return
      end
