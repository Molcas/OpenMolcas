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
      subroutine stat_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"


      if(ip(3).ge.1)then
        write(6,'(/,a,i16)')
     >    ' Total number of structure transformations :',n_applyt
        write(6,'(a,i17)')
     >    ' Total number of Hamiltonian applications :',n_applyh
        write(6,'(a,i11)')
     >    ' Total number of 2-electron density evaluations :',n_2el
        write(6,'(a,i21)')
     >    ' Total number of Hessian applications :',n_hess
        if(n_orbhess.gt.0)write(6,'(a,i8)')
     >    ' Total number of pure orbital Hessian applications :',
     >    n_orbhess
        if(n_cihess.gt.0)write(6,'(a,i13)')
     >    ' Total number of pure CI Hessian applications :',n_cihess
        write(6,'(a,i18,/)')
     >    ' Approximate memory usage (8-byte words) :',ibasemx
     >      -ibase0
        write(6,'(a,f10.3,a)')' CASVB at ',tim_cvb(cpu0),' CPU seconds'
        memused=0
        call stat1_cvb()
      endif
      return
      end
c  Format statement settings :
