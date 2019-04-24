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
      subroutine tuneinp_cvb()
      implicit real*8 (a-h,o-z)
      parameter(nstrin=37,ncmp=8)
      character*8 string(nstrin)
      character*1 true(1)
#include "tune_cvb.fh"
#include "tols_cvb.fh"
#include "trst_cvb.fh"
#include "davtune_cvb.fh"
#include "optze_cvb.fh"
#include "print_cvb.fh"
#include "idbl_cvb.fh"
      dimension iaux(1),daux(1)
      save string,true
      data string/'CNRMTOL ','SAFETY  ','SIGNTOL','ALFTOL  ',
     >            'DFXTOL  ','EXP12TOL','GRDWRNGT','EIGWRNG ',
     >            'SINGUL  ','DFX     ','SIGN    ','ZZMAX   ',
     >            'ZZMIN   ','DX      ','GRD     ','NOPTH1  ',
     >            'NOPTH2  ','DELOPTH1','DELOPTH2','HHREJFAC',
     >            'HHACCFAC','ZZACCLIM','HHTOL   ','DFXMIN  ',
     >            'ZZREJMIN','ZZREJMAX','SCALESMA','HHSTART ',
     >            'RESTHR  ','NORTITER','ORTHTHR ','FOLLOW  ',
     >            'MXDAV   ','LASTUPD ','ENDIFCLO','ENDTUNE ',
     >            'END     '/
      data true/'T'/

1000  call fstring_cvb(string,nstrin,istr,ncmp,2)
      if(istr.eq.36.or.istr.eq.37)then
c 'ENDTUNE' or 'END'
        istr=0
      endif
      if(istr.eq.1)then
c  'CNRMTOL'
        call real_cvb(daux,1,nread,1)
        cnrmtol=daux(1)
      elseif(istr.eq.2)then
c  'SAFETY'
        call real_cvb(daux,1,nread,1)
        safety=daux(1)
      elseif(istr.eq.3)then
c  'SIGNTOL'
        call real_cvb(daux,1,nread,1)
        signtol=daux(1)
      elseif(istr.eq.4)then
c  'ALFTOL'
        call real_cvb(daux,1,nread,1)
        alftol=daux(1)
      elseif(istr.eq.5)then
c  'DFXTOL'
        call real_cvb(daux,1,nread,1)
        dfxtol=daux(1)
      elseif(istr.eq.6)then
c  'EXP12TOL'
        call real_cvb(daux,1,nread,1)
        exp12tol=daux(1)
      elseif(istr.eq.7)then
c  'GRDWRNGT'
        call real_cvb(daux,1,nread,1)
        grdwrngtol=daux(1)
      elseif(istr.eq.8)then
c  'EIGWRNG '
        call real_cvb(daux,1,nread,1)
        eigwrntol=daux(1)
      elseif(istr.eq.9)then
c  'SINGUL'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.3)then
          write(6,*)' Illegal I index in SINGUL :',i
          call abend_cvb()
        endif
        call real_cvb(singul(i),1,nread,1)
      elseif(istr.eq.10)then
c  'DFX'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.6)then
          write(6,*)' Illegal I index in DFX :',i
          call abend_cvb()
        endif
        call real_cvb(dfx(i),1,nread,1)
      elseif(istr.eq.11)then
c  'SIGN'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.6)then
          write(6,*)' Illegal I index in SIGN :',i
          call abend_cvb()
        endif
        call real_cvb(sign(i),1,nread,1)
      elseif(istr.eq.12)then
c  'ZZMAX'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.6)then
          write(6,*)' Illegal I index in ZZMAX :',i
          call abend_cvb()
        endif
        call real_cvb(zzmax(i),1,nread,1)
      elseif(istr.eq.13)then
c  'ZZMIN'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.6)then
          write(6,*)' Illegal I index in ZZMIN :',i
          call abend_cvb()
        endif
        call real_cvb(zzmin(i),1,nread,1)
      elseif(istr.eq.14)then
c  'DX'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.3)then
          write(6,*)' Illegal I index in DX :',i
          call abend_cvb()
        endif
        call int_cvb(iaux,1,nread,1)
        j=iaux(1)
        if(j.lt.1.or.j.gt.6)then
          write(6,*)' Illegal J index in DX :',j
          call abend_cvb()
        endif
        call real_cvb(dx(i,j),1,nread,1)
      elseif(istr.eq.15)then
c  'GRD'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.3)then
          write(6,*)' Illegal I index in GRD :',i
          call abend_cvb()
        endif
        call int_cvb(iaux,1,nread,1)
        j=iaux(1)
        if(j.lt.1.or.j.gt.6)then
          write(6,*)' Illegal J index in GRD :',j
          call abend_cvb()
        endif
        call real_cvb(grd(i,j),1,nread,1)
      elseif(istr.eq.16)then
c  'NOPTH1'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in NOPTH1 :',i
          call abend_cvb()
        endif
        call int_cvb(nopth1(i),1,nread,1)
      elseif(istr.eq.17)then
c  'NOPTH2'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in NOPTH2 :',i
          call abend_cvb()
        endif
        call int_cvb(nopth2(i),1,nread,1)
      elseif(istr.eq.18)then
c  'DELOPTH1'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in DELOPTH1 :',i
          call abend_cvb()
        endif
        call real_cvb(delopth1(i),1,nread,1)
      elseif(istr.eq.19)then
c  'DELOPTH2'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in DELOPTH2 :',i
          call abend_cvb()
        endif
        call real_cvb(delopth2(i),1,nread,1)
      elseif(istr.eq.20)then
c  'HHREJFAC'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in HHREJFAC :',i
          call abend_cvb()
        endif
        call real_cvb(hhrejfac(i),1,nread,1)
      elseif(istr.eq.21)then
c  'HHACCFAC'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.5)then
          write(6,*)' Illegal I index in HHACCFAC :',i
          call abend_cvb()
        endif
        call int_cvb(iaux,1,nread,1)
        j=iaux(1)
        if(j.lt.1.or.j.gt.2)then
          write(6,*)' Illegal J index in HHACCFAC :',j
          call abend_cvb()
        endif
        call real_cvb(hhaccfac(i,j),1,nread,1)
      elseif(istr.eq.22)then
c  'ZZACCLIM'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.4)then
          write(6,*)' Illegal I index in ZZACCLIM :',i
          call abend_cvb()
        endif
        call int_cvb(iaux,1,nread,1)
        j=iaux(1)
        if(j.lt.1.or.j.gt.2)then
          write(6,*)' Illegal J index in ZZACCLIM :',j
          call abend_cvb()
        endif
        call real_cvb(zzacclim(i,j),1,nread,1)
      elseif(istr.eq.23)then
c  'HHTOL'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in HHTOL :',i
          call abend_cvb()
        endif
        call real_cvb(hhtol(i),1,nread,1)
      elseif(istr.eq.24)then
c  'DFXMIN'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in DFXMIN :',i
          call abend_cvb()
        endif
        call real_cvb(dfxmin(i),1,nread,1)
      elseif(istr.eq.25)then
c  'ZZREJMIN'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in ZZREJMIN :',i
          call abend_cvb()
        endif
        call real_cvb(zzrejmin(i),1,nread,1)
      elseif(istr.eq.26)then
c  'ZZREJMAX'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in ZZREJMAX :',i
          call abend_cvb()
        endif
        call real_cvb(zzrejmax(i),1,nread,1)
      elseif(istr.eq.27)then
c  'SCALESMALL'
        call int_cvb(iaux,1,nread,1)
        i=iaux(1)
        if(i.lt.1.or.i.gt.2)then
          write(6,*)' Illegal I index in SCALESMALL :',i
          call abend_cvb()
        endif
        call fstring_cvb(true,1,istr2,1,1)
        scalesmall(i)=(istr2.eq.1)
      elseif(istr.eq.28)then
c  'HHSTART'
        call real_cvb(daux,1,nread,1)
        hhstart=daux(1)
      elseif(istr.eq.29)then
c  'RESTHR'
        call real_cvb(daux,1,nread,1)
        resthr=daux(1)
      elseif(istr.eq.30)then
c  'NORTITER'
        call int_cvb(iaux,1,nread,1)
        nortiter=iaux(1)
      elseif(istr.eq.31)then
c  'ORTHTHR'
        call real_cvb(daux,1,nread,1)
        orththr=daux(1)
      elseif(istr.eq.32)then
c  'FOLLOW'
        call fstring_cvb(true,1,istr2,1,1)
        follow=(istr2.eq.1)
      elseif(istr.eq.33)then
c  'MXDAV'
        call int_cvb(iaux,1,nread,1)
        mxdav=iaux(1)
      elseif(istr.eq.34)then
c  'LASTUPD'
        call fstring_cvb(true,1,istr2,1,1)
        lastupd=(istr2.eq.1)
      elseif(istr.eq.35)then
c  'ENDIFCLOSE'
        call fstring_cvb(true,1,istr2,1,1)
        endifclose=(istr2.eq.1)
      endif
c 'END', 'ENDTUNE' or unrecognized keyword -- end of input :
      if(istr.ne.0)goto 1000
      return
      entry tuneprint_cvb()
      if(ip(3).lt.3)return
      if(imethod.ne.4)then
      write(6,'(/,2a)')' -------- Details of parameters used by 2nd-',
     >  'order optimizer: -------------'
      write(6,'(/,a,/)')' General parameters:'
      call fout_cvb(safety,'SAFETY',
     >  'Alpha safety in denominator, (H - alpha * I):')
      call fout_cvb(cnrmtol,'CNRMTOL',
     >  'Tolerance for size of update:')
      call fout_cvb(signtol,'SIGNTOL',
     >  'Tolerance for sign of Hessian eigenvalues:')
      call fout_cvb(alftol,'ALFTOL',
     >  'Convergence criterion on alpha:')
      call fout_cvb(dfxtol,'DFXTOL',
     >  'DFX tolerance for act/exp ratio:')
      call fout_cvb(exp12tol,'EXP12TOL',
     >  'Criterion on expected change of f(x):')
      call fout_cvb(grdwrngtol,'GRDWRNGTOL',
     >  'Gradient tol. for scaling small updates:')
      call fout_cvb(eigwrngtol,'EIGWRNGTOL',
     >  'Eigenvalue tol. for scaling small updates:')
      call lout_cvb(lastupd,'LASTUPD',
     >  'Perform update at convergence?')
      call lout_cvb(endifclose,'ENDIFCLOSE',
     >  'Exit if optimization close to convergence?')
      write(6,'(/,a)')' Convergence criteria:'
      write(6,'(/,a)')' Elements of arrays:'
      write(6,'(a)')' (1) ... Optimization is in global region.'
      write(6,'(a)')' (2) ... Optimization is in local region.'
      write(6,'(2a)')' (3) ... Optimization is close to wrong ',
     >  'stationary point.'
      call fouti_cvb(singul,3,'SINGUL(3)',
     >  'Thresholds for sing. Hessian (max abs eig):')
      write(6,'(/,a)')' Elements of arrays:'
      write(6,'(a)')' (*,1) ... Global region, non-singular Hessian.'
      write(6,'(a)')' (*,2) ... Global region, singular Hessian.'
      write(6,'(a)')' (*,3) ... Local region, non-singular Hessian.'
      write(6,'(a)')' (*,4) ... Local region, singular Hessian.'
      write(6,'(2a)')' (*,5) ... Wrong stationary point, ',
     >  'non-singular Hessian.'
      write(6,'(2a)')' (*,6) ... Wrong stationary point, ',
     >  'singular Hessian.'
      call fouti_cvb(sign,6,'SIGN(6)',
     >  'Threshold for sign of Hessian eigenvalues:')
      call fouti_cvb(zzmin,6,'ZZMIN(6)',
     >  'Mininum allowed act/exp ratio:')
      call fouti_cvb(zzmax,6,'ZZMAX(6)',
     >  'Maximum allowed act/exp ratio:')
      call fouti_cvb(dfx,6,'DFX(6)',
     >  'Maximum allowed change in f(x):')
      write(6,'(/,a)')' Elements of arrays:'
      write(6,'(a)')' (1,*) ... Use maximum absolute value in vector.'
      write(6,'(a)')' (2,*) ... Use norm of vector.'
      write(6,'(a)')' (3,*) ... Use RMS of elements in vector.'
      call foutij_cvb(dx,3,6,'DX(3,6)',
     >  'Maximum allowed change in variables:')
      call foutij_cvb(grd,3,6,'GRD(3,6)',
     >  'Maximum allowed gradient:')
      write(6,'(/,a,/)')' Trust region control:'
      call fout_cvb(hhstart,'HHSTART',
     >  'Initial trust region size:')
      call iout_cvb(nopth1(1),'NOPTH1(1)',
     >  'Number of steps (primary trust size opt):')
      call iout_cvb(nopth2(1),'NOPTH2(1)',
     >  'Number of steps (secondary trust size opt):')
      call iout_cvb(nopth1(2),'NOPTH1(2)',
     >  'Number of steps (primary trust size opt):')
      call iout_cvb(nopth2(2),'NOPTH2(2)',
     >  'Number of steps (secondary trust size opt):')
      call fouti_cvb(delopth1,2,'DELOPTH1(2)',
     >  'Primary change of trust region size:')
      call fouti_cvb(delopth2,2,'DELOPTH2(2)',
     >  'Secondary change of trust region size:')
      call fouti_cvb(hhmax,2,'HHMAX(2)',
     >  'Maximum allowed trust region size:')
      call fouti_cvb(zzrejmin,2,'ZZREJMIN(2)',
     >  'Minimum allowed act/exp ratio:')
      call fouti_cvb(zzrejmax,2,'ZZREJMAX(2)',
     >  'Maximum allowed act/exp ratio:')
      call fouti_cvb(dfxmin,2,'DFXMIN(2)',
     >  'Minimum allowed change in f(x):')
      call fouti_cvb(hhrejfac,2,'HHREJFAC(2)',
     >  'Trust region size scale factor for rejections:')
      call foutij_cvb(zzacclim,4,2,'ZZACCLIM(4,2)',
     >  'Act/exp regions for scaling accepted steps:')
      call foutij_cvb(hhaccfac,5,2,'HHACCFAC(5,2)',
     >  'Trust scale factors for accepted steps:')
      call fouti_cvb(hhtol,2,'HHTOL(2)',
     >  'Minimum allowed trust region size:')
      call lout_cvb(scalesmall(1),'SCALESMALL(1)',
     >  'Scale predicted steps smaller than trust?')
      call lout_cvb(scalesmall(2),'SCALESMALL(2)',
     >  'Scale predicted steps smaller than trust?')
      write(6,'(/,2a)')' -------------------------------------------',
     >  '------------------------------'
      elseif(imethod.eq.4)then
      write(6,'(/,2a,/)')' -------- Details of parameters used by Davi',
     >  'dson optimizer: --------------'
      call iout_cvb(mxdav,'MXDAV',
     >  'Maxium dimension of Davidson subspace:')
      call fout_cvb(resthr,'RESTHR',
     >  'Convergence criterion on residual norm:')
      call lout_cvb(follow,'FOLLOW',
     >  'Root following (for excited states):')
      call fout_cvb(orththr,'ORTHTHR',
     >  'Tolerance for orthogonality between vectors:')
      call iout_cvb(nortiter,'NORTITER',
     >  'Maximum number of orthogonalization attempts:')
      write(6,'(/,2a)')' -------------------------------------------',
     >  '------------------------------'
      endif
      return
      end
