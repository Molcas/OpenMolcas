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
      subroutine casinfo1_rvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "casinfo_cvb.fh"
      logical iphex,oldex

c  Information from molcas interface file 'JOBIPH' :
      write(6,'(2a)')' ------- Recover RASSCF-related information',
     >              ' --------------------------------------'
      call f_inquire('JOBIPH',iphex)
      call f_inquire('JOBOLD',oldex)
      if(iphex)then
        write(6,'(/,a)')' Using JOBIPH interface file.'
c        call systemf("cp -p JOBIPH JOBOLD")
       Call Copy_JobIph("JOBIPH","JOBOLD")
      elseif(oldex)then
        write(6,'(/,a)')' Using JOBOLD interface file.'
c        call systemf("cp -p JOBOLD JOBIPH")
       Call Copy_JobIph("JOBOLD","JOBIPH")
      else
        write(6,'(/,a)')' Error: need either JOBOLD or JOBIPH file!'
        call abend_cvb()
      endif

      call rdjobiph_cvb('JOBIPH')
      call setjobiph_cvb(iorcore_c,iorclos_c,iorocc_c,mxirrep,
     >  nstsym_c,weight_c,istnel_c,istsy_c,istms2_c,nstats_c,
     >  mxstt_ci,mxstsy_ci,nel_c,norb_c,i2s_c,isym_c,mcore_c,neltot_c)
      call rasscf(ireturn_rasscf)
      call clsfls_rasscf()
c  rasscf will have overwritten jobiph ...
c      call systemf("cp -p JOBOLD JOBIPH")
       Call Copy_JobIph("JOBOLD","JOBIPH")
      write(6,'(2a)')' ------- RASSCF-related information recovered',
     >              ' ------------------------------------'
      return
      end
