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
      subroutine makefile_cvb()

      call makeinit_cvb()
c  Memory allocation :
      call depend_cvb('MEM1','MEM0')
      call depend_cvb('MEM2','MEM1')
      call depend_cvb('MEM3','MEM2')
      call depend_cvb('MEM4','MEM3')
      call depend_cvb('MEM5','MEM4')
      call depend_cvb('MEM6','MEM5')
      call depend_cvb('MEM7','MEM6')
      call depend_cvb('MEMORY','MEM7')
c
      call decl_cvb('ORBS')
      call decl_cvb('CVB')
      call decl_cvb('ORBSTRY')
      call decl_cvb('CVBTRY')

      call depend_cvb('WFN','ORBS')
      call depend_cvb('WFN','CVB')
      call depend_cvb('WFNTRY','ORBSTRY')
      call depend_cvb('WFNTRY','CVBTRY')

      call depend_cvb('SVB','WFN')
      call depend_cvb('EVB','WFN')
      call depend_cvb('SVBTRY','WFNTRY')
      call depend_cvb('EVBTRY','WFNTRY')

      call depend_cvb('ORBFREE','MEM5')
      call depend_cvb('CIFREE','MEM5')
c
      call depend_cvb('ICONFS','MEM1')
c
      call depend_cvb('GENDET','ICONFS')
c
      call depend_cvb('SYMELM','MEM5')
      call decl_cvb('PRSYMELM')
c
      call depend_cvb('SYMINIT','MEM5')
      call depend_cvb('SYMINIT','SYMELM')
      call depend_cvb('ORBFREE','SYMINIT')
c
      call depend_cvb('CONSTRUC','MEM5')
      call depend_cvb('CONSTRUC','MEM3')
      call depend_cvb('CONSTRUC','SYMELM')
      call depend_cvb('CIFREE','CONSTRUC')
c
      call depend_cvb('RDINT','MEM0')
c
      call depend_cvb('RDCAS','MEM4')
c
      call depend_cvb('OOHESS','MEM6')
c
      call mkafter_cvb('SCALSTR','ICONFS')
c
c  Starting guess is chosen :
c  1) As input guess if applicable.
c  2) Otherwise from "start" fileid, if applicable.
c  3) Otherwise from previous opt. step if applicable.
c  4) Otherwise semi-random (orbs)/perfect-pairing (cvb).
      call touchdepend_cvb('GUESS','INPGS')
      call touchdepend_cvb('GUESS','STRTGS')
      call touchdepend_cvb('GUESS','RESTGS')
      call untouch_cvb('INPGS')
      call untouch_cvb('STRTGS')
      call untouch_cvb('RESTGS')
      call depend_cvb('ORBS','GUESS')
      call depend_cvb('CVB','GUESS')
      call depend_cvb('RESTGS','MEM2')
c
      call mkafter_cvb('GUESS','SCALSTR')
c
      call depend_cvb('TRWF','GUESS')
      call depend_cvb('SYMINIT','TRWF')
      call depend_cvb('CONSTRUC','TRWF')
c
      call depend_cvb('SYMORBS','TRWF')
      call depend_cvb('SYMORBS','SYMINIT')
      call depend_cvb('SYMORBS','ORBS')
c
      call depend_cvb('SYMCVB','TRWF')
      call depend_cvb('SYMCVB','SYMINIT')
c
      call depend_cvb('SYMCVB','CONSTRUC')
      call depend_cvb('SYMCVB','CVB')
c
      call depend_cvb('SYMWF','SYMORBS')
      call depend_cvb('SYMWF','SYMCVB')
c
      call touch_cvb('CASPRINT')
      call touch_cvb('CNFPRINT')
c
      call untouch_cvb('WRITEGS')
c
      call decl_cvb('CASCHECK')
c
      call mkafter_cvb('ICONFS','MEMORY')
      call mkafter_cvb('GENDET','ICONFS')
      call mkafter_cvb('SYMWF','GENDET')
      call mkafter_cvb('RDINT','SYMWF')
      call mkafter_cvb('RDCAS','RDINT')
c
      call depend_cvb('INIT','MEMORY')
      call depend_cvb('INIT','ICONFS')
      call depend_cvb('INIT','GENDET')
      call depend_cvb('INIT','SYMWF')
      call depend_cvb('INIT','RDINT')
      call depend_cvb('INIT','RDCAS')
c
      call untouch_cvb('ORBPERM')
      call depend_cvb('TRWF','ORBPERM')
      call untouch_cvb('TRNSPN')
      call depend_cvb('TRWF','TRNSPN')
      call mkafter_cvb('TRNSPN','GUESS')
c  Contents of CI vectors :
      call depend_cvb('CI-ORBS','ORBS')
      call depend_cvb('CI-CVB','CVB')
      call depend_cvb('CI-ALL','MEM4')

      call untouch_cvb('STAT')

      call decl_cvb('PRTSUM')

      return
      end
      subroutine rules_cvb(chr)
      character*(*) chr
      if(chr.eq.'MEM1')then
        call chop1_cvb()
      elseif(chr.eq.'MEM2')then
        call chop2_cvb()
      elseif(chr.eq.'MEM3')then
        call chop3_cvb()
      elseif(chr.eq.'MEM4')then
        call chop4_cvb()
      elseif(chr.eq.'MEM5')then
        call chop5_cvb()
      elseif(chr.eq.'MEM6')then
        call chop6_cvb()
      elseif(chr.eq.'MEM7')then
        call chop7_cvb()
      elseif(chr.eq.'ORBFREE')then
        call mkorbfree_cvb()
      elseif(chr.eq.'CIFREE')then
        call mkcifree_cvb()
      elseif(chr.eq.'ICONFS')then
        call mkiconfs_cvb()
      elseif(chr.eq.'GENDET')then
        call mkciinfo_cvb()
        call mkvbinfo_cvb()
      elseif(chr.eq.'SYMELM')then
        call mksymelm_cvb()
      elseif(chr.eq.'SYMINIT')then
        call mksyminit_cvb()
      elseif(chr.eq.'CONSTRUC')then
        call mkconstruc_cvb()
      elseif(chr.eq.'RDINT')then
        call mkrdint_cvb()
      elseif(chr.eq.'RDCAS')then
        call mkrdcas_cvb()
      elseif(chr.eq.'SYMORBS')then
        call mksymorbs_cvb()
      elseif(chr.eq.'SYMCVB')then
        call mksymcvb_cvb()
      elseif(chr.eq.'GUESS')then
        call mkguess_cvb()
      elseif(chr.eq.'ORBPERM')then
        call mkorbperm_cvb()
      elseif(chr.eq.'TRNSPN')then
        call mktrnspn_cvb()
      elseif(chr.eq.'STAT')then
        call stat_cvb()
      endif

      return
      end
      subroutine touchrules_cvb(chr)
      character*(*) chr
      if(chr.eq.'CI-ORBS')then
        call clearcnt_cvb(1)
      elseif(chr.eq.'CI-CVB')then
        call clearcnt_cvb(2)
      elseif(chr.eq.'CI-ALL')then
        call clearcnt_cvb(3)
      endif
      return
      end
