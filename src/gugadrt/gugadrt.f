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
      subroutine gugadrt(ireturn)
#include "gendrt.fh"
      data istep_occ/0,1,1,2/
      data mul_tab/1,2,3,4,5,6,7,8,
     *             2,1,4,3,6,5,8,7,
     *             3,4,1,2,7,8,5,6,
     *             4,3,2,1,8,7,6,5,
     *             5,6,7,8,1,2,3,4,
     *             6,5,8,7,2,1,4,3,
     *             7,8,5,6,3,4,1,2,
     *             8,7,6,5,4,3,2,1/
      sc0=seconds()

      call gugainit()

      call gugadrt_mole_inf()
      call gugadrt_paras_calculate()
      call arrange_orbital_molcas()

      call gugadrt_dbl_upwalk()       ! add by wyb 01.9.5
      call gugadrt_ext_downwalk()       ! add by wyb 01.9.5
      call gugadrt_active_drt()        ! add by wyb 01.9.5

      call add_info("CI_DIM",[dble(nci_dim)],1,1)
      call gugadrt_gugafinalize()
      ireturn=0

      sc1=seconds()
      sc=sc1-sc0
      write(6,910) sc
      return
910   format(/,1x,"end of job, takes ",f10.2,1x,"seconds"/)
      end
