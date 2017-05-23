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
      subroutine allocate_vplp_memory()
#include "drt_h.fh"
#include "pl_structure_h.fh"
      allocate(lp_coe(norb_dz+1:norb_inn+1,maxpl))
      allocate(lp_head(maxpl))
      allocate(lp_ltail(maxpl))
      allocate(lp_rtail(maxpl))
      allocate(lp_lwei(maxpl))
      allocate(lp_rwei(maxpl))
      allocate(vplp_w0(maxpl))
      allocate(vplp_w1(maxpl))
      allocate(logic_br(maxpl))
      allocate(lpnew_coe(norb_dz+1:norb_inn+1,maxpl))
      allocate(lpnew_head(maxpl))
      allocate(lpnew_ltail(maxpl))
      allocate(lpnew_rtail(maxpl))
      allocate(lpnew_lwei(maxpl))
      allocate(lpnew_rwei(maxpl))
      allocate(vplpnew_w0(maxpl))
      allocate(vplpnew_w1(maxpl))
      allocate(logic_newbr(maxpl))

      allocate(value_lpext(max_tmpvalue))
      allocate(value_lpext1(max_tmpvalue))
      allocate(value_lpext2(max_tmpvalue))
      allocate(index_lpext(max_tmpvalue))
      allocate(index_lpext1(max_tmpvalue))
      allocate(index_lpext2(max_tmpvalue))

      return
      end

      subroutine deallocate_vplp_memory()
#include "drt_h.fh"
#include "pl_structure_h.fh"

      deallocate(lp_coe)
      deallocate(lp_head)
      deallocate(lp_rtail)
      deallocate(lp_lwei)
      deallocate(lp_rwei)
      deallocate(vplp_w0)
      deallocate(vplp_w1)
      deallocate(logic_br)
      deallocate(lpnew_coe)
      deallocate(lpnew_head)
      deallocate(lpnew_ltail)
      deallocate(lpnew_rtail)
      deallocate(lpnew_lwei)
      deallocate(lpnew_rwei)
      deallocate(vplpnew_w0)
      deallocate(vplpnew_w1)
      deallocate(logic_newbr)

      deallocate(value_lpext)
      deallocate(value_lpext1)
      deallocate(value_lpext2)
      deallocate(index_lpext)
      deallocate(index_lpext1)
      deallocate(index_lpext2)

      return
      end

      subroutine change_vplp_pointer_arrays()
#include "drt_h.fh"
      REAL*8,pointer :: plptmp(:)
      integer, pointer :: iplptmp(:)
#include "pl_structure_h.fh"
      plptmp=>vplp_w0
      vplp_w0=>vplpnew_w0
      vplpnew_w0=>plptmp

      plptmp=>vplp_w1
      vplp_w1=>vplpnew_w1
      vplpnew_w1=>plptmp

      iplptmp=>lp_head
      lp_head=>lpnew_head
      lpnew_head=>iplptmp

      iplptmp=>lp_lwei
      lp_lwei=>lpnew_lwei
      lpnew_lwei=>iplptmp

      iplptmp=>lp_rwei
      lp_rwei=>lpnew_rwei
      lpnew_rwei=>iplptmp

      iplptmp=>lp_ltail
      lp_ltail=>lpnew_ltail
      lpnew_ltail=>iplptmp

      iplptmp=>lp_rtail
      lp_rtail=>lpnew_rtail
      lpnew_rtail=>iplptmp

      end

      subroutine change_coe_pointer_arrays()
#include "drt_h.fh"
      integer, dimension(:,:), pointer :: icoetmp
#include "pl_structure_h.fh"
      icoetmp=>lp_coe
      lp_coe=>lpnew_coe
      lpnew_coe=>icoetmp
      end

      subroutine change_br_pointer_arrays()
#include "drt_h.fh"
      logical, dimension(:), pointer :: logic_brtmp
#include "pl_structure_h.fh"
      logic_brtmp=>logic_br
      logic_br=>logic_newbr
      logic_newbr=>logic_brtmp
      end
