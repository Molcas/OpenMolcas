!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine allocate_vplp_memory()

use gugaci_global, only: index_lpext, index_lpext1, index_lpext2, logic_br, logic_newbr, lp_coe, lp_head, lp_ltail, lp_lwei, &
                         lp_rtail, lp_rwei, lpnew_coe, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, max_tmpvalue, &
                         maxpl, norb_dz, norb_inn, value_lpext, value_lpext1, value_lpext2, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use stdalloc, only: mma_allocate

implicit none

call mma_allocate(lp_coe,[norb_dz+1,norb_inn+1],[1,maxpl],label='lp_coe')
call mma_allocate(lp_head,maxpl,label='lp_head')
call mma_allocate(lp_ltail,maxpl,label='lp_ltail')
call mma_allocate(lp_rtail,maxpl,label='lp_rtail')
call mma_allocate(lp_lwei,maxpl,label='lp_lwei')
call mma_allocate(lp_rwei,maxpl,label='lp_rwei')
call mma_allocate(vplp_w0,maxpl,label='vplp_w0')
call mma_allocate(vplp_w1,maxpl,label='vplp_w1')
call mma_allocate(logic_br,maxpl,label='logic_br')
call mma_allocate(lpnew_coe,[norb_dz+1,norb_inn+1],[1,maxpl],label='lpnew_coe')
call mma_allocate(lpnew_head,maxpl,label='lpnew_head')
call mma_allocate(lpnew_ltail,maxpl,label='lpnew_ltail')
call mma_allocate(lpnew_rtail,maxpl,label='lpnew_rtail')
call mma_allocate(lpnew_lwei,maxpl,label='lpnew_lwei')
call mma_allocate(lpnew_rwei,maxpl,label='lpnew_rwei')
call mma_allocate(vplpnew_w0,maxpl,label='vplpnew_w0')
call mma_allocate(vplpnew_w1,maxpl,label='vplpnew_w1')
call mma_allocate(logic_newbr,maxpl,label='logic_newbr')

call mma_allocate(value_lpext,max_tmpvalue,label='value_lpext')
call mma_allocate(value_lpext1,max_tmpvalue,label='value_lpext1')
call mma_allocate(value_lpext2,max_tmpvalue,label='value_lpext2')
call mma_allocate(index_lpext,max_tmpvalue,label='index_lpext')
call mma_allocate(index_lpext1,max_tmpvalue,label='index_lpext1')
call mma_allocate(index_lpext2,max_tmpvalue,label='index_lpext2')

return

end subroutine allocate_vplp_memory

subroutine deallocate_vplp_memory()

use gugaci_global, only: index_lpext, index_lpext1, index_lpext2, logic_br, logic_newbr, lp_coe, lp_head, lp_ltail, lp_lwei, &
                         lp_rtail, lp_rwei, lpnew_coe, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, lpnew_rwei, value_lpext, &
                         value_lpext1, value_lpext2, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use stdalloc, only: mma_deallocate

implicit none

call mma_deallocate(lp_coe)
call mma_deallocate(lp_head)
call mma_deallocate(lp_ltail)
call mma_deallocate(lp_rtail)
call mma_deallocate(lp_lwei)
call mma_deallocate(lp_rwei)
call mma_deallocate(vplp_w0)
call mma_deallocate(vplp_w1)
call mma_deallocate(logic_br)
call mma_deallocate(lpnew_coe)
call mma_deallocate(lpnew_head)
call mma_deallocate(lpnew_ltail)
call mma_deallocate(lpnew_rtail)
call mma_deallocate(lpnew_lwei)
call mma_deallocate(lpnew_rwei)
call mma_deallocate(vplpnew_w0)
call mma_deallocate(vplpnew_w1)
call mma_deallocate(logic_newbr)

call mma_deallocate(value_lpext)
call mma_deallocate(value_lpext1)
call mma_deallocate(value_lpext2)
call mma_deallocate(index_lpext)
call mma_deallocate(index_lpext1)
call mma_deallocate(index_lpext2)

return

end subroutine deallocate_vplp_memory

subroutine change_vplp_pointer_arrays()

use gugaci_global, only: lp_head, lp_ltail, lp_lwei, lp_rtail, lp_rwei, lpnew_head, lpnew_ltail, lpnew_lwei, lpnew_rtail, &
                         lpnew_rwei, vplp_w0, vplp_w1, vplpnew_w0, vplpnew_w1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), allocatable :: iplptmp(:)
real(kind=wp), allocatable :: plptmp(:)

call move_alloc(vplp_w0,plptmp)
call move_alloc(vplpnew_w0,vplp_w0)
call move_alloc(plptmp,vplpnew_w0)

call move_alloc(vplp_w1,plptmp)
call move_alloc(vplpnew_w1,vplp_w1)
call move_alloc(plptmp,vplpnew_w1)

call move_alloc(lp_head,iplptmp)
call move_alloc(lpnew_head,lp_head)
call move_alloc(iplptmp,lpnew_head)

call move_alloc(lp_lwei,iplptmp)
call move_alloc(lpnew_lwei,lp_lwei)
call move_alloc(iplptmp,lpnew_lwei)

call move_alloc(lp_rwei,iplptmp)
call move_alloc(lpnew_rwei,lp_rwei)
call move_alloc(iplptmp,lpnew_rwei)

call move_alloc(lp_ltail,iplptmp)
call move_alloc(lpnew_ltail,lp_ltail)
call move_alloc(iplptmp,lpnew_ltail)

call move_alloc(lp_rtail,iplptmp)
call move_alloc(lpnew_rtail,lp_rtail)
call move_alloc(iplptmp,lpnew_rtail)

end subroutine change_vplp_pointer_arrays

subroutine change_coe_pointer_arrays()

use gugaci_global, only: lp_coe, lpnew_coe
use Definitions, only: iwp

implicit none
integer(kind=iwp), allocatable :: icoetmp(:,:)

call move_alloc(lp_coe,icoetmp)
call move_alloc(lpnew_coe,lp_coe)
call move_alloc(icoetmp,lpnew_coe)

end subroutine change_coe_pointer_arrays

subroutine change_br_pointer_arrays()

use gugaci_global, only: logic_br, logic_newbr
use Definitions, only: iwp

implicit none
logical(kind=iwp), allocatable :: logic_brtmp(:)

call move_alloc(logic_br,logic_brtmp)
call move_alloc(logic_newbr,logic_br)
call move_alloc(logic_brtmp,logic_newbr)

end subroutine change_br_pointer_arrays
