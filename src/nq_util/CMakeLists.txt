#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

set (sources
     allok2_funi.F90
     angular_grid.F90
     angular_prune.F90
     anmesh.F90
     aoadd_full.F90
     box_on_sphere.F90
     calcorboff.F90
     calcp2mocube.F90
     calc_pot1.F90
     calc_pot2.F90
     calc_pot2_inner.F90
     calcpuvxoff.F90
     comp_d.F90
     compute_d2mdx2.F90
     compute_d2odx2.F90
     compute_dmdx.F90
     compute_dodx.F90
     compute_grad.F90
     compute_m.F90
     compute_o.F90
     compute_rho.F90
     compute_tau.F90
     compute_t.F90
     converttabso.F90
     dede_funi.F90
     dft_grad.F90
     dft_int.F90
     do_batch.F90
     do_grid.F90
     do_index.F90
     do_nint_d.F90
     do_nintx.F90
     do_pi2.F90
     do_pi2grad.F90
     drvnq.F90
     drvnq_inner.F90
     dwdr.F90
     eval_rmax.F90
     eval_rmin.F90
     free_dede_funi.F90
     funi_init.F90
     funi_input.F90
     funi_print.F90
     gauleg.F90
     genradquad_b.F90
     genradquad_mhl.F90
     genradquad_mk.F90
     genradquad_pam.F90
     genradquad_ta.F90
     genvoronoi.F90
     get_subblock.F90
     g.F90
     ggl_grid.F90
     grid_on_disk.F90
     lebedev.F90
     lebedev_grid.F90
     libxc.F90
     libxc_interface.F90
     libxc_version.F90
     lobatto.F90
     lobatto_grid.F90
     mk_mos.F90
     mk_rho.F90
     mk_sos.F90
     modify_nq_grid.F90
     nbas_eff.F90
     nq_grid.F90
     nqgrid_init.F90
     nq_info.F90
     nq_init.F90
     nq_mo.F90
     nq_pdft.F90
     nq_structure.F90
     packpot1.F90
     pdftfock.F90
     pdftfock_inner.F90
     pdftmemalloc.F90
     pdftmemdealloc.F90
     phi_point.F90
     print_nq_info.F90
     process_coor.F90
     reset_nq_grid.F90
     resortd.F90
     rotgrd.F90
     setup_nq.F90
     subblock.F90
     symadp_full.F90
     transactmo2.F90
     transactmo.F90
     transfermo.F90
     translatedens.F90
     truncate_grid.F90
     unzipd1.F90
     unzipp2.F90
     w.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     libxc.F90
     nq_grid.F90
     nq_info.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
