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
     ampint.F90
     ampmem.F90
     ampr.F90
     assemble_dtdmu.F90
     assemble_dvdb.F90
     assemble_mgauss.F90
     assmbl.F90
     ass_px.F90
     ass_pxp.F90
     cassmbl.F90
     ccmbnmp.F90
     ccmbnve.F90
     ccrtcmp.F90
     cmbnke.F90
     cmbnke_giao.F90
     cmbnker.F90
     cmbnmp_giao.F90
     cmbnmv.F90
     cmbnve.F90
     cntint.F90
     cntmem.F90
     contact.F90
     cvelint.F90
     d1int.F90
     d1mem.F90
     darwin.F90
     dmsint.F90
     dmsmem.F90
     dtdmu_int.F90
     dtdmu_mem.F90
     dumint.F90
     dummem.F90
     efint.F90
     efmem.F90
     emfint.F90
     emfmem.F90
     epeint.F90
     epemem.F90
     kneint.F90
     kneint_giao.F90
     knemem.F90
     knemem_giao.F90
     knemmp.F90
     kntc.F90
     kntc_giao.F90
     m1int.F90
     m1mem.F90
     m2int.F90
     m2mem.F90
     mltint.F90
     mltint_giao.F90
     mltmem.F90
     mltmem_giao.F90
     mltmmp.F90
     mve.F90
     mveint.F90
     mvemem.F90
     naint.F90
     naint_giao.F90
     namem.F90
     namem_giao.F90
     oamint.F90
     oammem.F90
     omqint.F90
     omqmem.F90
     one2h5_crtmom.F90
     one2h5_fckint.F90
     one2h5_ovlmat.F90
     pam2.F90
     pam2int.F90
     pam2mem.F90
     p_int.F90
     p_mem.F90
     potint.F90
     pot_nuc.F90
     prjint.F90
     prjmem.F90
     pvint.F90
     pvmem.F90
     pxint.F90
     pxmem.F90
     pxpint.F90
     pxpmem.F90
     qpvint.F90
     qpvmem.F90
     sroint.F90
     sromem.F90
     util1.F90
     util2.F90
     util3.F90
     util4.F90
     util5.F90
     util8.F90
     veint.F90
     velint.F90
     vemem.F90
     vpint.F90
     vpmem.F90
     welint.F90
     welmem.F90
     xfdint.F90
     xfdmem.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list "")

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
