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
     atnear.F90
     atsymb.F90
     cav_hss.F90
     cavitation.F90
     cmbn_ef_dpnt.F90
     colchg.F90
     colour.F90
     coltss.F90
     datasol.F90
     dcdc.F90
     dcdr.F90
     derbet.F90
     derd.F90
     deriva.F90
     der_norm.F90
     derphi.F90
     dmat_cpcm.F90
     drdc.F90
     drdr.F90
     drv_ef_pcm.F90
     dsd.F90
     dver.F90
     fndbnd.F90
     fndsph.F90
     fndtess.F90
     gaubon.F90
     geoder.F90
     grdtr_alaska.F90
     gvwrite.F90
     hybnew.F90
     icolat.F90
     inprct.F90
     inter_pcm.F90
     ipbo.F90
     irowat.F90
     matpcm.F90
     mlt_pcm.F90
     nalpar.F90
     ncalph.F90
     numsolv.F90
     oneel_g_pcm.F90
     over.F90
     pcm_arrays.F90
     pcm_cav_grd.F90
     pcmdef.F90
     pcm_driver.F90
     pcm_ef_grd.F90
     pcmgrd1.F90
     pcmhss.F90
     pcm_init.F90
     pcmmmh.F90
     piat.F90
     polygen.F90
     potgrd.F90
     prgrad_pcm.F90
     prtcav.F90
     solvent_data.F90
     tessera.F90
     uatm.F90
     v_ef_pcm.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     pcm_arrays.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
