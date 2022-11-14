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
     cct3_add10.F90
     cct3_add21.F90
     cct3_add22.F90
     cct3_add32.F90
     cct3_add41.F90
     cct3_add42.F90
     cct3_add43.F90
     cct3_add44.F90
     cct3_add.F90
     cct3_divfok.F90
     cct3_expand0.F90
     cct3_expand1.F90
     cct3_expand2.F90
     cct3_expand3.F90
     cct3_expand40.F90
     cct3_expand41.F90
     cct3_expand.F90
     cct3.F90
     cct3_fokunpck1.F90
     cct3_fokunpck2.F90
     cct3_fokunpck3.F90
     cct3_fokunpck4.F90
     cct3_fokunpck5.F90
     cct3_getint.F90
     cct3_getmap.F90
     cct3_getmediate.F90
     cct3_global.F90
     cct3_grc0.F90
     cct3_grc22c.F90
     cct3_grc23c.F90
     cct3_grc24c.F90
     cct3_grc32c.F90
     cct3_grc32y.F90
     cct3_grc33c.F90
     cct3_grc34c.F90
     cct3_grc42c.F90
     cct3_grc42y.F90
     cct3_grc43c.F90
     cct3_grc43y.F90
     cct3_grc44c.F90
     cct3_map11.F90
     cct3_map21.F90
     cct3_map22.F90
     cct3_map31.F90
     cct3_map32.F90
     cct3_map41.F90
     cct3_map42.F90
     cct3_map.F90
     cct3_mc0c1a3b.F90
     cct3_mc0c1at3b.F90
     cct3_multc0.F90
     cct3_mult.F90
     cct3_multy0.F90
     cct3_mv0v1a3u.F90
     cct3_noperm.F90
     cct3_rea.F90
     cct3_t3grc0.F90
     defv.F90
     defvhlp1.F90
     defvhlp21.F90
     defvhlp22.F90
     defvhlp3.F90
     defvhlp4.F90
     defvhlp51.F90
     defvhlp52.F90
     defvhlp53.F90
     defvhlp54.F90
     defvhlp61.F90
     defvhlp62.F90
     defvhlp7.F90
     defvhlp81.F90
     defvhlp82.F90
     defvhlp9.F90
     ext.F90
     exth1.F90
     exth2.F90
     exth3.F90
     exth4.F90
     exth5.F90
     getintpos.F90
     minusa.F90
     setb.F90
     stz.F90
     t3addpck.F90
     t3aphlp1.F90
     t3aphlp2.F90
     t3aphlp3.F90
     t3aphlp4.F90
     t3aphlp5.F90
     t3aphlp6.F90
     t3aphlp7.F90
     t3aphlp8.F90
     t3dhlp1.F90
     t3dhlp2.F90
     t3dhlp3.F90
     t3dhlp4.F90
     t3div.F90
     t3initfiles.F90
     t3reaccsd.F90
     t3reainput.F90
     t3reaintsta.F90
     t3sgl.F90
     t3sglh11.F90
     t3sglh121.F90
     t3sglh122.F90
     t3sglh131.F90
     t3sglh132.F90
     t3sglh141.F90
     t3sglh142.F90
     t3sglh143.F90
     t3sglh211.F90
     t3sglh212.F90
     t3sglh221.F90
     t3sglh222.F90
     t3sglh223.F90
     t3sglh311.F90
     t3sglh312.F90
     t3sglh321.F90
     t3sglh322.F90
     t3sglh323.F90
     t3wresult.F90
     t3wrhead.F90
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     cct3_global.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
