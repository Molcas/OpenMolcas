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
     aixcheck.F90
     aixcls.F90
     aixerr.c
     aixfsz.F90
     aixmv.F90
     aixopn.F90
     aixprd.F90
     aixpwr.F90
     aixrd.F90
     aixrm.F90
     aixwr.F90
     allocdisk.c
     append_file.F90
     append_file_u.F90
     bdafile.F90
     cdafile.F90
     chdafile.F90
     cio.c
     daclos.F90
     daeras.F90
     dafile_checkarg.F90
     dafile.F90
     daname.F90
     daname_main.F90
     daname_mf.F90
     daname_mf_wa.F90
     daname_wa.F90
     ddafile.F90
     f_inquire.F90
     fast_io.F90
     fastio.F90
     fcopy.F90
     filesystem_wrapper.c
     find_lu.F90
     fioinit.F90
     fscb2unit.F90
     get_mbl_wa.F90
     handle2name.F90
     i1dafile.F90
     idafile.F90
     isfreeunit.F90
     is_opened.F90
     lu2desc.F90
     lu2handle.F90
     molcas_binaryopen_vanilla.F90
     molcas_open_ext2.F90
     molcas_open.F90
     mpdafile.F90
     prgm.F90
     prgminit.F90
     prgmtranslatec.F90
     prgmtranslate.F90
     prgmtranslate_master.F90
     text_file.F90
     zip.c
)

# Source files defining modules that should be available to other *_util directories
set (modfile_list
     prgm.F90
     text_file.F90
)

include (${PROJECT_SOURCE_DIR}/cmake/util_template.cmake)
