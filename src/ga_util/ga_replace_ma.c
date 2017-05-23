/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
#ifdef _GA_
#  include <stdlib.h>
#  include "molcastype.h"
#  include "mma.h"
#  include "ga.h"

#  ifdef _CAPITALS_
#    define ga_replace_ma GA_REPLACE_MA
#  else
#    ifndef ADD_
#      define ga_replace_ma ga_replace_ma_
#    endif
#  endif

void* replace_malloc (size_t bytes, int align, char* name) {
        return allomblck("GA", &bytes);
}

void replace_free (void* ptr) {
        freemblck(ptr);
}

void ga_replace_ma (void) {
        GA_Register_stack_memory(replace_malloc, replace_free);
}
#endif
