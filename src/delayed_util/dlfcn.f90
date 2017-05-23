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
MODULE ISO_C_UTILITIES
   USE ISO_C_BINDING ! Intrinsic module

   CHARACTER(KIND=C_CHAR), DIMENSION(1), SAVE, TARGET, PRIVATE :: dummy_string="?"

CONTAINS

   FUNCTION C_F_STRING(CPTR) RESULT(FPTR)
      ! Convert a null-terminated C string into a Fortran character array pointer
      TYPE(C_PTR), INTENT(IN) :: CPTR ! The C address
      CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR

      INTERFACE ! strlen is a standard C function from <string.h>
         ! int strlen(char *string)
         FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
            USE ISO_C_BINDING
            TYPE(C_PTR), VALUE :: string ! A C pointer
            INTEGER(KIND=C_INT) :: len
         END FUNCTION
      END INTERFACE

      IF (C_ASSOCIATED(CPTR)) THEN
         CALL C_F_POINTER(FPTR=FPTR, CPTR=CPTR, SHAPE=[strlen(CPTR)])
      ELSE
         ! To avoid segfaults, associate FPTR with a dummy target:
         FPTR=>dummy_string
      END IF

   END FUNCTION

END MODULE

MODULE DLFCN
   USE ISO_C_BINDING
   USE ISO_C_UTILITIES
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DLOpen, DLSym, DLClose, DLError, DLAddr ! DL API

   ! Valid modes for mode in DLOpen:
   INTEGER, PARAMETER, PUBLIC :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=256, RTLD_LOCAL=0
      ! Obtained from the output of the previously listed C program

   ! Struct for DLAddr
   TYPE, BIND(C), PUBLIC :: DL_info
      TYPE(C_PTR) :: dli_fname
      TYPE(C_FUNPTR) :: dli_fbase
      TYPE(C_PTR) :: dli_sname
      TYPE(C_FUNPTR) :: dli_saddr
   END TYPE DL_info

   INTERFACE ! All we need is interfaces for the prototypes in <dlfcn.h>
      FUNCTION DLOpen(file,mode) RESULT(handle) BIND(C,NAME="dlopen")
         ! void *dlopen(const char *file, int mode);
         USE ISO_C_BINDING
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: file
            ! C strings should be declared as character arrays
         INTEGER(C_INT), VALUE :: mode
         TYPE(C_PTR) :: handle
      END FUNCTION
      FUNCTION DLSym(handle,name) RESULT(funptr) BIND(C,NAME="dlsym")
         ! void *dlsym(void *handle, const char *name);
         USE ISO_C_BINDING
         TYPE(C_PTR), VALUE :: handle
         CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: name
         TYPE(C_FUNPTR) :: funptr ! A function pointer
      END FUNCTION
      FUNCTION DLClose(handle) RESULT(status) BIND(C,NAME="dlclose")
         ! int dlclose(void *handle);
         USE ISO_C_BINDING
         TYPE(C_PTR), VALUE :: handle
         INTEGER(C_INT) :: status
      END FUNCTION
      FUNCTION DLError() RESULT(error) BIND(C,NAME="dlerror")
         ! char *dlerror(void);
         USE ISO_C_BINDING
         TYPE(C_PTR) :: error
      END FUNCTION
      ! dladdr is a Glibc extension, not POSIX
      FUNCTION DLAddr(funptr,info) RESULT(output) BIND(C,NAME="dladdr")
         ! int dladdr(void *addr, Dl_info *info)
         USE ISO_C_BINDING
         TYPE(C_FUNPTR), VALUE :: funptr ! A function pointer
         TYPE(C_PTR), VALUE :: info
         INTEGER(C_INT) :: output
      END FUNCTION
   END INTERFACE

END MODULE
