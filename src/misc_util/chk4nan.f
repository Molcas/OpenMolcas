************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per Ake Malmqvist                                      *
*               Jesper Wisborg Krogh                                   *
************************************************************************
*  Chk4NAN
*
*> @brief
*>   Check whether an array contains NANs
*> @author Per &Aring;ke Malmqvist
*> @modified_by Jesper Wisborg Krogh
*>
*> @details
*> The routine checks the elements of the input array for NANs. If any are found,
*> up to 100 of these will be listed. Upon return \p iErr contains the number of
*> NANs found.
*>
*> @param[in]  nDim  Total dimension of array
*> @param[in]  Array Array to be checked
*> @param[out] Ierr  Return code
************************************************************************
      Subroutine Chk4NAN(nDim, Array, Ierr)
      Implicit None
      Real*8 Array, CheckSum
      Character*16 str16
      Integer nDim, iCount, I, Ierr
      Dimension Array(nDim)
*
      ICOUNT=0
      CHECKSUM=0.0D0
      DO I=1,nDim
         CHECKSUM=CHECKSUM+ARRAY(I)
      END DO
      WRITE(STR16,'(G16.7)') CHECKSUM
      CALL NORMAL(STR16)
      IF(STR16(1:1).eq.'N') THEN
         Write(6,*) '!!! WARNING !!!'
         Write(6,*) 'NANs encountered'
         Write(6,*)
         Write(6,*)' The numbers in the array will now be checked.'
         Write(6,*)' There are ',NDIM,' elements.'
         DO I=1,nDim
            WRITE(STR16,'(G16.7)') ARRAY(I)
            CALL NORMAL(STR16)
            IF(STR16(1:1).eq.'N') THEN
               ICOUNT=ICOUNT+1
               IF(iCount .le. 100) THEN
                  WRITE(6,*)' Element nr.', I,' is ',ARRAY(I)
               END IF
            END IF
         END DO
         IF(ICOUNT.GT.100) THEN
            WRITE(6,*)' ...too many. I give up here.'
         END IF
         Write(6,*) 'There were a total of ',iCount,' NANs'
      END IF
*
      iErr = iCount
      Return
      End
