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
      Subroutine prwf_vibrot(ndim,R)
      Implicit Real*8 (A-H,O-Z)
C Print the value of (unnormalized) vibrational wave functions.
CPAM97 For MOLCAS-4
C
#include "dimensions.fh"
      Parameter(nemax=nVib_Max,ndimmx=5000,lwork=nemax*ndimmx)
#include "intinp.fh"
      Dimension Vib(lwork)
      Dimension R(*)
*
      Call qEnter('PrWf_VibRot')

      write(6,*)
      Call CollapseOutput(1,'PRINTOUT OF VIBRATIONAL WAVE FUNCTIONS')
      write(6,*)

      ne1=nvib1+1
      ndim1=ndim+1
      nwork=ne1*ndim1
      ierr=0
      If(nwork.gt.lwork) then
        write(6,*)' PRWF_VibRot: Local work space too small.'
        write(6,*)' Please increase LWORK and recompile.'
        write(6,*)' Need at least:',nwork
        ierr=ierr+1
      End If
      If(ne1.gt.nemax) then
        write(6,*)' PRWF_VibRot: Local work space too small.'
        write(6,*)' Please increase NEMAX and recompile.'
        write(6,*)' Need at least:',ne1
        ierr=ierr+1
      End If
      If(ndim1.gt.ndimmx) then
        write(6,*)' PRWF_VibRot: Local work space too small.'
        write(6,*)' Please increase NDIMMX and recompile.'
        write(6,*)' Need at least:',ndim1
        ierr=ierr+1
      End If
      If(ierr.gt.0) Call Abend

C Loop over rotational quantum numbers
      Do J1=J1A,J2A
       Jad1=J1-J1A+1
       Write(6,'(1x,a,i3)')' Rotational quantum number J=',J1

C Read vibrational functions for this J-value.
       ist=1
       iadr1=iad12(Jad1)
       Do nv=1,ne1
        Call DDafile(Vibwvs,2,Vib(ist),ndim1,iadr1)
        ist=ist+ndim1
       End Do

C Write out the wave functions:
       do nvsta=0,ne1-1,5
         nvend=min(nvsta+4,ne1-1)
         write(6,*)
         write(6,'(5x,a12,8x,''v='',5(i2,13x))')
     &         'Radial dist.',(nv,nv=nvsta,nvend)
         do i=1,ndim
           write(6,'(1x,f12.6,5x,5f15.8)') r(i),
     &         (vib(i+ndim1*nv),nv=nvsta,nvend)
         end do
       end do

C End of loop over J1.
      end do
      Call CollapseOutput(0,'PRINTOUT OF VIBRATIONAL WAVE FUNCTIONS')
      write(6,*)
      Call qExit('PrWf_VibRot')
      return
      end
