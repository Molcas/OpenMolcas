      Integer function norder(icoord,int_code,lmax)

      Implicit None
      Integer lmax
      Integer icoord(lmax),int_code(lmax)
      Integer nb,isite

      nb=0
      Do isite=1,lmax
         nb=nb+icoord(isite)*int_code(isite)
      End Do

      norder=nb+1
      Return
      End
