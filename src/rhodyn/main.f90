program main
  use rhodyn_data
  implicit none
  character(8), parameter :: module_name = 'rhodyn'

  call start(module_name)
  call rhodyn()
  call finish(ireturn)
end
