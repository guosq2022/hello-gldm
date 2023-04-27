program main
  use hello_gldm, only: say_hello, gldm
  implicit none
  real(8) :: a0, z0, a1, z1, a2, z2
  real(8) :: qexp, t12

  call say_hello()

  a0 = 212.
  z0 = 86. !86. !84.
  a2 = 4.
  z2 = 2.
  a1 = a0 - a2
  z1 = z0 - z2
  qexp = 6.3847 !6.3847 !8.9542

  call gldm(a1, z1, a2, z2, qexp, t12)

end program main

