program main
  use hello_gldm, only: say_hello, gldm
  implicit none
  real(8) :: a0, z0, a1, z1, a2, z2
  real(8) :: qexp, t12

  call say_hello()

  a0 = 238.
  z0 = 92.
  a2 = 4.
  z2 = 2.
  a1 = a0 - a2
  z1 = z0 - z2
  qexp = 4.00

  call gldm(a1, z2, a2, z2, qexp, t12)

end program main
