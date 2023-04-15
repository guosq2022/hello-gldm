module hello_gldm
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, hello-gldm!"
    print *, " test github working!!!"
    print *, " test github working2!!!"
    print *, " test github working3!!!"
  end subroutine say_hello
end module hello_gldm
