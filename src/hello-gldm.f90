module hello_gldm
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, hello-gldm!"
  end subroutine say_hello
end module hello_gldm
