module firstSteps
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, firstSteps!"
  end subroutine say_hello
end module firstSteps
