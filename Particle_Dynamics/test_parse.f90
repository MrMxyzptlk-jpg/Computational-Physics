program test_parse
  use FoX_dom
  implicit none

  type(Node), pointer :: myDoc
  myDoc => parseFile("dummy.xml")
end program