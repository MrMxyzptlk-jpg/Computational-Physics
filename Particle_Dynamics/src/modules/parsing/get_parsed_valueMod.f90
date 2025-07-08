!############ Subroutines to inquire about and get the value of each kind of data ###############
MODULE get_parsed_valueMod
    use precisionMod
    use FoX_dom
    implicit none

    public get_parsed_value

    interface get_parsed_value
        module procedure get_parsed_string
        module procedure get_parsed_logical
        module procedure get_parsed_integerScalar
        module procedure get_parsed_integerVector
        module procedure get_parsed_realScalar
        module procedure get_parsed_realVector
    end interface
CONTAINS

subroutine get_parsed_string(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    character(len=*)                :: value

    if (hasAttribute(inputNode, name)) then
        value = getAttribute(inputNode, name)
    end if

end subroutine get_parsed_string

subroutine get_parsed_logical(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    logical                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_logical

subroutine get_parsed_integerScalar(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    integer                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_integerScalar

subroutine get_parsed_integerVector(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    integer                         :: value(3)
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
    end if

end subroutine get_parsed_integerVector

subroutine get_parsed_realScalar(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    real(pr)                        :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_realScalar

subroutine get_parsed_realVector(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    real(pr)                        :: value(3)
    character(len=88)               :: attr_string      ! Long enough for the total format_state specifier

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
    end if

end subroutine get_parsed_realVector

END MODULE get_parsed_valueMod