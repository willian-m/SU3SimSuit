!Willian Matioli Serenone
!Institution: Universidade de Sao Paulo
!             Instituto de Fisica de Sao Carlos
!e-mail: willian.serenone@usp.br
!######################################################################

!==============================
!MODULE: xml_parser
!Uses the FoX library to parse an input file and load simulation parameters 
!FoX library is an open-source library, licensed under the BSD license.
!You need to download and compile it before compiling this module.
!Download it at https://homepages.see.leeds.ac.uk/~earawa/FoX/
!Unpack the content under src/lib/Fox
!Execute configure script: ./configure
!Compile the library: make
!==============================
module xml_parser
    use FoX_sax
    use types_params
    implicit none
    logical :: in_nx,in_ny,in_nz,in_nt,in_beta !Used to load the parameters of the types_params_modules
    logical :: in_latt_init,in_number_MC_steps,in_measurement_step !Parameters used to get parameters of the Monte Carlo evolution
    logical :: in_mu,in_nu,in_rho,in_sigma !Parameters used on the computation of the correlation function
    integer :: mu,nu,rho,sigma

    private
    public :: read_xml
    public :: mu,nu,rho,sigma
    
    contains
    !Empty routine
    subroutine startDocument_handler
    end subroutine startDocument_handler

    subroutine startElement_handler(namespaceURI, localname, name, atts)
        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name
        type(dictionary_t), intent(in) :: atts
    
        select case (name)
            case('lattice')
                print*, "Entered lattice"
            case('nx')
                print *, "Enetered nx"
                in_nx = .true.
            case('ny')
                in_ny = .true.
            case('nz')
                in_nz = .true.
            case('nt')
                in_nt = .true.
            case('beta')
                in_beta = .true.
            case('mu')
                in_mu = .true.
            case('nu')
                in_nu = .true.
            case('rho')
                in_rho = .true.
            case('sigma')
                in_sigma = .true.
        end select
      end subroutine startElement_handler
    
      subroutine endElement_handler(namespaceURI, localname, name)
        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: localname
        character(len=*), intent(in) :: name
    
        select case (name)
            case('nx')
                in_nx = .false.
            case('ny')
                in_ny = .false.
            case('nz')
                in_nz = .false.
            case('nt')
                in_nt = .false.
            case('beta')
                in_beta = .false.
            case('mu')
                in_mu = .false.
            case('nu')
                in_nu = .false.
            case('rho')
                in_rho = .false.
            case('sigma')
                in_sigma = .false.
        end select
      end subroutine endElement_handler
    
      subroutine characters_handler(chars)
        
        character(len=*), intent(in) :: chars
    
        if (in_nx) then
            read(chars,*) nx
        else if (in_ny) then
            read(chars,*) ny
        else if (in_nz) then
            read(chars,*) nz
        else if (in_nt) then
            read(chars,*) nt
        else if (in_beta) then
            read(chars,*) beta
        else if (in_mu) then
            read(chars,*) mu
        else if (in_nu) then
            read(chars,*) nu
        else if (in_rho) then
            read(chars,*) rho
        else if (in_sigma) then
            read(chars,*) sigma
        end if
      end subroutine characters_handler
    
      subroutine endDocument_handler
      end subroutine endDocument_handler
    
      subroutine read_xml(input_xml)
        character(len=1024), intent(in) :: input_xml
        type(xml_t) :: xf
        integer :: iostat
        
        !Opens the file     
        call open_xml_file(xf, input_xml, iostat)
        if (iostat/=0) then
          print*, 'error opening file'
          stop
        endif

        call parse(xf, &
        startDocument_handler = startDocument_handler, &
        startElement_handler = startElement_handler, &
        endElement_handler = endElement_handler, &
        characters_handler = characters_handler, &
        endDocument_handler = endDocument_handler &
        )
    
        call close_xml_t(xf)

      end subroutine
end module