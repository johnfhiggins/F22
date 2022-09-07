program hello
    integer :: io 
    open(newunit=io, file="log.txt", status="new", action="write")
    write(io,*) 'Hello, Dean!'
    close(io)
end program hello