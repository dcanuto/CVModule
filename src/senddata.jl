function senddata(system::CVSystem,n::Int64,sender::Socket)
    # convert data to scientific notation to maintain precision
    a,b=sc(system.branches.Q[14][n+1,11]);
    c,d=sc(system.branches.P[14][n+1,11]);
    e,f=sc(system.liver.Q[n+1,1]);
    g,h=sc(system.liver.P[n+1,1]);

    # write data to IO buffer
    io = IOBuffer();
    write(io,string(a),string(" "),string(b),string(" "),string(c),string(" "),
        string(d),string(" "),string(e),string(" "),string(f),string(" "),
        string(g),string(" "),string(h),string(" "));

    # send data, wait for confirmation of receipt
    ZMQ.send(sender, Message(io))
    msg = ZMQ.recv(sender);

    # print confirmation every 100 time steps
    if n % 100 == 0
        str = unsafe_string(msg);
        print("Received ")
        print(str)
        print(", length: ")
        println(length(str))
        println("Last character: ")
        print("Converted message to float: ")
        println(parse(Float64, str[1:end-1]))
    end
end
