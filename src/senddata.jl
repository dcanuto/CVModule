function senddata(system::CVSystem,n::Int64,sender::Socket)
    ZMQ.send(sender, ZMQ.Message(string(system.branches.P[1][n+1,1])))
    msg = ZMQ.recv(sender);
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
