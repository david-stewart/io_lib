root -l <<EOF
    for (auto& name : vector<const char*>{"Class","OptMap","_fnc","_fmt"}) {
        cout << name << endl;
        .L Form("io_code/io%s.cc+",name)
    }
    .L 
   // .L io_code/ioClass.cc+
   // .L io_code/ioOptMap.cc+
   // .L io_code/io_fnc.cc+
   // .L io_code/io_fmt.cc+
EOF
