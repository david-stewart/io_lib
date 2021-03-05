root -l<<EOF
    .x io_loadlibs.C
    ioIntList badruns { "bad_run_v1.list", false };
    ioIntMap  id_time  { "run_id.list", 1, 2, false };

    ofstream fout { "good_run_id.list" };
    fout << Form("%10s %10s %10s", "#good-id","runId","dur(sec)") << endl;
    
    int cnt{0};
    for ( auto& id : id_time.keys() ) {
        if (!badruns(id))
            fout << Form("%10i %10i %10i", cnt++, id, id_time[id]) << endl;
    }
    cout << " done " << endl;
EOF
