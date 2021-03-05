root -l<<EOF
    .x io_loadlibs.C
    // the cut for 600 seconds doesn't appear to be used in these runs,
    // add the cut for now
    ioIntMap bad_runs { "bad_run.list", 0, 0, false };
    ioIntMap id_time  { "run_id.list",  1, 2, false };
    for (auto& key : id_time.keys()) {
        if (id_time[key] < 600 && !bad_runs.has(key)) {
            bad_runs[key] = 0;
            cout << "added bad run " << key << " with time " << id_time[key] << endl;
        }
    }
    ofstream fout {"bad_run_v1.list"};
    fout << "# Update to bad run list to fix including runs with < 600 second" << endl;
    for (auto& key : bad_runs.keys()) fout << key << endl;
EOF
