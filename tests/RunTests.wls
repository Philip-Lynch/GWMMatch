(* RunTests.wls *)
(* Run all GWMMatch unit tests from the command line or from a notebook. *)
(*                                                                        *)
(* Command line:                                                          *)
(*   wolframscript -file RunTests.wls                                     *)
(*                                                                        *)
(* Inside a notebook:                                                     *)
(*   Get["path/to/tests/RunTests.wls"]                                    *)

$testDir = DirectoryName[$InputFileName];

files = FileNames["Test*.wlt", $testDir];

Print["Running ", Length[files], " test suite(s)...\n"];

allPassed = True;

Do[
  Print["--- ", FileNameTake[f], " ---"];
  report = TestReport[f];
  Print["  Passed:  ", report["TestsSucceededCount"]];
  Print["  Failed:  ", report["TestsFailedCount"]];
  Print["  Errored: ", report["TestsFailedWithErrorsCount"]];
  If[report["TestsFailedCount"] > 0 || report["TestsFailedWithErrorsCount"] > 0,
    allPassed = False;
    Print["  FAILURES:"];
    Do[
      If[tr["Outcome"] =!= "Success",
        Print["    ", tr["TestID"], ": ", tr["Outcome"]];
        Print["      Expected: ", tr["ExpectedOutput"]];
        Print["      Actual:   ", tr["ActualOutput"]];
      ],
      {tr, report["TestResults"]}
    ]
  ];
  Print[],
  {f, files}
];

If[allPassed,
  Print["All tests passed."];
  Exit[0],
  Print["Some tests FAILED."];
  Exit[1]
]
