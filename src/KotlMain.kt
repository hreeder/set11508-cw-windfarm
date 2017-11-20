fun main(args : Array<String>) {
    var files: MutableList<String> = mutableListOf()

    (1..5).forEach { i -> files.add("./Scenarios/competition_$i.xml") }

    var practice_wfle = KusiakLayoutEvaluator()
    var practice_sc = WindScenario("./Scenarios/practice_0.xml")
    practice_wfle.initialize(practice_sc)

    var practice_solver = Solver(SolverSettings(
            evaluator = practice_wfle,
            max_hours_to_run = 1
    ))
    practice_solver.run_cw()

    files.forEach { filename ->
        var wfle = KusiakLayoutEvaluator()
        var sc = WindScenario(filename)
        wfle.initialize(sc)

        var solver = Solver(SolverSettings(
            evaluator = wfle,
            scenario_filename = filename,
            max_hours_to_run = 48 / files.size
        ))

        solver.run_cw()
    }
}