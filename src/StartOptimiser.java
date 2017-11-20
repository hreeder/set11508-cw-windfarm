public class StartOptimiser {
  public static int nSc;
  public static void main(String argv[]) throws Exception {
    String scenario_filename;
    try {
        scenario_filename = argv[0];
    } catch (ArrayIndexOutOfBoundsException ex) {
        // Use a default of scenario 0
        scenario_filename = "./Scenarios/practice_0.xml";
    }

	KusiakLayoutEvaluator eval = new KusiakLayoutEvaluator();
	WindScenario sc = new WindScenario(scenario_filename);
	eval.initialize(sc);
	Solver algorithm = new Solver(eval, scenario_filename);
	algorithm.run_cw();
  }
}
