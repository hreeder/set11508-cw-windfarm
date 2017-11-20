public class StartOptimiser {
  public static int nSc;
  public static void main(String argv[]) throws Exception {
    try {
        nSc = Integer.parseInt(argv[0]);
    } catch (ArrayIndexOutOfBoundsException ex) {
        // Use a default of scenario 0
        nSc = 0;
    }

	KusiakLayoutEvaluator eval = new KusiakLayoutEvaluator();
	WindScenario sc = new WindScenario("./Scenarios/practice_"+nSc+".xml");
	eval.initialize(sc);
	Solver algorithm = new Solver(eval, nSc);
	algorithm.run_cw();
  }
}
